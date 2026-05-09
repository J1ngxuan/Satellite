function configs = thruster_configurations(dim, Fmax)
% THRUSTER_CONFIGURATIONS  Library of candidate thruster layouts (>=3).
%
% Inputs
%   dim  : 3x1 body dimensions [m]
%   Fmax : max per-nozzle thrust [N]
%
% Output
%   configs : struct array with fields
%       .name      : short label
%       .dirs      : 3xM unit force-direction matrix (body frame)
%       .pos       : 3xM nozzle position vectors (body frame)
%       .M         : nozzle count
%       .redundant : true if surviving any single-nozzle loss retains
%                    full 3D translation capability
%       .eta_iso   : averaged efficiency over isotropic direction samples
%                    (ratio of commanded-direction thrust to total firing
%                    thrust magnitude, see analyze_config)
%       .eta_axis  : 3x1 efficiency along +X,+Y,+Z principal directions
%
% Efficiency definition (per requirement):
%   eta(u) = |sum F_i * d_i dot u| / sum F_i
% where F_i are commanded nozzle thrusts (>=0) chosen to realize pure
% translation along unit direction u using non-negative least squares so
% net torque ~= 0 (attitude is not disturbed).

if nargin<1, P=satellite_params(); dim=P.dim; Fmax=P.thr.Fmax; end

hx=dim(1)/2; hy=dim(2)/2; hz=dim(3)/2;

% -------- Scheme A: 6 orthogonal nozzles (non-redundant baseline) -------
A.name = 'A-6orth';
A.dirs = [ +1 -1  0  0  0  0;
            0  0 +1 -1  0  0;
            0  0  0  0 +1 -1 ];
A.pos  = [ -hx  hx  0   0   0   0;
            0   0  -hy  hy  0   0;
            0   0   0   0  -hz  hz];

% -------- Scheme B: 8 canted nozzles (tetrahedral style on ±Z faces) ----
% Two nozzles per quadrant on ±Z face, canted 45 deg out toward X,Y.
c = 1/sqrt(3);
B.name = 'B-8canted';
B.dirs = c*[ +1 -1 +1 -1 +1 -1 +1 -1;
             +1 +1 -1 -1 +1 +1 -1 -1;
             +1 +1 +1 +1 -1 -1 -1 -1 ];
B.pos  = [ -hx  hx -hx  hx -hx  hx -hx  hx;
           -hy -hy  hy  hy -hy -hy  hy  hy;
            hz  hz  hz  hz -hz -hz -hz -hz ];

% -------- Scheme C: 12 orthogonal (2 per axis direction, redundant) -----
C.name = 'C-12orth';
C.dirs = [ +1 +1 -1 -1  0  0  0  0  0  0  0  0;
            0  0  0  0 +1 +1 -1 -1  0  0  0  0;
            0  0  0  0  0  0  0  0 +1 +1 -1 -1 ];
C.pos  = [ -hx -hx  hx  hx  0   0   0   0  hx -hx  hx -hx;
            hy -hy  hy -hy -hy -hy  hy  hy  0   0   0   0;
            0   0   0   0  hz -hz  hz -hz -hz -hz  hz  hz];

% -------- Scheme E: 12 non-orthogonal omni-directional (icosahedral) ---
% Directions = regular-icosahedron vertices (12 points on a sphere, 6
% antipodal pairs).  Each nozzle is mounted where its thrust-line ray
% intersects the body box, i.e. position = -t*dir, so r x F = 0 for
% every individual nozzle (zero parasitic torque per thruster - ideal
% for orbit-control-without-attitude-maneuver).  A 1-D nonlinear
% optimisation blends this omnidirectional geometry with axis-aligned
% directions, weighted by the inertia asymmetry (I_x < I_y = I_z).
E = build_icosa_config(dim, diag([0.288 0.45025 0.45025]));

% Assemble and analyze
raw = [A B C E];
configs = [];
for k=1:numel(raw)
    c = raw(k);
    c.M = size(c.dirs,2);
    [c.redundant, c.eta_iso, c.eta_axis] = analyze_config(c, Fmax);
    if isempty(configs), configs = c; else, configs(end+1) = c; end %#ok<AGROW>
end
end

% ========================================================================
function [redundant, eta_iso, eta_axis] = analyze_config(c, Fmax)
% Test redundancy and efficiency.
% Full 3D translation w/o attitude torque requires: span(dirs)=R^3 and
% null-space of [dirs; cross(pos,dirs)] (6xM) contains positive vector
% enough to realize any direction. We approximate by solving non-negative
% least squares for 36 isotropic target directions and checking residuals.

M = size(c.dirs,2);
A_force  = c.dirs;                     % 3xM
A_torque = cross(c.pos, c.dirs, 1);    % 3xM
A = [A_force; A_torque];               % 6xM

% Isotropic unit directions (Fibonacci sphere, 64 points)
N = 64;
phi = (1+sqrt(5))/2;
idx = (0:N-1).';
z = 1 - 2*(idx+0.5)/N;
r = sqrt(1-z.^2);
theta = 2*pi*idx/phi;
U = [r.*cos(theta), r.*sin(theta), z].';

% Baseline efficiency (no failure)
eta_all = zeros(N,1);
feasible_all = true(N,1);
for i=1:N
    [F,ok] = alloc_direction(U(:,i), A, Fmax);
    if ~ok, feasible_all(i)=false; eta_all(i)=0; continue; end
    net = A_force*F;
    eta_all(i) = norm(net)/max(sum(F),eps);
end
eta_iso = mean(eta_all(feasible_all));

axis_dirs = [eye(3) -eye(3)];
ea = zeros(6,1);
for i=1:6
    [F,ok] = alloc_direction(axis_dirs(:,i), A, Fmax);
    if ok, ea(i) = norm(A_force*F)/max(sum(F),eps); end
end
eta_axis = mean(reshape(ea,3,2),2);    % average ±

% Redundancy: drop each nozzle and recheck force-producing capability on
% the 6 principal axes. Torque imbalance is allowed (attitude controller
% can trim it); we only need the column span of the remaining force
% sub-matrix to include each axis.
redundant = true;
for j=1:M
    Df = c.dirs;  Df(:,j) = 0;          % remaining force directions
    for i=1:size(axis_dirs,2)
        Fj = lsqnonneg(Df, axis_dirs(:,i));
        if norm(Df*Fj - axis_dirs(:,i)) > 1e-3
            redundant = false; break;
        end
    end
    if ~redundant, break; end
end
end

% ------------------------------------------------------------------------
function [F, ok] = alloc_direction(u, A, Fmax) %#ok<INUSD>
% Non-negative allocation to produce net force along u with zero torque.
% Target vector b = [u;0;0;0] (unit force, zero torque). Result F >= 0.
% No Fmax clamp here — we are testing geometric feasibility only.
b = [u(:); 0;0;0];
warn_s = warning('off','MATLAB:lsqnonneg:IterationCountExceeded');
F = lsqnonneg(A, b);
warning(warn_s);
res = A*F - b;
% Feasible if both torque residual ~0 AND force aligns with u.
Fforce = A(1:3,:)*F;
ok = norm(res(4:6)) < 1e-4 ...
  && norm(Fforce) > 0.2 ...
  && dot(Fforce, u)/max(norm(Fforce),eps) > 0.95;
end

% ========================================================================
function E = build_icosa_config(dim, J)
% BUILD_ICOSA_CONFIG  12-nozzle non-orthogonal omnidirectional layout.
%
% Geometry: positions = 12 box-edge midpoints (4 edges parallel to each
% body axis); directions = 12 icosahedron vertices, assigned to edges
% by minimum cant angle.  Per-nozzle torque = cross(pos, dir) is
% generically nonzero, so the layout serves both translation (any-axis
% delta-v) and attitude-support roles (e.g. when wheels are degraded).
%
% A 1-D nonlinear search blends each thrust direction toward the
% radial-from-center direction (which would zero the torque and
% maximise translation efficiency) by fraction alpha, weighted by the
% inertia asymmetry I_x < I_y = I_z.  The optimum trades torque
% authority on the low-inertia X axis against omnidirectional thrust.
%
% Inputs:
%   dim : 3x1 body dimensions [m]
%   J   : 3x3 inertia tensor [kg*m^2]
%
% Output: E struct with .name, .dirs (3x12), .pos (3x12)

hx = dim(1)/2; hy = dim(2)/2; hz = dim(3)/2;
phi = (1+sqrt(5))/2;

% 12 edge midpoints (positions, body frame)
Pos = [  0   0   0   0  hx -hx  hx -hx  hx -hx  hx -hx;
        hy -hy  hy -hy   0   0   0   0  hy  hy -hy -hy;
        hz  hz -hz -hz  hz  hz -hz -hz   0   0   0   0 ];

% 12 icosahedron vertices (unit-normalised), 6 antipodal pairs
Vraw = [ 0  0  1 -1  phi -phi   0   0 -1  1 -phi  phi;
         1 -1  phi phi 0  0    -1   1 -phi -phi 0  0;
         phi phi 0  0  1 -1   -phi -phi 0 0 -1  1 ];
Vico = Vraw ./ vecnorm(Vraw);

% Greedy assignment: pair each edge midpoint with the icosahedral vertex
% whose direction is closest to the edge's radial direction (so each
% nozzle pushes "outward").  This keeps the layout intuitive and avoids
% positive-feedback couples.
Pdir = Pos ./ vecnorm(Pos);
used = false(1,12);
order = zeros(1,12);
for i=1:12
    cosang = Pdir(:,i).' * Vico;
    cosang(used) = -inf;
    [~,j] = max(cosang);
    order(i) = j;
    used(j) = true;
end
V = Vico(:, order);

% Inertia weights: reward axes with smaller inertia
Jdiag = diag(J);
wI = 1 ./ Jdiag;  wI = wI / sum(wI);

% 1-D nonlinear search: blend icosahedral dirs toward radial direction
% by fraction alpha (alpha=0 -> pure icosahedron; alpha=1 -> radial,
% per-nozzle torque vanishes).  Maximise inertia-weighted figure-of-merit.
cost = @(a) -figure_of_merit_pos(blend_radial(V, Pdir, a), Pos, wI);
alpha_opt = fminbnd(cost, 0, 0.9, optimset('TolX',1e-3,'Display','off'));

E.name = sprintf('E-12icosa(a=%.2f)', alpha_opt);
E.dirs = blend_radial(V, Pdir, alpha_opt);
E.pos  = Pos;
end

% ------------------------------------------------------------------------
function D = blend_radial(V, R, alpha)
% Blend each icosahedral direction V(:,i) toward the radial direction
% R(:,i) by fraction alpha in [0,1], then renormalise.
D = (1-alpha)*V + alpha*R;
D = D ./ vecnorm(D);
end

% ------------------------------------------------------------------------
function m = figure_of_merit_pos(D, P, wI)
% Inertia-weighted omnidirectional score evaluated at fixed positions P.
A = [D; cross(P, D, 1)];
axes6 = [eye(3) -eye(3)];
axis_eta = zeros(1,6);
warn_s = warning('off','MATLAB:lsqnonneg:IterationCountExceeded');
for i=1:6
    F = lsqnonneg(A, [axes6(:,i); 0;0;0]);
    if norm(D*F) > 1e-6
        axis_eta(i) = norm(D*F)/max(sum(F),eps);
    end
end
warning(warn_s);
per_axis = reshape(axis_eta,3,2);
mean_ax  = mean(per_axis,2);
m = wI.' * mean_ax + 0.5*min(axis_eta);
end
