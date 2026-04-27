function [F, Fbody, Tbody, info] = thruster_ft_allocation(F_des, T_des, P, health, scale)
% THRUSTER_FT_ALLOCATION  Fault-tolerant non-negative allocation of M
% reaction thrusters to realize a commanded body force F_des and torque
% T_des, with per-thruster health flag and thrust-degradation scale.
%
% Inputs:
%   F_des  : 3x1 desired body force [N]
%   T_des  : 3x1 desired body torque [Nm] (typically 0 for "orbit control
%            without attitude maneuver" + small correction from attitude
%            controller for stabilisation)
%   P      : satellite param struct (uses P.thr.dirs, P.thr.pos, P.thr.Fmax)
%   health : 1xM logical (true = thruster operable)
%   scale  : 1xM in (0,1] thrust scale per thruster (1 = nominal,
%            <1 = degraded). Combined with health.
%
% Outputs:
%   F      : Mx1 commanded per-thruster thrust [N], 0..Fmax_eff
%   Fbody  : 3x1 actual delivered force on body [N]
%   Tbody  : 3x1 actual delivered torque on body [Nm]
%   info   : struct .feasible, .residual, .nfire
%
% Allocation: solve  min || A F - b ||^2 + lambda*||F||^2
%             s.t.   0 <= F_i <= Fmax_eff_i
% where A=[dirs; cross(pos,dirs)], b=[F_des; T_des].
% We use lsqnonneg by augmenting with a damping term.

M = P.thr.M;
if nargin<4 || isempty(health), health = true(1,M); end
if nargin<5 || isempty(scale),  scale  = ones(1,M);  end
Fmax_eff = P.thr.Fmax * scale .* double(health);

A_force  = P.thr.dirs;
A_torque = cross(P.thr.pos, P.thr.dirs, 1);
A = [A_force; A_torque];
% Zero-out failed columns
A(:, ~health | scale<=0) = 0;

b = [F_des(:); T_des(:)];

lambda = 1e-3;                          % control-effort regulariser
Aaug = [A; lambda*eye(M)];
baug = [b; zeros(M,1)];

warn_s = warning('off','MATLAB:lsqnonneg:IterationCountExceeded');
F = lsqnonneg(Aaug, baug);
warning(warn_s);

% Saturate at max effective thrust
F = min(F, Fmax_eff(:));

Fbody = A_force * F;
Tbody = A_torque * F;

info.residual = norm([Fbody;Tbody] - b);
info.feasible = info.residual < 0.05*max(norm(b),1e-3) || norm(b)<1e-9;
info.nfire    = nnz(F > 1e-6);
end
