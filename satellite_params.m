function P = satellite_params()
% SATELLITE_PARAMS  Returns struct of platform constants and interfaces.
%   P = SATELLITE_PARAMS()
%
% Output struct fields (units in SI unless noted):
%   P.mass        : platform mass [kg]
%   P.dim         : body dimensions [m] (X,Y,Z)
%   P.J           : 3x3 inertia tensor [kg*m^2] (body frame, principal)
%   P.Jinv        : inverse inertia tensor
%   P.wheel       : reaction wheel parameters
%       .Tmax     : max single-wheel torque [Nm]
%       .Hmax     : max single-wheel angular momentum [Nms]
%       .axes     : 3xN unit-vector matrix of wheel spin axes (body frame)
%       .N        : number of wheels
%   P.thr         : reaction thruster parameters
%       .Fmax     : max thrust per nozzle [N]
%       .Isp      : specific impulse [s]  (used to size mass loss term, optional)
%       .dirs     : 3xM unit-vector thrust direction matrix (body frame)
%       .pos      : 3xM nozzle position vectors (body frame, m)
%       .M        : number of nozzles
%   P.sensor      : sensor noise statistics
%       .star_sigma : 1-sigma star-tracker attitude noise [rad]
%       .gyro_sigma : 1-sigma gyro angular-rate noise [rad/s]
%       .gyro_bias  : 3x1 constant gyro bias [rad/s]
%   P.ctrl        : controller defaults (gains, deadbands)
%       .Kp_att     : proportional gain matrix attitude
%       .Kd_att     : derivative gain matrix attitude
%       .att_dead   : attitude pointing deadband [rad]
%   P.env         : environment
%       .mu       : Earth gravitational parameter [m^3/s^2]
%       .Re       : Earth equatorial radius [m]
%       .J2       : J2 zonal harmonic
%
% --- Mass and geometry --------------------------------------------------
P.mass = 30;                                    % [kg]
P.dim  = [0.350; 0.240; 0.240];                 % [m]
P.J    = diag([0.288, 0.45025, 0.45025]);       % [kg*m^2]
P.Jinv = inv(P.J);

% --- Reaction wheels (default 4 wheels in pyramid 35.26 deg from Z) -----
beta = atan(1/sqrt(2));                         % pyramid half-cone angle
P.wheel.Tmax = 0.02;                            % [Nm] max wheel torque
P.wheel.Hmax = 0.5;                             % [Nms] (representative)
P.wheel.axes = [ sin(beta)*[1 -1 -1  1];
                 sin(beta)*[1  1 -1 -1];
                 cos(beta)*[1  1  1  1] ];      % 3x4 unit vectors
P.wheel.N    = size(P.wheel.axes,2);

% --- Thrusters (default = 6 orthogonal pairs scheme C; overwritten by
%     thruster_configurations.m when running configuration trade study) --
P.thr.Fmax = 0.5;                               % [N]
P.thr.Isp  = 220;                               % [s] (cold-gas representative)
[P.thr.dirs, P.thr.pos] = local_default_thrusters(P.dim);
P.thr.M    = size(P.thr.dirs,2);

% --- Sensor noise (1-sigma) ---------------------------------------------
P.sensor.star_sigma = deg2rad(5/3600);          % 5 arcsec
P.sensor.gyro_sigma = deg2rad(0.01);            % 0.01 deg/s noise
P.sensor.gyro_bias  = deg2rad(0.005)*[1;-1;1];  % bias

% --- Controller defaults -------------------------------------------------
P.ctrl.Kp_att  = 0.06*eye(3);                   % proportional
P.ctrl.Kd_att  = 0.30*eye(3);                   % derivative
P.ctrl.att_dead = deg2rad(0.05);                % 0.05 deg deadband

% --- Environment ---------------------------------------------------------
P.env.mu = 3.986004418e14;
P.env.Re = 6378137;
P.env.J2 = 1.082626e-3;
end

% ------------------------------------------------------------------------
function [D, R] = local_default_thrusters(dim)
% Default scheme: 12 thrusters along +/-X, +/-Y, +/-Z body axes (2 each)
% positioned on opposite faces; gives full redundancy with single failure.
hx = dim(1)/2; hy = dim(2)/2; hz = dim(3)/2;
% direction vectors (body) – thrust direction is the direction the gas
% leaves the nozzle, so the FORCE on the body is opposite. We define dirs
% as the FORCE direction for clarity (so acceleration = sum F_i*dirs(:,i)).
D = [ +1 +1 -1 -1  0  0  0  0  0  0  0  0;
       0  0  0  0 +1 +1 -1 -1  0  0  0  0;
       0  0  0  0  0  0  0  0 +1 +1 -1 -1 ];
% positions: pair on each face so cross-product torques are small
R = [ -hx -hx  hx  hx  0   0   0   0   0   0   0   0;
       hy -hy  hy -hy -hy -hy  hy  hy  0   0   0   0;
       0   0   0   0  hz -hz  hz -hz -hz  hz -hz  hz];
end
