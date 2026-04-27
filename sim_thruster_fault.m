function out = sim_thruster_fault(P, mode, fault_list, dv_dir_body)
% SIM_THRUSTER_FAULT  Orbit-maneuver simulation while holding attitude,
% under thruster faults.  Hybrid actuation: reaction wheels run the
% attitude loop and absorb residual torques the thrusters cannot null;
% thrusters deliver the commanded delta-v with a small torque hint to
% off-load the wheels.  This reflects standard small-sat practice and is
% required when faulted nozzles produce uncancellable torque (e.g. a
% surviving +X nozzle off-axis from c.o.m.).
%
% Inputs:
%   P           : satellite parameter struct (thruster config selected)
%   mode        : 'degrade'  (multiple thrusters reduced-thrust)
%                 'failure'  (multiple thrusters completely failed)
%                 'nominal'  (no faults, baseline)
%   fault_list  : indices of thrusters affected (>=2 recommended)
%   dv_dir_body : 3x1 unit vector target dv direction in body frame
%
% Outputs struct fields:
%   t, q, w, err_deg, dv_achieved [m/s], dv_target [m/s], mode,
%   F_hist (MxN), Tw_hist (NwxN), residual_force, residual_torque

if nargin<2, mode = 'nominal'; end
if nargin<3, fault_list = []; end
if nargin<4, dv_dir_body = [1;0;0]; end

dt = 0.1;  T_burn = 40;  T_end = 60;
N  = round(T_end/dt);
t  = (0:N-1).'*dt;

% Target total delta-v and constant commanded thrust acceleration
dv_target_mag = 0.3;                           % [m/s]
acc_cmd_mag   = dv_target_mag / T_burn;        % [m/s^2]
F_cmd_mag     = P.mass * acc_cmd_mag;          % [N] body force

% Initial state
q  = [1;0;0;0];    w  = zeros(3,1);
q_cmd = q;
dv_acc = zeros(3,1);                           % accumulated body-frame dv

M  = P.thr.M;
Nw = P.wheel.N;
health = true(1,M);
scale  = ones(1,M);
switch mode
    case 'failure',  health(fault_list) = false;
    case 'degrade',  scale(fault_list)  = 0.4;   % 60% thrust reduction
    case 'nominal'   % no-op
    otherwise, error('Unknown mode');
end
wheel_health = true(1,Nw);                     % wheels assumed healthy

% Logs
q_log = zeros(4,N); w_log = zeros(3,N);
err_log = zeros(1,N);
F_hist  = zeros(M,N);
Tw_hist = zeros(Nw,N);
res_f   = zeros(1,N); res_t = zeros(1,N);
dv_log  = zeros(3,N);

for k = 1:N
    % Sensor measurements
    [qm, wm] = sensor_model(q, w, P.sensor);
    qe = qmult(qinv(q_cmd), qm);
    we = wm;

    % Commanded body force only during burn window
    if t(k) < T_burn
        F_body_cmd = F_cmd_mag * dv_dir_body(:);
    else
        F_body_cmd = [0;0;0];
    end

    % Thruster allocation runs first so its residual torque is known.
    % T_des thruster = 0; the NNLS allocator minimises both force and
    % torque error, but residual torque is generally non-zero under
    % fault (e.g. surviving off-axis +X nozzle generates parasitic Z
    % torque that no non-negative combination can cancel).
    [F, Fb, Tb_thr, info] = thruster_ft_allocation(F_body_cmd, [0;0;0], P, health, scale);

    % Wheel attitude loop with feedforward cancellation of Tb_thr.
    % Wheels supply (-Tb_thr) plus PD action, so the net body torque
    % from wheels+thrusters tracks the PD command.
    [Tw_body, Tw, ~] = wheel_attitude_controller(qe, we, P, wheel_health, Tb_thr);

    % Propagate attitude (control = wheels + thrusters; Td = environment)
    Td = 1e-5*[sin(0.01*t(k)); cos(0.01*t(k)); 0.2];
    [q, w] = attitude_dynamics(q, w, Tw_body + Tb_thr, Td, P.J, P.Jinv, dt);

    % Accumulate dv in body frame (small-attitude approx)
    dv_acc = dv_acc + (Fb/P.mass)*dt;

    q_log(:,k)=q; w_log(:,k)=w;
    err_log(k)=2*acos(min(1,abs(q(1))))*180/pi;
    F_hist(:,k)=F;
    Tw_hist(:,k)=Tw;
    res_f(k)=norm(Fb-F_body_cmd);
    res_t(k)=norm(Tb_thr);
    dv_log(:,k)=dv_acc;
end

out.t = t; out.q = q_log; out.w = w_log; out.err_deg = err_log;
out.F_hist = F_hist; out.Tw_hist = Tw_hist;
out.residual_force = res_f; out.residual_torque = res_t;
out.dv_achieved = dv_log; out.dv_target = dv_target_mag*dv_dir_body(:);
out.mode = mode; out.fault_list = fault_list;
end

% helpers
function qi = qinv(q), qi=[q(1);-q(2:4)]/(q.'*q); end
function q = qmult(a,b)
s1=a(1); v1=a(2:4); s2=b(1); v2=b(2:4);
q=[s1*s2-v1.'*v2; s1*v2+s2*v1+cross(v1,v2)];
end
