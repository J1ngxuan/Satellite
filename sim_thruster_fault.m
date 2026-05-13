function out = sim_thruster_fault(P, mode, fault_list, dv_dir_body, opts)
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
%   opts        : optional struct, .use_diagnosis (default true),
%                 .fault_time_s (default 1.5; set 0 for initial fault)
%
% Outputs struct fields:
%   t, q, w, err_deg, dv_achieved [m/s], dv_target [m/s], mode,
%   F_hist (MxN), Tw_hist (NwxN), residual_force, residual_torque

if nargin<2, mode = 'nominal'; end
if nargin<3, fault_list = []; end
if nargin<4, dv_dir_body = [1;0;0]; end
if nargin<5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'use_diagnosis'), opts.use_diagnosis = true; end
if ~isfield(opts, 'fault_time_s'), opts.fault_time_s = 1.5; end

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
true_scale_final  = ones(1,M);
switch mode
    case 'failure'
        true_scale_final(fault_list) = 0;
    case 'degrade'
        true_scale_final(fault_list)  = 0.4;   % 60% thrust reduction
    case 'nominal'   % no-op
    otherwise, error('Unknown mode');
end
wheel_health_est = true(1,Nw);                 % wheels assumed healthy
thruster_health_est = true(1,M);
thruster_scale_est = ones(1,M);
diag_state = [];
diag_out = local_default_diag(P);

% Logs
q_log = zeros(4,N); w_log = zeros(3,N);
err_log = zeros(1,N);
F_hist  = zeros(M,N);
Tw_hist = zeros(Nw,N);
res_f   = zeros(1,N); res_t = zeros(1,N);
dv_log  = zeros(3,N);
gamma_log = zeros(Nw + M,N);
thr_health_est_log = false(M,N);
thr_scale_est_log = zeros(M,N);
thr_scale_true_log = zeros(M,N);
fault_alarm_log = false(1,N);
diag_res_log = zeros(6,N);

for k = 1:N
    true_scale = ones(1,M);
    if t(k) >= opts.fault_time_s
        true_scale = true_scale_final;
    end
    true_health = true_scale > 0;

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
    if ~opts.use_diagnosis
        thruster_health_est = true_health;
        thruster_scale_est = true_scale;
    end
    [F, ~, Tb_pred, ~] = thruster_ft_allocation(F_body_cmd, [0;0;0], P, ...
        thruster_health_est, thruster_scale_est);
    F_actual = F .* true_scale(:);
    Fb = P.thr.dirs * F_actual;
    Tb_thr = cross(P.thr.pos, P.thr.dirs, 1) * F_actual;

    % Wheel attitude loop with feedforward cancellation of Tb_thr.
    % Wheels supply (-Tb_thr) plus PD action, so the net body torque
    % from wheels+thrusters tracks the PD command.
    [Tw_body, Tw, ~] = wheel_attitude_controller(qe, we, P, wheel_health_est, Tb_pred);

    % Propagate attitude (control = wheels + thrusters; Td = environment)
    Td = 1e-5*[sin(0.01*t(k)); cos(0.01*t(k)); 0.2];
    [q, w] = attitude_dynamics(q, w, Tw_body + Tb_thr, Td, P.J, P.Jinv, dt);

    % Accumulate dv in body frame (small-attitude approx)
    dv_acc = dv_acc + (Fb/P.mass)*dt;

    accel_meas = local_accel_measurement(Fb / P.mass, P.sensor);
    if opts.use_diagnosis
        u_diag = [Tw; F];
        actuator_meas = [Tw; F_actual];
        [diag_out, diag_state] = actuator_fault_diagnosis(u_diag, wm, accel_meas, ...
            t(k), P, diag_state, dt, actuator_meas);
        wheel_health_est = diag_out.wheel_health_est;
        thruster_health_est = diag_out.thruster_health_est;
        thruster_scale_est = diag_out.thruster_scale_est;
    else
        diag_out = local_truth_diag(P, wheel_health_est, true_health, true_scale);
    end

    q_log(:,k)=q; w_log(:,k)=w;
    err_log(k)=2*acos(min(1,abs(q(1))))*180/pi;
    F_hist(:,k)=F;
    Tw_hist(:,k)=Tw;
    res_f(k)=norm(Fb-F_body_cmd);
    res_t(k)=norm(Tb_thr);
    dv_log(:,k)=dv_acc;
    gamma_log(:,k)=diag_out.gamma_hat;
    thr_health_est_log(:,k)=diag_out.thruster_health_est(:);
    thr_scale_est_log(:,k)=diag_out.thruster_scale_est(:);
    thr_scale_true_log(:,k)=true_scale(:);
    fault_alarm_log(k)=diag_out.fault_alarm;
    diag_res_log(:,k)=diag_out.res_vec;
end

out.t = t; out.q = q_log; out.w = w_log; out.err_deg = err_log;
out.F_hist = F_hist; out.Tw_hist = Tw_hist;
out.residual_force = res_f; out.residual_torque = res_t;
out.dv_achieved = dv_log; out.dv_target = dv_target_mag*dv_dir_body(:);
out.mode = mode; out.fault_list = fault_list;
out.true_thruster_scale = thr_scale_true_log;
out.thruster_health_est = thr_health_est_log;
out.thruster_scale_est = thr_scale_est_log;
out.gamma_hat = gamma_log;
out.fault_alarm = fault_alarm_log;
out.detect_time = diag_out.detect_time;
out.diag_residual = diag_res_log;
out.use_diagnosis = opts.use_diagnosis;
out.fault_time_s = opts.fault_time_s;
end

% helpers
function qi = qinv(q), qi=[q(1);-q(2:4)]/(q.'*q); end
function q = qmult(a,b)
s1=a(1); v1=a(2:4); s2=b(1); v2=b(2:4);
q=[s1*s2-v1.'*v2; s1*v2+s2*v1+cross(v1,v2)];
end
function accel_meas = local_accel_measurement(accel_true, sensor)
accel_meas = accel_true(:);
if isfield(sensor, 'accel_bias')
    accel_meas = accel_meas + sensor.accel_bias(:);
end
if isfield(sensor, 'accel_sigma')
    accel_meas = accel_meas + sensor.accel_sigma * randn(3, 1);
end
end
function diag_out = local_default_diag(P)
diag_out = local_truth_diag(P, true(1, P.wheel.N), true(1, P.thr.M), ones(1, P.thr.M));
diag_out.detect_time = NaN;
end
function diag_out = local_truth_diag(~, wheel_health, thruster_health, thruster_scale)
gamma = [double(wheel_health(:)); thruster_scale(:)];
diag_out.gamma_hat = gamma;
diag_out.fault_alarm = any(~wheel_health) || any(~thruster_health) || any(thruster_scale < 1);
diag_out.detect_time = NaN;
diag_out.res_vec = zeros(6, 1);
diag_out.wheel_health_est = logical(wheel_health);
diag_out.thruster_health_est = logical(thruster_health);
diag_out.thruster_scale_est = thruster_scale;
end
