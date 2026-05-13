function out = sim_wheel_failure(P, fail_case, opts)
% SIM_WHEEL_FAILURE  Attitude-hold simulation with reaction-wheel faults.
%
% fail_case:
%   'none'   - all 4 wheels healthy
%   'one'    - wheel #1 failed
%   'two'    - wheels #1 and #2 failed (triggers wheel+thruster co-control)
%
% Optional opts fields:
%   .health          - 1xN logical wheel health override
%   .use_diagnosis   - true enables online diagnosis (default true)
%   .fault_time_s    - time when true faults become active (default 1.5 s)
%   .double_strategy - 'assist_immediate' (default), 'two_wheel_only',
%                      or 'wheel_first'
%   .err_threshold_deg - pointing threshold for wheel_first [deg]
%   .window_s        - rolling evaluation window for wheel_first [s]
%
% Output struct fields: t, q, w, err_deg, Twheels, mode

if nargin<2, fail_case = 'none'; end
if nargin<3 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'double_strategy'), opts.double_strategy = 'assist_immediate'; end
if ~isfield(opts, 'err_threshold_deg'), opts.err_threshold_deg = 0.1; end
if ~isfield(opts, 'window_s'), opts.window_s = 30; end
if ~isfield(opts, 'use_diagnosis'), opts.use_diagnosis = true; end
if ~isfield(opts, 'fault_time_s'), opts.fault_time_s = 1.5; end

dt = 0.1;  T_end = 200;
N  = round(T_end/dt);
t  = (0:N-1).'*dt;

% Initial condition: 5 deg offset about Z, 0.3 deg/s rates
q  = [cos(deg2rad(2.5)); 0; 0; sin(deg2rad(2.5))];
w  = deg2rad([0.3; -0.2; 0.1]);
q_cmd = [1;0;0;0];

% True health vector per case. The controller receives estimated health
% when diagnosis is enabled; this true vector only affects actual dynamics.
switch fail_case
    case 'none', true_health_final = true(1,P.wheel.N);
    case 'one',  true_health_final = [false true true true];
    case 'two',  true_health_final = [false false true true];
    case 'custom', true_health_final = true(1,P.wheel.N);
    otherwise, error('Unknown fail_case');
end
if isfield(opts, 'health') && ~isempty(opts.health)
    true_health_final = logical(opts.health);
end
if numel(true_health_final) ~= P.wheel.N
    error('Wheel health vector must have length P.wheel.N');
end
if sum(true_health_final) > 2
    % The strategy choice is only relevant after loss of three-axis wheel
    % authority. Keep nominal and one-wheel-failure behavior unchanged.
    opts.double_strategy = 'assist_immediate';
end

% Log buffers
q_log = zeros(4,N);  w_log = zeros(3,N);
err_log = zeros(1,N); Tw_log = zeros(P.wheel.N,N);
mode_log = zeros(1,N);
gamma_log = zeros(P.wheel.N + P.thr.M, N);
wheel_health_est_log = false(P.wheel.N, N);
wheel_health_true_log = false(P.wheel.N, N);
fault_alarm_log = false(1, N);
diag_res_log = zeros(6, N);
assist_enabled = false;
assist_time = NaN;
window_N = max(1, round(opts.window_s/dt));
diag_state = [];
diag_out = local_default_diag(P);
wheel_health_est = true(1, P.wheel.N);
thruster_health_est = true(1, P.thr.M);
thruster_scale_est = ones(1, P.thr.M);

for k = 1:N
    true_health = true(1, P.wheel.N);
    if t(k) >= opts.fault_time_s
        true_health = true_health_final;
    end
    wheel_scale_true = double(true_health(:));

    % Sensor measurements
    [qm, wm] = sensor_model(q, w, P.sensor);
    % Error quaternion (inv(cmd) * meas) -> rotation FROM cmd TO meas
    qe = qmult(qinv(q_cmd), qm);
    we = wm;        % target rate = 0

    % Wheel allocation
    if ~opts.use_diagnosis
        wheel_health_est = true_health;
    end
    [~, Tw, info] = wheel_attitude_controller(qe, we, P, wheel_health_est);
    Tw_actual = Tw .* wheel_scale_true;
    Tw_body = P.wheel.axes * Tw_actual;

    % Co-control: if the wheel subsystem cannot reach 3 axes, either keep
    % using the two-wheel projection for evaluation or add thruster couples.
    T_total = Tw_body;
    mode = 1;                 % 1 = wheels only
    if ~info.feasible
        switch opts.double_strategy
            case 'assist_immediate'
                assist_enabled = true;
                if isnan(assist_time), assist_time = t(k); end
            case 'two_wheel_only'
                assist_enabled = false;
            case 'wheel_first'
                if ~assist_enabled && k > window_N
                    recent_err = err_log((k-window_N):(k-1));
                    if mean(recent_err) > opts.err_threshold_deg
                        assist_enabled = true;
                        assist_time = t(k);
                    end
                end
            otherwise
                error('Unknown double_strategy');
        end
    end

    if ~info.feasible && assist_enabled
        T_missing = -P.ctrl.Kp_att*sign_q(qe) - P.ctrl.Kd_att*we - Tw_body;
        T_missing = max(min(T_missing, 0.05), -0.05);
        [F_thr, Fb, ~, ~] = thruster_ft_allocation([0;0;0], T_missing, P, ...
                                                thruster_health_est, thruster_scale_est);
        F_thr_actual = F_thr;
        Tb = cross(P.thr.pos, P.thr.dirs, 1) * F_thr_actual;
        T_total = Tw_body + Tb;
        % Fb is near-zero for torque-only couples; ignore its orbit effect
        mode = 2;             % 2 = wheel + thruster co-control
        Fb_unused = Fb; %#ok<NASGU>
    end

    % Disturbance torque (small constant + sinusoid)
    Td = 1e-5*[sin(0.01*t(k)); cos(0.01*t(k)); 0.5];

    % Propagate
    [q, w] = attitude_dynamics(q, w, T_total, Td, P.J, P.Jinv, dt);

    if mode ~= 2
        F_thr = zeros(P.thr.M, 1);
        F_thr_actual = zeros(P.thr.M, 1);
        Fb = zeros(3, 1);
    end
    accel_meas = local_accel_measurement(Fb / P.mass, P.sensor);
    if opts.use_diagnosis
        u_diag = [Tw; F_thr];
        actuator_meas = [Tw_actual; F_thr_actual];
        [diag_out, diag_state] = actuator_fault_diagnosis(u_diag, wm, accel_meas, ...
            t(k), P, diag_state, dt, actuator_meas);
        wheel_health_est = diag_out.wheel_health_est;
        thruster_health_est = diag_out.thruster_health_est;
        thruster_scale_est = diag_out.thruster_scale_est;
    else
        diag_out = local_truth_diag(P, true_health, true(1, P.thr.M), ones(1, P.thr.M));
    end

    % Log
    q_log(:,k)  = q;
    w_log(:,k)  = w;
    err_log(k)  = 2*acos(min(1,abs(q(1))))*180/pi;
    Tw_log(:,k) = Tw;
    mode_log(k) = mode;
    gamma_log(:,k) = diag_out.gamma_hat;
    wheel_health_est_log(:,k) = diag_out.wheel_health_est(:);
    wheel_health_true_log(:,k) = true_health(:);
    fault_alarm_log(k) = diag_out.fault_alarm;
    diag_res_log(:,k) = diag_out.res_vec;
end

out.t = t; out.q = q_log; out.w = w_log;
out.err_deg = err_log; out.Twheels = Tw_log; out.mode = mode_log;
out.fail_case = fail_case;
out.strategy = opts.double_strategy;
out.assist_time = assist_time;
out.assist_used = any(mode_log == 2);
out.health = true_health_final;
out.true_wheel_health = wheel_health_true_log;
out.wheel_health_est = wheel_health_est_log;
out.gamma_hat = gamma_log;
out.fault_alarm = fault_alarm_log;
out.detect_time = diag_out.detect_time;
out.diag_residual = diag_res_log;
out.use_diagnosis = opts.use_diagnosis;
out.fault_time_s = opts.fault_time_s;
idx = t >= (t(end)-20);
out.steady_mean_deg = mean(err_log(idx));
out.steady_max_deg = max(err_log(idx));
end

% --- helpers -----------------------------------------------------------
function qi = qinv(q)
qi = [q(1); -q(2:4)] / (q.'*q);
end
function q = qmult(a,b)
s1=a(1); v1=a(2:4); s2=b(1); v2=b(2:4);
q = [s1*s2 - v1.'*v2; s1*v2 + s2*v1 + cross(v1,v2)];
end
function v = sign_q(qe)
v = qe(2:4); if qe(1)<0, v=-v; end
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
