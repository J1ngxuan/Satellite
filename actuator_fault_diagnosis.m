function [diag_out, state] = actuator_fault_diagnosis(u_cmd, omega_meas, accel_meas, t, P, state, dt, actuator_meas)
% ACTUATOR_FAULT_DIAGNOSIS  Online RLS diagnosis for wheels and thrusters.
%
% The estimator identifies a continuous actuator health vector gamma:
%   gamma(1:P.wheel.N)       : reaction-wheel health
%   gamma(P.wheel.N+1:end)   : thruster thrust scale
%
% Inputs use physical actuator commands, u_cmd = [wheel_torque; thrust].
% The regression maps those commands into measured body torque and force.
% Optional actuator_meas supplies direct actuator telemetry in the same
% units as u_cmd; when present it improves wheel fault observability.

if nargin < 6 || isempty(state)
    state = local_init_state(P);
end
if nargin < 7 || isempty(dt)
    if isfield(state, 'dt')
        dt = state.dt;
    else
        dt = 0.1;
    end
end
state.dt = dt;
if nargin < 8
    actuator_meas = [];
end

u_cmd = u_cmd(:);
omega_meas = omega_meas(:);
accel_meas = accel_meas(:);

num_wheel = P.wheel.N;
num_thr = P.thr.M;
num_act = num_wheel + num_thr;
if numel(u_cmd) ~= num_act
    error('actuator_fault_diagnosis:CommandSize', ...
        'u_cmd must have length %d.', num_act);
end

[u_filt, omega_filt, accel_filt, dot_omega_est, state.filter] = ...
    local_filter_signals(u_cmd, omega_meas, accel_meas, dt, P, state.filter);

B_rw = P.wheel.axes;
B_thr = cross(P.thr.pos, P.thr.dirs, 1);
A_thr = P.thr.dirs;

M_meas = P.J * dot_omega_est + cross(omega_filt, P.J * omega_filt);
F_meas = P.mass * accel_filt;
Y = [M_meas; F_meas];
H_mat = [B_rw, B_thr; zeros(3, num_wheel), A_thr] * diag(u_filt);

res_for_alarm = Y - H_mat * ones(num_act, 1);
if t > local_diag_param(P, 'start_ignore_s', 0.5) && ~state.is_faulted
    torque_alarm = norm(res_for_alarm(1:3)) > local_diag_param(P, 'res_torque_threshold', 0.01);
    force_alarm = norm(res_for_alarm(4:6)) > local_diag_param(P, 'res_force_threshold', 0.01);
    state.is_faulted = torque_alarm || force_alarm;
    if state.is_faulted && isnan(state.detect_time)
        state.detect_time = t;
    end
end

Y_comp = Y;
lambda = 0.998;
if state.is_faulted
    if state.timer == 0
        delta_theta = pinv(H_mat) * res_for_alarm;
        state.theta_hat = ones(num_act, 1) + delta_theta;
        state.P = eye(num_act);
        state.P(1:num_wheel, 1:num_wheel) = eye(num_wheel) * 2000;
        state.P(num_wheel+1:end, num_wheel+1:end) = eye(num_thr) * 1000;
        lambda = 0.80;
        state.timer = 1;
    else
        res_est = Y - H_mat * state.theta_hat;
        res_norm_est = norm(res_est);
        Ki = 0.3;
        state.res_integral = state.res_integral + res_est * dt;
        Y_comp = Y + Ki * state.res_integral;
        lambda = 0.82 + (0.999 - 0.82) * exp(-150 * res_norm_est^2);
        lambda = max(min(lambda, 0.999), 0.80);
        state.timer = state.timer + 1;
    end
end

if ~isempty(actuator_meas)
    actuator_meas = actuator_meas(:);
    active_cmd = abs(u_cmd(1:num_wheel)) > max(1e-4, 0.02 * P.wheel.Tmax);
    if any(active_cmd)
        ratio = abs(actuator_meas(1:num_wheel)) ./ max(abs(u_cmd(1:num_wheel)), eps);
        ratio = max(min(ratio, 1.2), 0);
        blend = 0.35;
        idx = find(active_cmd);
        state.wheel_gamma_direct(idx) = ...
            (1 - blend) * state.wheel_gamma_direct(idx) + blend * ratio(idx);
        state.theta_hat(1:num_wheel) = state.wheel_gamma_direct;
        failed_direct = any(ratio(active_cmd) < local_diag_param(P, 'wheel_binary_threshold', 0.5));
        if t > local_diag_param(P, 'start_ignore_s', 0.5) && failed_direct && ~state.is_faulted
            state.is_faulted = true;
            state.detect_time = t;
            state.timer = 1;
        end
    end
end

theta_internal = state.theta_hat;
if state.timer > max(1, round(0.25 / dt))
    theta_internal(1:num_wheel) = double(state.theta_hat(1:num_wheel) > ...
        local_diag_param(P, 'wheel_binary_threshold', 0.5));
end

error_val = Y_comp - H_mat * theta_internal;
W_mat = diag([1; 1; 1; 50; 50; 50]);
PHt = state.P * H_mat';
S = lambda * eye(6) + H_mat * PHt * W_mat + 1e-9 * eye(6);
K = (PHt * W_mat) * pinv(S);
state.theta_hat = state.theta_hat + K * error_val;
state.P = (state.P - K * H_mat * state.P) / lambda;
state.theta_hat = max(min(state.theta_hat, local_diag_param(P, 'gamma_max', 1.2)), 0);
state.theta_hat(1:num_wheel) = state.wheel_gamma_direct;

gamma_hat = state.theta_hat;
wheel_health_est = gamma_hat(1:num_wheel).' >= local_diag_param(P, 'wheel_binary_threshold', 0.5);
thruster_scale_est = max(min(gamma_hat(num_wheel+1:end).', 1), 0);
thruster_health_est = thruster_scale_est > local_diag_param(P, 'thruster_health_threshold', 0.1);

diag_out.gamma_hat = gamma_hat;
diag_out.fault_alarm = state.is_faulted;
diag_out.detect_time = state.detect_time;
diag_out.P_diag = diag(state.P);
diag_out.res_vec = res_for_alarm;
diag_out.error_val = error_val;
diag_out.lambda = lambda;
diag_out.u_filt = u_filt;
diag_out.omega_filt = omega_filt;
diag_out.accel_filt = accel_filt;
diag_out.dot_omega_est = dot_omega_est;
diag_out.wheel_health_est = wheel_health_est;
diag_out.thruster_health_est = thruster_health_est;
diag_out.thruster_scale_est = thruster_scale_est;
end

% -------------------------------------------------------------------------
function state = local_init_state(P)
num_wheel = P.wheel.N;
num_thr = P.thr.M;
num_act = num_wheel + num_thr;
state.P = eye(num_act);
state.P(1:num_wheel, 1:num_wheel) = eye(num_wheel) * 0.1;
state.P(num_wheel+1:end, num_wheel+1:end) = eye(num_thr) * 0.01;
state.theta_hat = ones(num_act, 1);
state.wheel_gamma_direct = ones(num_wheel, 1);
state.is_faulted = false;
state.timer = 0;
state.detect_time = NaN;
state.res_integral = zeros(6, 1);
state.filter = struct();
state.dt = 0.1;
end

% -------------------------------------------------------------------------
function [u_filt, omega_filt, accel_filt, dot_omega_est, filter] = ...
    local_filter_signals(u_cmd, omega_meas, accel_meas, dt, P, filter)
if isempty(fieldnames(filter))
    filter.om_x1 = omega_meas;
    filter.om_x2 = zeros(3, 1);
    filter.u_x1 = u_cmd;
    filter.u_x2 = zeros(numel(u_cmd), 1);
    filter.ac_x1 = accel_meas;
    filter.ac_x2 = zeros(3, 1);
end

omega_n = local_diag_param(P, 'filter_omega_n', 8);
zeta = local_diag_param(P, 'filter_zeta', 0.9);

[filter.om_x1, filter.om_x2] = local_second_order_step( ...
    filter.om_x1, filter.om_x2, omega_meas, dt, omega_n, zeta);
[filter.u_x1, filter.u_x2] = local_second_order_step( ...
    filter.u_x1, filter.u_x2, u_cmd, dt, omega_n, zeta);
[filter.ac_x1, filter.ac_x2] = local_second_order_step( ...
    filter.ac_x1, filter.ac_x2, accel_meas, dt, omega_n, zeta);

omega_filt = filter.om_x1;
dot_omega_est = filter.om_x2;
u_filt = filter.u_x1;
accel_filt = filter.ac_x1;
end

% -------------------------------------------------------------------------
function [x1, x2] = local_second_order_step(x1, x2, signal, dt, omega_n, zeta)
err = x1 - signal;
dx1 = x2;
dx2 = -omega_n^2 * err - 2 * zeta * omega_n * x2;
x1 = x1 + dx1 * dt;
x2 = x2 + dx2 * dt;
end

% -------------------------------------------------------------------------
function value = local_diag_param(P, name, default_value)
if isfield(P, 'diag') && isfield(P.diag, name) && ~isempty(P.diag.(name))
    value = P.diag.(name);
else
    value = default_value;
end
end
