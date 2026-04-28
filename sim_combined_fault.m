function out = sim_combined_fault(P, opts)
% SIM_COMBINED_FAULT  Coupled attitude/orbit maneuver under actuator faults.
%
% The simulation reuses the existing wheel attitude controller and
% thruster force/torque allocator. Fault states are injected through health
% and thrust-scale vectors; no control-law internals are modified.
%
% Required input:
%   P : satellite parameter struct with selected thruster layout.
%
% Optional opts fields:
%   .wheel_health      1xN logical vector, true = wheel healthy
%   .thruster_health   1xM logical vector, true = thruster healthy
%   .thruster_scale    1xM thrust scale, 1 nominal, 0.4 degraded
%   .dv_dir_body       3x1 target delta-v direction in body frame
%   .dt                integration step [s], default 0.1
%   .T_burn            maneuver burn duration [s], default 40
%   .T_end             simulation duration [s], default 80
%   .dv_target_mag     target delta-v magnitude [m/s], default 0.3
%   .att_threshold_deg attitude failure threshold, default 0.1
%   .dv_threshold_rel  delta-v relative error threshold, default 0.05
%
% Output includes time histories and fault-tree summary fields.

if nargin < 2 || isempty(opts)
    opts = struct();
end

M = P.thr.M;
Nw = P.wheel.N;

opts = local_default(opts, 'wheel_health', true(1, Nw));
opts = local_default(opts, 'thruster_health', true(1, M));
opts = local_default(opts, 'thruster_scale', ones(1, M));
opts = local_default(opts, 'dv_dir_body', [1; 0; 0]);
opts = local_default(opts, 'dt', 0.1);
opts = local_default(opts, 'T_burn', 40);
opts = local_default(opts, 'T_end', 80);
opts = local_default(opts, 'dv_target_mag', 0.3);
opts = local_default(opts, 'att_threshold_deg', 0.1);
opts = local_default(opts, 'dv_threshold_rel', 0.05);
opts = local_default(opts, 'sat_window_s', 20);

wheel_health = logical(opts.wheel_health);
thruster_health = logical(opts.thruster_health);
thruster_scale = opts.thruster_scale;
dv_dir_body = opts.dv_dir_body(:);
if norm(dv_dir_body) < eps
    error('dv_dir_body must be nonzero.');
end
dv_dir_body = dv_dir_body / norm(dv_dir_body);

if numel(wheel_health) ~= Nw
    error('wheel_health must have length P.wheel.N.');
end
if numel(thruster_health) ~= M || numel(thruster_scale) ~= M
    error('thruster_health and thruster_scale must have length P.thr.M.');
end

dt = opts.dt;
T_burn = opts.T_burn;
T_end = opts.T_end;
num_steps = round(T_end / dt);
t = (0:num_steps-1).' * dt;

dv_target_mag = opts.dv_target_mag;
acc_cmd_mag = dv_target_mag / T_burn;
F_cmd_mag = P.mass * acc_cmd_mag;

% Use the wheel-failure initial condition so the combined case evaluates
% attitude recovery as well as burn-time disturbance rejection.
q = [cos(deg2rad(2.5)); 0; 0; sin(deg2rad(2.5))];
w = deg2rad([0.3; -0.2; 0.1]);
q_cmd = [1; 0; 0; 0];
dv_acc = zeros(3, 1);
wheel_axes = P.wheel.axes;
wheel_axes(:, ~wheel_health) = 0;
wheel_rank_loss = rank(wheel_axes, 1e-6) < 3;

q_log = zeros(4, num_steps);
w_log = zeros(3, num_steps);
err_log = zeros(1, num_steps);
F_hist = zeros(M, num_steps);
Tw_hist = zeros(Nw, num_steps);
dv_log = zeros(3, num_steps);
mode_log = ones(1, num_steps);
wheel_sat_log = false(1, num_steps);
thr_sat_log = false(1, num_steps);
thr_feasible_log = false(1, num_steps);
res_force_log = zeros(1, num_steps);
res_torque_log = zeros(1, num_steps);

for k = 1:num_steps
    [qm, wm] = sensor_model(q, w, P.sensor);
    qe = qmult(qinv(q_cmd), qm);
    we = wm;

    if t(k) < T_burn
        F_body_cmd = F_cmd_mag * dv_dir_body;
    else
        F_body_cmd = zeros(3, 1);
    end

    assist_active = wheel_rank_loss;
    if assist_active
        % First allocate the requested maneuver force. The projected wheel
        % response identifies the torque component that needs thruster assist.
        [~, ~, Tb_force, ~] = thruster_ft_allocation(F_body_cmd, zeros(3, 1), ...
            P, thruster_health, thruster_scale);
        [Tw_body, ~] = wheel_attitude_controller(qe, we, P, ...
            wheel_health, Tb_force);

        T_pd = -P.ctrl.Kp_att * sign_q(qe) - P.ctrl.Kd_att * we;
        T_thr_des = T_pd - (Tw_body + Tb_force);
        T_thr_des = max(min(T_thr_des, 0.05), -0.05);
        [F, Fb, Tb_thr, info_thr] = thruster_ft_allocation(F_body_cmd, ...
            T_thr_des, P, thruster_health, thruster_scale);
    else
        T_thr_des = zeros(3, 1);
        [F, Fb, Tb_thr, info_thr] = thruster_ft_allocation(F_body_cmd, ...
            T_thr_des, P, thruster_health, thruster_scale);
    end
    [Tw_body, Tw, info_wheel] = wheel_attitude_controller(qe, we, P, ...
        wheel_health, Tb_thr);

    Td = 1e-5 * [sin(0.01 * t(k)); cos(0.01 * t(k)); 0.2];
    [q, w] = attitude_dynamics(q, w, Tw_body + Tb_thr, Td, P.J, P.Jinv, dt);
    dv_acc = dv_acc + (Fb / P.mass) * dt;

    Fmax_eff = P.thr.Fmax * thruster_scale(:) .* double(thruster_health(:));
    q_log(:, k) = q;
    w_log(:, k) = w;
    err_log(k) = 2 * acos(min(1, abs(q(1)))) * 180 / pi;
    F_hist(:, k) = F;
    Tw_hist(:, k) = Tw;
    dv_log(:, k) = dv_acc;
    mode_log(k) = 1 + double(assist_active);
    wheel_sat_log(k) = info_wheel.saturated;
    active_thrusters = Fmax_eff > 0;
    thr_sat_log(k) = any(active_thrusters & F >= (Fmax_eff - 1e-7));
    thr_feasible_log(k) = info_thr.feasible;
    res_force_log(k) = norm(Fb - F_body_cmd);
    res_torque_log(k) = norm(Tb_thr - T_thr_des);
end

idx_final = t >= (t(end) - 20);
idx_prev = t >= max(0, t(end) - 40) & t < (t(end) - 20);
if ~any(idx_prev)
    idx_prev = idx_final;
end

steady_mean_deg = mean(err_log(idx_final));
steady_max_deg = max(err_log(idx_final));
prev_mean_deg = mean(err_log(idx_prev));
dv_error_rel = norm(dv_acc - dv_target_mag * dv_dir_body) / max(dv_target_mag, eps);
wheel_sat_frac = mean(wheel_sat_log(idx_final));
thruster_sat_frac = mean(thr_sat_log(idx_final));
err_not_converged = steady_mean_deg > 0.95 * max(prev_mean_deg, eps);

attitude_fail = steady_mean_deg > opts.att_threshold_deg;
orbit_fail = dv_error_rel > opts.dv_threshold_rel;
co_control_fail = (wheel_sat_frac > 0.5 || thruster_sat_frac > 0.5) && ...
    (err_not_converged || attitude_fail);

out.t = t;
out.q = q_log;
out.w = w_log;
out.err_deg = err_log;
out.F_hist = F_hist;
out.Tw_hist = Tw_hist;
out.dv_achieved = dv_log;
out.dv_target = dv_target_mag * dv_dir_body;
out.mode = mode_log;
out.wheel_health = wheel_health;
out.thruster_health = thruster_health;
out.thruster_scale = thruster_scale;
out.assist_used = any(mode_log == 2);
out.steady_mean_deg = steady_mean_deg;
out.steady_max_deg = steady_max_deg;
out.dv_error_rel = dv_error_rel;
out.wheel_sat_frac = wheel_sat_frac;
out.thruster_sat_frac = thruster_sat_frac;
out.thruster_feasible_frac = mean(thr_feasible_log(idx_final));
out.residual_force_final = mean(res_force_log(idx_final));
out.residual_torque_final = mean(res_torque_log(idx_final));
out.attitude_fail = attitude_fail;
out.orbit_fail = orbit_fail;
out.co_control_fail = co_control_fail;
out.top_event = attitude_fail || orbit_fail || co_control_fail;
out.branch = local_branch(attitude_fail, orbit_fail, co_control_fail);
end

% -------------------------------------------------------------------------
function opts = local_default(opts, name, value)
if ~isfield(opts, name) || isempty(opts.(name))
    opts.(name) = value;
end
end

% -------------------------------------------------------------------------
function qi = qinv(q)
qi = [q(1); -q(2:4)] / (q.' * q);
end

% -------------------------------------------------------------------------
function q = qmult(a, b)
s1 = a(1);
v1 = a(2:4);
s2 = b(1);
v2 = b(2:4);
q = [s1 * s2 - v1.' * v2;
     s1 * v2 + s2 * v1 + cross(v1, v2)];
end

% -------------------------------------------------------------------------
function v = sign_q(qe)
v = qe(2:4);
if qe(1) < 0
    v = -v;
end
end

% -------------------------------------------------------------------------
function branch = local_branch(attitude_fail, orbit_fail, co_control_fail)
parts = cell(1, 3);
num_parts = 0;
if attitude_fail
    num_parts = num_parts + 1;
    parts{num_parts} = 'attitude_hold';
end
if orbit_fail
    num_parts = num_parts + 1;
    parts{num_parts} = 'orbit_maneuver';
end
if co_control_fail
    num_parts = num_parts + 1;
    parts{num_parts} = 'co_control';
end
if num_parts == 0
    branch = 'none';
else
    branch = strjoin(parts(1:num_parts), '+');
end
end
