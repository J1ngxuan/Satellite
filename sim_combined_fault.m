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
%   .attitude_first_orbit true enables attitude retargeting when the
%                       requested body direction is unavailable after
%                       thruster faults, default true
%   .attitude_maneuver_s fixed pre-burn attitude maneuver time [s],
%                       default 20
%   .sat_window_s      steady-state evaluation window after burn [s],
%                       default 20
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
opts = local_default(opts, 'attitude_first_orbit', true);
opts = local_default(opts, 'attitude_maneuver_s', 20);
opts = local_default(opts, 'retarget_force_tol', 0.05);
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

dv_target_mag = opts.dv_target_mag;
acc_cmd_mag = dv_target_mag / T_burn;
F_cmd_mag = P.mass * acc_cmd_mag;
dv_dir_inertial = dv_dir_body;

burn_dir_body = dv_dir_body;
attitude_first_used = false;
retarget_feasible = true;
retarget_force_rel = 0;
if opts.attitude_first_orbit
    [burn_dir_body, attitude_first_used, retarget_feasible, retarget_force_rel] = ...
        local_select_burn_direction(P, thruster_health, thruster_scale, ...
        F_cmd_mag, dv_dir_body, opts.retarget_force_tol);
end
q_cmd_burn = local_quat_from_two_vectors(burn_dir_body, dv_dir_inertial);
burn_start_s = double(attitude_first_used) * opts.attitude_maneuver_s;
burn_stop_s = burn_start_s + T_burn;
T_end = max(T_end, burn_stop_s + opts.sat_window_s);
num_steps = round(T_end / dt);
t = (0:num_steps-1).' * dt;

% Use the wheel-failure initial condition so the combined case evaluates
% attitude recovery as well as burn-time disturbance rejection.
q = [cos(deg2rad(2.5)); 0; 0; sin(deg2rad(2.5))];
w = deg2rad([0.3; -0.2; 0.1]);
dv_acc = zeros(3, 1);          % accumulated inertial-frame delta-v
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
    q_cmd = q_cmd_burn;
    [qm, wm] = sensor_model(q, w, P.sensor);
    qe = qmult(qinv(q_cmd), qm);
    we = wm;

    if t(k) >= burn_start_s && t(k) < burn_stop_s
        F_body_cmd = F_cmd_mag * burn_dir_body;
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
        T_thr_delta = T_pd - (Tw_body + Tb_force);
        T_thr_delta = max(min(T_thr_delta, 0.05), -0.05);
        T_thr_des = Tb_force + T_thr_delta;
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
    dv_acc = dv_acc + (local_quat_to_dcm(q) * Fb / P.mass) * dt;

    Fmax_eff = P.thr.Fmax * thruster_scale(:) .* double(thruster_health(:));
    q_log(:, k) = q;
    w_log(:, k) = w;
    qe_true = qmult(qinv(q_cmd), q);
    err_log(k) = 2 * acos(min(1, abs(qe_true(1)))) * 180 / pi;
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

window_s = opts.sat_window_s;
final_start_s = max(burn_stop_s, t(end) - window_s);
idx_final = t >= final_start_s;
idx_prev = t >= max(0, final_start_s - window_s) & t < final_start_s;
if ~any(idx_prev)
    idx_prev = idx_final;
end

steady_mean_deg = mean(err_log(idx_final));
steady_max_deg = max(err_log(idx_final));
prev_mean_deg = mean(err_log(idx_prev));
dv_error_rel = norm(dv_acc - dv_target_mag * dv_dir_inertial) / max(dv_target_mag, eps);
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
out.dv_target = dv_target_mag * dv_dir_inertial;
out.dv_dir_inertial = dv_dir_inertial;
out.burn_dir_body = burn_dir_body;
out.q_cmd_burn = q_cmd_burn;
out.burn_start_s = burn_start_s;
out.burn_stop_s = burn_stop_s;
out.final_window_start_s = final_start_s;
out.attitude_first_used = attitude_first_used;
out.retarget_feasible = retarget_feasible;
out.retarget_force_rel = retarget_force_rel;
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
function [burn_dir_body, retargeted, feasible, force_rel] = ...
    local_select_burn_direction(P, health, scale, F_cmd_mag, target_dir_body, force_tol)
target_dir_body = target_dir_body(:) / norm(target_dir_body);

[~, Fb_direct] = thruster_ft_allocation(F_cmd_mag * target_dir_body, ...
    zeros(3, 1), P, health, scale);
direct_force_rel = norm(Fb_direct - F_cmd_mag * target_dir_body) / ...
    max(F_cmd_mag, eps);
if direct_force_rel <= force_tol
    burn_dir_body = target_dir_body;
    retargeted = false;
    feasible = true;
    force_rel = direct_force_rel;
    return;
end

active = logical(health(:).') & scale(:).' > 0;
candidates = target_dir_body;
dirs = P.thr.dirs(:, active);
for i = 1:size(dirs, 2)
    u = dirs(:, i);
    if norm(u) > eps
        candidates(:, end+1) = u / norm(u); %#ok<AGROW>
    end
end
candidates = local_unique_dirs(candidates);

best_idx = 1;
best_score = [Inf Inf Inf];
force_rel_values = Inf(1, size(candidates, 2));
for i = 1:size(candidates, 2)
    u = candidates(:, i);
    [~, Fb, Tb] = thruster_ft_allocation(F_cmd_mag * u, zeros(3, 1), ...
        P, health, scale);
    f_rel = norm(Fb - F_cmd_mag * u) / max(F_cmd_mag, eps);
    force_rel_values(i) = f_rel;
    angle_cost = acos(max(-1, min(1, dot(u, target_dir_body))));
    score = [double(f_rel > force_tol), angle_cost, norm(Tb)];
    if local_score_lt(score, best_score)
        best_score = score;
        best_idx = i;
    end
end

burn_dir_body = candidates(:, best_idx);
force_rel = force_rel_values(best_idx);
feasible = force_rel <= force_tol;
retargeted = feasible && norm(cross(burn_dir_body, target_dir_body)) > 1e-6;
if ~feasible
    burn_dir_body = target_dir_body;
    retargeted = false;
end
end

% -------------------------------------------------------------------------
function dirs = local_unique_dirs(dirs)
keep = true(1, size(dirs, 2));
for i = 1:size(dirs, 2)
    if ~keep(i)
        continue;
    end
    dirs(:, i) = dirs(:, i) / norm(dirs(:, i));
    for j = i+1:size(dirs, 2)
        if norm(dirs(:, i) - dirs(:, j) / norm(dirs(:, j))) < 1e-9
            keep(j) = false;
        end
    end
end
dirs = dirs(:, keep);
end

% -------------------------------------------------------------------------
function tf = local_score_lt(a, b)
tol = 1e-12;
tf = false;
for i = 1:numel(a)
    if a(i) < b(i) - tol
        tf = true;
        return;
    elseif a(i) > b(i) + tol
        return;
    end
end
end

% -------------------------------------------------------------------------
function q = local_quat_from_two_vectors(a, b)
a = a(:) / norm(a);
b = b(:) / norm(b);
c = max(-1, min(1, dot(a, b)));
if c > 1 - 1e-12
    q = [1; 0; 0; 0];
elseif c < -1 + 1e-12
    axis = cross(a, [1; 0; 0]);
    if norm(axis) < 1e-9
        axis = cross(a, [0; 1; 0]);
    end
    axis = axis / norm(axis);
    q = [0; axis];
else
    axis = cross(a, b);
    q = [1 + c; axis];
    q = q / norm(q);
end
end

% -------------------------------------------------------------------------
function R = local_quat_to_dcm(q)
q = q(:) / norm(q);
s = q(1);
x = q(2);
y = q(3);
z = q(4);
R = [1 - 2*(y*y + z*z), 2*(x*y - s*z),     2*(x*z + s*y);
     2*(x*y + s*z),     1 - 2*(x*x + z*z), 2*(y*z - s*x);
     2*(x*z - s*y),     2*(y*z + s*x),     1 - 2*(x*x + y*y)];
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
