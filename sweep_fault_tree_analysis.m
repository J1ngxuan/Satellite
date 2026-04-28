function results = sweep_fault_tree_analysis(P, opts)
% SWEEP_FAULT_TREE_ANALYSIS  Enumerate combined actuator faults.
%
% Top event: attitude/orbit-control mission failure.
% Intermediate events:
%   1) attitude hold failure: final 20 s mean pointing error > 0.1 deg
%   2) orbit maneuver failure: final delta-v relative error > 5%
%   3) unrecovered co-control: long saturation and non-converging attitude
%
% The default sweep covers:
%   - wheel states: no failure, any one wheel failed, any two wheels failed
%   - thruster states: nominal, any two degraded to 40%, any two failed
%   - burn directions: +X, +Y, +Z

if nargin < 2 || isempty(opts)
    opts = struct();
end

opts = local_default(opts, 'dt', 0.5);
opts = local_default(opts, 'T_end', 60);
opts = local_default(opts, 'T_burn', 40);
opts = local_default(opts, 'dv_dirs', eye(3));
opts = local_default(opts, 'keep_timeseries', false);
opts = local_default(opts, 'seed', 20260427);
opts = local_default(opts, 'max_thruster_pairs_per_type', Inf);

wheel_cases = local_wheel_cases(P.wheel.N);
thruster_cases = local_thruster_cases(P.thr.M, opts.max_thruster_pairs_per_type);
dv_dirs = opts.dv_dirs;

num_rows = numel(wheel_cases) * numel(thruster_cases) * size(dv_dirs, 2);
results = repmat(local_empty_result(), num_rows, 1);

row = 0;
for iw = 1:numel(wheel_cases)
    wc = wheel_cases(iw);
    for it = 1:numel(thruster_cases)
        tc = thruster_cases(it);
        for id = 1:size(dv_dirs, 2)
            row = row + 1;
            rng(opts.seed + row - 1);

            sim_opts = struct();
            sim_opts.wheel_health = wc.health;
            sim_opts.thruster_health = tc.health;
            sim_opts.thruster_scale = tc.scale;
            sim_opts.dv_dir_body = dv_dirs(:, id);
            sim_opts.dt = opts.dt;
            sim_opts.T_end = opts.T_end;
            sim_opts.T_burn = opts.T_burn;

            sim = sim_combined_fault(P, sim_opts);

            results(row).case_id = row;
            results(row).wheel_case = wc.name;
            results(row).wheel_failed = wc.failed;
            results(row).wheel_fault_order = numel(wc.failed);
            results(row).thruster_case = tc.name;
            results(row).thruster_fault_type = tc.fault_type;
            results(row).thruster_faulted = tc.faulted;
            results(row).dv_dir_label = local_dir_label(dv_dirs(:, id));
            results(row).top_event = sim.top_event;
            results(row).branch = sim.branch;
            results(row).attitude_fail = sim.attitude_fail;
            results(row).orbit_fail = sim.orbit_fail;
            results(row).co_control_fail = sim.co_control_fail;
            results(row).steady_mean_deg = sim.steady_mean_deg;
            results(row).steady_max_deg = sim.steady_max_deg;
            results(row).dv_error_rel = sim.dv_error_rel;
            results(row).assist_used = sim.assist_used;
            results(row).wheel_sat_frac = sim.wheel_sat_frac;
            results(row).thruster_sat_frac = sim.thruster_sat_frac;
            results(row).thruster_feasible_frac = sim.thruster_feasible_frac;
            if opts.keep_timeseries
                results(row).out = sim;
            end
        end
    end
end
end

% -------------------------------------------------------------------------
function opts = local_default(opts, name, value)
if ~isfield(opts, name) || isempty(opts.(name))
    opts.(name) = value;
end
end

% -------------------------------------------------------------------------
function cases = local_wheel_cases(num_wheels)
defs = {};
defs{end+1} = [];
for i = 1:num_wheels
    defs{end+1} = i; %#ok<AGROW>
end
pairs = nchoosek(1:num_wheels, 2);
for i = 1:size(pairs, 1)
    defs{end+1} = pairs(i, :); %#ok<AGROW>
end

cases = repmat(struct('name', '', 'failed', [], 'health', []), numel(defs), 1);
for i = 1:numel(defs)
    failed = defs{i};
    health = true(1, num_wheels);
    health(failed) = false;
    cases(i).failed = failed;
    cases(i).health = health;
    if isempty(failed)
        cases(i).name = 'W0';
    else
        cases(i).name = sprintf('W%d_%s', numel(failed), local_join_ints(failed, '_'));
    end
end
end

% -------------------------------------------------------------------------
function cases = local_thruster_cases(num_thrusters, max_pairs_per_type)
pairs = nchoosek(1:num_thrusters, 2);
if isfinite(max_pairs_per_type)
    pairs = pairs(1:min(size(pairs, 1), max_pairs_per_type), :);
end

num_cases = 1 + 2 * size(pairs, 1);
cases = repmat(struct('name', '', 'fault_type', '', 'faulted', [], ...
    'health', true(1, num_thrusters), 'scale', ones(1, num_thrusters)), ...
    num_cases, 1);

cases(1).name = 'T0';
cases(1).fault_type = 'nominal';

row = 1;
for i = 1:size(pairs, 1)
    row = row + 1;
    faulted = pairs(i, :);
    scale = ones(1, num_thrusters);
    scale(faulted) = 0.4;
    cases(row).name = sprintf('TD_%s', local_join_ints(faulted, '_'));
    cases(row).fault_type = 'degrade';
    cases(row).faulted = faulted;
    cases(row).scale = scale;
end

for i = 1:size(pairs, 1)
    row = row + 1;
    faulted = pairs(i, :);
    health = true(1, num_thrusters);
    health(faulted) = false;
    cases(row).name = sprintf('TF_%s', local_join_ints(faulted, '_'));
    cases(row).fault_type = 'failure';
    cases(row).faulted = faulted;
    cases(row).health = health;
end
end

% -------------------------------------------------------------------------
function label = local_dir_label(v)
[~, idx] = max(abs(v));
axis_name = 'XYZ';
if v(idx) >= 0
    label = ['+' axis_name(idx)];
else
    label = ['-' axis_name(idx)];
end
end

% -------------------------------------------------------------------------
function s = local_join_ints(values, sep)
if isempty(values)
    s = '';
    return;
end
parts = arrayfun(@(x) sprintf('%d', x), values, 'UniformOutput', false);
s = strjoin(parts, sep);
end

% -------------------------------------------------------------------------
function r = local_empty_result()
r = struct('case_id', NaN, ...
           'wheel_case', '', ...
           'wheel_failed', [], ...
           'wheel_fault_order', NaN, ...
           'thruster_case', '', ...
           'thruster_fault_type', '', ...
           'thruster_faulted', [], ...
           'dv_dir_label', '', ...
           'top_event', false, ...
           'branch', '', ...
           'attitude_fail', false, ...
           'orbit_fail', false, ...
           'co_control_fail', false, ...
           'steady_mean_deg', NaN, ...
           'steady_max_deg', NaN, ...
           'dv_error_rel', NaN, ...
           'assist_used', false, ...
           'wheel_sat_frac', NaN, ...
           'thruster_sat_frac', NaN, ...
           'thruster_feasible_frac', NaN, ...
           'out', []);
end
