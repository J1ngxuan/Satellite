function results = sweep_double_wheel_failure(P, strategies)
% SWEEP_DOUBLE_WHEEL_FAILURE  Compare double-wheel fault strategies.
%
% The sweep stays within the research target: wheels are either healthy or
% completely failed. It evaluates whether the two surviving wheels can hold
% attitude before thruster assist is required.

if nargin < 2 || isempty(strategies)
    strategies = {'two_wheel_only', 'wheel_first', 'assist_immediate'};
end

pairs = nchoosek(1:P.wheel.N, 2);
num_cases = size(pairs, 1);
num_strat = numel(strategies);
results = repmat(local_empty_result(), num_cases * num_strat, 1);

row = 0;
for i = 1:num_cases
    failed = pairs(i, :);
    health = true(1, P.wheel.N);
    health(failed) = false;

    for j = 1:num_strat
        row = row + 1;
        opts = struct();
        opts.health = health;
        opts.double_strategy = strategies{j};
        opts.err_threshold_deg = 0.1;
        opts.window_s = 30;

        out = sim_wheel_failure(P, 'custom', opts);

        results(row).failed = failed;
        results(row).survivors = find(health);
        results(row).strategy = strategies{j};
        results(row).steady_mean_deg = out.steady_mean_deg;
        results(row).steady_max_deg = out.steady_max_deg;
        results(row).assist_used = out.assist_used;
        results(row).assist_time = out.assist_time;
        results(row).meets_0p1deg = out.steady_mean_deg <= 0.1;
        results(row).out = out;
    end
end
end

function r = local_empty_result()
r = struct('failed', [], ...
           'survivors', [], ...
           'strategy', '', ...
           'steady_mean_deg', NaN, ...
           'steady_max_deg', NaN, ...
           'assist_used', false, ...
           'assist_time', NaN, ...
           'meets_0p1deg', false, ...
           'out', []);
end
