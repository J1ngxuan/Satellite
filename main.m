function main()
% MAIN  Top-level driver for the satellite fault-tolerant control study.
%
% Pipeline:
%   1. Print platform parameters.
%   2. Build thruster configuration library, evaluate redundancy and
%      omnidirectional efficiency, recommend a configuration.
%   3. Run reaction-wheel fault simulations (no fail / 1 fail / 2 fail).
%   4. Run thruster degradation simulation (>=2 nozzles at reduced thrust).
%   5. Run thruster total-failure simulation (>=2 nozzles dead).
%   6. Plot results and print numerical summaries.

clc; close all;
fprintf('==============================================================\n');
fprintf('  Satellite Fault-Tolerant Attitude/Orbit Control Simulation\n');
fprintf('==============================================================\n');

P = satellite_params();
fprintf('Mass = %.1f kg | Inertia diag = [%.4f %.4f %.4f] kg*m^2\n', ...
        P.mass, P.J(1,1), P.J(2,2), P.J(3,3));
fprintf('Reaction wheels: %d (pyramid)  | Wheel Tmax = %.3f Nm\n', ...
        P.wheel.N, P.wheel.Tmax);

% ---------- (1) Thruster configuration trade study ---------------------
fprintf('\n--- Thruster Configuration Trade Study ---\n');
configs = thruster_configurations(P.dim, P.thr.Fmax, P.J);
fprintf('%-12s %4s %10s %12s %12s\n','Name','M','Redundant','eta_iso','eta_axis_avg');
for k=1:numel(configs)
    c = configs(k);
    fprintf('%-12s %4d %10s %12.4f %12.4f\n', ...
        c.name, c.M, ternary(c.redundant,'YES','NO'), c.eta_iso, mean(c.eta_axis));
end

% Selection rule: must be redundant; among those pick highest eta_iso,
% break ties by fewest nozzles.
ok = [configs.redundant];
if any(ok)
    cand = configs(ok);
    [~,bestIdx] = max(arrayfun(@(c) c.eta_iso - 1e-3*c.M, cand));
    chosen = cand(bestIdx);
else
    [~,bestIdx] = max([configs.eta_iso]); chosen = configs(bestIdx);
end
fprintf('\n>> Selected configuration: %s (M=%d, eta_iso=%.4f)\n', ...
        chosen.name, chosen.M, chosen.eta_iso);

% Adopt the chosen layout
P.thr.dirs = chosen.dirs; P.thr.pos = chosen.pos; P.thr.M = chosen.M;

% ---------- (2) Reaction-wheel failure simulations ---------------------
fprintf('\n--- Reaction-Wheel Fault Simulations ---\n');
rng(42);                                     % reproducibility
res_w0 = sim_wheel_failure(P,'none');
res_w1 = sim_wheel_failure(P,'one');
res_w2 = sim_wheel_failure(P,'two');

fprintf('Steady-state pointing error (final 20 s mean):\n');
print_steady('  no-fail ', res_w0);
print_steady('  1-wheel ', res_w1);
print_steady('  2-wheel ', res_w2);

figure('Name','Reaction Wheel Faults','Color','w');
subplot(3,1,1);
plot(res_w0.t, res_w0.err_deg, res_w1.t, res_w1.err_deg, ...
     res_w2.t, res_w2.err_deg,'LineWidth',1.2);
grid on; ylabel('Pointing err [deg]');
legend('no fail','1 wheel fail','2 wheels fail (co-control)','Location','best');
title('Attitude pointing error under reaction-wheel faults');

subplot(3,1,2);
plot(res_w2.t, res_w2.w'*180/pi,'LineWidth',1.2);
grid on; ylabel('\omega [deg/s]'); legend('\omega_x','\omega_y','\omega_z');
title('Body rates (2-wheel-fail case)');

subplot(3,1,3);
plot(res_w0.t, res_w0.mode, res_w1.t, res_w1.mode, ...
     res_w2.t, res_w2.mode,'LineWidth',1.4);
ylim([0.5 2.5]); yticks([1 2]); yticklabels({'wheels','wheels+thr'});
grid on; ylabel('Active actuators'); xlabel('Time [s]');
legend('no fail','1 wheel','2 wheels','Location','east');
title('Co-control activation (rank<3 triggers thruster assist)');

fprintf('\n--- Double-Wheel Failure Strategy Sweep ---\n');
strategies = {'two_wheel_only', 'wheel_first', 'assist_immediate'};
res_sweep = sweep_double_wheel_failure(P, strategies);
print_double_sweep(res_sweep);
plot_double_sweep(res_sweep, strategies);

% ---------- (3) Thruster degradation simulation ------------------------
fprintf('\n--- Thruster Thrust-Degradation Simulation ---\n');
deg_list = [1 5];                            % >=2 thrusters degraded
res_deg  = sim_thruster_fault(P,'degrade', deg_list, [1;0;0]);
res_nom  = sim_thruster_fault(P,'nominal',[],         [1;0;0]);
fprintf('  Nominal :  achieved dV = [%.3f %.3f %.3f] m/s\n',  res_nom.dv_achieved(:,end));
fprintf('  Degraded:  achieved dV = [%.3f %.3f %.3f] m/s (faulty: %s)\n', ...
        res_deg.dv_achieved(:,end), num2str(deg_list));

% ---------- (4) Thruster total-failure simulation ---------------------
fprintf('\n--- Thruster Total-Failure Simulation ---\n');
fail_list = [2 6];                            % >=2 thrusters dead
res_fail  = sim_thruster_fault(P,'failure', fail_list, [1;0;0]);
fprintf('  Failure :  achieved dV = [%.3f %.3f %.3f] m/s (failed: %s)\n', ...
        res_fail.dv_achieved(:,end), num2str(fail_list));

figure('Name','Thruster Faults','Color','w');
subplot(3,1,1); hold on;
plot(res_nom.t, vecnorm(res_nom.dv_achieved),'LineWidth',1.3,'DisplayName','nominal');
plot(res_deg.t, vecnorm(res_deg.dv_achieved),'LineWidth',1.3,'DisplayName','degrade');
plot(res_fail.t,vecnorm(res_fail.dv_achieved),'LineWidth',1.3,'DisplayName','failure');
yline(norm(res_nom.dv_target),'k--','target');
grid on; ylabel('|dV| [m/s]'); legend('Location','southeast');
title('Achieved dV vs. time');

subplot(3,1,2);
plot(res_deg.t, res_deg.err_deg, res_fail.t, res_fail.err_deg,'LineWidth',1.2);
grid on; ylabel('Att err [deg]'); legend('degrade','failure');
title('Attitude error during burn (orbit ctrl w/o attitude maneuver)');

subplot(3,1,3);
plot(res_fail.t, res_fail.F_hist','LineWidth',1.0);
grid on; ylabel('Thrust per nozzle [N]'); xlabel('Time [s]');
title(sprintf('Thrust commands (failure case, dead = %s)', mat2str(fail_list)));

% ---------- (5) Combined fault-tree sweep ------------------------------
fprintf('\n--- Fault Tree / Combined Fault Sweep ---\n');
res_fault_tree = sweep_fault_tree_analysis(P);
print_fault_tree_summary(res_fault_tree);
plot_fault_tree_summary(res_fault_tree);

fprintf('\nAll simulations finished. Figures opened.\n');
end

% ----------------------------------------------------------------------
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
function print_steady(label, r)
idx = r.t >= (r.t(end)-20);
fprintf('%s : mean = %.4f deg, max = %.4f deg\n', ...
        label, mean(r.err_deg(idx)), max(r.err_deg(idx)));
end

function print_double_sweep(results)
fprintf('%-9s %-17s %10s %10s %8s %10s\n', ...
    'Failed', 'Strategy', 'Mean[deg]', 'Max[deg]', '<0.1?', 'Assist[s]');
for i = 1:numel(results)
    r = results(i);
    assist_s = 'none';
    if r.assist_used
        assist_s = sprintf('%.1f', r.assist_time);
    end
    fprintf('[%d %d]     %-17s %10.4f %10.4f %8s %10s\n', ...
        r.failed(1), r.failed(2), r.strategy, r.steady_mean_deg, ...
        r.steady_max_deg, ternary(r.meets_0p1deg, 'YES', 'NO'), assist_s);
end
end

function plot_double_sweep(results, strategies)
num_strat = numel(strategies);
num_cases = numel(results) / num_strat;
means = reshape([results.steady_mean_deg], num_strat, num_cases).';
assist = reshape([results.assist_used], num_strat, num_cases).';

labels = cell(num_cases, 1);
for i = 1:num_cases
    failed = results((i-1)*num_strat + 1).failed;
    labels{i} = sprintf('[%d,%d]', failed(1), failed(2));
end

figure('Name','Double Wheel Failure Strategy Sweep','Color','w');
bar(means, 'grouped');
hold on; yline(0.1, 'k--', '0.1 deg criterion');
grid on; ylabel('Final 20 s mean pointing error [deg]');
set(gca, 'XTickLabel', labels);
xlabel('Failed wheel pair');
legend(strrep(strategies, '_', '\_'), 'Location', 'best');
title('Two-wheel-only control before thruster fallback');

for i = 1:num_cases
    for j = 1:num_strat
        if assist(i,j)
            text(i + (j-2)*0.22, means(i,j), ' *', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
end
end

function print_fault_tree_summary(results)
num_cases = numel(results);
top = [results.top_event];
att = [results.attitude_fail];
orb = [results.orbit_fail];
coop = [results.co_control_fail];
assist = [results.assist_used];
att_first = [results.attitude_first_used];

fprintf('Fault-tree cases: %d | top events: %d (%.1f%%)\n', ...
    num_cases, nnz(top), 100*mean(top));
fprintf('  branches: attitude=%d, orbit=%d, co-control=%d | thruster assist used=%d | attitude-first burns=%d\n', ...
    nnz(att), nnz(orb), nnz(coop), nnz(assist), nnz(att_first));

fprintf('%-7s %-8s %8s %8s %9s %8s %s\n', ...
    'WheelN', 'ThrMode', 'Cases', 'Top[%]', 'AttMax', 'DvMax[%]', 'Assist[%]');
orders = unique([results.wheel_fault_order]);
types = {'nominal', 'degrade', 'failure'};
for i = 1:numel(orders)
    for j = 1:numel(types)
        mask = [results.wheel_fault_order] == orders(i) & ...
            strcmp({results.thruster_fault_type}, types{j});
        if ~any(mask)
            continue;
        end
        subset = results(mask);
        fprintf('%-7d %-8s %8d %8.1f %9.4f %8.2f %8.1f\n', ...
            orders(i), types{j}, numel(subset), ...
            100*mean([subset.top_event]), max([subset.steady_mean_deg]), ...
            100*max([subset.dv_error_rel]), 100*mean([subset.assist_used]));
    end
end
end

function plot_fault_tree_summary(results)
[rate, att_max, row_labels, col_labels] = fault_tree_matrix(results);

figure('Name','Fault Tree Combined Fault Sweep','Color','w');
subplot(1,2,1);
imagesc(rate);
axis tight;
colorbar;
clim([0 100]);
set(gca, 'XTick', 1:numel(col_labels), 'XTickLabel', col_labels, ...
    'YTick', 1:numel(row_labels), 'YTickLabel', row_labels);
xlabel('Thruster fault mode');
ylabel('Wheel fault case');
title('Top-event rate [%]');

subplot(1,2,2);
imagesc(att_max);
axis tight;
colorbar;
set(gca, 'XTick', 1:numel(col_labels), 'XTickLabel', col_labels, ...
    'YTick', 1:numel(row_labels), 'YTickLabel', row_labels);
xlabel('Thruster fault mode');
title('Max final-20-s mean attitude error [deg]');
end

function [rate, att_max, row_labels, col_labels] = fault_tree_matrix(results)
row_labels = unique({results.wheel_case}, 'stable');
col_labels = {'nominal', 'degrade', 'failure'};
rate = zeros(numel(row_labels), numel(col_labels));
att_max = zeros(numel(row_labels), numel(col_labels));

for i = 1:numel(row_labels)
    for j = 1:numel(col_labels)
        mask = strcmp({results.wheel_case}, row_labels{i}) & ...
            strcmp({results.thruster_fault_type}, col_labels{j});
        subset = results(mask);
        if isempty(subset)
            rate(i,j) = NaN;
            att_max(i,j) = NaN;
        else
            rate(i,j) = 100 * mean([subset.top_event]);
            att_max(i,j) = max([subset.steady_mean_deg]);
        end
    end
end
end
