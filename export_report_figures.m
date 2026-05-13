function export_report_figures()
% EXPORT_REPORT_FIGURES  Export layout and simulation figures for the
% project delivery report.
%
% Outputs:
%   ./figures/layout_a_6orth.png
%   ./figures/layout_b_8canted.png
%   ./figures/layout_c_12orth.png
%   ./figures/layout_e_12icosa.png
%   ./figures/reaction_wheel_faults.png
%   ./figures/double_wheel_strategy_sweep.png
%   ./figures/thruster_faults.png
%   ./figures/thruster_failure_feedforward_compare.png
%   ./figures/fault_tree_combined_sweep.png

clc;
close all;

outdir = fullfile(pwd, 'figures');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

P = satellite_params();
configs = thruster_configurations(P.dim, P.thr.Fmax, P.J);

layout_files = { ...
    'layout_a_6orth.png', ...
    'layout_b_8canted.png', ...
    'layout_c_12orth.png', ...
    'layout_e_12icosa.png'};

for k = 1:numel(configs)
    export_layout_figure(configs(k), P.dim, fullfile(outdir, layout_files{k}));
end

chosen = select_configuration(configs);
P.thr.dirs = chosen.dirs;
P.thr.pos = chosen.pos;
P.thr.M = chosen.M;

% Match the main simulation sequence so the exported figures and reported
% numbers stay consistent.
rng(42);
res_w0 = sim_wheel_failure(P, 'none');
res_w1 = sim_wheel_failure(P, 'one');
res_w2 = sim_wheel_failure(P, 'two');
export_reaction_wheel_figure(res_w0, res_w1, res_w2, ...
    fullfile(outdir, 'reaction_wheel_faults.png'));

strategies = {'two_wheel_only', 'wheel_first', 'assist_immediate'};
res_sweep = sweep_double_wheel_failure(P, strategies);
export_double_wheel_sweep_figure(res_sweep, strategies, ...
    fullfile(outdir, 'double_wheel_strategy_sweep.png'));

deg_list = [1 5];
fail_list = [2 6];
res_nom = sim_thruster_fault(P, 'nominal', [], [1; 0; 0]);
res_deg = sim_thruster_fault(P, 'degrade', deg_list, [1; 0; 0]);
res_fail = sim_thruster_fault(P, 'failure', fail_list, [1; 0; 0]);
export_thruster_fault_figure(res_nom, res_deg, res_fail, fail_list, ...
    fullfile(outdir, 'thruster_faults.png'));

% Use the same noise realization for the feedforward comparison.
rng(1234);
res_fail_ff = sim_thruster_fault_variant(P, 'failure', fail_list, [1; 0; 0], true);
rng(1234);
res_fail_noff = sim_thruster_fault_variant(P, 'failure', fail_list, [1; 0; 0], false);
export_feedforward_compare_figure(res_fail_ff, res_fail_noff, ...
    fullfile(outdir, 'thruster_failure_feedforward_compare.png'));

res_fault_tree = sweep_fault_tree_analysis(P);
export_fault_tree_figure(res_fault_tree, ...
    fullfile(outdir, 'fault_tree_combined_sweep.png'));

fprintf('Exported report figures to %s\n', outdir);
end

% -------------------------------------------------------------------------
function chosen = select_configuration(configs)
ok = [configs.redundant];
if any(ok)
    cand = configs(ok);
    [~, bestIdx] = max(arrayfun(@(c) c.eta_iso - 1e-3*c.M, cand));
    chosen = cand(bestIdx);
else
    [~, bestIdx] = max([configs.eta_iso]);
    chosen = configs(bestIdx);
end
end

% -------------------------------------------------------------------------
function export_layout_figure(config, dim, outfile)
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 700]);
ax = axes(fig);
hold(ax, 'on');
draw_body_box(ax, dim);

M = size(config.dirs, 2);
clr = lines(max(M, 7));
arrow_scale = 0.18 * max(dim);

for i = 1:M
    p = config.pos(:, i);
    d = config.dirs(:, i);
    c = clr(mod(i-1, size(clr, 1)) + 1, :);
    quiver3(ax, p(1), p(2), p(3), ...
        arrow_scale * d(1), arrow_scale * d(2), arrow_scale * d(3), ...
        0, 'LineWidth', 2.0, 'MaxHeadSize', 0.55, 'Color', c);
    text(ax, p(1), p(2), p(3), sprintf('  #%d', i), ...
        'Color', c, 'FontSize', 10, 'FontWeight', 'bold');
end

title(ax, sprintf('%s | M=%d | eta_{iso}=%.4f | eta_{axis}=%.4f', ...
    config.name, config.M, config.eta_iso, mean(config.eta_axis)), ...
    'Interpreter', 'tex');
subtitle(ax, 'Body box and thrust-force directions in body frame');
xlabel(ax, 'X [m]');
ylabel(ax, 'Y [m]');
zlabel(ax, 'Z [m]');
grid(ax, 'on');
axis(ax, 'equal');
padding = 0.22;
xlim(ax, [-dim(1)/2-padding, dim(1)/2+padding]);
ylim(ax, [-dim(2)/2-padding, dim(2)/2+padding]);
zlim(ax, [-dim(3)/2-padding, dim(3)/2+padding]);
view(ax, 35, 24);
box(ax, 'on');

exportgraphics(fig, outfile, 'Resolution', 300);
close(fig);
end

% -------------------------------------------------------------------------
function draw_body_box(ax, dim)
hx = dim(1) / 2;
hy = dim(2) / 2;
hz = dim(3) / 2;

verts = [ ...
    -hx -hy -hz;
     hx -hy -hz;
     hx  hy -hz;
    -hx  hy -hz;
    -hx -hy  hz;
     hx -hy  hz;
     hx  hy  hz;
    -hx  hy  hz];

faces = [ ...
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8];

patch(ax, 'Vertices', verts, 'Faces', faces, ...
    'FaceColor', [0.88 0.93 1.00], ...
    'FaceAlpha', 0.18, ...
    'EdgeColor', [0.35 0.45 0.60], ...
    'LineWidth', 1.0);

plot3(ax, [0 0.15], [0 0], [0 0], 'r-', 'LineWidth', 2);
plot3(ax, [0 0], [0 0.15], [0 0], 'g-', 'LineWidth', 2);
plot3(ax, [0 0], [0 0], [0 0.15], 'b-', 'LineWidth', 2);
text(ax, 0.16, 0, 0, 'X', 'Color', 'r', 'FontWeight', 'bold');
text(ax, 0, 0.16, 0, 'Y', 'Color', 'g', 'FontWeight', 'bold');
text(ax, 0, 0, 0.16, 'Z', 'Color', 'b', 'FontWeight', 'bold');
end

% -------------------------------------------------------------------------
function export_reaction_wheel_figure(res_w0, res_w1, res_w2, outfile)
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 1050]);

subplot(4, 1, 1);
plot(res_w0.t, res_w0.err_deg, res_w1.t, res_w1.err_deg, ...
     res_w2.t, res_w2.err_deg, 'LineWidth', 1.2);
grid on;
ylabel('Pointing err [deg]');
legend('no fail', '1 wheel fail', '2 wheels fail', 'Location', 'best');
title('Attitude pointing error under reaction-wheel faults');

subplot(4, 1, 2);
plot(res_w2.t, res_w2.w' * 180 / pi, 'LineWidth', 1.2);
grid on;
ylabel('\omega [deg/s]');
legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'best');
title('Body rates in the 2-wheel-failure case');

subplot(4, 1, 3);
plot(res_w0.t, res_w0.mode, res_w1.t, res_w1.mode, ...
     res_w2.t, res_w2.mode, 'LineWidth', 1.4);
ylim([0.5 2.5]);
yticks([1 2]);
yticklabels({'wheels', 'wheels+thr'});
grid on;
ylabel('Active actuators');
xlabel('Time [s]');
legend('no fail', '1 wheel', '2 wheels', 'Location', 'east');
title('Automatic co-control activation after rank loss');

subplot(4, 1, 4);
plot(res_w2.t, res_w2.gamma_hat(1:4, :)', 'LineWidth', 1.2);
hold on;
yline(0.5, 'k--', 'health threshold');
if any(res_w2.fault_alarm)
    xline(res_w2.detect_time, 'r--', 'alarm');
end
grid on;
ylim([-0.05 1.15]);
ylabel('\gamma wheel');
xlabel('Time [s]');
legend('W1', 'W2', 'W3', 'W4', 'Location', 'best');
title('Online wheel health estimates');

exportgraphics(fig, outfile, 'Resolution', 300);
close(fig);
end

% -------------------------------------------------------------------------
function export_double_wheel_sweep_figure(results, strategies, outfile)
num_strat = numel(strategies);
num_cases = numel(results) / num_strat;
means = reshape([results.steady_mean_deg], num_strat, num_cases).';
assist = reshape([results.assist_used], num_strat, num_cases).';

labels = cell(num_cases, 1);
for i = 1:num_cases
    failed = results((i-1) * num_strat + 1).failed;
    labels{i} = sprintf('[%d,%d]', failed(1), failed(2));
end

fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 950 650]);
bar(means, 'grouped');
hold on;
yline(0.1, 'k--', '0.1 deg criterion');
grid on;
ylabel('Final 20 s mean pointing error [deg]');
xlabel('Failed wheel pair');
set(gca, 'XTickLabel', labels);
legend(strrep(strategies, '_', '\_'), 'Location', 'best');
title('Double-wheel failure strategies: two-wheel hold before thruster fallback');

for i = 1:num_cases
    for j = 1:num_strat
        if assist(i, j)
            text(i + (j - 2) * 0.22, means(i, j), ' *', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
end

exportgraphics(fig, outfile, 'Resolution', 300);
close(fig);
end

% -------------------------------------------------------------------------
function export_thruster_fault_figure(res_nom, res_deg, res_fail, fail_list, outfile)
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 1050]);

subplot(4, 1, 1);
hold on;
plot(res_nom.t, vecnorm(res_nom.dv_achieved), 'LineWidth', 1.3, 'DisplayName', 'nominal');
plot(res_deg.t, vecnorm(res_deg.dv_achieved), 'LineWidth', 1.3, 'DisplayName', 'degrade');
plot(res_fail.t, vecnorm(res_fail.dv_achieved), 'LineWidth', 1.3, 'DisplayName', 'failure');
yline(norm(res_nom.dv_target), 'k--', 'target');
grid on;
ylabel('|dV| [m/s]');
legend('Location', 'southeast');
title('Achieved delta-v versus time');

subplot(4, 1, 2);
plot(res_deg.t, res_deg.err_deg, res_fail.t, res_fail.err_deg, 'LineWidth', 1.2);
grid on;
ylabel('Att err [deg]');
legend('degrade', 'failure', 'Location', 'best');
title('Attitude error during the burn');

subplot(4, 1, 3);
plot(res_fail.t, res_fail.F_hist', 'LineWidth', 1.0);
grid on;
ylabel('Thrust per nozzle [N]');
xlabel('Time [s]');
title(sprintf('Nozzle thrust commands in the failure case (dead = %s)', mat2str(fail_list)));

subplot(4, 1, 4);
hold on;
num_wheel = size(res_deg.Tw_hist, 1);
deg_list = res_deg.fault_list(res_deg.fault_list <= size(res_deg.F_hist, 1));
fail_list = fail_list(fail_list <= size(res_fail.F_hist, 1));
legend_entries = {};
if ~isempty(deg_list)
    plot(res_deg.t, res_deg.gamma_hat(num_wheel + deg_list, :)', 'LineWidth', 1.2);
    legend_entries = [legend_entries, arrayfun(@(i) sprintf('degrade T%d', i), ...
        deg_list, 'UniformOutput', false)];
end
if ~isempty(fail_list)
    plot(res_fail.t, res_fail.gamma_hat(num_wheel + fail_list, :)', '--', 'LineWidth', 1.2);
    legend_entries = [legend_entries, arrayfun(@(i) sprintf('failed T%d', i), ...
        fail_list, 'UniformOutput', false)];
end
yline(0.1, 'k--', 'usable threshold');
if any(res_deg.fault_alarm)
    xline(res_deg.detect_time, 'Color', [0.85 0.1 0.1], 'LineStyle', ':');
end
if any(res_fail.fault_alarm)
    xline(res_fail.detect_time, 'Color', [0.85 0.1 0.1], 'LineStyle', '--');
end
grid on;
ylim([-0.05 1.15]);
ylabel('\gamma thr');
xlabel('Time [s]');
if ~isempty(legend_entries)
    legend(legend_entries, 'Location', 'best');
end
title('Online thruster effectiveness estimates');

exportgraphics(fig, outfile, 'Resolution', 300);
close(fig);
end

% -------------------------------------------------------------------------
function export_feedforward_compare_figure(res_ff, res_noff, outfile)
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 800]);

subplot(2, 1, 1);
semilogy(res_ff.t, max(res_ff.err_deg, 1e-3), 'LineWidth', 1.4);
hold on;
semilogy(res_noff.t, max(res_noff.err_deg, 1e-3), 'LineWidth', 1.4);
grid on;
ylabel('Pointing err [deg]');
legend('with feedforward', 'without feedforward', 'Location', 'northwest');
title('Attitude error under total thruster failure');

subplot(2, 1, 2);
plot(res_ff.t, vecnorm(res_ff.w) * 180 / pi, 'LineWidth', 1.4);
hold on;
plot(res_noff.t, vecnorm(res_noff.w) * 180 / pi, 'LineWidth', 1.4);
grid on;
xlabel('Time [s]');
ylabel('|\omega| [deg/s]');
legend('with feedforward', 'without feedforward', 'Location', 'northwest');
title('Body-rate response under total thruster failure');

exportgraphics(fig, outfile, 'Resolution', 300);
close(fig);
end

% -------------------------------------------------------------------------
function export_fault_tree_figure(results, outfile)
[rate, att_max, row_labels, col_labels] = fault_tree_matrix(results);

fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 760]);

subplot(1, 2, 1);
imagesc(rate);
axis tight;
colorbar;
clim([0 100]);
set(gca, 'XTick', 1:numel(col_labels), 'XTickLabel', col_labels, ...
    'YTick', 1:numel(row_labels), 'YTickLabel', row_labels, ...
    'TickLabelInterpreter', 'none');
xlabel('Thruster fault mode');
ylabel('Wheel fault case');
title('Fault-tree top-event rate [%]');

subplot(1, 2, 2);
imagesc(att_max);
axis tight;
colorbar;
set(gca, 'XTick', 1:numel(col_labels), 'XTickLabel', col_labels, ...
    'YTick', 1:numel(row_labels), 'YTickLabel', row_labels, ...
    'TickLabelInterpreter', 'none');
xlabel('Thruster fault mode');
title('Max final-20-s mean attitude error [deg]');

exportgraphics(fig, outfile, 'Resolution', 300);
close(fig);
end

% -------------------------------------------------------------------------
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
            rate(i, j) = NaN;
            att_max(i, j) = NaN;
        else
            rate(i, j) = 100 * mean([subset.top_event]);
            att_max(i, j) = max([subset.steady_mean_deg]);
        end
    end
end
end

% -------------------------------------------------------------------------
function out = sim_thruster_fault_variant(P, mode, fault_list, dv_dir_body, use_feedforward)
% Variant of the mission burn simulation with optional wheel feedforward.

if nargin < 2, mode = 'nominal'; end
if nargin < 3, fault_list = []; end
if nargin < 4, dv_dir_body = [1; 0; 0]; end
if nargin < 5, use_feedforward = true; end

dt = 0.1;
T_burn = 40;
T_end = 60;
N = round(T_end / dt);
t = (0:N-1).' * dt;

dv_target_mag = 0.3;
acc_cmd_mag = dv_target_mag / T_burn;
F_cmd_mag = P.mass * acc_cmd_mag;

q = [1; 0; 0; 0];
w = zeros(3, 1);
q_cmd = q;
dv_acc = zeros(3, 1);

M = P.thr.M;
Nw = P.wheel.N;
health = true(1, M);
scale = ones(1, M);

switch mode
    case 'failure'
        health(fault_list) = false;
    case 'degrade'
        scale(fault_list) = 0.4;
    case 'nominal'
        % no-op
    otherwise
        error('Unknown mode');
end

wheel_health = true(1, Nw);

q_log = zeros(4, N);
w_log = zeros(3, N);
err_log = zeros(1, N);
F_hist = zeros(M, N);
dv_log = zeros(3, N);

for k = 1:N
    [qm, wm] = sensor_model(q, w, P.sensor);
    qe = qmult(qinv(q_cmd), qm);
    we = wm;

    if t(k) < T_burn
        F_body_cmd = F_cmd_mag * dv_dir_body(:);
    else
        F_body_cmd = [0; 0; 0];
    end

    [F, Fb, Tb_thr] = thruster_ft_allocation(F_body_cmd, [0; 0; 0], P, health, scale);

    if use_feedforward
        T_ff = Tb_thr;
    else
        T_ff = zeros(3, 1);
    end

    [Tw_body, ~] = wheel_attitude_controller(qe, we, P, wheel_health, T_ff);

    Td = 1e-5 * [sin(0.01 * t(k)); cos(0.01 * t(k)); 0.2];
    [q, w] = attitude_dynamics(q, w, Tw_body + Tb_thr, Td, P.J, P.Jinv, dt);

    dv_acc = dv_acc + (Fb / P.mass) * dt;

    q_log(:, k) = q;
    w_log(:, k) = w;
    err_log(k) = 2 * acos(min(1, abs(q(1)))) * 180 / pi;
    F_hist(:, k) = F;
    dv_log(:, k) = dv_acc;
end

out.t = t;
out.q = q_log;
out.w = w_log;
out.err_deg = err_log;
out.F_hist = F_hist;
out.dv_achieved = dv_log;
out.dv_target = dv_target_mag * dv_dir_body(:);
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
