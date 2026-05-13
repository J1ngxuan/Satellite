function export_architecture_diagrams()
% EXPORT_ARCHITECTURE_DIAGRAMS  Generate two report figures:
%   figures/architecture_block_diagram.png  - closed-loop data flow
%   figures/cocontrol_state_machine.png     - cooperative-control modes
%
% Layouts use orthogonal polyline routing through pre-allocated corridors
% so cables never cross blocks. Labels are placed on long straight segments
% to avoid arrowhead overlap.

clc; close all;

outdir = fullfile(pwd, 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

draw_architecture(fullfile(outdir, 'architecture_block_diagram.png'));
draw_state_machine(fullfile(outdir, 'cocontrol_state_machine.png'));

fprintf('Exported architecture diagrams to %s\n', outdir);
end

% =========================================================================
function draw_architecture(file_path)

fig = figure('Name', 'System Architecture', 'Color', 'w', ...
    'Position', [60 60 1320 780]);
ax = axes('Parent', fig, 'Position', [0.02 0.02 0.96 0.92]); %#ok<LAXES>
hold(ax, 'on'); axis(ax, 'off'); axis(ax, [0 100 0 60]);
ax.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax);

c_ref   = [0.86 0.92 1.00];
c_ctrl  = [0.78 0.90 0.78];
c_diag  = [1.00 0.86 0.66];
c_plant = [0.85 0.85 0.85];
c_sens  = [0.95 0.85 0.95];
c_logic = [1.00 0.94 0.72];

% ---- Blocks (cx, cy, w, h) ------------------------------------------
blk(ax, 10, 51, 16, 6, c_ref,   {'\bf 指令 / Reference', ...
    'q_{cmd},  \omega_{cmd}=0', '\Deltav 方向 dv\_dir\_body'});
blk(ax, 32, 51, 14, 6, c_ref,   {'\bf 姿态误差', ...
    'q_e = q_{cmd}^{-1} \otimes q_{meas}', '\omega_e = \omega_{meas}'});
blk(ax, 58, 51, 18, 6, c_ctrl,  {'\bf 飞轮控制器', ...
    'wheel\_attitude\_controller', 'PD + 加权伪逆'});
blk(ax, 85, 39, 14, 24, c_plant,{'\bf 动力学 / Plant', ...
    'attitude\_dynamics', '(RK4, 四元数 + Euler)', ...
    'orbit\_dynamics', '(RK4, 二体 + J_2 + 推力)'});
blk(ax, 58, 36, 18, 6, c_ctrl,  {'\bf 推力器分配', ...
    'thruster\_ft\_allocation (NNLS)', 'health/scale 闸门 + 饱和'});
blk(ax, 32, 36, 14, 6, c_logic, {'\bf 协同控制决策', ...
    'rank(A_w)<3 \Rightarrow assist', ...
    '方向不可达 \Rightarrow 姿态优先', ...
    '|T_w|>0.7T_{max} \Rightarrow 节流'});
blk(ax, 58, 20, 18, 6, c_sens,  {'\bf 传感器 sensor\_model', ...
    '星敏 q_{meas},  陀螺 \omega_{meas}', '加速度计 a_{meas}'});
blk(ax, 24, 20, 22, 6, c_diag,  {'\bf RLS 故障诊断', ...
    'actuator\_fault\_diagnosis', ...
    '输入: [T_w;F]_{cmd}, \omega_{meas}, a_{meas}'});
blk(ax, 32,  9, 50, 4, c_logic, {'\bf 健康向量总线  \gamma_{hat}', ...
    'wheel\_health\_est (N_w) | thruster\_health\_est, scale\_est (M)'});

% ---- arrows (polyline routing) --------------------------------------
% A1: Reference -> Error
arrow_poly(ax, [18 25], [51 51]);
% A2: F_body_cmd -> Thruster allocator (top entry)
arrow_poly(ax, [10 10 58 58], [48 45 45 39], 'F_{body,cmd}', 2);
% A3: Error -> Wheel controller
arrow_poly(ax, [39 49], [51 51]);
% A4: Error -> Logic (state output q_e, w_e for rank/throttle check)
arrow_poly(ax, [32 32], [48 39], 'q_e, \omega_e', 1);
% A5: Logic -> Thruster allocator (assist / T_thr,des)
arrow_poly(ax, [39 49], [36 36], 'assist / T_{thr,des}', 1);
% A6: Logic -> Wheel controller (assist flag, scale)
%     horizontal at y=42 to avoid overlap with A12 horizontal at y=47.
arrow_poly(ax, [36 36 53 53], [39 42 42 48], 'assist / scale', 2);
% A7: T_b,thr feedforward Thruster -> Wheel controller
arrow_poly(ax, [63 63], [39 48], 'T_{b,thr} 前馈', 1);
% A8: Wheel controller -> Plant
arrow_poly(ax, [67 78], [51 47], 'T_w', 1);
% A9: Thruster allocator -> Plant
arrow_poly(ax, [67 78], [36 41], 'F, T_{b,thr}', 1);
% A10: Plant -> Sensor
arrow_poly(ax, [85 85 67], [27 20 20], 'q, \omega, a', 2);
% A11: Sensor -> RLS
arrow_poly(ax, [49 35], [20 20], '\omega_{meas}, a_{meas}', 1);
% A12: q_meas, omega_meas feedback Sensor -> Error
arrow_poly(ax, [49 42 42 32 32], [23 23 47 47 48], 'q_{meas}, \omega_{meas}', 2);
% A13: Command tap [T_w; F]_cmd -> RLS
arrow_poly(ax, [49 49 32 32], [33 28 28 23], '[T_w; F]_{cmd}', 2);
% A14: RLS -> Health bus
arrow_poly(ax, [24 24], [17 11]);
% A15: Bus -> Wheel controller (main vertical at x=72, label on long seg)
arrow_poly(ax, [57 72 72 67], [9 9 51 51], '\gamma_w  飞轮健康', 2);
% A16: Bus -> Thruster allocator (branch at x=72, y=36)
branch_dot(ax, 72, 36);
arrow_poly(ax, [72 67], [36 36], '\gamma_t  推力器 health/scale', 1);

% Title
text(ax, 50, 58, ['系统架构  —  sensors \rightarrow RLS 诊断 \rightarrow ' ...
    '飞轮 / 推力器分配 \rightarrow 动力学'], ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

% Legend chips along bottom
legend_chip(ax,  2, 1.5, c_ref,   '指令 / 参考');
legend_chip(ax, 19, 1.5, c_ctrl,  '控制器 / 分配');
legend_chip(ax, 36, 1.5, c_diag,  'RLS 诊断');
legend_chip(ax, 53, 1.5, c_plant, '动力学');
legend_chip(ax, 68, 1.5, c_sens,  '传感器');
legend_chip(ax, 83, 1.5, c_logic, '决策 / 健康总线');

exportgraphics(fig, file_path, 'Resolution', 200);
end

% =========================================================================
function draw_state_machine(file_path)

fig = figure('Name', 'Co-control State Machine', 'Color', 'w', ...
    'Position', [60 60 1320 780]);
ax = axes('Parent', fig, 'Position', [0.02 0.02 0.96 0.92]); %#ok<LAXES>
hold(ax, 'on'); axis(ax, 'off'); axis(ax, [0 100 0 60]);
ax.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax);

c_nom   = [0.78 0.92 0.78];
c_wheel = [0.85 0.90 1.00];
c_coop  = [1.00 0.85 0.72];
c_burn  = [0.95 0.85 0.95];
c_att1  = [1.00 0.94 0.72];
c_thr   = [0.98 0.78 0.78];
c_done  = [0.82 0.94 0.94];

% Top row (top-left = nominal, then 1-wheel, then co-control)
sblk(ax, 18, 48, 18, 8, c_nom, ...
    {'\bf S_0  标称', 'mode = 1, 仅飞轮控姿', 'rank(A_w) = 3, 无机动'});
sblk(ax, 50, 48, 18, 8, c_wheel, ...
    {'\bf S_1  单飞轮故障', '飞轮加权伪逆容错', 'feasible = true'});
sblk(ax, 82, 48, 18, 8, c_coop, ...
    {'\bf S_2  飞轮+推力器协同', 'mode = 2, assist\_active', 'T_{missing} 力偶补偿'});

% Middle row (left = attitude-first, center = burn, right = throttle)
sblk(ax, 18, 30, 18, 8, c_att1, ...
    {'\bf S_3  姿态优先重定向', 'attitude\_first\_used = true', ...
    '调姿 attitude\_maneuver\_s'});
sblk(ax, 50, 30, 18, 8, c_burn, ...
    {'\bf S_4  轨控点火', '推力器执行 \Deltav', '飞轮控姿 + T_{b,thr} 前馈'});
sblk(ax, 82, 30, 18, 8, c_thr, ...
    {'\bf S_5  力矩安全节流', 'force\_scale \in grid', '|T_w| \leq 0.7 T_{max}'});

% Bottom row (burn complete)
sblk(ax, 50, 12, 18, 8, c_done, ...
    {'\bf S_6  机动完成', 'dv\_acc \geq target(1-tol)', '回到姿态保持'});

% Initial entry marker
plot(ax, 5, 48, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
arrow_poly(ax, [6 9], [48 48]);

% --- Transitions -----------------------------------------------------
arrow_poly(ax, [27 41], [48 48], '1 飞轮失效');
arrow_poly(ax, [59 73], [48 48], 'rank(A_w)<3');
% S0 -> S2 direct (top arc, peak well below title at y=58)
arrow_arc(ax, 18, 52, 82, 52, +3, '2 飞轮同时失效');

% Burn requests
arrow_poly(ax, [50 50], [44 34], sprintf('\\Deltav 请求\n方向可达'));
arrow_poly(ax, [18 18], [44 34], sprintf('\\Deltav 请求\n方向不可达'));

% S3 -> S4 after slew
arrow_poly(ax, [27 41], [30 30], '调姿完成');

% S4 -> S5 throttle trigger
arrow_poly(ax, [59 73], [30 30], '|T_w|>0.7T_{max}');
% S5 -> S4 throttle resolved (small downward arc to avoid colliding with the above arrow)
arrow_arc(ax, 73, 27, 59, 27, -3, '减小 scale');

% S4 -> S6 burn complete
arrow_poly(ax, [50 50], [26 16], 'dv 达标');

% S2 -> S4 (co-control engaged during burn) — diagonal in open corridor
arrow_poly(ax, [74 58], [44 34], 'rank<3 期间');

% S6 -> S0 return loop, routed around the left of the chart
arrow_poly(ax, [41 5 5 9], [12 12 48 48], '机动结束', 2);

% Title
text(ax, 50, 58, ['协同控制状态机  ' ...
    '(sim\_wheel\_failure / sim\_combined\_fault)'], ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

% Footnote
text(ax, 2, 2.5, ...
    ['\rm 转移由在线诊断 \gamma_{hat} 触发: rank(A_w) 由 wheel\_health\_est 计算; ' ...
     '方向可行性由 thruster\_ft\_allocation 残差检验。'], ...
    'FontSize', 9, 'Color', [0.2 0.2 0.2]);

exportgraphics(fig, file_path, 'Resolution', 200);
end

% =========================================================================
% --- drawing helpers -----------------------------------------------------
function blk(ax, cx, cy, w, h, col, lines)
x = cx - w/2; y = cy - h/2;
rectangle(ax, 'Position', [x y w h], 'Curvature', [0.06 0.14], ...
    'FaceColor', col, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.0);
text(ax, cx, cy, lines, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 9);
end

function sblk(ax, cx, cy, w, h, col, lines)
x = cx - w/2; y = cy - h/2;
rectangle(ax, 'Position', [x y w h], 'Curvature', [0.35 0.55], ...
    'FaceColor', col, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);
text(ax, cx, cy, lines, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 9);
end

function arrow_poly(ax, xs, ys, label, label_seg)
% Draw polyline through (xs, ys) with a single arrowhead at the last vertex.
% Optional label is placed at the midpoint of segment index `label_seg`
% (1-based). If `label_seg` is omitted, the longest segment is used.
xs = xs(:).'; ys = ys(:).';
plot(ax, xs, ys, 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
ang = atan2(ys(end) - ys(end-1), xs(end) - xs(end-1));
L = 1.0;
x2 = xs(end); y2 = ys(end);
xa = x2 - L*cos(ang - pi/9);
ya = y2 - L*sin(ang - pi/9);
xb = x2 - L*cos(ang + pi/9);
yb = y2 - L*sin(ang + pi/9);
patch(ax, [x2 xa xb], [y2 ya yb], [0.15 0.15 0.15], 'EdgeColor', 'none');
if nargin >= 4 && ~isempty(label)
    if nargin < 5 || isempty(label_seg)
        seg_len = hypot(diff(xs), diff(ys));
        [~, label_seg] = max(seg_len);
    end
    mx = (xs(label_seg) + xs(label_seg + 1)) / 2;
    my = (ys(label_seg) + ys(label_seg + 1)) / 2;
    % offset label perpendicular to the segment so it doesn't sit on the line
    dx = xs(label_seg + 1) - xs(label_seg);
    dy = ys(label_seg + 1) - ys(label_seg);
    segL = hypot(dx, dy);
    if segL < eps
        ox = 0; oy = 1;
    else
        ox = -dy / segL; oy = dx / segL;
    end
    text(ax, mx + 0.8*ox, my + 0.8*oy + 0.5, label, ...
        'HorizontalAlignment', 'center', 'FontSize', 8, ...
        'Color', [0.0 0.2 0.5], 'Interpreter', 'tex');
end
end

function arrow_arc(ax, x1, y1, x2, y2, bulge, label)
% Quadratic bezier with explicit signed bulge perpendicular to chord.
t = linspace(0, 1, 30).';
dx = x2 - x1; dy = y2 - y1;
L = hypot(dx, dy);
if L < eps
    px = ones(size(t))*x1; py = ones(size(t))*y1;
else
    nx = -dy / L; ny = dx / L;
    cx = (x1 + x2)/2 + bulge * nx;
    cy = (y1 + y2)/2 + bulge * ny;
    px = (1-t).^2 * x1 + 2*(1-t).*t * cx + t.^2 * x2;
    py = (1-t).^2 * y1 + 2*(1-t).*t * cy + t.^2 * y2;
end
plot(ax, px, py, 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
ang = atan2(py(end) - py(end-1), px(end) - px(end-1));
L2 = 1.0;
xa = px(end) - L2*cos(ang - pi/9);
ya = py(end) - L2*sin(ang - pi/9);
xb = px(end) - L2*cos(ang + pi/9);
yb = py(end) - L2*sin(ang + pi/9);
patch(ax, [px(end) xa xb], [py(end) ya yb], [0.15 0.15 0.15], 'EdgeColor', 'none');
if nargin >= 7 && ~isempty(label)
    mid = round(numel(t)/2);
    text(ax, px(mid), py(mid) + sign(bulge)*1.2, label, ...
        'HorizontalAlignment', 'center', 'FontSize', 8, ...
        'Color', [0.0 0.2 0.5], 'Interpreter', 'tex');
end
end

function branch_dot(ax, x, y)
plot(ax, x, y, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerEdgeColor', 'none');
end

function legend_chip(ax, x, y, col, label)
rectangle(ax, 'Position', [x y 2.5 2.0], 'Curvature', [0.4 0.4], ...
    'FaceColor', col, 'EdgeColor', [0.25 0.25 0.25]);
text(ax, x + 3.2, y + 1.0, label, 'FontSize', 9, ...
    'VerticalAlignment', 'middle');
end
