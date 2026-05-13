%% 平台基本信息 (保持不变)
sat_mass = 30;
sat_J = [0.288, 0, 0; 0, 0.45025, 0; 0, 0, 0.45025];
L_x = 0.350/2; L_y = 0.240/2; L_z = 0.240/2;

%% 12推进器配置方案
% 1-4号: X轴推力; 5-8号: Y轴推力; 9-12号: Z轴推力
% 奇数号和偶数号成对安装在面上的对称位置

% 推力方向矩阵 D_thr (3x12)
D_thr = [ 1, 1,-1,-1,  0, 0, 0, 0,  0, 0, 0, 0;  % X轴
          0, 0, 0, 0,  1, 1,-1,-1,  0, 0, 0, 0;  % Y轴
          0, 0, 0, 0,  0, 0, 0, 0,  1, 1,-1,-1]; % Z轴

% 力矩分配矩阵 B_thr (3x12)
% 以X轴推进器为例：1,2号在+X面，1号偏+Y产生-Z力矩，2号偏-Y产生+Z力矩
B_thr = [ 0,    0,    0,    0,   L_z, -L_z,  L_z, -L_z, -L_y,  L_y, -L_y,  L_y; % Roll (X)
         L_z, -L_z,  L_z, -L_z,   0,    0,    0,    0,    L_x, -L_x,  L_x, -L_x; % Pitch (Y)
        -L_y,  L_y, -L_y,  L_y,  L_x, -L_x,  L_x, -L_x,   0,    0,    0,    0]; % Yaw (Z)

% 飞轮布局 (保持金字塔构型)
alpha_rw = deg2rad(54.74);
B_rw = [cos(alpha_rw), 0, -cos(alpha_rw), 0;
        0, cos(alpha_rw), 0, -cos(alpha_rw);
        sin(alpha_rw), sin(alpha_rw), sin(alpha_rw), sin(alpha_rw)];
%% 故障识别算法参数预设
fault_threshold_low = 0.1;  % 失效判定阈值
rls_forgetting_factor = 0.98; % 最小二乘法遗忘因子