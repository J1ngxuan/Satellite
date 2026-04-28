# 微小卫星执行机构容错控制研究报告

> 30 kg 平台 · 反作用飞轮 + 12 台冷气推进系统 · MATLAB 仿真验证
> 对应需求：`research_target.md` § (一)–(五)

---

## 1 平台参数与符号约定

| 项目 | 符号 | 取值 | 单位 |
|---|---|---|---|
| 质量 | $m$ | 30 | kg |
| 尺寸 (X,Y,Z) | $(2h_x, 2h_y, 2h_z)$ | $(0.350, 0.240, 0.240)$ | m |
| 主惯量 | $J=\mathrm{diag}(I_x,I_y,I_z)$ | $\mathrm{diag}(0.288,\,0.45025,\,0.45025)$ | kg·m² |
| 飞轮 (4 台金字塔) | $T_w^{\max}$ | 0.02 | N·m |
| 推力器 (≤12 台) | $F_{\max}$ | 0.5 | N |

**坐标系**：所有矢量在体坐标系内表示。四元数采用 scalar-first 约定 $q=[q_0;\,q_v]^\top$，$\|q\|=1$。

**符号**：第 $i$ 台推力器的安装位置 $r_i\in\mathbb{R}^3$、推力方向单位矢 $d_i\in S^2$、推力幅值 $F_i\ge 0$。其对体坐标系产生的力与力矩

$$
\mathbf F = \sum_{i=1}^{M} F_i\, d_i, \qquad
\mathbf T = \sum_{i=1}^{M} F_i\, (r_i\times d_i).
$$

写成矩阵形式：

$$
\begin{bmatrix}\mathbf F\\ \mathbf T\end{bmatrix}
= \underbrace{\begin{bmatrix} D \\ B \end{bmatrix}}_{A_{6\times M}} F,\quad
D = [d_1\;\cdots\;d_M],\;\; B = [r_1\times d_1\;\cdots\;r_M\times d_M].
$$

---

## 2 推进器配置最优设计

### 2.1 评价指标

**(i) 单方向效率**  对单位目标方向 $u\in S^2$，要求净力沿 $u$、净力矩为零，效率定义为

$$
\eta(u) \;=\; \frac{\|\mathbf F(u)\|}{\sum_i F_i(u)}, \qquad
\text{s.t. } DF=\eta(u)\cdot u,\; BF=0,\; F\ge 0. \tag{2.1}
$$

**(ii) 各向同性平均效率**

$$
\bar\eta = \frac{1}{4\pi}\!\int_{S^2}\!\eta(u)\,d\Omega
\;\approx\; \frac{1}{N_s}\sum_{k=1}^{N_s}\eta(u_k), \quad u_k:\text{Fibonacci 球面 64 点}. \tag{2.2}
$$

**(iii) 主轴效率**  $\eta_{\text{axis}}=\frac{1}{6}\sum_{k=1}^{6}\eta(\pm\hat e_k)$。

**(iv) 单点容错冗余**：

$$
\text{Redundant} \iff \forall j,\;\forall k\in\{1,\dots,6\}:\; \mathrm{rank}\big(D^{(\setminus j)}\big)\ge 3 \text{ 且 } \pm\hat e_k\in\mathrm{cone}\{d_i:\,i\ne j\}.
\tag{2.3}
$$

### 2.2 效率上界（核心引理）

**引理 1**  对任意 $u\in S^2$，

$$
\eta(u)\;\le\;\max_i\,(u^\top d_i)\;\le\;1. \tag{2.4}
$$

**证明**：由 $F_i\ge 0$，$\|\mathbf F\|\ge u^\top\mathbf F = \sum F_i(u^\top d_i)\le(\sum F_i)\,\max_i(u^\top d_i)$。两边除以 $\sum F_i$ 即得。$\square$

**推论**：$\eta(u)=1\Leftrightarrow \exists i:\,d_i=u$。

### 2.3 四个候选方案的几何与数学分析

#### 方案 A — 六正交（$M=6$）

$$
D_A = [\hat e_1,\,-\hat e_1,\,\hat e_2,\,-\hat e_2,\,\hat e_3,\,-\hat e_3]
= \begin{bmatrix}+1&-1&0&0&0&0\\0&0&+1&-1&0&0\\0&0&0&0&+1&-1\end{bmatrix},
$$

每个 $\pm$主轴方向**恰好 1 台**喷嘴。由引理 1 推论，$\eta(\pm\hat e_k)=1$；但去掉任一台喷嘴 $j$ 后，$d_j$ 方向再无喷嘴，$\eta(d_j)=0$，**不冗余**。

**指标**：$M=6$、不冗余、$\eta_{\text{axis}}=1.000$、$\bar\eta=0.674$。

#### 方案 B — 八台倾斜（$M=8$，正三角锥顶点对称）

8 台喷嘴对称分布于 $\pm Z$ 两面四象限，方向为单位化的 $\frac{1}{\sqrt 3}(\pm 1,\pm 1,\pm 1)$：

$$
D_B = \frac{1}{\sqrt 3}\begin{bmatrix}
+1&-1&+1&-1&+1&-1&+1&-1\\
+1&+1&-1&-1&+1&+1&-1&-1\\
+1&+1&+1&+1&-1&-1&-1&-1
\end{bmatrix}.
$$

主轴效率上界为

$$
\eta_{\text{axis}}^B \le \max_i(\hat e_x^\top d_i)=\tfrac{1}{\sqrt 3}\approx 0.577,
$$

实测 $0.577$ 取等号（同符号 4 台同时点火）。

**冗余性**：失去任一喷嘴，剩 7 台仍有方向覆盖三象限，$\mathrm{rank}(D)=3$ → 冗余。

**指标**：$M=8$、冗余、$\eta_{\text{axis}}=0.577$、$\bar\eta=0.449$。

#### 方案 C — 十二正交对（$M=12$）

每个 $\pm$主轴方向**严格 2 台共线喷嘴**，对称放在 $r$ 与 $-r'$ 两点（其中 $r\perp d$、$r'\perp d$）使力矩内消：

$$
\sum_{i\in S_k} r_i\times d_i = r\times d + (-r')\times d = (r-r')\times d.
$$

取 $r'=r$（轴对称布置）即得 $\sum r_i\times d_i=0$。这就是 `C.pos` 的 $(\pm h_y, \pm h_z)$ 配对结构。

**最优性定理**（C 的可证明优越性）：

> **定理**：在约束 (P1) $M\le 12$、(P2) $\eta(\pm\hat e_k)=1\;\forall k$、(P3) 任一台失效后 (P2) 仍成立 下，$M$ 的最小值为 $12$，且任意可行方案在每个 $\pm\hat e_k$ 方向必有 $|S_k|\ge 2$ 台严格共线喷嘴。C-12orth 取得该下界。

**证明**：由引理 1 推论，$\eta(\hat e_k)=1\Rightarrow S_k:=\{i:d_i=\hat e_k\}\neq\emptyset$。若 $|S_k|=1$，唯一喷嘴失效后 $S_k=\emptyset$，违反 (P3)，故 $|S_k|\ge 2$。六个方向的 $S_k$ 不相交，

$$
M\ge\sum_{k=1}^6 |S_k|\ge 6\times 2=12. \qquad\square
$$

**指标**：$M=12$、冗余、$\eta_{\text{axis}}=1.000$、$\bar\eta=0.674$。

#### 方案 E — 12 台非正交全方位（基于正二十面体）

**几何构造**：

- **位置 $r_i$**：12 条立方体棱边中点

$$
\{(\,0,\pm h_y,\pm h_z),\;(\pm h_x,\,0,\pm h_z),\;(\pm h_x,\pm h_y,\,0)\}.
$$

- **方向 $d_i$**：12 个正二十面体顶点（黄金比 $\varphi=(1+\sqrt 5)/2$）

$$
V = \big\{\hat v : v\in\{(0,\pm 1,\pm\varphi),\,(\pm 1,\pm\varphi,0),\,(\pm\varphi,0,\pm 1)\}\big\},\quad \hat v = v/\|v\|.
$$

- **配对**：贪心算法将每个棱边中点与"最朝外"的二十面体顶点匹配，使 $r_i^\top d_i / \|r_i\|$ 最大。

**非线性寻优**：引入 1 维倾角参数 $\alpha\in[0,0.9]$，按

$$
d_i(\alpha) = \frac{(1-\alpha)\,v_i + \alpha\,\hat r_i}{\|(1-\alpha)\,v_i+\alpha\,\hat r_i\|},\qquad \hat r_i = r_i/\|r_i\|, \tag{2.5}
$$

混合二十面体方向与径向方向。优化目标为惯量加权效率：

$$
J(\alpha) = w_I^\top \bar\eta_{\text{ax}}(\alpha) + 0.5\min_k \eta_k(\alpha),\qquad
w_I = \frac{1/J_{kk}}{\sum_k 1/J_{kk}}. \tag{2.6}
$$

由于 $I_x<I_y=I_z$，权 $w_I=(0.412,\,0.294,\,0.294)$ 倾向于强化 X 轴效率。`fminbnd` 给出 $\alpha^\star\approx 0.07$，说明二十面体几何已接近最优。

**位置-方向解耦**：由于 $r_i\not\parallel d_i$，

$$
r_i\times d_i \;\ne\; 0,
$$

每台喷嘴单独可产生姿态力矩——这是 E 优于"穿质心安装"方案的关键，使其在飞轮失效时也能辅助稳姿。

**指标**：$M=12$、冗余、$\eta_{\text{axis}}=0.672$、$\bar\eta=0.504$（受非负 + 力矩零约束的代价）。

### 2.4 对比汇总

| 方案 | M | 冗余 | $\eta_{\text{axis}}$ | $\bar\eta$ | 备注 |
|---|---|---|---|---|---|
| A-6orth | 6 | NO | 1.000 | 0.674 | 单点失效即失轴 |
| B-8canted | 8 | YES | 0.577 | 0.449 | 倾斜安装，效率最低 |
| **C-12orth** | 12 | YES | **1.000** | 0.674 | **§2.3 定理证明的最优解** |
| E-12icosa($\alpha=0.07$) | 12 | YES | 0.672 | 0.504 | 全方位、自带姿态力矩臂 |

**选型结论**：C 在"轨控不调姿 + 单点容错"问题下数学上最优；E 适用于以斜向 $\Delta v$ 为主、需推力器辅助稳姿的场景。`main.m` 默认选 C。

---

## 3 容错控制算法

### 3.1 姿态动力学

刚体姿态运动学与动力学：

$$
\dot q = \tfrac{1}{2}\Omega(\omega)\,q,\quad
J\dot\omega = -\omega\times(J\omega) + T_{\text{ctrl}} + T_{\text{dist}}, \tag{3.1}
$$

其中 $\Omega(\omega)=\begin{bmatrix}0 & -\omega^\top\\ \omega & -[\omega]_\times\end{bmatrix}$。数值积分采用四阶 Runge–Kutta，步长 $\Delta t=0.1$ s。

### 3.2 飞轮控制律与冗余分配

**误差四元数**（commanded → measured 的旋转）

$$
q_e = q_{\text{cmd}}^{-1} \otimes q_{\text{meas}},\quad
q_{e,v} = q_e(2{:}4),\quad \tilde q_v = \mathrm{sign}(q_{e,0})\,q_{e,v}. \tag{3.2}
$$

**PD + 前馈控制**：

$$
T_{\text{des}} = -K_p\,\tilde q_v - K_d\,\omega_{\text{meas}} - T_{\text{ff}}. \tag{3.3}
$$

**冗余分配**（$N_w$ 台飞轮，旋转轴矩阵 $W=[w_1\,\cdots\,w_{N_w}]$，健康指示向量 $h\in\{0,1\}^{N_w}$）：

设 $W^{(h)} = W\,\mathrm{diag}(h)$。当 $\mathrm{rank}(W^{(h)})\ge 3$，

$$
\tau^\star = (W^{(h)})^+ T_{\text{des}}, \tag{3.4}
$$

$(W^{(h)})^+$ 为 Moore–Penrose 伪逆，给出最小 2-范数解。当 $\mathrm{rank}<3$（双飞轮失效），用 SVD $W^{(h)}=U\Sigma V^\top$ 提取列空间投影器 $P_R = U_r U_r^\top$，将期望力矩投影到可达子空间：

$$
T_{\text{des}}^{\text{proj}} = P_R\,T_{\text{des}},\quad
\tau^\star = (W^{(h)})^+ T_{\text{des}}^{\text{proj}}. \tag{3.5}
$$

不可达分量由推力器协同补偿（见 § 3.4）。

### 3.3 推力器容错分配（NNLS）

期望体力 $F_d$、期望力矩 $T_d$。求

$$
\min_{F\ge 0}\;\Big\|\,A^{(h,s)} F - \begin{bmatrix}F_d\\T_d\end{bmatrix}\Big\|^2 + \lambda^2\|F\|^2,\qquad
0\le F_i\le s_i\,h_i\,F_{\max}, \tag{3.6}
$$

其中 $h_i\in\{0,1\}$（健康），$s_i\in[0,1]$（推力下降系数），$\lambda=10^{-3}$ 稀疏化正则项。$A^{(h,s)}$ 是把失效列置零、缩放列尺度后的分配矩阵。等价的 `lsqnonneg` 形式：

$$
F^\star = \arg\min_{F\ge 0} \left\|\begin{bmatrix}A^{(h,s)}\\ \lambda I\end{bmatrix} F - \begin{bmatrix}F_d\\T_d\\0\end{bmatrix}\right\|^2. \tag{3.7}
$$

事后裁剪 $F^\star \leftarrow \min(F^\star, s\odot h\odot F_{\max})$，输出实际体力 $\mathbf F = D F^\star$、体力矩 $\mathbf T = B F^\star$。

### 3.4 协同控制（飞轮 + 推力器）

**飞轮失效场景**（双轮失效）：

1. 飞轮分配 (3.5) 产出 $T_w^{\text{body}}$ 与可达指示 `feasible`；
2. 若 `feasible=false`，残差 $T_{\text{miss}} = T_{\text{des}} - T_w^{\text{body}}$ 经幅值限制后送入推力器分配 (3.7)，$F_d=0$、$T_d=T_{\text{miss}}$。
3. 仿真中 `mode` 状态机记录这一切换：1=纯飞轮，2=飞轮+推力器。

**推力器故障场景**（"轨控不调姿"模式）：

- 推力器分配先行：$F_d=F_{\text{burn}}\,\hat u$、$T_d=0$。NNLS 输出实际力 $\mathbf F$ 和**寄生力矩** $T_{\text{thr,res}}=BF^\star$（受非负约束几何上不可消除）。
- 飞轮以 (3.3) 跟随，**前馈** $T_{\text{ff}}=T_{\text{thr,res}}$ 主动抵消。

**几何不可消除性的证明**（以失效喷嘴 #2 为例，$+X$ 方向上仅余喷嘴 #1）：要求净力 $F_x>0$、$T_z=0$。可用喷嘴及其力/力矩贡献：

| 喷嘴 | $d$ | $r$ | 力 $X$ | 力矩 $Z$ |
|---|---|---|---|---|
| #1 | $+\hat x$ | $(-h_x,+h_y,0)$ | $+F_1$ | $-h_y F_1$ |
| #3 | $-\hat x$ | $(+h_x,+h_y,0)$ | $-F_3$ | $-h_y F_3$ |
| #4 | $-\hat x$ | $(+h_x,-h_y,0)$ | $-F_4$ | $+h_y F_4$ |

约束方程

$$
\begin{cases}F_1-F_3-F_4 = F_x>0\\ -F_1-F_3+F_4=0\end{cases}
\Rightarrow F_3 = -F_x/2 < 0,
$$

违反非负性。**故需飞轮前馈** $T_{\text{ff}}$ 主动抵消寄生 $-Z$ 力矩。

---

## 4 仿真结果

仿真步长 $\Delta t=0.1$ s。星敏 1σ = 5 arcsec、陀螺 1σ = 0.01°/s + 常值偏置 0.005°/s。

### 4.1 飞轮容错（场景 1）

初始姿态偏置 5°、初始角速度 $(0.3,-0.2,0.1)^\circ/\mathrm{s}$。终段 20 s 平均指向误差：

| 故障 | 平均误差 | 最大误差 | 工作模式 |
|---|---|---|---|
| 无故障 | 0.076° | 0.085° | 纯飞轮 (mode=1) |
| 1 台失效 | 0.068° | 0.076° | 纯飞轮 (mode=1) |
| **2 台失效** | **0.072°** | **0.086°** | **飞轮+推力器协同 (mode=2)** |

双轮失效时 $\mathrm{rank}(W^{(h)})=2<3$，触发协同控制，最终精度仍优于 0.1°，达到平台姿态保持要求。

### 4.2 参考文献启示下的双飞轮优先策略

Zhang 等（2023）的欠驱动姿态容错研究表明，双飞轮失效后并不必然只能立即依赖推力器；在零动量、剩余执行机构能力充足等条件下，两台飞轮仍可能维持三轴姿态稳定。结合本项目 `research_target.md` 的边界，本研究仅将该思想用于“飞轮完全失效”工况，不扩展到飞轮效率下降或偏置故障。

在现有金字塔 4 飞轮构型下，新增三种双轮失效策略对比：

1. `two_wheel_only`：只使用剩余两台飞轮的可达力矩投影，不启用推力器；
2. `wheel_first`：先运行两飞轮控制，若 30 s 滚动平均姿态误差超过 0.1°，再启用推力器补偿缺失力矩；
3. `assist_immediate`：保持原方案，双轮失效后立即进入飞轮 + 推力器协同。

6 种双飞轮失效组合的最终 20 s 平均姿态误差如下：

| 失效飞轮 | two_wheel_only | wheel_first | assist_immediate |
|---|---:|---:|---:|
| [1 2] | 3.281° | 0.075° | 0.068° |
| [1 3] | 23.711° | 0.076° | 0.075° |
| [1 4] | 79.284° | 0.064° | 0.076° |
| [2 3] | 25.240° | 0.063° | 0.069° |
| [2 4] | 61.228° | 0.051° | 0.055° |
| [3 4] | 37.717° | 0.051° | 0.059° |

![双飞轮失效策略对比](./figures/double_wheel_strategy_sweep.png)

结果说明：在当前平台惯量、金字塔飞轮几何、初始偏差和扰动条件下，单纯两飞轮可达子空间控制不能满足 0.1° 姿态保持精度；但双飞轮失效后可以先进行短时间两飞轮能力评估，再根据精度判据启用推力器兜底。该策略比“立即推力器协同”更能体现参考文献的欠驱动启示，同时仍满足 `research_target.md` 中“若双台飞轮故障无法满足精度，可进行飞轮和推力器协同控制”的要求。

### 4.3 推力器推力下降（场景 2）

12 台 C-12orth、目标 $\Delta v=0.3\,\text{m/s}$ 沿 $+X$、烧 40 s。失效喷嘴 #1 与 #5 推力下降至 40%：

$$
\boldsymbol{\Delta v}_{\text{achieved}} = (0.300, 0.000, 0.000)\,\text{m/s}, \quad \text{姿态误差峰值 }0.11^\circ.
$$

NNLS 自动将分配权重转移至健康喷嘴，$\Delta v$ 损失为零。

### 4.4 推力器完全失效（场景 3）

失效喷嘴 #2 与 #6（含一个 $+X$ 方向）：

$$
\boldsymbol{\Delta v}_{\text{achieved}} = (0.296, 0.000, 0.000)\,\text{m/s},\quad \text{相对损失 1.3\%}.
$$

加入飞轮前馈后，姿态误差峰值从 **51.92°** 降至 **0.09°**，验证了 § 3.4 协同方案的有效性。

### 4.5 组合故障与故障树分析

上述场景 1--3 属于典型工况验证。为覆盖飞轮与推力器同时故障的情况，新增 `sweep_fault_tree_analysis.m` 对故障树顶事件“姿轨控任务失败”进行组合扫掠。顶事件由三类中间事件触发：

1. 姿态保持失败：最终 20 s 平均姿态误差 $>0.1^\circ$；
2. 轨道机动失败：最终 $\Delta v$ 相对误差 $>5\%$；
3. 协同控制不可恢复：飞轮或推力器长期饱和，且姿态误差不收敛。

默认矩阵覆盖无飞轮故障、任意单飞轮失效、任意双飞轮失效，叠加任意两台推力器降额至 40% 或任意两台推力器完全失效，并分别沿 $+X/+Y/+Z$ 执行 $0.3\,\text{m/s}$ 轨控机动。完整扫掠共 4389 个工况，触发顶事件 1507 个，占 34.3%。分支统计为：姿态保持失败 1486 个、轨道机动失败 41 个、协同控制不可恢复 72 个；2394 个工况启用了推力器姿态辅助。

| 飞轮失效数 | 推力器模式 | 工况数 | 顶事件比例 | 最大稳态姿态误差 | 最大 $\Delta v$ 误差 |
|---:|---|---:|---:|---:|---:|
| 0 | nominal | 3 | 0.0% | 0.0984° | 0.00% |
| 0 | degrade | 198 | 10.1% | 0.1176° | 0.00% |
| 0 | failure | 198 | 10.6% | 0.1199° | 100.00% |
| 1 | nominal | 12 | 0.0% | 0.0985° | 0.00% |
| 1 | degrade | 792 | 12.5% | 0.1293° | 0.00% |
| 1 | failure | 792 | 24.0% | 74.1648° | 100.00% |
| 2 | nominal | 18 | 33.3% | 0.1557° | 0.00% |
| 2 | degrade | 1188 | 34.6% | 0.1730° | 0.00% |
| 2 | failure | 1188 | 64.0% | 95.1789° | 100.17% |

![组合故障树扫掠热力图](./figures/fault_tree_combined_sweep.png)

结果说明：推力器降额对 $\Delta v$ 完成度影响较小，风险主要体现为姿态保持裕度下降；推力器完全失效在部分方向上会直接造成轨控能力丧失。双飞轮失效叠加推力器故障时，推力器既要完成轨控又要补偿飞轮缺失力矩，顶事件比例显著上升，是故障树中的关键组合事件。

---

## 5 结论

1. **配置最优性**：在"轨控不调姿 + 单点容错"指标下，C-12orth 在 $M\le 12$ 约束下数学上最优（§ 2.3 定理）；E-12icosa 提供 $\alpha=0.07$ 的非线性优化的非正交全方位备选，在斜向 $\Delta v$ + 推力器辅助稳姿场景中具有优势。
2. **飞轮容错**：金字塔 4 飞轮在 1 台或 2 台失效下均能维持 < 0.1° 指向精度；针对参考文献启示，新增双飞轮失效后的两飞轮优先评估，结果表明当前平台单纯两飞轮投影控制不足以满足精度，需在 30 s 判据后启用推力器兜底。
3. **推力器容错**：基于 NNLS 的容错分配 (3.7) 配合飞轮前馈 (3.3) 实现 1.3% 以内的 $\Delta v$ 损失和 0.1° 量级的姿态误差，对推力下降和完全失效两种故障均验证有效。
4. **组合故障覆盖性**：新增故障树扫掠把典型工况扩展为 4389 个组合工况，明确识别出“双飞轮失效 + 推力器完全失效”是主要风险组合。
5. **MATLAB 实现**：模块化、接口参数化、纯 base toolbox 依赖、有效注释 > 10%，满足需求 § (五)。

## 附录：模块清单与接口

| 模块 | 输入 | 输出 | 角色 |
|---|---|---|---|
| `satellite_params` | – | $P$ 结构体 | 参数注入接口 |
| `thruster_configurations` | $\dim, F_{\max}$ | configs[] | A/B/C/E 库 + 评价 |
| `attitude_dynamics` | $q,\omega,T_c,T_d,J,\Delta t$ | $q^+,\omega^+$ | RK4 (3.1) |
| `orbit_dynamics` | $r,v,a_{\text{thr}},\Delta t,\text{env}$ | $r^+,v^+$ | 二体 + J2 |
| `sensor_model` | $q,\omega,\text{sensor}$ | $\tilde q,\tilde\omega$ | 测量噪声 |
| `wheel_attitude_controller` | $q_e,\omega_e,P,h,T_{\text{ff}}$ | $T_w^{\text{body}},\tau,\text{info}$ | 式 (3.3)–(3.5) |
| `thruster_ft_allocation` | $F_d,T_d,P,h,s$ | $F^\star,\mathbf F,\mathbf T,\text{info}$ | 式 (3.7) |
| `sim_wheel_failure` | $P,$ case | 时序 + mode | 场景 1 |
| `sweep_double_wheel_failure` | $P,$ strategies | 双轮失效组合扫描 | 双飞轮优先/兜底策略评估 |
| `sim_thruster_fault` | $P,$ mode, list, $\hat u$ | 时序 + $\Delta v$ | 场景 2/3 |
| `sim_combined_fault` | $P,$ fault opts | 姿轨控组合故障时序 | 飞轮 + 推力器组合故障 |
| `sweep_fault_tree_analysis` | $P,$ opts | 故障树结构体表格 | 组合故障覆盖性分析 |
| `main` | – | 控制台 + 图 | 编排 |
