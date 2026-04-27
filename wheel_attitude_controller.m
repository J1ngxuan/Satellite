function [Tcmd, Twheels, info] = wheel_attitude_controller(qe, we, P, health, T_ff)
% WHEEL_ATTITUDE_CONTROLLER  PD attitude controller with fault-tolerant
% control allocation across N reaction wheels.
%
% Inputs:
%   qe     : 4x1 quaternion error (commanded * inv(measured))
%   we     : 3x1 body-rate error [rad/s]
%   P      : satellite parameter struct (uses P.wheel.*, P.ctrl.*, P.J)
%   health : 1xN logical vector (true = wheel healthy). [] -> all healthy.
%   T_ff   : 3x1 feedforward body torque to cancel a known external
%            torque (e.g. uncancellable thruster residual during burn).
%            Defaults to 0.  Wheels will produce -T_ff plus PD action.
%
% Outputs:
%   Tcmd    : 3x1 commanded body torque produced by surviving wheels [Nm]
%   Twheels : Nx1 individual wheel torque commands [Nm]
%   info    : struct .feasible, .saturated, .rank
%
% Fault-tolerant allocation:
%   We seek wheel torques tau (Nx1) such that A*tau = T_des, where the
%   columns of A corresponding to failed wheels are zeroed. The minimum
%   2-norm solution that respects per-wheel torque limits is obtained via
%   weighted pseudo-inverse with subsequent saturation. If the surviving
%   axes are < 3, the controller projects T_des onto the achievable
%   subspace (range of A).

N = P.wheel.N;
if nargin<4 || isempty(health), health = true(1,N); end
if nargin<5 || isempty(T_ff),   T_ff   = zeros(3,1); end

% --- Outer-loop PD law in body frame + feedforward ---------------------
% Use vector-part of error quaternion (sign-correct shortest path)
qv = qe(2:4);
if qe(1) < 0, qv = -qv; end
T_des = -P.ctrl.Kp_att*qv - P.ctrl.Kd_att*we - T_ff(:);

% Saturate desired torque magnitude per axis to avoid runaway
T_max_axis = P.wheel.Tmax * sum(abs(P.wheel.axes),2);
T_des = max(min(T_des, T_max_axis), -T_max_axis);

% --- Allocation matrix with health gating ------------------------------
A = P.wheel.axes;            % 3xN
A(:, ~health) = 0;
r = rank(A, 1e-6);
info.rank = r;

if r >= 3
    % Standard weighted pseudo-inverse (minimum-norm wheel torque)
    Apinv = pinv(A);
    tau = Apinv * T_des;
    info.feasible = true;
else
    % Project desired torque onto reachable subspace (range of A)
    [U,~,~] = svd(A,'econ');
    Up = U(:,1:r);
    T_des_proj = Up*(Up.'*T_des);
    tau = pinv(A)*T_des_proj;
    info.feasible = false;
    T_des = T_des_proj;
end

% --- Per-wheel saturation -----------------------------------------------
sat = abs(tau) > P.wheel.Tmax;
info.saturated = any(sat);
tau = max(min(tau, P.wheel.Tmax), -P.wheel.Tmax);
tau(~health) = 0;

Twheels = tau;
Tcmd    = A*tau;
end
