function [q1, w1] = attitude_dynamics(q0, w0, Tctrl, Tdist, J, Jinv, dt)
% ATTITUDE_DYNAMICS  One-step RK4 propagation of rigid-body attitude.
%
% Inputs:
%   q0    : 4x1 quaternion [q0(scalar); qv] at t (unit norm)
%   w0    : 3x1 body angular velocity [rad/s]
%   Tctrl : 3x1 control torque on body [Nm]
%   Tdist : 3x1 disturbance torque [Nm]
%   J     : 3x3 inertia tensor [kg*m^2]
%   Jinv  : 3x3 inverse inertia
%   dt    : time step [s]
%
% Outputs:
%   q1    : updated quaternion (unit-normalised)
%   w1    : updated angular velocity
%
% Kinematics: qdot = 0.5 * Omega(w) * q
% Dynamics : Jdotw = T - w x (J w)
%
% RK4 on the augmented state [q;w].

y0 = [q0(:); w0(:)];
f  = @(y) att_ode(y, Tctrl, Tdist, J, Jinv);
k1 = f(y0);
k2 = f(y0 + 0.5*dt*k1);
k3 = f(y0 + 0.5*dt*k2);
k4 = f(y0 +     dt*k3);
y1 = y0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

q1 = y1(1:4)/norm(y1(1:4));
w1 = y1(5:7);
end

function dy = att_ode(y, Tctrl, Tdist, J, Jinv)
q = y(1:4); w = y(5:7);
qs = q(1); qv = q(2:4);
qdot = 0.5*[ -qv.'*w;
             qs*w + cross(qv,w) ];
wdot = Jinv*(Tctrl + Tdist - cross(w, J*w));
dy = [qdot; wdot];
end
