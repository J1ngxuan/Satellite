function [r1, v1] = orbit_dynamics(r0, v0, a_thr, dt, env)
% ORBIT_DYNAMICS  RK4 propagate Keplerian + J2 + thruster acceleration.
%
% Inputs:
%   r0,v0 : 3x1 ECI position [m] and velocity [m/s]
%   a_thr : 3x1 thruster acceleration in ECI [m/s^2] (assumed constant
%           across dt; rotated from body frame upstream)
%   dt    : time step [s]
%   env   : struct with mu, Re, J2
%
% Outputs:
%   r1,v1 : updated ECI state
%
% Two-body + J2 perturbation:
%   a = -mu r / |r|^3 + a_J2 + a_thr

y0 = [r0(:); v0(:)];
f  = @(y) orb_ode(y, a_thr, env);
k1 = f(y0);
k2 = f(y0 + 0.5*dt*k1);
k3 = f(y0 + 0.5*dt*k2);
k4 = f(y0 +     dt*k3);
y1 = y0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
r1 = y1(1:3);
v1 = y1(4:6);
end

function dy = orb_ode(y, a_thr, env)
r = y(1:3); v = y(4:6);
rn  = norm(r);
ag  = -env.mu*r/rn^3;
% J2 perturbation in ECI (assuming z = polar axis)
zr2 = (r(3)/rn)^2;
fac = -1.5*env.J2*env.mu*env.Re^2/rn^5;
aJ2 = fac * [ r(1)*(1-5*zr2);
              r(2)*(1-5*zr2);
              r(3)*(3-5*zr2) ];
dy = [v; ag + aJ2 + a_thr];
end
