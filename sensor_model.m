function [q_meas, w_meas] = sensor_model(q_true, w_true, sensor, rng_state)
% SENSOR_MODEL  Synthesize noisy star-tracker + gyro measurements.
%
% Inputs:
%   q_true    : 4x1 true attitude quaternion
%   w_true    : 3x1 true body rate [rad/s]
%   sensor    : struct with fields star_sigma, gyro_sigma, gyro_bias
%   rng_state : optional persistent RNG state (unused; randn directly)
%
% Outputs:
%   q_meas    : measured quaternion (small-angle perturbed)
%   w_meas    : measured body rate (with bias + white noise)

if nargin < 4, rng_state = []; end %#ok<NASGU>

% Star tracker: per-axis small-angle error
theta = sensor.star_sigma * randn(3,1);
dq    = [1; 0.5*theta];
dq    = dq/norm(dq);
q_meas = qmult(q_true, dq);
q_meas = q_meas/norm(q_meas);

% Gyro: bias + white noise
w_meas = w_true + sensor.gyro_bias + sensor.gyro_sigma*randn(3,1);
end

function q = qmult(qa, qb)
% Hamilton product with scalar-first convention.
s1 = qa(1); v1 = qa(2:4);
s2 = qb(1); v2 = qb(2:4);
q  = [ s1*s2 - v1.'*v2;
       s1*v2 + s2*v1 + cross(v1,v2) ];
end
