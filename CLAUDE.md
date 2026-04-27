# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

MATLAB implementation of a fault-tolerant attitude/orbit control study for a 30 kg microsatellite. The full requirement spec is in `research_target.md` (Chinese). The deliverables are: (i) a thruster-configuration trade study with ≥3 schemes, (ii) FT control under reaction-wheel and thruster faults, (iii) numerical simulations for each fault mode.

## How to run

There is no build step. From MATLAB with this folder on the path:

```matlab
main           % runs the full pipeline: config trade study + 3 fault sims + plots
```

Run a single module in isolation, e.g.:

```matlab
P = satellite_params();
configs = thruster_configurations(P.dim, P.thr.Fmax);  % just the trade study
out = sim_wheel_failure(P, 'two');                     % just one scenario
```

Use the matlab MCP tools (`run_matlab_file`, `evaluate_matlab_code`, `check_matlab_code`) to execute and lint without leaving the agent. The MATLAB desktop is visible to the user; figures appear in MATLAB's UI.

## Architecture

The code is structured as a chain of pure-ish modules driven by `main.m`. Data flows through a single parameter struct `P` (created by `satellite_params`) that carries inertia, sensor noise stats, controller gains, environment, and **interfaces for actuator configuration** — the reason for this design is requirement §(五): swapping in a new wheel or thruster layout must be a single struct edit, not a code change.

Pipeline (top-down):

1. **Configuration layer** — `thruster_configurations.m` builds candidate layouts (A: 6-orth, B: 8-canted, C: 12-orth) and scores each on (a) **redundancy** (drop each nozzle, check if the remaining columns still span all 6 axis force directions via `lsqnonneg`) and (b) **efficiency η = |net force|/Σ|F_i|** averaged over a 64-point Fibonacci sphere. `main.m` selects the redundant scheme with highest η, then writes its `dirs/pos` back into `P.thr` before any simulation runs.

2. **Dynamics** — `attitude_dynamics.m` (RK4 on `[q;ω]` with Euler's eqn) and `orbit_dynamics.m` (RK4, two-body + J2 + thruster acceleration). Quaternion convention is **scalar-first**, unit-normalised after every step.

3. **Sensing** — `sensor_model.m` adds star-tracker small-angle attitude noise and gyro bias + ARW.

4. **Controllers** — two parallel implementations of fault-tolerant control allocation:
   - `wheel_attitude_controller.m`: PD on the **inverse-quaternion error** (`qe = qmult(qinv(q_cmd), q_meas)` — sign matters; see history below) → 3×N weighted pseudo-inverse over wheel spin axes. When `rank(A) < 3` (≥2 wheels failed), it projects the desired torque onto the achievable subspace and signals `info.feasible = false` so the caller can engage thruster co-control.
   - `thruster_ft_allocation.m`: NNLS (`lsqnonneg`) with column gating (health) and per-column saturation (degradation scale), regularised by a small λI block to keep solutions sparse and bounded. Solves force + torque jointly so a single call serves both attitude-hold and Δv-on-arbitrary-direction.

5. **Scenarios** — `sim_wheel_failure.m` covers the three wheel-fault cases ('none' / 'one' / 'two'); the 'two' case automatically engages wheel+thruster co-control when the wheel allocator reports infeasible. `sim_thruster_fault.m` covers both 'degrade' (per-thruster scale < 1) and 'failure' (per-thruster health = false) modes during a translational burn while holding attitude — this is the "轨控不调姿" requirement. **Under thruster fault the wheels run the attitude loop**, with the thruster allocator's residual torque (`Tb_thr`) fed forward into the wheel controller as a known disturbance to cancel. Without this feedforward the PD loop would only generate the cancellation torque after the angle grew to ~50°, since uncancellable nozzle-asymmetry torques can be too large for the unaugmented PD steady-state.

## Conventions and gotchas

- **Quaternion error sign**: use `qe = qmult(qinv(q_cmd), q_meas)`, not `qmult(q_cmd, qinv(q_meas))`. The latter inverts the control sign and the spacecraft tumbles to ~180°.
- **Redundancy test is intentionally relaxed**: it only checks that the remaining force directions still span ±X/±Y/±Z, not that pure-translation-with-zero-torque is achievable from a non-negative combination. The strict test rejects scheme C even though it is operationally redundant — the attitude loop closes the residual torque. Don't tighten this without also re-architecting the allocator.
- **Δv burn target is sized to be achievable under fault** (0.3 m/s over 40 s ⇒ 0.225 N axial command vs. 0.5 N nozzle limit). Raising it past the surviving +X capacity will saturate the allocator and the failure case will under-deliver — this is a hardware limit, not a solver bug.
- **Comment ratio ≥10%** is a hard requirement from §(五). Each module's header documents inputs, outputs, units, and frame.
- All modules use `lsqnonneg` from base MATLAB — no toolbox dependency. `check_matlab_code` reports clean.
