# Repository Guidelines

## Project Structure & Module Organization

This repository contains a MATLAB fault-tolerant satellite attitude/orbit control study. Top-level `.m` files are the source modules and scripts:

- `main.m` runs the full simulation pipeline and produces interactive MATLAB plots.
- `satellite_params.m` defines shared platform, actuator, controller, sensor, and environment parameters.
- `thruster_configurations.m`, `thruster_ft_allocation.m`, and `wheel_attitude_controller.m` implement actuator layout and control allocation logic.
- `attitude_dynamics.m`, `orbit_dynamics.m`, and `sensor_model.m` provide reusable simulation primitives.
- `sim_wheel_failure.m` and `sim_thruster_fault.m` define fault scenarios.
- `export_report_figures.m` regenerates PNG assets under `figures/`.
- `research_target.md`, `research_report.md`, and `research_report_full.md` are project documentation; `reference/` stores source papers.

Keep new simulation modules at the repository root unless the project is reorganized consistently.

## Build, Test, and Development Commands

There is no build step. Run from MATLAB with this folder as the working directory or on the MATLAB path:

```matlab
main                    % full trade study, fault simulations, and plots
export_report_figures   % regenerate report figures in ./figures
```

For focused work, call modules directly:

```matlab
P = satellite_params();
out = sim_wheel_failure(P, 'two');
```

Use MATLAB Code Analyzer before committing edits:

```matlab
checkcode main.m
```

## Coding Style & Naming Conventions

Use MATLAB function files with one primary function per file. Name files and functions in lower snake case, matching the existing pattern (`sim_thruster_fault.m`, `orbit_dynamics.m`). Use 4-space indentation, descriptive variable names, and explicit units in comments or struct fields. Preserve scalar-first quaternion conventions and pass shared configuration through the `P` struct rather than hard-coding parameters in scenario code.

## Testing Guidelines

No formal test suite is currently present. Validate changes by running `main` end-to-end and any affected scenario directly. For numerical changes, compare key printed summaries and regenerated figures against prior behavior. Add focused MATLAB test scripts only when introducing reusable logic with clear expected outputs.

## Commit & Pull Request Guidelines

This checkout has no Git history available, so use concise imperative commit messages such as `Fix thruster allocation saturation` or `Regenerate report figures`. Pull requests should describe the changed scenario or model assumption, list MATLAB commands run, and include updated figures when outputs under `figures/` change.

## Agent-Specific Instructions

Do not tighten the relaxed thruster redundancy test or change quaternion error sign without documenting the control impact. Treat files in `figures/` as generated artifacts from `export_report_figures.m`.
