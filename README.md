# Satellite Fault-Tolerant Control Study

MATLAB simulation study for fault-tolerant satellite attitude and orbit control.
The repository contains reaction-wheel and thruster-fault scenarios, controller
and allocation logic, shared platform parameters, report text, references, and
generated report figures.

## Quick Start

Run MATLAB with this folder as the working directory or on the MATLAB path:

```matlab
main
```

Regenerate report figures:

```matlab
export_report_figures
```

Run a focused scenario:

```matlab
P = satellite_params();
out = sim_wheel_failure(P, 'two');
```

## Repository Notes

- Source MATLAB modules live at the repository root.
- `figures/` contains generated PNG report assets from `export_report_figures.m`.
- `reference/` stores source papers used by the report.
- Local tool settings such as `.claude/settings.local.json` are intentionally ignored.
