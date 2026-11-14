# Numerical reproduction

The primary entry is `run_supplement.m`. It resolves dependencies relative to its own location and restores the previous MATLAB path and working directory when execution ends. Use a writable working copy. Mechanical force and length units are N and m unless a column name specifies a display conversion. Input and reference arrays retain their stored precision.

## Environment

The code targets MATLAB R2024b and has been executed in MATLAB 24.2.0.2712019, Windows PCWIN64. The main mechanical, prediction and plotting paths use Base MATLAB and its Java virtual machine. PDF, SVG and native FIG export use MATLAB. EMF export is a Windows graphics operation. Other operating systems require a supported export substitute for EMF.

Optional Python implementation comparisons use SciPy and OpenPile. The recorded dependency set is `research_base/research/pisa_model/requirements_reference.txt`. The primary MATLAB analyses do not import these Python packages. Their separate installation and comparison commands are described in DEPENDENCIES.md.

## Commands

| Command | Work performed | Output |
|---|---|---|
| `run_supplement('checks')` | Verifies ten auxiliary input files and eleven mechanical sources, resolves four result pointers, checks 120 applicability records and solves two prescribed displacement states | `reproduction_runs/checks_*` and `research_base/research/QA/entrypoint_checks_*` |
| `run_supplement('smoke')` | Executes small PISA design, traditional clay and sand baseline, and field model/specimen flows | New `design_smoke_*` or `smoke_*` folders under the corresponding results/runs directories |
| `run_supplement('audit')` | Refits and recalibrates predictors, repeats group bootstrap calculations and all 120 policy runs, and performs 25 consistency and information-isolation checks | New `research_base/auxiliary_evaluation/runs/audit_*` |
| `run_supplement('figures')` | Regenerates quantitative figures from full-precision numerical reference results | `research_base/research/figures`, `research_base/mesh_verification/figures`, and `research_base/auxiliary_evaluation/plotting/*/figures` |
| `run_supplement('formal')` | Recalculates the main mechanical design, geometry mechanism, conventional baseline and field-comparison analyses and their figures | New main-model result folders and updated local current-result pointers |
| `run_supplement('generate_auxiliary')` | Generates the 260 development/calibration/test geometries and 128 candidate geometries, solves paired loads, refits predictors and reruns the policies | `auxiliary_generation/runs/generated_*` and newly generated solver cache entries |

The checks, smoke, audit and figures commands have been executed with the distributed MATLAB sources. All 25 auxiliary replay checks passed. The smoke flow completed one PISA design condition, two conventional benchmark conditions and four field model/specimen conditions without failures. The figure entry generated 33 quantitative group figures. Full-precision results of the formal mechanical campaign are supplied as numerical reference outputs.

The audit reuses the recorded measured development cost ratio and evaluates the same response arrays. It is an internal reproduction, not a new independent evaluation population. Auxiliary regeneration measures new costs. Changes in wall time can consequently change policy paths even when the mechanical responses agree.

## Data and paths

`research_base/research/pisa_model` contains the PISA reaction laws, beam interpolation, assembly, nonlinear solver and module checks. `research_base/research/code` contains the design comparison, geometric partial and total derivatives, adjoint channels and plotting routines. Traditional and field analyses occupy `design_baselines` and `field_validation`.

The four `current_*.json` result pointers use local relative paths. Their target directories contain full-precision reference results. Files under `sources/design_benchmarks`, `published_data` and `soil_profiles` contain numerical-curve extraction and interpreted-profile data. The soil-profile directory retains the vector-coordinate and physical-unit data used for extraction.

`research_base/auxiliary_evaluation/inputs/formal` contains the ten inputs listed in `input_manifest.csv`. The auxiliary source code operates on one row per geometry after the maximum over two prescribed load cases is taken. The 80 development, 60 calibration and 120 test geometries remain grouped. The separate finite candidate pool contains 128 geometries.

Full-precision response arrays are supplied as reference results. Auxiliary generation creates solver outputs and cache keys from the mechanical sources, configuration, geometry, treatment and loads. Cache reuse across MATLAB releases requires a separate numerical equivalence check.

## Targeted numerical checks

The source files `run_pisa_validation.m`, `run_pisa_crosschecks.m`, `run_cowden_profile_verification.m`, `compare_openpile_curves.m`, `verify_design_baselines.m`, `verify_field_numerics.m` and `run_targeted_mesh_review.m` provide focused checks. Add their containing directories to the MATLAB path and inspect their function signatures before use. They write fresh output in their stated local result directories.

The strict nonlinear residual criterion is 10^-9. The 80-element Clay D1 updated and frozen branches failed that criterion. The 40-element calculations retained the original criterion. Its 3×10^-9 retry is a separate diagnostic. Graphical reading envelopes, numerical mesh differences and physical comparison errors have distinct interpretations.

The main PISA derivatives assume fixed wall thickness, embedded length, load height and external load. Soil-law geometry is changed separately from actual beam geometry. The finite diameter comparison uses full nonlinear re-equilibration. Local adjoint sensitivity is evaluated at the nominal state. Calibration-domain and displacement-domain flags remain part of the result tables.

Reference MAT files retain the computed arrays and model parameters. Stored anonymous callbacks are represented by their mathematical expressions and captured parameters without local executable bindings. Executable callbacks are constructed by the supplied configuration functions during a new calculation.
