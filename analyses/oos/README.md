# analyses/oos — out-of-sample (70/30) validation

Trains the model on a random **70%** of MRDR patients, predicts the held-out **30%**, and compares the predictions to those patients' **observed** outcomes. This tests *generalisation* (not just in-sample fit) and refreshes the 70/30 validation published in the 2024 PLOS ONE paper for the current model.

The split is made **before imputation** and the two folds are imputed **separately**, so the held-out 30% never informs the fit (no leakage).

This is the model's **mainstay validation, re-run for each version** — it replaced the in-sample registry acceptance test (the old `validation/` folder was retired to `scratch/`; its comparison engine now lives here as `validate_outcomes.do`).

## Layout
```
analyses/oos/
├── README.md
├── run.do                   analysis runbook: prep -> simulate -> validate (deterministic + bootstrap record)
├── simulate.do              simulation dispatcher (70%-trained coeffs on the 30% cohort)
├── validate_oos.do          point-estimate: compares simulated 30% vs observed targets
├── validate_outcomes.do     shared comparison engine: sim .dta vs target CSVs  [moved from validation/]
├── bootstrap_validation.do  HPC: 95% bootstrap prediction-interval coverage vs held-out (percentile method)
├── prep/
│   ├── oos_split.do         70/30 patient split (fixed seed) -> ${data_path}/oos/oos_split.dta
│   ├── oos_targets.do       observed outcomes for the 30%  -> targets/*.csv  (reuses generate_benchmarks.do)
│   └── oos_cohort.do        30% baseline covariates as a sim cohort -> patients/oos_cohort.dta
├── outcomes/txr_oos.do      regimen definitions (identical to base_model)
├── coefficients/            70%-trained coefficients (+ bootstrap/)   [built by run.do step 2]
├── patients/                oos_cohort.dta                            [built by run.do step 4]
├── targets/                 observed 30% target CSVs                  [built by run.do step 3]
├── simulated/               simulation output (+ bootstrap/)          [built by simulate.do]
└── results/                 bootstrap PI coverage report + csv        [built by bootstrap_validation.do]
```

## How it reuses the main pipeline
The main prep/validation scripts take optional arguments so they serve both the main model and OOS, unchanged when the argument is empty:

- `prep/multiple_imputation.do` — 5th arg `$sample` (`train`/`test`): filters to that fold (via the split crosswalk) and writes under `${data_path}/oos/`.
- `prep/risk_equations.do` — 8th arg `$sample`: reads the OOS-fold imputed data; output is keyed by `$analysis` (so `oos` lands here in `coefficients/`).
- `prep/generate_benchmarks.do` — `$bench_in` / `$bench_out`: `oos_targets.do` points it at the test fold and `targets/`.
- `validate_outcomes.do` (the shared comparison engine, now in `analyses/oos/`) — `$val_targets` / `$val_simfile`: `validate_oos.do` points it at `targets/` and the OOS simulation.

## Workflow
`run.do` is the runbook — the full sequence in order, split into a **deterministic** track (runs locally, top to bottom) and a **bootstrap** track (heavy; the MI/risk-equation/simulation steps run on the HPC and are recorded but left commented). Run from the repository root.

**Deterministic (point estimate):**
1. `prep/oos_split.do` — split once (fixed seed).
2. `multiple_imputation.do … train` and `… test`.
3. `risk_equations.do oos oos 1995 2040 0 . . train`.
4. `prep/oos_targets.do` — observed 30% targets.
5. `prep/oos_cohort.do` — 30% simulation cohort.
6. `simulate.do` — simulate the 30% (point estimate).
7. `validate_oos.do` — compare the point estimate to the observed targets (fixed tolerances).

**Bootstrap prediction intervals (HPC; the headline metric):** train bootstrap MI + risk equations, then `simulate.do 1 1 500` (500 simulations of the held-out 30%), then `bootstrap_validation.do` — 95% percentile interval per outcome, coverage vs the held-out observed. See `run.do` bootstrap steps (a)–(c) for the commands, and the `/* … */` shell block at the end of `run.do` for the cluster transfer/submit — run those lines in a VS Code terminal after `source hpc/env.sh` (git-ignored; supplies the machine-specific paths).

## Open items before a publication-grade run
- **`oos_cohort.do`** mirrors `population_1995_2040.do`'s cohort schema — verify the column set still matches `core/load_patients.do` / `core/mata_setup.do` against live data.
- **Prediction-interval calibration** (the headline OOS metric) is implemented in `bootstrap_validation.do` (percentile method: does each held-out observed value fall in the bootstrap 95% interval?). It needs the 500 bootstrap simulations (`simulate.do 1 1 500`); designed to run on the HPC, which keeps the bootstrap outputs — pull back `results/oos_bootstrap_validation.md` + `.csv` for review. Smoke-tested locally (logic verified on a single sim replicated 5×); not yet run over genuine bootstraps.
- **Match the 2024 methods**: confirm the split (simple vs stratified) and per-fold imputation choices against the PLOS ONE paper for comparability.
- All scripts are **untested against live MRDR data** (built without drive access) — expect a debugging pass on first run.
