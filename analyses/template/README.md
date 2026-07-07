# Analysis template

A copy-me skeleton for a **new analysis**. It follows the repo's [analysis layout convention](../../README.md#analysis-layout-convention): a `simulate.do` dispatcher, a `run.do` runbook, and per-line regimen + override files under `outcomes/`.

## How to use

1. **Copy** this folder to `analyses/<your_analysis>/`.
2. **Rename** `outcomes/txr_template.do` → `outcomes/txr_<coeffs>.do` to match the `$coeffs` you'll set.
3. **Edit** `simulate.do`'s Configuration block (`$analysis`, `$int`, `$line`, `$coeffs`, `$data`, cohort…) and `run.do`'s steps.
4. **Fit** coefficients: `do "prep/risk_equations.do" <analysis> <coeffs> 1995 2040 0` → `coefficients/coefficients_<coeffs>.mmat`. It loads your `outcomes/txr_<coeffs>.do` regimen lists.
5. **Provide** the simulation cohort under `patients/`:
   - *Whole-population analysis* (like `base_model`) — point `$data` at a `patients/population_*.dta` cohort; nothing to build.
   - *Line-specific decision analysis* (like `transport_dvd`, assessing a drug at line L) — build a cohort that reaches line L: `patients/cohort_pool.do` (build the entry pool once) → `patients/draw_cohort.do` (draw a fixed-size, seeded cohort → `patients_<analysis>_<line>.dta`).
6. **Run**: `do "analyses/<your_analysis>/run.do"` (or `simulate.do` alone for one point-estimate run).

Keep the two `outcomes/sim_*_override.do` files only if you need them (see below); otherwise delete them — the engine uses the standard models when they're absent.

## Files

| File / folder | Role |
|---|---|
| `run.do` | Runbook — the full pipeline (fit → simulate → validate) + the bootstrap HPC plumbing. |
| `simulate.do` | Dispatcher — one run, configured by the globals block; optional positionals `boot min_bs max_bs [scenario]`. |
| `outcomes/txr_<coeffs>.do` | Per-line regimen code lists (`$TXR_L1..L9`); `gen_txr` in `prep/risk_equations.do` builds the dummies. |
| `outcomes/sim_bcr_override.do` | *(optional)* Replace the response (BCR) draw at `$line` — auto-run by `core/outcomes/sim_bcr.do`. |
| `outcomes/sim_txr_override.do` | *(optional)* Force the regimen at `$line` — auto-run by `core/outcomes/sim_txr.do`. |
| `coefficients/` | Fitted `coefficients_<coeffs>.mmat` (+ `bootstrap/`). |
| `patients/cohort_pool.do` | *(optional, line-specific)* Build once — simulate populations, keep line-L reachers → `cohort_pool_<line>.dta`. |
| `patients/draw_cohort.do` | *(optional, line-specific)* Draw a fixed-size, seeded cohort from the pool → `patients_<analysis>_<line>.dta`. |
| `patients/` | The simulation cohort `.dta` (built by the two scripts above, or a `population_*.dta` you point `$data` at). |
| `simulated/` | Run outputs (`bootstrap/` for resamples; a `<scenario>/` subfolder if you use scenarios). |
| `results/` | Analysis-level CSVs / figures / `results.md`. |

## The outcome-override mechanism

`core/outcomes/sim_bcr.do` and `sim_txr.do` run the **default** draw for every line, then — **only when the current `Line` equals `$line`** — auto-execute `outcomes/sim_bcr_override.do` / `sim_txr_override.do` **if the file exists**. Nothing to register: just set `$line` and drop the file in `outcomes/`. Both templates here are documented, inert-by-default skeletons; the full worked versions (trial-calibrated response for second-line DVd, scenario- and bootstrap-aware) are in `analyses/transport_dvd/outcomes/`.

## Line-specific cohorts (optional)

If your analysis assesses an intervention at a specific line L (e.g. a new second-line drug), you need a synthetic cohort of patients who *reach* line L with a realistic case-mix — not the whole incident population. Two scripts under `patients/` build it, mirroring `transport_dvd`:

1. **`cohort_pool.do`** — the expensive build-once step. Simulates each independent incident population from diagnosis, keeps everyone who reaches line L inside a fixed **case-mix window**, blanks the line-L regimen (so the pool is *arm-agnostic* — the dispatcher/override assigns intervention vs comparator), and rebases the diagnosis clock → `cohort_pool_<line>.dta`. The case-mix is a property of the *window*; don't widen it to change size.
2. **`draw_cohort.do`** — draws a fixed-size, fixed-seed sample from the pool → `patients_<analysis>_<line>.dta`. Cohort **size** is a precision choice, independent of case-mix.

**Choosing N:** `transport_dvd` sizes N by a Monte-Carlo-precision rule (MC SD of the incremental QALY ≤ a fraction of the parameter SD), via `patients/sample_size/ce_precision.do` (computes σ_pp) + `ce_sample_size.do` (O'Hagan ANOVA variance decomposition). Those are specific to a two-arm cost-effectiveness/PSA setup, so they're **not** in this template — copy them from `transport_dvd` if you need formal sizing; otherwise pick N large and check stability.

## Reference analyses

- `analyses/base_model/` — the simplest single-run analysis (no scenarios, no overrides).
- `analyses/transport_dvd/` — scenarios (`A_trial`/`B_transport`/`C_mrdr`), a coefficient-generation step, and both overrides.
- `analyses/oos/run.do` — the complete HPC bootstrap shell block (rsync + `sbatch` array jobs) to copy from.
