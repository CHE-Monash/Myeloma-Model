# Default analysis

The reference analysis. One model, run two ways, selected by `$scenario` in `simulate.do`:

| `$scenario` | `$coeffs` (fit on) | `$data` (simulate) | Purpose |
|---|---|---|---|
| `""` (blank) | `full` — the full (100%) MRDR registry | `synthetic` — the synthetic incidence population | **Projection**: current-practice projections + costs + PDF report |
| `outsample` | `train` — the 70% training fold | `test` — the held-out real 30% | **Out-of-sample validation**: compare predictions to those patients' observed outcomes |

`$coeffs` picks which patients the risk equations were fit on; `$data` picks which cohort is simulated. In-sample vs out-of-sample is just the fit×cohort combination — the mainstay validation is `train` fit predicting the `test` fold it never saw.

## Run it

```stata
// Projection (default):
do "analyses/default/simulate.do"

// Out-of-sample validation (point estimate):
do "analyses/default/simulate.do" 0 . . outsample
do "analyses/default/validate_outsample.do"
```

`run.do` is the full ordered runbook (both tracks, deterministic + the bootstrap HPC plumbing).

## Layout

```
analyses/default/
├── run.do                    Runbook: projection + out-of-sample, local + HPC
├── simulate.do               Dispatcher ($scenario selects projection / outsample)
├── validate_outsample.do     Point-estimate: simulated 30% vs observed targets
├── validate_outcomes.do      Shared comparison engine (targets vs simfile)
├── bootstrap_validation.do   HPC: 95% bootstrap prediction-interval coverage + OS curve
├── plot_os_curve.do          Redraw the OS validation curve locally from the CSV
├── outcomes/
│   ├── txr_full.do           Canonical per-line regimen list
│   └── txr_train.do          sources txr_full.do (identical regimens, no drift)
├── coefficients/             coefficients_full.mmat (100%) + coefficients_train.mmat (70%)
├── patients/                 patients_test.dta (held-out 30%; git-ignored, MRDR-restricted)
├── prep/                     split.do (70/30), test_cohort.do, test_targets.do  [need the MRDR drive]
├── targets/                  observed 30% target CSVs (aggregate; committed)
├── simulated/                run outputs (all_0_synthetic.dta; outsample/all_0_test.dta)
└── results/                  bootstrap coverage + OS curve
```

## Reproducibility / data

The **projection** side is reproducible from the committed code + `coefficients_full.mmat` + the synthetic
population (`patients/synthetic_1995_2040_*.dta`, regenerated from the MRDR-derived risk equations).

The **out-of-sample** side is *inspect-only*: the 70%/30% folds are real MRDR patients and cannot be shared
(git-ignored, restricted). What is public — the code, the 70%-trained `coefficients_train.mmat`, and the
aggregate `targets/*.csv` — lets the validation be audited, but not re-run without the restricted registry.
