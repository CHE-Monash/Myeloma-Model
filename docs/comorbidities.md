# Comorbidities in the model — revision log

**Purpose.** Track an in-progress change to how patient comorbidities enter the model. Point a future
session here to understand what changed, why, and what is still outstanding.

## What is changing

Historically the only comorbidity term in the model was **`CMc`** — an ordinal 0/1/2/3+ *count* of
comorbidities (`CM_CKD + CM_CRD + CM_PLM + CM_DBT + CM_MLG`, capped at 3+) — and it appeared in **only
the two ASCT logits**. We are moving to the **four individual comorbidity flags** as covariates:

- `CM_CKD` — chronic kidney disease (eGFR ≤ 59)
- `CM_CRD` — moderate-to-severe cardiac disease
- `CM_PLM` — moderate-to-severe chronic lung disease
- `CM_DBT` — diabetes (insulin / oral hypoglycaemics)

**Why the four flags over the `CMc` composite.** (1) They reveal per-condition heterogeneity the count
hides — e.g. in the ASCT decision, lung disease weighs ~2× kidney disease, which `CMc`'s equal-weight
assumption flattens. (2) Internal consistency: the same four drive OS and treatment assignment. (3) It
lets `CMc` be retired.

## Where they now enter (this pass — three equations)

Fitted in `prep/risk_equations.do`:
- **OS** (`streg … b0.OS#b5.BCR CM_CKD CM_CRD CM_PLM CM_DBT`) — the four **added** (OS never used `CMc`).
  In-sample, all four are significant and correctly signed (HR ≈ 1.19–1.32), and `CM_CKD` survives
  conditional on `RISS`. (An earlier variant also carried a `prev_dur` frailty term ahead of the four
  flags; `prev_dur` was later removed — see `os_line_specific_experiment.md` §7 — so it is **not** in
  the current OS equation and the `bOS` layout below reflects its absence.)
- **ASCT at diagnosis** (`DN_SCT` logit) — `i.CMc` **replaced** by the four flags.
- **ASCT at L1 end** (`L1_SCT` logit) — `i.CMc` **replaced** by the four flags.

Engine plumbing (so the simulation runs with the four flags):
- `core/mata_setup.do` — loads `vCRD`/`vPLM`/`vDBT` alongside the existing `vCKD`.
- `core/outcomes/sim_os.do` — `coefCols`/design matrix: `CM_CKD`=58, `CM_CRD`=59, `CM_PLM`=60,
  `CM_DBT`=61, `_cons`=62 (`bOS` now 1×63 incl. `ln_p` at 63; `aux = bOS[1, cols(bOS)]`).
- `core/outcomes/sim_asct_dn.do`, `sim_asct_l1.do` — design matrix swaps the four `CMc` dummies for the
  four flags (coefficient-vector column counts unchanged: 16 and 21, so `extreme_value.do`'s SCT column
  is unchanged; the OS `_cons` crank is column **62** — see `core/tests/extreme_value.do`).
- `core/process_data.do` — carries `CM_CRD/PLM/DBT` into the simulated output (enables
  comorbidity-stratified validation).

No cohort-build changes were needed: `patients/population_1995_2040_*.dta` and the OOS cohort both
already carry all four flags, and `core/load_patients.do` loads all columns.

## Imputation fix — derive `CM_CKD` from imputed eGFR (do NOT impute the flag)

Adding `CM_CKD` to the OS `mi estimate: streg` first threw **r(459) "estimation sample varies between
m=1 and m=2"**: `CM_CKD` was passively derived from eGFR with no completeness floor, so one patient
(ID_BS 7356) had imputation-varying missingness. The fix mirrors how RISS is handled:
- `prep/data_extraction.do` creates observed `CM_CKD` from observed eGFR (`eGFR<=59` → 1, else 0) and
  adds it to the three global variable lists.
- `prep/multiple_imputation.do` registers **eGFR as *regular* (not imputed)**; drops `CM_CKD` up front;
  **imputes continuous eGFR** (regress) inside the chained model; then **derives `CM_CKD` from the
  imputed eGFR with a floor** (`CM_CKD=1 if eGFR<=59`, `=0 if eGFR>59`, `=0` for any residual missing on
  the estimation sample). `CM_CKD` is **kept OUT of the chained imputation regression** — including it
  caused *"logit failed to converge on observed data"* (quasi-separation: `CM_CKD ≈ f(SerumCreatinine)`
  via eGFR). Impute the continuous driver, derive the flag — never impute both.

## Validation outcome (OOS 70/30)

**Result: whole-population OS is ~neutral, but the flags calibrate per-stratum — which is the point
(transportability).** Adding the four flags redistributes OS risk without materially moving the
whole-population curve (the flags net out at the population mean), and the comorbidity-stratified OOS
target (`os_wholepop_cm.csv`) confirms the per-stratum calibration the composite `CMc` could not give.
The bootstrap OOS block (`bootstrap_validation.do`, guarded on all four flags) confirmed it. Full
numbers and the wider OOS state are in `os_line_specific_experiment.md` §7.

## Outstanding

- **Remove `CMc`.** It is now unused by any equation but is still *derived* in the pipeline
  (`prep/multiple_imputation.do`, `prep/population_1995_2040.do`) and carried through the engine. Strip
  it out in a later pass (MI derivation + carryforward + reshape lists + `mata_setup` `vCM*` +
  `process_data`). Deferred so this pass stays reviewable.
- **May extend comorbidities to other risk equations.** This pass covers OS + the two ASCT logits only.
  The same four flags could plausibly inform TXD/TFI/TXR or maintenance — to be assessed, not assumed.
- **Validation done (adopted).** The OOS harness was run: whole-population OS did not degrade, and the
  new **comorbidity-stratified OS target** (`os_wholepop_cm.csv`) confirmed the per-stratum calibration
  the flags were added for. See "Validation outcome" above and `os_line_specific_experiment.md` §7.
  (Still to firm up: re-read FMI/DF at **10 imputations** — earlier reads were at 2.)
- **`CM_LVR`/`CM_PNR`/`CM_MLG` not included.** LVR/PNR are rare and not in the registry (estimation) MI;
  MLG (other malignancy) is available but weaker and was dropped from the ASCT model by moving off
  `CMc`. Revisit only if a specific equation motivates them.
