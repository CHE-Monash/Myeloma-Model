# Overall survival — the modelling problem and how we solved it (per-line OS)

**Purpose.** Decision record for the model's overall-survival (OS) component: the fit problem, what we
ruled out and **why** (so the dead ends aren't re-tried), the per-line rebuild that was **adopted**,
and the validation. Supersedes the earlier "line-specific experiment (not adopted)" framing — per-line
OS is now the production model on `main`. Companion: `docs/comorbidities.md`.

---

## 1. The problem, and how it evolved

The original OS model was a **single Weibull clocked from diagnosis** (`streg … b0.OS#b5.BCR`), resampled
at each pathway stage *conditional on survival to the accumulated time* (`mTSD`) so the resamples
telescope to the marginal. Two successive symptoms, one root cause:

1. **Die-too-fast (original).** Simulated OS ran *below* the registry KM, worse for poorer responders and
   at longer horizons (e.g. L3 CR 5-yr ~33% sim vs 53% registry). Concentrated in weak responders (MR/SD),
   growing with horizon.
2. **Over-predict weak responders (after the Jul-2026 calibration pass).** Once comorbidities, log-normal
   TFI, the TXD fix and conditional pathway targets landed (§4), the whole-population OS *matched* but the
   per-BCR breakdown **flipped**: weak responders were now over-predicted (L1 No-ASCT VGPR/PR/SD/PD
   +11…+22 pp), masked in the aggregate.

Both signs are the **same mechanism**: the engine's stage-conditioning on accumulated time. A heavy TFI →
large `mTSD` at later stages → OS resampled longer, telescoping back onto the line-start curves. This is
what the per-line rebuild (§5) removes.

## 2. Data-quality caveat — the registry tail is partly artefactual

Some "under-prediction" is **not** model error. Early-era patients with indolent smouldering myeloma
(SMM) misdiagnosed as active MM survive far longer, inflating the registry tail and some weak-responder
survival. Fingerprint: registry **MR/SD** survival is anomalously high — e.g. L1 MR 5-yr (49.8) *above*
PR (44.4), clinically backwards, exactly where indolent cases collect. **If the tail is an artefact, a
model that under-predicts it is closer to true-MM survival, and fitting it would fit the artefact.** So a
heavier-tailed family is not an obvious fix, and the residual MR/SD misses (§6) are most plausibly this.

## 3. What we ruled out (and why — don't re-try)

- **Heavier-tailed OS family** — risks fitting the §2 SMM tail. Not pursued.
- **Gamma-frailty Weibull** (`scratch/os_frailty_check.do`) — θ significant in only ~7/20 cells; its
  effect is at 10 yr where it grossly over-predicts the (SMM-inflated) tail. Leverage is per-BCR *shape*,
  not frailty. Not adopted.
- **Per-BCR Weibull shape** (`ancillary(i.BCR)`) — the tempting one, ruled out **twice**:
  - *In-sample it wins* (AIC/BIC beat shared shape by 212/177; per-BCR 5-yr MAD 2.4→1.8) by giving weak
    responders a heavier tail — BUT it **degrades the out-of-sample whole-population OS** (the aggregate
    that LY/QALY projections depend on) by ~3–5 pp at every horizon, pushing all horizons outside their
    bootstrap PI and dropping OS coverage 63.8%→53.2%. An OOS decomposition (shared / +prev_dur /
    +ancillary / both) pinned the loss entirely on the **ancillary** (good-responder majority gets p>1 →
    lighter tail → lower aggregate survival); `prev_dur` was benign. Classic over-fit.
  - *Re-tried on the per-line engine* (branch history, 2026-07-06) to lift the residual L1 CR miss — **a
    wash** (OOS unchanged 143/172, CR unmoved, CR<VGPR inversion persisting). Confirmed by the
    FAMILY/TRANSPORT/ENGINE decomposition (§5): the CR miss is transportability + plateau fragmentation,
    not shape. **Retired for good.** If CR ever proves miscalibrated *in-sample*, the lever is an AFT
    (log-normal) family on the L1 stages, not `ancillary`.
- **`prev_dur` frailty** (prior-line duration as a robustness proxy) — OOS-neutral, and after the TFI/TXD
  tail fixes it became redundant/counter-productive (removing it barely moved OS). Removed.
- **`prev_BCR`, `prev_prev_dur`** (trajectory proxies) — both washed out against the per-BCR benchmarks
  (redistribute hazard *within* a cell, average out). Not adopted.
- **Time clock as a covariate** (`os_clock_check.do`) — line-start vs from-diagnosis moves cell survival
  ≤3 pp and in the wrong direction for the gap; the clock is not the residual.

## 4. What we adopted — the calibration fixes (Jul 2026)

- **Comorbidities** — four flags (`CM_CKD/CRD/PLM/DBT`) added to OS and the two ASCT logits (replacing the
  `CMc` composite). Whole-pop-neutral but calibrates per-stratum (transportability). Detail + the CM_CKD
  imputation fix in `docs/comorbidities.md`.
- **TFI → log-normal** (`$dTFI`) — the thin Weibull TFI tail progressed everyone too soon (over-stated
  reach). Log-normal / log-logistic both hug the KM tail (`scratch/tfi_family_check.do`); kept log-normal
  (finite variance). Engine gained `llogistic`/`lnormal` branches in `calcSurvTime`/`calcSurvProb`
  (verified vs `streg predict`). This fixed the per-transition reach.
- **TXD** — dropped the L1 `Duration<730` curtailment (right-censor instead). *Caveat:* the observed TXD
  target is inflated by **missing treatment-END dates** (patients read as "still on treatment"), so the
  target over-states on-treatment fractions and the residual TXD misses are partly not model error. TFI is
  clocked from the end date, so those patients drop out → TFI validates cleanly, TXD doesn't. `$dTXD` = Weibull.
- **Pathway targets → conditional per-transition** `P(reach L | reached L−1)` (AJ CIF, origin = prior-line
  reach date), so line-to-line errors no longer compound.
- **Per-line OS** — the structural fix; see §5.

## 5. The diagnosis, and the per-line OS rebuild

**Diagnosis (`scratch/os_weakresp_check.do`).** The per-BCR over-prediction decomposes additively into
`SIM − KM = FAMILY + TRANSPORT + ENGINE` (reconciles to ±0.2 pp):
- **FAMILY ≈ 0** — a shared-shape Weibull refit *from line start* already fits the weak-responder KM, and
  freeing the per-BCR shape doesn't improve it. **The shape is not the constraint.**
- **ENGINE dominant (~+13 pp at L1)** — the engine's OS sits well above the line-start fit it samples: the
  `mTSD`/heavy-TFI telescoping. Where the fit already matches KM, that lift becomes pure over-prediction.
- **TRANSPORT secondary** — the rare MR/PD cells, small-N 70→30 sampling noise (partly §2 SMM).

So the fix is **how the engine samples OS**, not the equation shape.

**The rebuild (adopted, on `main`).** OS is now **one Weibull per pathway stage, each clocked from that
stage's own entry event** — 13 models: `OS_DN`, `OS_L1`, `OS_L1_NoASCT`, `OS_L1_ASCT`, `OS_L2`/`_End`,
`OS_L3`/`_End`, `OS_L4`/`_End`, `OS_L5`/`_End`, `OS_L6plus`. Each is `origin()`/`exit()` window-censored to
its stage and carries the four comorbidity flags. The draw is **unconditional at a line's start**
(`vElapsed = 0`), so there is no accumulated-time lift; the result is stored back on the diagnosis clock
for `sim_mort`. (`L6plus` is a single conditional model for the sparse deep tail.)
- Fit: `prep/risk_equations.do` (OS block, saved via `save_coefs` → `$Coeffs` → matsave).
- Engine: `core/outcomes/sim_os.do` (per-stage firing map; `core/tests/extreme_value.do` OS crank).
- Downstream (`process_data`, `sim_mort`, validators) unchanged — OS still lands on the diagnosis clock.

## 6. Validation — the result

| | pass rate | whole-pop OS 3/5/10 yr (obs 74.7/58.4/30.6) |
|---|---|---|
| OOS 70/30 (held-out) | **143/172 (83.1%)** — up from 132/172 | 74.6 / 59.0 / 30.1 |
| In-sample (full data, base_model) | **163/172 (94.8%)** | 68.7 / 51.6 / 25.8 (cohort mix — see below) |

The weak-responder over-prediction is resolved (L1 No-ASCT PR/MR/SD from +15…+17 pp to within tolerance).
**Every remaining fail is a known data artefact, not a per-line-OS defect:**

- **CR under-prediction was an OOS sampling artefact.** In OOS, L1 No-ASCT CR under-predicted (−8.0/−11.8)
  with a CR<VGPR inversion. **In-sample it vanishes** — CR is calibrated at every line (L1 −0.4/+6.4) and
  CR>VGPR ordering is correct. It was the 70→30 split (train CR ran below the held-out test CR). Nothing to
  fix.
- **MR/SD tail misses = §2 SMM artefact.** The only in-sample OS fails are L1 MR and L2/L3 SD at 5 yr —
  the SMM-inflated categories; under-predicting them is arguably closer to true MM.
- **base_model whole-pop/ASCT run low = cohort mix.** `base_model` simulates the `population_1995_2040`
  incidence cohort, which is **older** than the registry (median age at diagnosis 72.9 vs 68.4; 44% aged
  75+ vs 27.5%). ASCT is ~0% above 75, so the incidence cohort transplants less (30.9% vs registry 44.7%);
  the *same* ASCT model on real registry patients (OOS) gives 39.7%, and reweighting the population's
  age-specific ASCT rates to the registry age mix reproduces 40.1%. The OS shortfall (~−6 to −8 pp) is the
  same effect; the OOS (real registry patients) matched whole-pop OS to within a point.
- **TXD fails = missing-end-date target corruption** (§4), pre-existing.

## 7. Open caveats and future directions

1. **Risk-equation transportability (registry → incidence population).** The equations are fit on the
   *younger* MRDR registry but applied to the *older* incidence population in `base_model`. The ASCT gap
   is the visible symptom; the age extrapolation of other equations (OS, TXD/TFI) is worth checking, and
   argues for age-standardising or reweighting when the population cohort drives a headline output.
2. **Clean the SMM benchmark.** Flag/exclude likely-SMM early-era long survivors and re-benchmark; then
   the residual weak-responder (MR/SD) tail is trustworthy and tells us whether real model error remains.
3. **CR plateau fragmentation (minor).** Per-line splits a CR patient's multi-year plateau across stages;
   in-sample this is immaterial (CR calibrated). If it ever matters, an AFT (log-normal) family on the L1
   OS stages is the lever — not `ancillary` (§3).

## 8. Technical reference

- **Production OS (current):** per-line, 13 stage models, each `streg Age Age2 Male i.ECOGcc i.RISS
  CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L<n>, d(weibull)`, window-censored (`origin()`/`exit()`), clocked from
  line start, stored on the diagnosis clock. No ancillary, no `prev_dur`. Fit in `prep/risk_equations.do`;
  sampled in `core/outcomes/sim_os.do`. Age enters as **age at the stage** (engine-updated), matching the fit.
- **Benchmarks:** `scratch/benchmarks/*.csv` (registry KM/horizon/CIF by line × BCR), from
  `prep/generate_benchmarks.do` (no args → in-sample; OOS targets in `analyses/oos/targets/`). Regenerate
  after any structural change — coefficient files and benchmarks both go stale.
- **Diagnostics (`scratch/`, git-ignored):** `os_weakresp_check.do` (family/transport/engine decomposition),
  `os_frailty_check.do`, `os_ancillary_check.do`, `os_clock_check.do`, `tfi_family_check.do`, `prevdur_check.do`.
