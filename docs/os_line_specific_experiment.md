# Overall-survival fit problem — briefing & the line-specific experiment (not adopted)

**Purpose.** A self-contained brief on the open problem with the model's overall-survival (OS) fit
to the MRDR registry, the experiment tried to fix it (per-line OS equations, on branch
`os-line-specific`, not merged), and directions worth exploring. Point a future session here to
understand the issue and brainstorm.

---

## 1. The problem in one paragraph

Simulated OS runs **below** the registry Kaplan–Meier — the model makes patients **die too fast**.
The shortfall is small for good responders and **grows for poorer responders and at longer
horizons**. This holds for **both** the production single-equation OS model **and** the
line-specific model tried on the branch; the per-line model is in fact modestly **worse**. So the
gap is not fixed by decomposing OS per line — it is a distributional / data-quality issue.

Direction of the error (important): the models predict **shorter** OS, not longer. "Tail" below =
the long-horizon part of the survival curve (5-yr+ survivors); the model's curve falls off too
fast, so it under-counts long survivors.

## 2. Validation output — both models vs registry

Simulated OS **from line start**, % surviving, by line × best clinical response (BCR:
1=CR 2=VGPR 3=PR 4=MR 5=SD 6=PD). `bench` = registry KM (target); `old` = single-equation model
(production); `new` = per-line model (branch). `o−b` / `n−b` = model minus benchmark (negative =
below registry = too much death).

**3-year survival**
```
              bench  old  new | o-b  n-b        bench  old  new | o-b  n-b
L1(noASCT)                          L2
  CR    77.8   81   73 |  +3   -5     CR   77.0   75   71 |  -2   -6
  VGPR  67.9   75   65 |  +7   -3     VGPR 68.2   66   59 |  -2   -9
  PR    66.9   66   53 |  -1  -14     PR   59.3   54   47 |  -5  -12
  MR    65.7   53   44 | -12  -22     MR   53.3   43   35 | -10  -18
  SD    45.7   49   37 |  +4   -8     SD   59.2   44   36 | -15  -23
  PD    39.0   43   35 |  +4   -4     PD   25.0   30   19 |  +5   -6
ASCT (BCR_SCT 1-4)                  L3
  CR    90.1   93   92 |  +2   +2     CR   69.5   57   51 | -12  -19
  VGPR  89.4   90   89 |   0   -1     VGPR 55.6   47   43 |  -8  -13
  PR    89.2   88   87 |  -1   -3     PR   46.0   40   34 |  -6  -12
  MR    82.2   84   79 |  +1   -3     MR   39.8   30   24 |  -9  -16
                                      SD   41.8   29   22 | -12  -19
                                      PD   18.8   16   10 |  -3   -9
```

**5-year survival**
```
              bench  old  new | o-b  n-b        bench  old  new | o-b  n-b
L1(noASCT)                          L2
  CR    55.8   64   56 |  +8   -0     CR   65.7   54   50 | -11  -16
  VGPR  48.7   52   44 |  +3   -5     VGPR 50.4   43   38 |  -7  -13
  PR    44.4   39   31 |  -5  -14     PR   43.7   31   26 | -13  -18
  MR    49.8   30   24 | -20  -25     MR   30.4   21   17 | -10  -13
  SD    31.2   26   19 |  -5  -12     SD   42.5   22   18 | -21  -25
  PD    33.6   23   18 | -11  -16     PD   14.7   12    8 |  -2   -7
ASCT                                L3
  CR    79.7   85   81 |  +5   +1     CR   53.3   33   28 | -20  -25
  VGPR  77.7   78   74 |   0   -4     VGPR 33.4   22   20 | -11  -13
  PR    76.8   74   69 |  -2   -7     PR   31.5   18   15 | -13  -16
  MR    72.1   66   59 |  -6  -13     MR   20.6   13   10 |  -8  -11
                                      SD   28.0   12    9 | -16  -19
                                      PD   12.9    6    4 |  -7   -9
```

Reading it:
- **Both models sit below the registry**, worse as response worsens (CR ≈ ok → MR/SD large gaps)
  and worse at 5-yr than 3-yr. ASCT is the closest (both near the registry).
- **`new` is ~2–8 pp further below than `old`** in nearly every cell (L2, L3: every cell). So the
  per-line rewrite did not help — it is a net regression.
- Provenance: run `scratch/os_sim_vs_fit.do` (bench vs sim) and `scratch/os_stage_decomp.do`
  (per-stage self-consistency) against a fresh simulate; `old` reproduced by checking out the
  single-equation `sim_os.do` + coefficients and re-simulating.

## 3. Data-quality caveat — the registry tail may be partly artefactual

Part of the "under-prediction" is likely **not** model error: early-era patients with indolent
smouldering myeloma (SMM) misdiagnosed as active MM and treated would survive far longer, inflating
the **registry tail** (and some weak-responder survival). Fingerprint in the table: registry L2
stable-disease (SD) 3-yr survival (~59%) is anomalously tied with PR and above MR — clinically
backwards, exactly where indolent cases would collect. If the heavy tail is an artefact, a model
that under-predicts it may be closer to true-MM survival, and **fitting the tail would fit the
artefact.** So a heavier-tailed distribution family is *not* an obvious fix.

## 4. What was tried, and what is ruled out

- **Per-line / stage-specific OS** (branch `os-line-specific`): separate parametric model per
  pathway point (DN, L1 start, L1 end split by ASCT, L2–L5 start/end, L6+), each `exit()`-censored
  at the next event. Window-censoring roughly halved an initial per-line penalty, but the composed
  result is still worse than the single equation (composition loses survival a single coherent
  conditioned model preserves). **Verdict: not adopted.**
- **Engine is not the problem.** A per-stage decomposition (model evaluated with age-at-stage, as
  the engine uses) gives `sim ≈ model` at every stage — the engine faithfully samples whatever OS
  model it is given; there is no sampling bug. (An earlier apparent "end-stage leak" was a
  diagnostic error: the check used age-at-diagnosis while the engine uses age-at-stage.)
- **Heavier-tailed family**: deferred — see §3; would risk fitting the SMM artefact.
- **Gamma-frailty Weibull** (single-equation, unshared frailty; `scratch/os_frailty_check.do`,
  checked against the registry benchmarks by line × BCR). **Verdict: not adopted.** Two findings:
  (a) *The family is adequate at the horizons that matter.* An intercept-only Weibull fit **within**
  each BCR group tracks KM at 3–5 yr with gaps mostly ±0–3 pp — e.g. L1 MR at 5 yr fits to **−1**
  vs the production engine's **−20** on that cell. So the engine's shortfall is **structural** — it
  shares **one Weibull shape across all BCR and lines** and a monotone BCR scale that cannot
  reproduce the registry's clinically-backwards MR>PR ordering — not a tail-family problem.
  (b) *Frailty mostly fattens the tail, and over-fits it.* θ is significant (LR p<0.05) in only
  ~7/20 cells; its main effect is at 10 yr, where in the one well-populated cell with a trustworthy
  10-yr benchmark (L3 PR: bench 4.7) it predicts **15.2** — a gross over-prediction, i.e. it fits
  the §3 SMM-artefact tail. It does close a couple of poor-responder 5-yr cells (L2/L3 PD, −5→0) but
  net trades under- for over-prediction. The leverage is per-BCR / per-line **shape**, not frailty.
- **Shape-relaxed Weibull** (`scratch/os_ancillary_check.do`): the follow-up to (b) — free the
  Weibull *shape* on a stacked line-start model with scale already free per line×BCR cell, so only
  the shape varies across specs. Results (AIC / BIC / mean|gap| at 5 yr, 22 cells):
    - M0 one shared shape — 18736 / 18950 / 2.38
    - **MB `ancillary(i.BCR)` — 18524 / 18773 / 1.81**  ← wins
    - MO `ancillary(i.OS)`   — 18619 / 18855 / 2.25

  **`ancillary(i.BCR)` wins decisively on AIC *and* BIC (Δ vs M0 = 212 / 177) and gives the best
  5-yr accuracy**, fixing the weak-responder cells that died too fast (L2/L3 SD & PD 5-yr gaps
  −3…−7 → ≈0). **`ancillary(i.OS)` (shape by line) improves likelihood but not benchmark accuracy —
  not adopted.** Note: even shared-shape M0 sits within ~2 pp once scale is free per cell (vs the
  engine's −20 on L1 MR), so per-BCR shape is a real but *secondary* refinement; the larger residual
  is the pooled OS-indicator structure, **not** the time clock (measured innocuous — next bullet).
  **Candidate fix: a per-BCR Weibull shape in the OS equation** — cheap to fit (`ancillary(i.BCR)`)
  and to simulate (shape vector indexed by BCR, same inversion in `core/outcomes/sim_os.do`).
- **Clock: line-start vs from-diagnosis** (`scratch/os_clock_check.do`). A controlled contrast —
  same patients/response/covariates, only the time origin differs — measuring what an earlier
  *inference* had attributed the residual to. **Result: the clock is nearly innocuous.** L1 (TSD≈0)
  is a passing sanity check (DIAG≈LS); at L2/L3 the from-diagnosis clock moves 3–5 yr survival only
  ~0–3 pp and, where it does, *lifts* it (slightly heavier tail) — the wrong sign for the die-too-fast
  gap. So the diagnosis clock the engine uses is **not** the residual. Time-since-diagnosis-at-line-entry
  (TSD) as a covariate is significant (L3 p<0.001, ΔAIC −12; negative coef = slower-to-reach-line →
  better OS, the indolent signal) but shifts the cell-level benchmarks ≤1 pp — it captures individual
  heterogeneity, not aggregate accuracy. Useful later for lifetime/composition coherence, not for the
  KM gaps. **Net: the leverage is per-BCR shape; clock and TSD are not where the accuracy is.**

### 4a. Baseline for the `ancillary(i.BCR)` test — production shared-shape OS, re-simulated 2026-07-02

The **actual simulated** OS vs registry KM from the current (shared-shape) engine — the number to
beat once `ancillary(i.BCR)` is wired in. Unlike the fit-only diagnostics above, this is the full
telescoped simulation (`scratch/os_sim_baseline.do`, reads `analyses/base_model/simulated/` +
`scratch/benchmarks/`; no drive needed). `s-b` = Sim − Bench (pp); negative = simulated OS below
registry (dies too fast).

```
                 3yr                5yr                10yr
  L1 (No ASCT)   bench sim  s-b     bench sim  s-b     bench sim  s-b
    1 CR          77.8 81.1 +3       55.8 64.3 +9       22.3 32.4 +10
    2 VGPR        67.9 75.1 +7       48.7 53.0 +4       17.0 20.0 +3
    3 PR          66.9 66.3 -1       44.4 41.3 -3       17.8 11.4 -6
    4 MR          65.7 54.2 -12      49.8 31.9 -18      23.8  8.5 -15
    5 SD          45.7 50.1 +4       31.2 27.7 -4       16.0  6.1 -10
    6 PD          39.0 44.2 +5       33.6 25.6 -8          .  6.2  .
  L1 ASCT
    1 CR          90.1 92.5 +2       79.7 84.7 +5       58.7 58.3 +0
    2 VGPR        89.4 89.5 +0       77.7 77.9 +0       39.9 47.6 +8
    3 PR          89.2 88.2 -1       76.8 75.2 -2       41.1 41.8 +1
    4 MR          82.2 83.3 +1       72.1 66.1 -6       38.8 30.0 -9
  L2
    1 CR          77.0 76.5 -1       65.7 57.3 -8       41.1 24.6 -17
    2 VGPR        68.2 67.6 -1       50.4 46.4 -4       21.0 16.2 -5
    3 PR          59.3 56.8 -3       43.7 34.0 -10      18.5  8.9 -10
    4 MR          53.3 46.0 -7       30.4 24.7 -6          .  4.9  .
    5 SD          59.2 45.7 -13      42.5 24.2 -18      17.7  5.5 -12
    6 PD          25.0 32.3 +7       14.7 15.1 +0          .  2.4  .
  L3
    1 CR          69.5 64.9 -5       53.3 41.7 -12         . 10.6  .
    2 VGPR        55.6 53.3 -2       33.4 28.2 -5       10.6  5.5 -5
    3 PR          46.0 44.8 -1       31.5 22.6 -9        4.7  4.1 -1
    4 MR          39.8 34.3 -5       20.6 15.8 -5          .  2.4  .
    5 SD          41.8 33.4 -8       28.0 15.0 -13         .  2.5  .
    6 PD          18.8 18.7 +0       12.9  7.7 -5          .  1.3  .
```

Headline: mean |s-b| ≈ **4.0 pp at 3 yr, 7.0 pp at 5 yr** (22 cells). Die-too-fast is concentrated
in the weak responders (MR/SD) and grows with horizon — the classic signature. **Crucially, this
simulated 5-yr MAD (~7 pp) is far larger than the line-start *fit* MAD (~2 pp; §4 shape-relaxed
work), confirming the big residual is the telescoped simulation (progression through rising-hazard
`OS` cells), not the per-line fit.** So `ancillary(i.BCR)` — which improves the line-start fit — is
expected to help but may only *partially* close this simulated gap; whatever remains is the
telescoping/collapsing (or the §3 tail artefact), to be read off by re-running this table.

### 4b. After: `ancillary(i.BCR)` wired in, re-simulated (branch `os-ancillary`)

Per-BCR Weibull shapes fitted (p): CR 1.47, VGPR 1.38, PR 1.28, SD 1.05, MR 0.86, PD 0.40 — good
responders get p>1 (lighter tail), weak responders p<1 (heavier tail), as intended. Re-simulated
`s-b`:

```
                 3yr                5yr                10yr
  L1 (No ASCT)   bench sim  s-b     bench sim  s-b     bench sim  s-b
    1 CR          77.8 81.3 +3       55.8 63.2 +7       22.3 28.2 +6
    2 VGPR        67.9 75.0 +7       48.7 52.3 +4       17.0 17.3 +0
    3 PR          66.9 66.9 +0       44.4 41.9 -3       17.8 10.9 -7
    4 MR          65.7 55.7 -10      49.8 34.1 -16      23.8  8.7 -15
    5 SD          45.7 51.1 +5       31.2 28.7 -3       16.0  6.3 -10
    6 PD          39.0 52.9 +14      33.6 31.5 -2          .  7.3  .
  L1 ASCT
    1 CR          90.1 93.3 +3       79.7 84.6 +5       58.7 55.2 -4
    2 VGPR        89.4 90.2 +1       77.7 78.4 +1       39.9 45.7 +6
    3 PR          89.2 88.5 -1       76.8 75.4 -1       41.1 40.3 -1
    4 MR          82.2 83.1 +1       72.1 66.8 -5       38.8 31.2 -8
  L2
    1 CR          77.0 75.1 -2       65.7 55.0 -11      41.1 21.6 -19
    2 VGPR        68.2 67.8 +0       50.4 46.0 -4       21.0 14.9 -6
    3 PR          59.3 58.1 -1       43.7 34.9 -9       18.5  8.6 -10
    4 MR          53.3 48.4 -5       30.4 26.7 -4          .  5.4  .
    5 SD          59.2 47.2 -12      42.5 25.7 -17      17.7  5.7 -12
    6 PD          25.0 34.8 +10      14.7 16.8 +2          .  2.8  .
  L3
    1 CR          69.5 64.5 -5       53.3 40.7 -13         .  9.5  .
    2 VGPR        55.6 54.1 -1       33.4 29.0 -4       10.6  5.5 -5
    3 PR          46.0 45.8 +0       31.5 23.3 -8        4.7  4.2 +0
    4 MR          39.8 36.2 -4       20.6 17.2 -3          .  2.7  .
    5 SD          41.8 36.1 -6       28.0 16.5 -11         .  2.6  .
    6 PD          18.8 22.3 +3       12.9  9.4 -4          .  1.5  .
```

**Mean |s-b|: 3 yr 4.05 → 4.27 (slightly worse); 5 yr 7.00 → 6.23 (~11% better).** Verdict: a
**modest, mixed** improvement, exactly as §4a predicted.
- *Intended direction confirmed:* every MR/SD cell improved 1–2 pp at 3 and 5 yr — but the gaps stay
  large (MR/SD still −11…−17 at 5 yr).
- *Side effects:* PD's extreme shape (p=0.40) now **over**-predicts at 3 yr (L1 PD +5→+14); good-responder
  CR at L2/L3 worsens slightly at 5 yr (p>1 → lighter tail → lower long survival).
- *Core problem untouched:* the simulated MAD (~6 pp) stays ~3× the line-start fit MAD (~2 pp), so the
  dominant residual is the **telescoping**, not the fit. Per-BCR shape shaves only the fit-sized part.

Bottom line: `ancillary(i.BCR)` is a small, defensible refinement in the right direction for weak
responders, but not the fix for die-too-fast. The remaining lever is structural (the telescoped
progression) and/or cleaning the §3 SMM-artefact benchmark. If adopted, consider taming the PD
shape (its p=0.40 over-lifts PD at 3 yr).

### 4d. Tried and NOT adopted: previous-line BCR (`i.prev_BCR`)

On top of §4c, added the previous line's best clinical response as a disease-biology / frailty proxy
(`i.b5.prev_BCR` main effect; passive, since BCR is imputed and BCR_L*/BCR_SCT are passive). Fitted
log-hazard effects (base SD): CR **+0.50**, VGPR +0.38, PR +0.21, MR −0.14, PD +0.02 — i.e. a *good*
prior response predicts a *higher* subsequent hazard. That sign is not wrong: current BCR is already
in the model, so conditional on current response, a good prior response means the patient *declined*
to get here (CR→PR = aggressive relapse) while a poor prior response means they *improved*. So
`prev_BCR` captures response **trajectory (declining vs improving)**, a real but different signal.

Validation (mean |s-b|): 3 yr 3.91 → **4.18** (worse), 5 yr 5.86 → **5.82** (flat, noise). **It does
not help.** Cause: `prev_BCR` conditional on current BCR only *redistributes* hazard **within** each
(line, current-BCR) cell; since the sim's within-cell prior-response mix ≈ the registry's, it
**averages out** against the current-BCR benchmarks (leaving a small Jensen residual that slightly
hurt the good-responder cells §4c had fixed). `prev_dur` escaped this because it is continuous and
shifts the tail strongly; a categorical trajectory term conditional on current response does not.
A line-varying version (`i.OS#i.prev_BCR`) would face the same wash-out per line, so it is not
pursued. **Verdict: not adopted; §4c (ancillary + prev_dur) remains the best model.** The stubborn
MR/SD residual is now most plausibly the §3 SMM-artefact benchmark, not a missing OS covariate.

*Fit-level head-to-head* (`scratch/os_predictor_compare.do`, AIC on the diagnosis-clock estimation
data): base 8868; **+prev_dur −52 (p<0.001)**; +prev_BCR **+1 (p=0.11, NOT a predictor alone)**;
+both −73. So `prev_dur` and `prev_BCR` are **not substitutes** — `prev_dur` is the dominant
standalone predictor and `prev_BCR` carries no OS signal on its own (a suppression: alone, its
protective-duration and harmful-decline components cancel; it only turns significant *on top of*
`prev_dur`, as the residual decline-direction). `both` has the best AIC, but that `prev_BCR`
increment does not help the sim validation (above). Confirms: **`prev_dur` alone is the right single
covariate; duration is the fundamental OS signal, response-alone is not.**

### 4e. Tried and NOT adopted: second duration lag (`prev_prev_dur`)

Tested whether *more* of the trajectory helps by adding a second lag — the duration of the line two
back (`(TXD_{L-2}+TFI_{L-2})/30.4375`; engine `mTSD[.,2·seg-4] - mTSD[.,2·seg-6]`, seg≥4). This is the
identifiable way to add trajectory beyond one lag, avoiding the clock-collinearity of raw
time-to-line (which on the diagnosis clock **is** `_t0`). Fitted coefficient **−0.0049**, ~3.5× weaker
than `prev_dur` (−0.017) — recency dominates, as expected.

Validation (mean |s-b|): 3 yr 3.91 → **4.00** (worse), 5 yr 5.86 → **5.73** (better) — a wash, both
within noise. Cell-level it's reshuffling, not real gain: the second lag helped the L3 good responders
it can reach (L3 CR 5-yr −13→−9, VGPR −5→−2), but the joint refit drifted L1 MR worse (−16→−19; L1 has
no `prev_prev_dur`). **Verdict: not adopted — one lag (`prev_dur`) was enough; the most-recent line's
pace is where the trajectory signal lives.** §4c remains the committed best.

### 4f. Out-of-sample whole-population OS (pre-§4c baseline)

§4a–4e are in-sample (registry KM = the fit data). The **OOS 70/30 harness** (`analyses/oos`) trains
on 70%, predicts the held-out 30%. It previously had only per-line×BCR OS targets, so a
**whole-population OS** target (from diagnosis, all patients, 3/5/10 yr) was added
(`prep/generate_benchmarks.do` → `os_wholepop.csv`; `analyses/oos/bootstrap_validation.do` →
`OS/ALL/{3,5,10}yr`; point-estimate block in `validate_outcomes.do`).

Pre-§4c OOS model (whatever was last trained there — **not** ancillary/prev_dur), held-out 30% (N=1888):

| horizon | observed | predicted (median) | bootstrap 95% PI | inside? |
|---|---|---|---|---|
| 3 yr  | 74.7 | 72.6 | [70.9, 74.2] | just outside (+0.6) |
| 5 yr  | 58.4 | 57.3 | [55.3, 59.0] | inside |
| 10 yr | 30.6 | 28.3 | [26.2, 30.2] | just outside (+0.4) |

Point estimate is within **1–2 pp** at every horizon (all PASS at the 10% tolerance). But observed
sits *fractionally above* the (very tight, ~3–4 pp wide) aggregate 95% PI at 3 and 10 yr — same
direction as 5 yr — i.e. a small **systematic ~1.5–2 pp under-prediction**: the die-too-fast,
quantified at the whole-population level. **So the whole population predicts well in aggregate; the
wobble is the BCR breakdown** (16% of point tests fail; overall PI coverage 59.5%, OS 63.8%), errors
both ways that cancel. The tight aggregate PI also suggests the bootstrap PIs are a touch narrow
(parameter uncertainty only).

**Re-training with §4c made it WORSE, not better** (point estimate matches the bootstrap median, so
the effect is real, not an artefact):

| horizon | observed | pre-§4c median [PI] | §4c median [PI] |
|---|---|---|---|
| 3 yr  | 74.7 | 72.6 [70.9, 74.2] | **68.8 [66.3, 71.1]** |
| 5 yr  | 58.4 | 57.3 [55.4, 59.0] | **53.7 [51.2, 56.0]** |
| 10 yr | 30.6 | 28.3 [26.2, 30.2] | **25.8 [23.4, 28.2]** |

§4c pulled whole-population OS ~4 pp **lower** at every horizon — under-prediction tripled (~2 pp →
~5–6 pp; observed now 3–8 pp *above* the §4c PI) — and coverage fell (overall 59.5% → 55.6%, OS
63.8% → 53.2%). So **§4c improves the in-sample per-BCR fit but degrades out-of-sample generalisation**
— an over-fitting / systematic-bias signature. **Mechanism: the ancillary per-BCR shape, not `prev_dur`.**
Whole-population OS is L1-dominated (`prev_dur = 0` at L1, so only the shape moves it); the ancillary
gives the good-responder majority p>1 (lighter tail → lower long survival), pulling the aggregate down.
The shared shape (~1.0) gave a near-perfect aggregate (73.0 vs obs 74.7); spreading the shapes to fit
each cell broke it. Since population projections (LY/QALY) depend on the aggregate, **§4c may be a net
negative for the model's primary purpose.** Decomposition underway: **prev_dur-only** (L2+ only,
expected benign for the aggregate) then **ancillary-only** (the suspected culprit).

**Decomposition — prev_dur-only (shared shape + prev_dur, no ancillary).** Re-trained + re-simulated
OOS. Whole-population OS, bootstrap median [PI]:

| horizon | observed | pre-§4c | §4c | **prev_dur-only** |
|---|---|---|---|---|
| 3 yr  | 74.7 | 72.6 [70.9,74.2] ✗ | 68.8 [66.3,71.1] ✗ | **72.5 [70.9,73.9] ✗**(+0.8) |
| 5 yr  | 58.4 | 57.3 [55.4,59.0] ✓ | 53.7 [51.2,56.0] ✗ | **57.0 [55.2,58.8] ✓** |
| 10 yr | 30.6 | 28.3 [26.2,30.2] ✗ | 25.8 [23.4,28.2] ✗ | **29.0 [27.0,31.1] ✓** |

(Deterministic point estimate 72.4/56.7/30.1, matches the median.) **Removing the ancillary restored
the aggregate** — prev_dur-only is the best of the three (5-yr & 10-yr inside the PI; 3-yr out by
0.8 pp), and OS coverage recovers to 61.7% (§4c 53.2%, pre-§4c 63.8%; overall 56.3%). So the ancillary
shape is confirmed as the OOS aggregate-killer, and **`prev_dur` alone is OOS-neutral-to-slightly-positive**
(≈ pre-§4c at 3/5 yr, slightly better at 10 yr; a defensible keep).

**Decomposition — ancillary-only (per-BCR shape `ancillary(b5.BCR)`, no prev_dur).** Re-trained on the
70% (`bOS` 1×65, `_cons` col 58, no prev_dur covariate; per-BCR shapes p = CR 1.46, VGPR 1.33, PR 1.31,
MR 0.97, SD 1.04, PD 0.35) and re-simulated on the held-out 30% (N=1888). Whole-population OS,
bootstrap median [PI] (498 resamples):

| horizon | observed | pre-§4c | §4c | prev_dur-only | **ancillary-only** |
|---|---|---|---|---|---|
| 3 yr  | 74.7 | 72.6 [70.9,74.2] ✗ | 68.8 [66.3,71.1] ✗ | 72.5 [70.9,73.9] ✗ | **69.9 [67.7,72.2] ✗** |
| 5 yr  | 58.4 | 57.3 [55.4,59.0] ✓ | 53.7 [51.2,56.0] ✗ | 57.0 [55.2,58.8] ✓ | **54.8 [52.3,57.1] ✗** |
| 10 yr | 30.6 | 28.3 [26.2,30.2] ✗ | 25.8 [23.4,28.2] ✗ | 29.0 [27.0,31.1] ✓ | **25.7 [23.2,28.2] ✗** |

(Deterministic point estimate 69.7/54.3/25.7 matches the median.) **Ancillary-only reproduces the §4c
degradation almost exactly** — 69.9/54.8/25.7 vs §4c's 68.8/53.7/25.8, within ~1 pp at every horizon and
~3–5 pp *below* both pre-§4c and prev_dur-only, with the observed value *outside* (above) the interval at
all three horizons. So the ancillary carries essentially the *entire* §4c whole-population loss on its
own; adding `prev_dur` on top (→ §4c) barely moves the aggregate (prev_dur is a no-op at L1, which
dominates the from-diagnosis curve). **This rules out an ancillary×prev_dur interaction** — the two
covariates act independently on the aggregate: ancillary pulls it down, prev_dur leaves it (at 3/5 yr)
or nudges it up (10 yr). OS-family coverage falls to **51.1%** (24/47; overall 54.0%) — the lowest of
the four, confirming the ancillary as the coverage-killer too.

### 4f-conclusion. OOS decomposition — verdict

The four OOS whole-population variants isolate each covariate's effect on the metric that drives the
model's primary output (LY/QALY projections depend on the aggregate OS curve, not the per-BCR cells).
Each variant differs only in the OS equation's shape and frailty terms:

| variant (OS-equation spec) | 3/5/10-yr median [PI] | horizons inside PI | OS-family coverage |
|---|---|---|---|
| **shared shape** (one Weibull shape for all responses, no frailty — *current production*) | 72.6 / 57.3 / 28.3 | 5-yr only | 63.8% |
| **shared shape + prior-line-duration frailty** (adds `prev_dur`: slow progressors get lower hazard) | 72.5 / 57.0 / 29.0 | 5- & 10-yr | 61.7% |
| **per-BCR shape** (each response group its own Weibull shape, no frailty) | 69.9 / 54.8 / 25.7 | none | 51.1% |
| **per-BCR shape + prior-line-duration frailty** (both; the committed §4c model) | 68.8 / 53.7 / 25.8 | none | 53.2% |

(Observed 74.7 / 58.4 / 30.6; all four rows are bootstrap medians.)

**Conclusion: the per-BCR (ancillary) Weibull shape is the OOS aggregate-killer; the prior-line-duration
frailty (`prev_dur`) is benign.** Compare the pairs: adding the per-BCR shape (rows 1→3, or 2→4) drops
whole-population OS ~3–5 pp at every horizon and pushes every horizon *outside* its PI; adding the
frailty term (rows 1→2, or 3→4) barely moves the aggregate. So the loss is the shape, not the frailty.
The mechanism (confirmed in §4f above): the whole-population OS-from-diagnosis curve is L1-dominated,
where `prev_dur = 0`, so only the shape moves it. The shared shape (~1.0) gives a near-perfect
aggregate; spreading the shapes to fit each in-sample BCR cell gives the good-responder *majority* p>1
(lighter tail → lower long-horizon survival), which drags the population curve ~4–5 pp below observed
at every horizon. This is a textbook over-fit: the ancillary buys a ~1 pp improvement in the in-sample
per-BCR 5-yr MAD (§4b/§4c) at the cost of ~3–5 pp on out-of-sample aggregate OS and a 10 pp drop in OS
PI coverage (63.8% → 53.2%).

**Recommendation:**
1. **Drop the ancillary.** Its only benefit (per-BCR in-sample shape) does not generalise and it
   materially harms the aggregate that population projections depend on. Revert `sim_os.do` /
   `risk_equations.do` / `extreme_value.do` to the shared-shape OS.
2. **Keep `prev_dur` (or at least do not fear it).** Alone it is OOS-neutral-to-slightly-positive
   (≈ pre-§4c at 3/5 yr, best-of-all at 10 yr, coverage ~unchanged) while attacking the telescoping
   die-too-fast in the later-line good-responder cells in-sample (§4c). It is a defensible keep; if
   parsimony is preferred, pre-§4c (shared shape, no prev_dur) is the safe fallback and is within
   1–2 pp of observed at every horizon.
3. The residual ~2 pp systematic under-prediction that survives in *every* variant (pre-§4c included)
   is the die-too-fast quantified at the population level; per §3/§5 it is most plausibly the SMM
   tail artefact in the benchmark, not a missing OS covariate — a benchmark-cleaning problem, not a
   model-specification one.

### 4c. After: `ancillary(i.BCR)` + `prev_dur` (progression-pace/frailty covariate)

Attacks the telescoping directly: adds the **previous line's start-to-start duration** (`prev_dur =
TXD_{L-1}+TFI_{L-1}`, months) to the OS regression as an observed-frailty proxy — a patient who
progressed slowly (long prior line) is robust and gets a lower forward hazard. Wired on branch
`os-ancillary` (fit `prep/risk_equations.do` by `Line`; engine `core/outcomes/sim_os.do` by OS-level
segment; both months). Fitted coefficient **−0.0169 per month** (negative = longer prior line lowers
hazard, as `os_clock_check` predicted; e.g. a 35-mo prior line ⇒ ~45% hazard cut).

Mean |s-b| across 22 cells, vs the earlier steps:

| model | 3-yr | 5-yr |
|---|---|---|
| §4a shared shape (production) | 4.05 | 7.00 |
| §4b + `ancillary(i.BCR)` | 4.27 | 6.23 |
| **§4c + `ancillary` + `prev_dur`** | **3.91** | **5.86** |

**`prev_dur` is the first change to beat production on BOTH horizons** — it recovers the 3-yr
regression `ancillary` introduced and pushes 5-yr further down (~16% vs production). The gains land
exactly where `ancillary` had hurt — **later-line good responders** (long prior line ⇒ robust):
L2 CR 5-yr −11→−8, L3 CR −13→−11, L2 VGPR −4→−2. The mechanism differentiates *within* each response
group by pace (fast-progressing PD gets little lift, correctly). L1 is ~unchanged (prev_dur=0 there;
small drift is the refit). **Not fixed:** weak-responder MR/SD die-too-fast (L2 SD −17→−16, L1 MR
−16→−18) — fast progressors get little frailty lift, and that residual is most likely the §3 SMM
tail artefact. Verdict: real, correctly-signed pull that validates the telescoping hypothesis;
the best combination so far, modest in absolute terms.

## 5. Open questions / directions to brainstorm

1. **Clean the OS target.** Flag/exclude likely-SMM early-era long survivors (diagnosis era +
   indolent course) and re-benchmark. Then the residual gap — especially the **3-yr weak-responder
   gap in the best-populated data** (e.g. L2/L3 MR/SD) — is trustworthy and tells us whether real
   model error remains.
2. **Restrict the analysis horizon** (e.g. 5–6 yr) if the far tail stays untrustworthy — but note
   this truncates life-years/QALYs for genuinely long-surviving MM, a modelling decision to document.
3. **Adopt a per-BCR Weibull shape.** §4 settles this: `ancillary(i.BCR)` wins on AIC/BIC and 5-yr
   accuracy; `ancillary(i.OS)` (shape by line) and gamma-frailty do not. **Wired on branch
   `os-ancillary`**: `prep/risk_equations.do` adds `ancillary(i.BCR)` to the OS `streg` and stores a
   1×6 per-BCR log-shape `bOS_p` (extracted by name from `e(b_mi)`); `core/outcomes/sim_os.do` uses
   `bOS_p[vBCR]` as the per-patient Weibull shape; `core/mata_functions.do` `calcSurvTime`/`calcSurvProb`
   accept a colvector `aux`. Pending: regenerate coefficients → re-simulate → re-run
   `scratch/os_sim_baseline.do` and compare against the §4a baseline.
4. **Is the weak-responder gap real or artefact?** Disentangling (1) from genuine under-fit is the
   crux; both point at the same cells, so a cleaned benchmark is the enabler.

## 6. Technical context (for whoever picks this up)

- **Production (single-equation) OS:** `streg Age Age2 Male i.ECOGcc i.RISS b0.OS#b5.BCR`, one
  Weibull, clocked from diagnosis, resampled each stage in `core/outcomes/sim_os.do` with
  conditioning (`vRN = U·S(TSD)`) so the resamples telescope to the marginal. `OS` is a collapsed
  line indicator (0=Pre … 7=L6+); spec in `prep/risk_equations.do` (main). Age enters as **age at
  the stage** (the engine updates it), matching the fit.
- **Per-line model:** branch `os-line-specific`, `prep/risk_equations.do` + `core/outcomes/sim_os.do`
  there; coefficients in that branch's `coefficients_base_model.mmat`.
- **Diagnostics (local, `scratch/` is gitignored):** `os_sim_vs_fit.do` (bench vs sim vs a fit
  column — the fit column is only meaningful for the single-equation model; for censored per-line
  models it extrapolates a window hazard and should be ignored) and `os_stage_decomp.do`
  (per-stage sim-vs-model self-consistency, using age-at-stage).
- Benchmarks: `scratch/benchmarks/os_*.csv` (registry KM by line × BCR), from
  `prep/generate_benchmarks.do`.

## 7. Duration-tail refit → pathway recalibration → prev_dur removal, and the weak-responder OS flip (2026-07-06)

A second OOS push, on top of §4 (comorbidities now in OS + ASCT; see `docs/comorbidities.md`). The
chain of changes and what each did:

**TXD `Duration < 730` truncation removed.** The L1 TXD `streg`s (ASCT splines + NoASCT) were fit only
on patients with duration <2yr, starving the tail → the model stopped L1 treatment far too early
(24-mo on-treatment ≈0% vs observed 8–23%). Dropping the filter (right-censoring the still-on-treatment
patients instead of excluding them) helped the 12-mo cells but not the 24-mo tail. **Key caveat (data
quality): the observed TXD target is inflated by missing treatment-END dates** — patients whose
treatment really ended but whose end date is blank appear "still on treatment", so the target itself is
unreliable and the model's shorter durations may be closer to truth. TFI is clocked *from* the end date
so those patients drop out of TFI entirely → TFI is on the clean subset (this is why TFI validates well
and TXD doesn't). Tried `$dTXD` = Gompertz (worse, shorter); reverted to **Weibull**.

**TFI → log-normal (was Weibull).** Diagnosed the over-progression (per-transition reach too high) as a
**thin Weibull TFI tail**: it matched 12/24-mo but progressed everyone too soon. `scratch/tfi_family_check.do`
(refit the TFI quantity per family, fitted-S(t) vs KM at horizons) → **log-logistic and log-normal both
hug the KM tail** (tail|err| ≈0.005 vs Weibull 0.037) while still matching 12/24-mo; Gompertz plateaus
and misses the near-term. Loglogistic and lognormal came out **empirically identical** on reach and OS
(A/B'd), so kept **log-normal** (finite variance; loglogistic γ≈0.56 → infinite variance). Engine:
added `llogistic` + `lnormal` branches to `calcSurvTime`/`calcSurvProb` in `core/mata_functions.do`
(both AFT; verified vs `streg predict` to machine precision). `$dTFI = "lognormal"`. This is the change
that **fixed the reach** (see next).

**Pathway targets made CONDITIONAL per-transition.** Old reach targets were cumulative-from-L1
Aalen–Johansen CIF (`P(reach L | started L1)`), so errors compounded line-by-line and the sim's
lifetime reach vs the observed's follow-up-truncated reach inflated the gap. Rewrote to
**`P(reach L | reached L-1)`** — observed AJ CIF with origin = the previous line's reach date, sim =
`reached L / reached L-1` — in `prep/generate_benchmarks.do`, `analyses/oos/bootstrap_validation.do`,
`analyses/oos/validate_outcomes.do`. (ASCT rate unchanged — already conditional on L1-end reachers.)

**`prev_dur` removed entirely.** After the tails changed, OS went from die-too-fast to over-predicting
the 5–10yr tail (whole-pop +3). Chased `prev_dur` (the (TXD+TFI) frailty covariate): `prevdur_check.do`
showed the *simulated* `prev_dur` runs 2–3× the observed at every percentile (median 23 vs 12, p90 138
vs 44 mo), driven by the TFI component — but the observed is **censoring-selected short** (only
completers-who-progressed are seen), so it's partly a comparison artefact. Tried a p99 cap on
`vPrevDur` (no effect), then removed `prev_dur` from OS altogether (fit + engine + `extreme_value`
column). **OS barely moved** → `prev_dur` was NOT the OS lever, and it is now redundant/counter-
productive (the lift is over-predicting, and it would lift more). **The actual OS-lift channel is the
engine's stage-conditioning**: OS is resampled conditional on survival to `mTSD` (cumulative time), so
a heavier TFI → larger `mTSD` at later stages → OS resampled longer. So the OS-tail-over and the
residual reach-over are **the same coupled effect** (heavy TFI keeps too many patients alive-and-
progressing late), not a `prev_dur` problem. `prev_dur` confirmed fully redundant; removed.

### 7a. Current model (2026-07-06)
OS: shared Weibull, `Age Age2 Male i.ECOGcc i.RISS b0.OS#b5.BCR CM_CKD CM_CRD CM_PLM CM_DBT` (NO
ancillary, NO prev_dur). TFI: **log-normal**. TXD: **Weibull**, no `<730`. ASCT logits: four comorbidity
flags (not `CMc`). Conditional pathway targets. bOS 1×63 (CKD=58,CRD=59,PLM=60,DBT=61,_cons=62,ln_p=63).
Deterministic OOS: whole-pop OS 74.4/61.6/34.0 (obs 74.7/58.4/30.6), reach L1→L2/L2→L3/L3→L4/L4→L5 =
77.9/75.7/74.8/69.3 (obs 71.2/74.4/67.9/61.1; L2→L3 passes). **132/172 tests pass (76.7%).**

### 7b. OPEN ISSUE — weak-responder OS over-predicted per-BCR (the whole-pop masks it)
The whole-pop OS looks fine (+3), but **OS-by-BCR is badly over-predicted for weak responders**,
worst at **L1 non-ASCT**: VGPR +14.8/+22.2, PR +16.4/+15.9, SD +17.2/+13.9, PD +17.5/+15.9 (sim−obs,
3yr/5yr %); only CR is calibrated. L2/L3 show the same in places. This is the **flip of the original
§4a die-too-fast weak-responder problem** — same cells, opposite sign. It is **not** a `prev_dur`
artefact (`prev_dur`=0 at L1). This is the **flip of the original §4a die-too-fast weak-responder
problem** — same cells, opposite sign.

### 7c. RESOLVED — it is the engine coupling, NOT the OS shape (`scratch/os_weakresp_check.do`, 2026-07-06)
`os_weakresp_check.do` decomposes the per-BCR over-prediction (line-start clock, 3 & 5 yr) into three
additive layers — **FAMILY** (a shared-shape Weibull refit *in-sample on the test fold* vs the test KM),
**TRANSPORT** (train-fitted vs test-fitted, both on test covariates → 70→30 generalisation), **ENGINE**
(the engine's simulated OS vs the line-start fit it samples → the from-diagnosis `mTSD`-conditioning /
heavy-TFI lift). `SIM − KM = FAMILY + TRANSPORT + ENGINE` (reconciles to ±0.2pp). L1 No-ASCT, 3yr,
weak-responder averages:

- **FAMILY ≈ 0** (VGPR −1.3, PR −2.6, SD −0.7, MR +0.1). A shared-shape line-start Weibull **already
  fits the weak-responder KM**, and the **`ancillary(i.bcr)` shape probe (MB) does not improve it**
  (−1.0/+1.9/−1.6/−1.3). **→ the shared shape is NOT the constraint; the per-BCR/per-line `ancillary`
  route is a confirmed dead end** — it was never going to buy accuracy, which is why it only damaged the
  aggregate (§4). Retire it.
- **ENGINE is the dominant driver** for common weak responders — VGPR/PR/SD ≈ **+13pp** (the engine's
  OS sits ~13–16pp above the line-start fit). The lift is roughly BCR-independent in absolute pp: it's
  the **heavy-TFI telescoping** — patients who progress through long TFIs get OS resampled upward at
  deeper stages (large `mTSD`), and that lands back on the line-start curves. Where the fit already
  matches KM (weak responders) the lift becomes pure over-prediction; where the fit under-predicts (CR,
  −11pp) the same lift compensates it back to KM.
- **TRANSPORT is a real secondary driver** for the rare/noisy cells — **MR/PD ≈ +10–13pp**, with ~0
  engine contribution. MR/PD are the smallest cells and the §3 SMM-artefact suspects (their held-out KM
  is likely inflated), so this is partly not-model-error.

**Implication for the fix.** The line-start refit (`FIT_test`) *is* a per-line-from-line-start OS model
and it tracks KM well → the fix is in **how the engine samples OS**, not the OS equation's shape:
  - **(a) per-line OS clocked from line start** — the old v4.0 rebuild, now re-motivated to *remove* the
    over-prediction; principled but a real `sim_os.do` + `risk_equations.do` rework, and whole-pop OS
    must be re-verified to recompose;
  - **(b) attenuate/cap the `mTSD` conditioning coupling in `sim_os.do`** — cheap to prototype, isolates
    the engine lift with one parameter and bounds what (a) would buy. **Recommended next step: prototype
    (b) first.**

Related unrelated lever (for residual over-progression, not this issue): a per-line **"no further
treatment" / discontinuation** off-ramp — `tfi_plateau_check.do` found **no clean cure plateau**, so a
finite-mixture is not identifiable; a calibrated discontinuation probability is the remaining option
(not implemented).

### 7d. Branch `os-per-line` — per-line OS on the new engine (2026-07-06, pending validation)
Acting on §7c (the over-prediction is the engine's from-diagnosis `mTSD` conditioning, not the OS
shape), branch `os-per-line` (off the post-calibration `main`) replaces the single from-diagnosis
`OS#BCR` model with **one fitted OS model per pathway stage, each clocked from that stage's own entry
event** — so the engine draws survival on the line clock (`vElapsed = 0` at a line's first stage →
**unconditional**, no accumulated-time lift). This is the §2/§4 per-line design **re-applied on the new
engine**, now re-motivated to *remove* over-prediction (not the old under-prediction), and it keeps the
four comorbidity flags.
- `prep/risk_equations.do`: the OS block is now 13 stage models — `OS_DN`, `OS_L1`, `OS_L1_NoASCT`,
  `OS_L1_ASCT`, `OS_L2`, `OS_L2_End`, `OS_L3`, `OS_L3_End`, `OS_L4`, `OS_L4_End`, `OS_L5`, `OS_L5_End`,
  `OS_L6plus` — each `origin()`/`exit()` window-censored to its stage, `streg … CM_CKD CM_CRD CM_PLM
  CM_DBT i.BCR_L<n>, d($dOS)`, saved via `save_coefs` (→ `$Coeffs` → matsave). L6+ is a single
  conditional model from L6 start (sparse deep tail).
- `core/outcomes/sim_os.do`: per-stage firing map; design matrix = Age, Age2, Male, ECOG(0,1,2),
  RISS(1,2,3), the four comorbidity flags, BCR block (width read as `cols(vCoef) − 15`), cons; draws on
  the line clock and stores on the diagnosis clock (`originTSD + residual`) for `sim_mort`.
- `core/tests/extreme_value.do`: OS crank routed through `xv_bump_os(d)` (bumps every stage intercept
  at `cols−1`) since there is no single `bOS` anymore.
- Unchanged: `mata_setup.do` already loads the full ECOG/RISS dummies + comorbidity vectors; downstream
  (`process_data`, `sim_mort`, validators) unaffected (OS still lands on the diagnosis clock).

**Outcome — ADOPTED (merged to `main`).** Deterministic OOS **143/172 (76.7%→83.1%, +11 tests)**.
Whole-population OS preserved (74.6/59.0/30.1 vs obs 74.7/58.4/30.6, all pass), and the weak-responder
over-prediction is largely resolved: L1 No ASCT PR/MR/SD move from +15…+17pp over to within tolerance;
L2/L3 similarly. The over-prediction was the engine's `mTSD`/heavy-TFI coupling (§7c), and clocking each
stage from its own line start removes it. **Residual:** L1 No ASCT **CR now UNDER-predicts** (−8.0/−11.8),
and L2 SD swings under (−12/−17, partly the §3 SMM artefact) — the mirror of the old problem, the CR lift
having been removed. Per §7c the CR miss is ~half a 70→30 **transport** artefact (train CR ran below test
CR — should shrink in the full-data fit) and ~half **plateau fragmentation** (no single stage sees a CR
patient's multi-year plateau). Non-OS fails (TXD, TFI, pathways) are unchanged and pre-existing.

**Stale coefficients note:** a structural OS change invalidates the tracked coefficient files — after
merging, `analyses/*/coefficients/coefficients_*.mmat` must be **regenerated** (`risk_equations.do`) on
the drive before any sim runs; the committed single-`bOS` files will error in the per-line `sim_os.do`.

### 7e. `ancillary(i.BCR)` on L1-L3 — TRIED, no gain, NOT adopted (2026-07-06)
To lift the L1 No ASCT CR under-prediction (§7d), the shape was freed by response via `ancillary(i.BCR)`
on the L1-L3 OS models (`sim_os.do` split into a main equation + a per-patient log-shape block). **It was
a wash: OOS unchanged at 143/172, CR unmoved (−7.7/−11.4), the CR<VGPR inversion persisting** (sim CR
74.7 < VGPR 77.6), while L1 VGPR tipped FAIL and L3 CR tipped PASS — net zero. Confirms the §7c `SHAPE(MB)`
probe: the L1 CR miss is **not** a per-stage shape/PH problem — it is transportability + plateau
fragmentation, which a shape parameter cannot fix. **Reverted** (branch reset to the plain per-line).
Don't re-try per-BCR shape for this; if CR proves miscalibrated *in-sample* (validate `base_model`), the
lever is an AFT family (log-normal) on the L1 stages, not `ancillary`.
