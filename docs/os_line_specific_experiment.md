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

## 5. Open questions / directions to brainstorm

1. **Clean the OS target.** Flag/exclude likely-SMM early-era long survivors (diagnosis era +
   indolent course) and re-benchmark. Then the residual gap — especially the **3-yr weak-responder
   gap in the best-populated data** (e.g. L2/L3 MR/SD) — is trustworthy and tells us whether real
   model error remains.
2. **Restrict the analysis horizon** (e.g. 5–6 yr) if the far tail stays untrustworthy — but note
   this truncates life-years/QALYs for genuinely long-surviving MM, a modelling decision to document.
3. **Adopt a per-BCR Weibull shape.** §4 settles this: `ancillary(i.BCR)` wins on AIC/BIC and 5-yr
   accuracy; `ancillary(i.OS)` (shape by line) and gamma-frailty do not. Remaining work is the
   engine wiring — fit `ln_p` as a function of BCR in `prep/risk_equations.do` and carry a shape
   vector indexed by BCR in `core/outcomes/sim_os.do` (currently one shared scalar `aux`).
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
