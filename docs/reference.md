# Monash Myeloma Model — Reference (codes, structures, glossary)

A quick-reference for the load-bearing conventions used throughout the engine: the pathway-point (OMC) and
event coding, regimen and response codes, the Mata data structures, the characteristic vectors, and the
common-random-number (CRN) slot registry. For the architecture narrative see
[`technical_review.md`](technical_review.md); for costs/utilities see [`economic_inputs.md`](economic_inputs.md).

## Glossary

| Term | Meaning |
|---|---|
| **OMC** | Outcome-milestone checkpoint — a discrete pathway point (1–19; see below) |
| **DN** | Diagnosis |
| **LoT / L** | Line of therapy (1–9) |
| **L`x`S / L`x`E** | Start / end of line `x` |
| **BCR** | Best clinical response |
| **TXR** | Treatment regimen (drug combination) |
| **TXD** | Treatment duration (months on a line) |
| **TFI** | Treatment-free interval (months between finishing a line and the next event) |
| **OS** | Overall survival |
| **TNE** | Time to next event (duration of the current OMC) |
| **TSD** | Time since diagnosis (months) |
| **MOR** | Mortality flag |
| **OC** | Final outcome (death/censor time + mortality flag) |
| **ASCT / SCT** | Autologous stem-cell transplant |
| **MNT** | Maintenance therapy |
| **ISS / R-ISS** | (Revised) International Staging System |
| **ECOG** | ECOG performance status |
| **CRN** | Common random numbers |
| **CM_\*** | Comorbidity flags |
| **MRDR** | Australia and New Zealand Myeloma and Related Diseases Registry |

## Pathway points (OMC 1–19)

The patient journey is discretised into **19 checkpoints**. The engine walks them in order; outcome
matrices are `Obs × 19` to match.

| OMC | Point | OMC | Point | OMC | Point |
|---|---|---|---|---|---|
| 1 | DN (diagnosis) | 8 | L4S | 14 | L7S |
| 2 | L1S | 9 | L4E | 15 | L7E |
| 3 | L1E | 10 | L5S | 16 | L8S |
| 4 | L2S | 11 | L5E | 17 | L8E |
| 5 | L2E | 12 | L6S | 18 | L9S |
| 6 | L3S | 13 | L6E | 19 | L9E |
| 7 | L3E | | | | |

**Rule:** OMC 1 = diagnosis; for line `L`, start = `2L`, end = `2L+1`.

## Event codes (`Event0` / `Event1`)

The MRDR long event-history codes each row by its bounding events. `Event1` is the event at that row's
time point; `Event0` is the previous event. Risk equations condition on these (e.g. ASCT intent is fit
where `Event0 == 3`; ASCT receipt where `Event1 == 11`).

| Code | Event |
|---|---|
| `3` | Diagnosis |
| `L×10` (10, 20, … 90) | Start of line L |
| `L×10+1` (11, 21, … 91) | End of line L |
| `100` | Post-ASCT response point (L1-end transplant branch; feeds `BCR_SCT`) |
| `104` | Death |

## State

`State` (column 1 of `mState`) is the cohort's **entry pathway point**, equal to the OMC at which a patient
enters the simulation:

- **From diagnosis:** `State = 1` (OMC 1 = DN).
- **Line-`L` entry cohort** (built by `cohort_pool.do`): `State = 2L` (e.g. an L2-entry cohort sets
  `State = 4` = OMC 4 = L2S).

Outcome modules gate on `mState[.,1] <= OMC` so a patient becomes active once the pathway reaches their
entry point.

## Regimen codes (TXR / MRDR `Regimen`)

Per-line regimen lists are declared in each analysis's `outcomes/txr_<coeffs>.do`; any code not listed for
a line is pooled into `0 = other`.

| Code | Regimen | | Code | Regimen |
|---|---|---|---|---|
| 0 | Other (pooled) | | 9 | VTd — Bort/Thal/Dexa |
| 2 | TCd — Thal/Cyclo/Dexa | | 31 | VRd — Bort/Lena/Dexa |
| 4 | VCd — Bort/Cyclo/Dexa | | 49 | Kd — Carf/Dexa |
| 5 | Vd — Bort/Dexa | | 56 | Pd — Poma/Dexa |
| 7 | Rd — Lena/Dexa | | 80 | DVd — Dara/Bort/Dexa |

## BCR categories (1–6, best → worst)

| Code | Response |
|---|---|
| 1 | CR — complete response |
| 2 | VGPR — very good partial response |
| 3 | PR — partial response |
| 4 | MR — minimal response |
| 5 | SD — stable disease |
| 6 | PD — progressive disease |

## Comorbidity flags

Four binary flags enter the OS and ASCT-eligibility equations:

| Variable | Comorbidity |
|---|---|
| `CM_CKD` | Chronic kidney disease (renal impairment) |
| `CM_CRD` | Moderate-to-severe cardiac disease |
| `CM_PLM` | Moderate-to-severe chronic lung disease |
| `CM_DBT` | Diabetes (insulin / oral hypoglycaemics) |

## Mata data structures

Built in `core/mata_setup.do`. Line-varying outcomes are matrices (`Obs × 19` for pathway-point series,
`Obs × 9` for per-line series); fixed characteristics are vectors.

| Matrix | Shape | Contents |
|---|---|---|
| `mState` | `Obs×2` | entry State, diagnosis date |
| `mAge` | `Obs×19` | age at each pathway point |
| `mOS` | `Obs×19` | simulated overall-survival time from diagnosis (months) |
| `mTNE` | `Obs×19` | time to next event (duration of the current OMC) |
| `mTSD` | `Obs×19` | cumulative time since diagnosis (`TSD_DN = 0`) |
| `mMOR` | `Obs×19` | mortality flag (1 dead, 0 alive, `.` not reached) |
| `mOC` | `Obs×2` | final outcome: death/censor time, mortality flag |
| `mTXR` | `Obs×9` | regimen code per line |
| `mTXD` | `Obs×9` | treatment duration (months) per line |
| `mBCR` | `Obs×10` | best clinical response per line, + post-ASCT response (`BCR_SCT`) |
| `mTFI` | `Obs×9` | treatment-free interval (`TFI_DN` first; TFI for line L is column L+1) |
| `mRN` | `Obs×74` | pre-drawn CRN uniform matrix (freed before `process_data`) |

Each matrix has paired `rm*`/`cm*` row/column-label matrices.

### Characteristic vectors

`vID`, `vAge`, `vAge2` (age²), `vMale`; ECOG as `vECOG` and dummies `vECOG0/1/2`;
staging as `vRISS`/`vRISS1..3` and `vISS`/`vISS1..3`; age-threshold indicators `vAge70`, `vAge75`;
comorbidity flags `vCKD`, `vCRD`, `vPLM`, `vDBT`; the constant `vCons`; and treatment-receipt vectors
`vSCT_DN` (ASCT intent at diagnosis), `vSCT_L1` (ASCT receipt at L1), `vMNT` (maintenance receipt).

## Common-random-number (CRN) slot registry

`core/rng_slots.do` is the single source of truth for the columns of `mRN` (built once per cohort in
`mata_setup.do` as `runiform(Obs, rn_K())`, seed `$crn_seed_base + $b`). Every stochastic event reads its
**fixed** column for its **fixed** patient row, so both arms of a comparison see identical randomness —
the basis for variance-reduced cost-effectiveness. Total columns **K = 74**.

| Cols | Event | Keyed by | Accessor |
|---|---|---|---|
| 1–19 | Overall survival | OMC 1..19 | `rn_os(omc)` |
| 20–28 | Best clinical response | line 1..9 | `rn_bcr(line)` |
| 29 | BCR after ASCT | — | `rn_bcr_asct()` |
| 30–38 | Treatment regimen | line 1..9 | `rn_txr(line)` |
| 39–43 | TXD at L1 (5 sub-draws) | sub 1..5 | `rn_txd_l1(sub)` |
| 44–51 | TXD at L2+ | line 2..9 | `rn_txd(line)` |
| 52 | TFI at diagnosis | — | `rn_tfi_dn()` |
| 53–54 | TFI at L1 (ASCT / no-ASCT) | branch 1,2 | `rn_tfi_l1(branch)` |
| 55–62 | TFI at L2+ | line 2..9 | `rn_tfi(line)` |
| 63 | ASCT at diagnosis | — | `rn_asct_dn()` |
| 64 | ASCT at L1 | — | `rn_asct_l1()` |
| 65–66 | Maintenance (ASCT / no-ASCT) | branch 1,2 | `rn_mnt(branch)` |
| 67–74 | Reserved for analysis overrides | i 1..8 | `rn_override(i)` |

TXD-L1 sub-indices: 1–3 = ASCT splines, 4 = no-ASCT, 5 = continuous therapy. `rnDraw(idx, slot)` returns
`mRN[idx, slot]`. An override that **replaces** a core event reuses that event's slot (e.g. a BCR override
calls `rn_bcr($line)`); one that introduces a **new** stochastic event takes a reserved `rn_override(i)`
column.

## Risk equations & coefficients

The **50 risk equations** (per-line OS, BCR ordered logits, TXR multinomials, TXD/TFI parametric survival,
ASCT and maintenance logits) are fitted by `prep/risk_equations.do` into a per-analysis Mata coefficient
set `analyses/<analysis>/coefficients/coefficients_<coeffs>.mmat`, loaded at run time via
`mata matuse`. See [`prep/README.md`](../prep/README.md) and [`technical_review.md`](technical_review.md).
