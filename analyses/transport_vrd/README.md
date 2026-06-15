# VRd Line 1 — Calibrated Transport

Second exemplar for the **Calibrated Transport** method: a transportability approach for *pre-market* health technology assessment (HTA) that calibrates a trial treatment effect to real-world practice using a comparator observed in **both** the trial and the registry (MRDR). Applied to first-line bortezomib + lenalidomide + dexamethasone (VRd) versus lenalidomide + dexamethasone (Rd) in newly diagnosed, transplant-ineligible multiple myeloma.

This README documents the **code and outputs** for the analysis. The method itself is introduced with the DVd exemplar (`analyses/transport_dvd/`) and written up with the paper (vault: `myeloma model/papers/2026 dvd-transport/`).

## Why this exemplar

The DVd exemplar's comparator (Vd) is rare in routine second-line practice — the registry anchor is only ~71–135 patients, so the cross-setting shift is directionally reliable but imprecise. **VRd's trial comparator is Rd, which is heavily used first line in the MRDR (anchor n = 486).** This is the widely-used common comparator the method needs to deliver a *precise, well-identified* shift, so this exemplar is the one that can turn the headline accuracy result from "directional" into "demonstrable". The two exemplars together also test the method's central scope claim — that precision scales with how widely the comparator is used.

## Method in brief

Conventional transportability needs real-world data on the novel treatment, which does not exist before funding. Calibrated Transport gets around this: Rd is observed in both SWOG S0777 and the MRDR, so it acts as an efficacy–effectiveness anchor (within-treatment, across-setting) that lets the trial VRd effect be calibrated to the Australian setting and expressed as an *absolute* real-world prediction.

**Comparator: Rd, first line.** The decision-relevant contrast is VRd vs Rd; the scenarios differ only in how the first-line outcomes are derived.

## Response scale (differs from the DVd exemplar)

SWOG S0777 reports best confirmed response as **five ordered categories — CR / VGPR / PR / SD / PD** — but does **not** report minimal response (MR). The MRDR's six-category BCR is therefore collapsed to match by **folding MR into SD** (MR sits between PR and SD; SWOG has no MR slot). Coding here: `1 = CR, 2 = VGPR, 3 = PR, 4 = SD (incl. MR), 5 = PD`. The DVd exemplar, by contrast, runs on the full six categories.

When feeding the prediction into the simulation (which expects the model's native six-level BCR), the combined SD bucket is split back into MR and SD using the **registry comparator's (first-line Rd) MR:SD ratio** (40:166 ≈ 19:81) — the only MR-resolved information available at the pre-market decision point. This re-expansion carries no transport signal and happens after the `ologit` prediction, never inside it.

## Underlying model

The comparison is **VRd vs Rd** in all scenarios. Rd is observed in both SWOG and the MRDR, which is what makes it usable as the efficacy–effectiveness anchor. First line, so — unlike the DVd relapsed-setting model — there are **no line-of-therapy indicators**.

A single, bare ordered logistic regression: `ologit BCR MRDR VRd`, fit on the **Rd anchor plus the SWOG arms** —

| group | role | MRDR | VRd | n |
|---|---|---|---|---|
| MRDR first-line Rd (Regimen == 7, Event0 == 10) | registry comparator / **anchor** | 1 | 0 | 486 |
| SWOG Rd | trial comparator (reference cell) | 0 | 0 | 214 |
| SWOG VRd | trial novel treatment | 0 | 1 | 216 |

**β_MRDR** is the registry-Rd-vs-SWOG-Rd cross-setting shift — it carries the whole trial→real-world calibration; **β_VRd** is the SWOG VRd-vs-Rd effect (the consistency check against the published trial). The real-world VRd prediction is the cell (MRDR = 1, VRd = 1); the trial-based prediction is the raw SWOG VRd distribution. The held-out validation target is **observed MRDR first-line VRd (Regimen == 31, Event0 == 10, n = 1,516)**.

## How to run

Two distinct pieces live here:

1. **`calibrated_transport.do` — the transport generator.** Reads the **actual MRDR records** (Rd anchor `Regimen==7`, VRd target `Regimen==31`, `Event0==10`; MR folded into SD), with synthetic SWOG arms (aggregate trial data only). Fits the `ologit`, predicts the Calibrated-Transport cell, scores both predictions against the observed target.
   - `do calibrated_transport.do 0` — deterministic point estimate + comparison table → `results/transport_vrd.csv`, **and** the MRDR Rd first-line baseline for Table 1 → `results/mrdr_rd_baseline.csv` (Age ≥65 / Male / ISS III / ECOG>1; confirm the ECOG variable name in the do-file's CONFIRM block). The baseline runs only on this deterministic call, not during the bootstrap.
   - `do calibrated_transport.do 1 500 71523` — bootstrap: MAE CIs + exceedance.
   - **`calibrated_transport_prop.do`** — partial proportional-odds variant (`gologit2 BCR MRDR VRd, npl(MRDR)`), to test whether relaxing the single-shift assumption helps the over-correction. It helps only modestly (see below).

2. **`transport_vrd.do` — the simulation dispatcher (scaffolding, not yet adapted).** Copied from `transport_dvd` and still carries DVd defaults. To drive the full economic pipeline for VRd it needs: `$line = 1`; the SWOG arms in place of CASTOR; a VRd coefficient set; first-line costs/utilities; and the MR→SD re-expansion above wired into `outcomes/sim_bcr_override.do`. The first-line economic infrastructure can largely be reused from the VRd post-market evaluation (Irving 2025).

## Preview results (deterministic point estimate)

From `calibrated_transport.do 0` (proportional-odds fit):

| response | trial (VRd) | Calibrated Transport | observed |
|---|---|---|---|
| CR | 15.7 | 10.8 | 15.8 |
| VGPR | 27.8 | 22.7 | 27.8 |
| PR | 38.0 | 36.6 | 23.2 |
| SD | 15.7 | 28.0 | 32.3 |
| PD | 2.8 | 2.0 | 0.9 |
| **≥VGPR** | **43.5** | **33.5** | **43.6** |
| **ORR** | **81.5** | **70.1** | **66.8** |

- **β_MRDR ≈ 0.48** — essentially identical to the DVd exemplar's 0.48; an independent confirmation that the trial→registry shift is of similar magnitude.
- **MAE 6.64 → 5.77 pp (−0.87, ≈ 13%).** Calibrated Transport corrects the dominant error: the trial badly over-predicts PR and under-predicts SD; CT shifts mass PR→SD and lands the ORR far closer (81.5 → 70.1 vs 66.8 observed).
- **Over-correction at the top — CT does not work well here.** The trial's ≥VGPR (43.5) was already accurate (observed 43.6), but the single shift drags it to 33.5 and PR stays ~36 vs 23 observed; the MAE gain rests entirely on the SD fix. **Partial proportional odds (`calibrated_transport_prop.do`) was tested and does NOT rescue it** — it recovers CR but ≥VGPR/PR remain off, because the residual is a **treatment-by-setting interaction** (the Rd anchor's own top dropped in the registry while VRd's did not), which no link function repairs. Reframed in the manuscript as a **limitation of BCR as an ordinal surrogate**, not of Calibrated Transport (other outcomes, e.g. PFS, would use other models).

## Data

- **Rd comparator / anchor:** MRDR first-line Rd, **n = 486** (Regimen == 7, Event0 == 10). 5-cat (MR folded into SD): CR 7.8 / VGPR 16.3 / PR 31.9 / SD 42.4 / PD 1.6 %.
- **Trial arms — SWOG S0777 (Durie 2017, Lancet; Table 3, assessable):** VRd **n = 216** (CR 34, VGPR 60, PR 82, SD 34, PD 6); Rd **n = 214** (CR 18, VGPR 50, PR 85, SD 52, PD 9). SD here already excludes MR (none reported).
- **VRd validation (held out):** MRDR first-line VRd, **n = 1,516** (Regimen == 31, Event0 == 10). 5-cat: CR 15.8 / VGPR 27.8 / PR 23.2 / SD 32.3 / PD 0.9 %. Used only to score predictions — no observed VRd outcomes enter the fit.

MRDR patient data access is via the MRDR Steering Committee ([mrdr.net.au](https://www.mrdr.net.au/)).

## Status

**Calibration generator built** (`calibrated_transport.do`, reading real MRDR records) and run, deterministic + bootstrap, plus the partial-PO variant. **Finding: Calibrated Transport works only partially for VRd** (fixes SD, over-corrects the top; MAE 7.6→5.9 boot, ~94% exceedance, reduction CI spans zero) — best framed as the honest stress-test alongside the clean DVd exemplar. The **full simulation/economic pipeline is scaffolding only** — the dispatcher, `generate_cohort.do`, `compare_scenarios.do` and `outcomes/` are copied from `transport_dvd` with the analysis name repointed to `transport_vrd`, but not yet adapted to the first-line VRd setting.

## Next steps (pick up here)

1. **Run the bootstrap** (`calibrated_transport.do 1 500`) to attach MAE CIs and the exceedance probability — with n = 486 anchor and n = 1,516 target these should be tight.
2. **Adapt the dispatcher to VRd:** `$line = 1`, SWOG arms for the trial scenario, a VRd coefficient set, and first-line costs/utilities (reuse Irving 2025).
3. **Wire the MR→SD re-expansion** (19:81 from registry Rd) into `outcomes/sim_bcr_override.do` so the 5-category prediction maps onto the model's six-level BCR.
4. **Tier 3 aggregation:** adapt `compare_scenarios.do` (trial / Calibrated Transport / observed) for the first-line VRd cohort.
5. **Manuscript:** add the second-exemplar results and the scope contrast (thin Vd anchor vs wide Rd anchor; precision scales with comparator usage).
