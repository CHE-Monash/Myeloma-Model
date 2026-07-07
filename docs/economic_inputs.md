# Monash Myeloma Model — Economic Inputs

Reference for the cost and utility parameters the model uses to turn a simulated treatment pathway into
costs and quality-adjusted life years (QALYs). All of the economic logic lives in **`core/process_data.do`**
(the `process_data` program), which runs at the end of every simulation pass. This document lists the
parameter values, how they are applied, and how discounting works.

- **Currency:** Australian dollars (AUD).
- **Perspective / price basis:** the unit costs below are fixed values embedded in `process_data.do`.
  The dispatcher's `$cost_year` global (`2020` for `transport_dvd`, `2025` for `base_model`) is a **label
  only** — it does **not** reprice the hardcoded values. Record the actual price year with each analysis and
  cross-check the source citations in the relevant manuscript/supplement (drug acquisition costs are
  PBS-based; hospital/community/emergency and ASCT costs from Australian costing sources).
- **Discount rate:** `$drate` (default **0.05** = 5% per year, the PBAC reference rate).

> ⚠️ These parameters are currently **hardcoded as locals** in `process_data.do` (no external cost file).
> To reprice or run a cost sensitivity analysis, edit the `local c*`/`local u*` values there. Keep this
> document in step with those edits.

## Discounting

Costs and QALYs are discounted with a **continuous (uniform-over-interval) approximation**. A stream that
accrues a total undiscounted amount `X` uniformly over an interval `[t0, t1]` (times in months) has
discounted value

```
X · [ (1+r)^(−t0/12) − (1+r)^(−t1/12) ] / ( ln(1+r) · (t1−t0)/12 )
```

where `r = $drate`. One-off costs (ASCT) are discounted at the event time: `X · (1+r)^(−t/12)`. Discounting
uses **relative-time markers** (`TSD_*_ref`): for a from-diagnosis run (`$line == 1`) time 0 is diagnosis;
for a line-`L` entry run (`$line > 1`) time 0 is the start of line `L`.

## Treatment (drug) costs

Each regimen has a **unit cost per cycle** (or per course) applied over the simulated treatment duration
`TXD_L*` (months). `TXD × 30.4375 / cycle_days` converts months to number of cycles. The bortezomib-based
induction regimens are **capped** at a fixed number of cycles (fixed-duration induction); the others accrue
per cycle for the whole `TXD`.

| Regimen | TXR code | Unit cost (AUD) | Cycle | Cap | Applied as |
|---|---|---|---|---|---|
| VCd (Bort/Cyclo/Dexa) | 4 | 902 | 21 d | 4 cycles | `902 · min(4, TXD·30.4375/21)` |
| VRd (Bort/Lena/Dexa) | 31 | 1,776 | 21 d | 5 cycles | `1776 · min(5, TXD·30.4375/21)` |
| Vd (Bort/Dexa) | 5 | 724 | 21 d | 8 cycles | `724 · min(8, TXD·30.4375/21)` |
| Rd (Lena/Dexa) | 7 | 1,608 | 28 d | — | `1608 · TXD·30.4375/28` |
| Kd (Carf/Dexa) | 49 | 15,025 | 28 d | — | `15025 · TXD·30.4375/28` |
| DVd (Dara/Bort/Dexa) | 80 | 12,110 | 28 d | — | `12110 · TXD·30.4375/28` |
| Pd (Poma/Dexa) | 56 | 2,291 | 28 d | — | `2291 · TXD·30.4375/28` |
| Other (VTd / TCd / Td / …) | 0 | 1,612 | 28 d | — | `1612 · TXD·30.4375/28` |

## ASCT and maintenance (L1 only)

| Item | Local | Cost (AUD) | Applied when | Applied as |
|---|---|---|---|---|
| Autologous stem-cell transplant | `cASCT` | 41,723 | `SCT_L1 == 1` | one-off at L1 end |
| Maintenance therapy | `cMNT` | 1,329 | `MNT == 1` | per 28-day cycle over `TFI_L1`: `1329 · TFI_L1·30.4375/28` |

## Non-treatment (health-state) costs

Applied as an **annual rate** over the survival time from the starting line (`OC_TIME_L`, months):
`(cHosp + cComm + cEmer) · OC_TIME_L/12`.

| Component | Local | Annual cost (AUD) |
|---|---|---|
| Hospitalisation | `cHosp` | 38,743 |
| Community care | `cComm` | 10,928 |
| Emergency | `cEmer` | 2,476 |
| **Total non-treatment** | | **52,147 / year** |

## Utilities (QALYs)

QALYs are `Σ (time_in_state / 12) · utility_state`. The model uses a small set of health-state utilities by
pathway phase:

| Health state | Local | Utility |
|---|---|---|
| Treatment-free interval (off treatment / remission) | `uTFI` | 0.72 |
| On first-line treatment (`TXD_L1`) | `uTXD_L1` | 0.63 |
| On second-line treatment (`TXD_L2`) | `uTXD_L2` | 0.67 |
| Post-second-line (everything after `TXD_L2`, and all of L3+) | `uPostL2` | 0.63 |

How they map onto the pathway depends on the starting line (`$line`):

- **`$line == 1` (from diagnosis):** pre-L1 TFI (`uTFI`) → on-L1 (`uTXD_L1`) → L1 TFI (`uTFI`) → on-L2
  (`uTXD_L2`) → everything after L2 (`uPostL2`).
- **`$line == 2`:** on-L2 (`uTXD_L2`) → post-L2 (`uPostL2`).
- **`$line >= 3`:** the whole run from line `L` uses `uPostL2` (0.63).

## Output variables

`process_data` writes cost and QALY variables to the per-patient dataset (all AUD / QALYs):

- **Costs:** `cost_tx_L<line>`, `cost_tx_asct`, `cost_tx_mnt`, `cost_tx` (total treatment), `cost_nt`
  (non-treatment), `cost_total`. A **`_d` suffix marks the discounted** counterpart (e.g. `cost_total_d`).
- **QALYs:** `qaly_tfi_DN`, `qaly_txd_L1`, `qaly_tfi_L1`, `qaly_txd_L2`, `qaly_post_L2`, `qaly_total`
  (and their `_d` discounted forms), collapsing to `qaly_txd_L2 + qaly_post_L2` at `$line == 2` and a
  single `qaly_total` at `$line >= 3`.

The ICER machinery (e.g. `transport_dvd`) collapses `cost_total_d` and `qaly_total_d` by arm.

## Caveats & maintenance notes

- **No time-varying or health-state-specific costs beyond the annual non-treatment bundle** — disease
  management cost is a flat annual rate over survival, not resolved by line or response.
- **Utilities are phase-based, not response-based** — BCR (response) does not directly modulate utility.
- **The `Other` (code 0) drug cost is a single blended value** covering all pooled regimens.
- **Reprice/sensitivity:** edit the `local c*` / `local u*` values in `core/process_data.do`; there is no
  parameter file to override them at run time. Confirm the price year and source citations against the
  analysis manuscript before publication.

See also: [`reference.md`](reference.md) for the pathway-point / regimen / duration variables these
calculations consume, and [`technical_review.md`](technical_review.md) for the simulation architecture.
