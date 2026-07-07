# Base Model Analysis

## Overview

This is the base Monash Myeloma Model configuration that includes **all available treatment regimens** in the risk equations. This represents the full treatment landscape as observed in the MRDR registry data.

## What's Different About This Model

Unlike the intervention analyses, the base model keeps **all** treatment regimens in play. At each line the most common regimens are modelled as distinct categories and everything else is pooled into an "other" category. The per-line regimen lists are defined in `outcomes/txr_base_model.do`:

| Line | Distinctly modelled regimens (code) | Pooled |
|---|---|---|
| **L1** | VCd (4), **VRd** (31) | all other regimens → "other" |
| **L2** | Rd (7), DVd (80) | " |
| **L3** | Kd (49), Rd (7) | " |
| **L4** | Kd (49), Pd (56) | " |
| **L5–L9** | — (none listed) | all regimens → "other" |

VCd = Bortezomib/Cyclophosphamide/Dexamethasone, VRd = Bortezomib/Lenalidomide/Dexamethasone, Rd = Lenalidomide/Dexamethasone, DVd = Daratumumab/Bortezomib/Dexamethasone, Kd = Carfilzomib/Dexamethasone, Pd = Pomalidomide/Dexamethasone. **VRd is included at L1** — the key contrast with the VRd post-market analysis, which excludes it.

## When to Use This Model

Use the base model for:
- Standard population projections with current treatment practices
- Comparing effectiveness across all available regimens
- Baseline estimates for health economic models

## Key Files

- **Runbook**: `run.do` (risk equations → simulate; deterministic track + bootstrap/HPC record)
- **Dispatcher**: `simulate.do` (configure via globals; bootstrap via `do simulate.do 1 1 500`)
- **Coefficients**: `coefficients/coefficients_base_model.mmat` (all regimen effects)
- **Population**: a cohort from the repo-root `patients/` folder (e.g. `population_1995_2040_1.dta`)

## Usage

Configure the run by editing the globals at the top of `simulate.do` (`$int`, `$line`, `$data`, `$min_id`/`$max_id`, …), then run it from the repository root:

```stata
do "analyses/base_model/simulate.do"
```

## Comparison with Other Analyses

| Analysis | VRd Status | Purpose |
|----------|------------|---------|
| **Base Model** | ✅ Included | Current practice projections |
| **VRd Post-Market** | ❌ Excluded | VRd impact assessment |

The VRd post-market analysis excludes VRd to estimate what outcomes would have been without VRd availability.
