# Base Model Analysis

## Overview

This is the base Monash Myeloma Model configuration that includes **all available treatment regimens** in the risk equations. This represents the full treatment landscape as observed in the MRDR registry data.

## What's Different About This Model

The base model includes **ALL** treatment regimens in the risk equations:

### Line 1 (LoT 1) Regimens Included
- **VCd** (Bortezomib, Cyclophosphamide, Dexamethasone) - 58%
- **VRd** (Bortezomib, Lenalidomide, Dexamethasone) - 15% ✅ **Included**
- **VTd** (Bortezomib, Thalidomide, Dexamethasone) - 8% 
- **Other regimens** - 19%

### Line 2 (LoT 2) Regimens Included  
- **Rd** (Lenalidomide, Dexamethasone) - 16%
- **DVd** (Daratumumab, Bortezomib, Dexamethasone) - 11%
- **Other regimens** - 73%

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
