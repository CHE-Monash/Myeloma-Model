# EpiMAP Myeloma Validation System

## Overview

This system extracts validation benchmarks from MRDR training data and compares simulation results against them to ensure model accuracy.

## Files

1. **generate_benchmarks.do** - Extract benchmarks from MRDR training data (run once)
2. **validate_simulation.do** - Compare simulation to benchmarks (run after each simulation)
3. **validation/load_benchmarks.do** - Helper script to load benchmarks into Mata

## Workflow

### Step 1: Generate Benchmarks (One-time)

```stata
// Run on MRDR training data
use "data/mrdr_training.dta", clear
do "generate_benchmarks.do"
```

**Outputs:**
- validation/benchmarks/OS_L1.csv
- validation/benchmarks/OS_ASCT.csv
- validation/benchmarks/OS_L2.csv
- validation/benchmarks/TFI_L1_NoASCT.csv
- validation/benchmarks/TFI_L1_ASCT.csv
- validation/benchmarks/TFI_L2.csv
- validation/benchmarks/TFI_L3.csv
- validation/benchmarks/TXD_L1.csv
- validation/benchmarks/TXD_L2.csv
- validation/benchmarks/BCR_distributions.csv
- validation/benchmarks/Pathways.csv

### Step 2: Run Simulation

```stata
do "main_simulation.do"
```

### Step 3: Validate Results

```stata
do "validate_simulation.do"
```

**Tests performed:**
1. BCR distributions by line (±5% tolerance)
2. TFI medians by BCR (±20% tolerance)
3. OS at 3-year, 5-year by BCR (±10% tolerance)
4. TXD medians by BCR (±20% tolerance)
5. Treatment pathways (±5% tolerance)

## Validation Tolerances

| Outcome | Metric | Tolerance | Reason |
|---------|--------|-----------|--------|
| BCR | Percentages | ±5% | Stable distributions |
| TFI | Median | ±20% | Censoring, parametric underestimation |
| OS | 3/5-year survival | ±10% | Well-measured timepoints |
| TXD | Median | ±20% | Similar to TFI |
| Pathways | Percentages | ±5% | Cumulative proportions |

## Expected Discrepancies

**TFI underestimation (10-20%):**
- Parametric models underestimate due to censored training data
- Weibull tail behavior differs from Kaplan-Meier
- Acceptable if ordering (CR > VGPR > PR) preserved

**Poor responder categories (MR, SD):**
- Wider tolerance due to sparse data (N<200)
- Higher uncertainty in parameter estimates

**Late-line outcomes (L5+):**
- Limited training data
- Higher acceptable deviation

## Interpreting Results

**Pass rate >90%:** Model accurately reproduces training data
**Pass rate 75-90%:** Minor discrepancies, review failed tests
**Pass rate <75%:** Significant model issues, investigate code

## When to Re-generate Benchmarks

- MRDR training data updated
- Imputation strategy changed
- Variable definitions changed
- Model structure changed (e.g., BCR categories)

## Integration with Git

Commit benchmark CSV files to version control to track:
- Changes in MRDR data over time
- Model performance across versions
- Regression detection

## Notes

- Benchmarks use mi extract 1 (first imputation) for deterministic results
- Validation uses simulation output, not live MRDR data
- Variable names in validate_simulation.do must match simulation output
