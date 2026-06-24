# VRd LoT 1 Post-Market Analysis

## Research Question

**"What is the clinical impact of VRd at LoT 1 compared to historical practice?"**

## Key Methodology

This analysis **excludes VRd from the risk equations** to simulate what would have happened to VRd patients under historical treatment patterns.

## Treatment Regimens in Risk Equations

### What's Excluded
- **VRd** (Bortezomib, Lenalidomide, Dexamethasone) ❌ **EXCLUDED**

### What's Included
- **VCd** (Bortezomib, Cyclophosphamide, Dexamethasone) ✅
- **VTd** (Bortezomib, Thalidomide, Dexamethasone) ✅
- **Other historical regimens** ✅

## How It Works

1. **Identify VRd patients**: Patients who would receive VRd based on their characteristics
2. **Simulate without VRd**: What regimen would they have gotten historically?
3. **Compare outcomes**: VRd effectiveness vs historical alternatives

## Analysis Scenarios

The scenario is selected by the `$int` global in `vrd_post.do`; run the dispatcher from the repository root.

### Scenario 1: Historical Practice (SoC)

Set `$int "SoC"`, then:

```stata
do "analyses/vrd_post/vrd_post.do"
```

- VRd patients get alternative regimens (VCd, VTd, Other)
- Shows what outcomes would have been without VRd

### Scenario 2: VRd Available

Set `$int "VRd"`, then:

```stata
do "analyses/vrd_post/vrd_post.do"
```

- VRd patients receive VRd
- Shows actual VRd impact

## Key Files

- **Dispatcher**: `vrd_post.do` (configure via globals; `$int` toggles SoC / VRd)
- **Coefficients**: `coefficients/coefficients_vrd_post.mmat` (VRd excluded)
- **Patient cohort**: `patients/patients_vrd_l1_post.dta` (VRd-eligible patients; loaded via `$cohort_file`)
- **Bootstrap coefficients**: samples in `coefficients/bootstrap/`

## Expected Results

Compare survival, response rates, and treatment pathways between scenarios to quantify VRd's clinical benefit.

## Files Structure

```
analyses/vrd_post/
├── README.md                       # This file
├── vrd_post.do                     # Dispatcher (configure via globals)
├── coefficients/
│   ├── coefficients_vrd_post.mmat  # VRd-excluded coefficients
│   └── bootstrap/                  # bootstrap coefficient samples
├── patients/
│   └── patients_vrd_l1_post.dta    # VRd-eligible cohort
└── simulated/                      # run outputs
```
