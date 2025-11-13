# DVd Line 2 Methodological Analysis

**Bridging Trial Evidence and Real-World Data: A Novel Prediction Method for Pre-Funding Health Technology Assessment**

## Overview

This analysis demonstrates a novel hybrid approach for predicting treatment outcomes when clinical trial comparators have limited representation in real-world practice. Using daratumumab, bortezomib and dexamethasone (DVd) at second-line therapy as a case study, we compare three prediction methods to evaluate their accuracy in forecasting real-world outcomes.

## Research Context

### The Problem
When health technology assessment (HTA) agencies evaluate new therapies, they face a common challenge:
- Clinical trials compare new treatments against comparators that may be uncommon in real-world practice
- Trial populations often differ from patients who will actually receive the treatment
- Standard economic evaluations rely heavily on trial data, which may not reflect real-world effectiveness

### The DVd Case Study
- **DVd funded**: July 2020 in Australia (PBS listing)
- **Trial comparator**: Vd (bortezomib + dexamethasone) from CASTOR trial
- **Real-world challenge**: Only n=52 Vd patients at LoT 2 in MRDR (2013-2019)
- **Question**: Can we accurately predict DVd outcomes using limited Vd registry data?

## Three Prediction Scenarios

### Scenario A: Traditional Trial-Based Approach
**Method**: Uses published trial BCR distribution directly  
**Data source**: CASTOR trial published results (DVd vs Vd)  
**BCR distribution**: CR=19%, VGPR=40%, PR=24%, MR=4%, SD=10%, PD=2%  
**Represents**: Standard HTA practice of using trial data without adjustment

### Scenario B: Novel Hybrid Method ⭐
**Method**: Mixed-level data approach combining:
- Individual-level Vd patient data from MRDR (2013-2019, n=52)
- Aggregate DVd trial arm characteristics from CASTOR trial
- Multinomial logit regression with bootstrapped uncertainty

**Key Innovation**: Leverages registry heterogeneity to predict outcomes for new therapy  
**Uncertainty**: 500 bootstrap iterations reflecting n=52 sample limitation

### Scenario C: Post-Market Reality (Gold Standard)
**Method**: Uses actual DVd patient data from MRDR  
**Data source**: MRDR 2020-2024 (post-funding period)  
**Purpose**: Validation target - "what actually happened"  
**Represents**: True real-world outcomes

## Hypothesis

**Scenario B predictions will be closer to Scenario C (reality) than Scenario A predictions**, demonstrating that registry-based hybrid methods can improve pre-funding prediction accuracy even with limited comparator data.

## Key Methodological Contributions

1. **Mixed-level data integration**: Combines individual and aggregate data sources
2. **Small sample handling**: Demonstrates method works with n=52 comparator patients
3. **Uncertainty quantification**: Bootstrap approach captures sampling uncertainty
4. **Policy relevance**: Directly applicable to HTA decision-making process
5. **Generalizable framework**: Can be applied to other therapies with similar challenges

## Repository Structure
```
analyses/dvd_l2_method/
├── README.md                           # This file
├── dvd_l2_method.do                    # Main analysis dispatcher
│
├── data/
│   ├── coefficients/
│   │   ├── scenario_a_trial/           # Trial-based coefficients
│   │   │   ├── coefficients_trial.mmat
│   │   │   └── bootstrap/              # (500 iterations)
│   │   │
│   │   ├── scenario_b_hybrid/          # Hybrid method coefficients
│   │   │   ├── coefficients_hybrid.mmat
│   │   │   ├── bcr_predictions_bootstrap.dta
│   │   │   └── bootstrap/              # (500 iterations)
│   │   │
│   │   └── scenario_c_reality/         # Post-market coefficients
│   │       ├── coefficients_reality.mmat
│   │       └── bootstrap/              # (500 iterations)
│   │
│   ├── patients/
│   │   ├── dvd_predicted_patients.dta  # Patients predicted for DVd
│   │   └── vd_comparator_2013_2019.dta # Vd patients (n=52)
│   │
│   └── simulated/
│       ├── scenario_a/                  # Trial-based results
│       ├── scenario_b/                  # Hybrid method results
│       └── scenario_c/                  # Reality results
│
└── scripts/
    ├── 00_prepare_vd_data.do            # Extract Vd patients 2013-2019
    ├── 01_prepare_coefficients_scenario_a.do
    ├── 02_prepare_coefficients_scenario_b.do
    ├── 03_prepare_coefficients_scenario_c.do
    ├── 04_run_all_scenarios.do          # Master execution script
    └── 05_compare_scenarios.do          # Generate comparison tables/figures
```

## Data Requirements

### Input Data (Pre-2020)
- **Vd patients**: MRDR registry data, LoT 2, 2013-2019 (n=52)
  - Patient characteristics: age, sex, ECOG, R-ISS, previous BCR
  - Observed BCR to Vd treatment
  
- **Trial data**: CASTOR trial aggregate characteristics
  - DVd arm mean age, sex distribution, ECOG, R-ISS
  - Published BCR distribution for DVd and Vd arms

### Validation Data (Post-2020)
- **DVd patients**: MRDR registry data, LoT 2, 2020-2024
  - Observed BCR distribution
  - Used only for validation (Scenario C)

## How to Run the Analysis

### Prerequisites
```stata
// Ensure you're on the dvd-l2-method branch
// Required Stata version: 15.0 or higher
// Required packages: Standard EpiMAP Myeloma environment
```

### Step 1: Prepare Data and Coefficients
```stata
// Extract Vd patient data (2013-2019)
do "analyses/dvd_l2_method/scripts/00_prepare_vd_data.do"

// Generate coefficients for all three scenarios
do "analyses/dvd_l2_method/scripts/01_prepare_coefficients_scenario_a.do"
do "analyses/dvd_l2_method/scripts/02_prepare_coefficients_scenario_b.do"
do "analyses/dvd_l2_method/scripts/03_prepare_coefficients_scenario_c.do"
```

### Step 2: Run Simulations

#### Scenario A: Trial-Based (No Bootstrap Needed)
```stata
do "EpiMAP_Myeloma.do" dvd_l2_method DVd 2 trial Population 1 1000 0
```

#### Scenario B: Hybrid Method (With Bootstrap)
```stata
// Single run (uses mean predictions)
do "EpiMAP_Myeloma.do" dvd_l2_method DVd 2 hybrid Population 1 1000 0

// Full bootstrap (500 iterations) - for uncertainty quantification
do "EpiMAP_Myeloma.do" dvd_l2_method DVd 2 hybrid Population 1 1000 1 1 500
```

#### Scenario C: Reality (For Validation)
```stata
do "EpiMAP_Myeloma.do" dvd_l2_method DVd 2 reality Population 1 1000 0
```

#### Run All Scenarios at Once
```stata
do "analyses/dvd_l2_method/scripts/04_run_all_scenarios.do"
```

### Step 3: Compare Results
```stata
// Generate comparison tables and figures
do "analyses/dvd_l2_method/scripts/05_compare_scenarios.do"
```

## Key Outputs

### Primary Outcomes
1. **BCR Distribution Comparison**
   - Predicted vs observed BCR for each scenario
   - Absolute differences from reality (Scenario C)

2. **Cost-Effectiveness Analysis**
   - ICER for each scenario
   - Probability cost-effective at £50K/QALY threshold
   - Distance from "true" ICER (Scenario C)

3. **Uncertainty Analysis**
   - 95% confidence intervals for Scenario B predictions
   - Distribution of ICERs across 500 bootstrap iterations

### Figures
- **Figure 1**: BCR distribution comparison (all scenarios)
- **Figure 2**: Cost-effectiveness plane (all scenarios)
- **Figure 3**: Bootstrap distribution of Scenario B predictions
- **Figure 4**: Prediction accuracy by scenario

### Tables
- **Table 1**: Patient characteristics (Vd 2013-2019 vs DVd trial)
- **Table 2**: Predicted vs observed BCR by scenario
- **Table 3**: Cost-effectiveness results by scenario
- **Table 4**: Prediction accuracy metrics

## Technical Implementation Details

### BCR Override Mechanism
The analysis uses conditional logic in `core/outcomes/sim_bcr.do` to override BCR predictions at LoT 2 for DVd patients based on the scenario:
```stata
// Scenario A: Use trial BCR directly
if ("$Coeffs" == "trial") {
    // Apply CASTOR trial BCR distribution (19% CR, 40% VGPR, etc.)
}

// Scenario B: Use hybrid method predictions
else if ("$Coeffs" == "hybrid") {
    // Apply predicted BCR from multinomial logit on Vd patients
    // With bootstrap iteration-specific predictions if $Boot == 1
}

// Scenario C: No override needed
// Uses standard risk equations that include DVd
```

### Bootstrap Uncertainty Implementation
For Scenario B, uncertainty is captured by:
1. Resampling n=52 Vd patients with replacement (500 iterations)
2. Re-estimating multinomial logit on each bootstrap sample
3. Predicting DVd BCR using aggregate trial characteristics
4. Storing 500 different BCR probability vectors
5. Each simulation iteration uses one bootstrap prediction

This reflects sampling uncertainty from the limited Vd comparator data.

## Validation Strategy

### Internal Validation
- Out-of-sample prediction: Training on 2013-2019, validating on 2020-2024
- Bootstrap validation: 500 iterations for Scenario B
- Temporal separation: Clean pre/post funding periods

### External Validation
- Compare predictions to published DVd trial results
- Assess generalizability across patient subgroups (age, R-ISS, ECOG)

## Limitations

1. **Small comparator sample**: Only n=52 Vd patients (addressed via bootstrap)
2. **Temporal changes**: Treatment patterns may have evolved 2013-2019 to 2020-2024
3. **Trial-registry differences**: DVd trial population may differ from MRDR
4. **Single therapy**: Case study approach - generalizability to be established
5. **Registry data quality**: Depends on MRDR completeness and accuracy

## Comparison to Standard Approaches

| Feature | Traditional (A) | Our Method (B) | Post-Market (C) |
|---------|----------------|----------------|-----------------|
| **Data source** | Trial only | Trial + Registry | Registry only |
| **Population** | Trial patients | Real-world patients | Real-world patients |
| **Timing** | Pre-funding ✓ | Pre-funding ✓ | Post-funding only |
| **Heterogeneity** | Limited | High | High |
| **Sample size** | Trial N | Registry n=52 | Registry N |
| **Uncertainty** | Trial CI | Bootstrap | Observed |
| **Policy relevance** | Standard HTA | Enhanced HTA | Validation only |

## Associated Publication

**Status**: In preparation

**Proposed title**: "Bridging Trial Evidence and Real-World Data for Pre-Funding Economic Evaluation: A Mixed-Level Data Approach for Predicting Treatment Outcomes in Multiple Myeloma"

**Target journals**: 
- Medical Decision Making
- Value in Health
- PharmacoEconomics
- Health Economics

**Authors**: [To be determined]

## Citation

When using this analysis, please cite:
```
Irving A, Petrie D, Harris A, et al. (2025). 
DVd Line 2 Methodological Analysis. 
EpiMAP Myeloma Repository. 
https://github.com/[your-repo]/analyses/dvd_l2_method
```

## Version History

- **v0.1** (2025-01): Initial structure and Scenario A implementation
- **v0.2** (2025-01): Added Scenario B with bootstrap uncertainty
- **v0.3** (2025-02): Added Scenario C and comparison scripts
- **v1.0** (TBD): Final version for publication

## Contact

For questions about this analysis:
- **Lead analyst**: Adam Irving ([adam.irving@monash.edu](mailto:adam.irving@monash.edu))
- **Project lead**: [Name]
- **GitHub issues**: [Link to repository issues]

## Acknowledgments

This analysis uses data from the Australia & New Zealand Myeloma and Related Diseases Registry (MRDR). We thank the patients, clinicians, and research staff at participating centres for their contributions to the registry.

The EpiMAP Myeloma model is supported by Monash University's MASSIVE High Performance Computing facility.

## License

This analysis is part of the EpiMAP Myeloma project, released under GPL v3.0.

---

**Branch**: `dvd-l2-method`  
**Last updated**: 3/11/25  
**Status**: 🚧 In development
