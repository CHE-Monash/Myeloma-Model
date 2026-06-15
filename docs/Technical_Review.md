---

editor_options: 
  markdown: 
    wrap: 72
---

# Monash Myeloma Model - Technical Review

## Executive Summary

The Monash Myeloma Model is a sophisticated **discrete-event simulation (DES)** model for multiple myeloma disease progression and treatment outcomes. Version 2.1 represents a major architectural transformation from loop-based to **vectorised matrix operations** using Stata's Mata language, achieving \~700× performance improvements whilst maintaining identical outputs.

------------------------------------------------------------------------

## Architecture Overview

### Core Design Philosophy

The simulation follows a **modular, state-based architecture** where: - Patient cohorts flow through discrete outcome milestone checkpoints (OMCs) - Each checkpoint represents a specific disease state (diagnosis, line start, line end) - Outcomes at each checkpoint are predicted using evidence-based risk equations - All computations operate on entire patient cohorts simultaneously (vectorised)

### Key Architectural Components

```         
repo/
├── core/                          # Simulation engine
│   ├── simulation_engine.do       # Main orchestration logic
│   ├── mata_setup.do              # Matrix/vector initialisation  
│   ├── load_patients.do           # Patient cohort loading
│   ├── mata_functions.do          # Utility functions
│   ├── process_data.do            # Post-simulation data assembly
│   └── outcomes/                  # Outcome prediction modules
│       ├── sim_os.do             # Overall survival
│       ├── sim_bcr.do            # Best clinical response
│       ├── sim_txr.do            # Treatment regimen
│       ├── sim_txd*.do           # Treatment duration
│       ├── sim_tfi*.do           # Treatment-free intervals
│       ├── sim_asct*.do          # ASCT eligibility/outcomes
│       ├── sim_mnt.do            # Maintenance therapy
│       ├── sim_mort.do           # Mortality tracking
│       └── sim_age.do            # Age progression
├── analyses/                      # Analysis configurations
│   ├── base_model/               # Base model (all regimens)
│   ├── vrd_l1_post/              # VRd post-market analysis
│   └── dvd_l2_method/            # DVd L2 methodology
├── patients/                      # Patient data
│   ├── population/               # Population projections
│   └── predicted/                # Trial-based cohorts
└── tests/                        # Validation suite
    ├── validate_vectors.do
    └── test_*.do
```

------------------------------------------------------------------------

## Simulation Flow

### High-Level Execution Pipeline

``` stata
1. Load Coefficients    → mata matuse "coefficients.mmat"
2. Load Patients        → load_patients (from .dta file)
3. Setup Vectors        → mata_setup (initialise all matrices)
4. Run Simulation       → simulation (30 outcome equations)
5. Process Results      → process_data (assemble output dataset)
6. Validate & Report    → validation.do, generate_report.do
```

### Outcome Milestone Checkpoints (OMC)

The simulation progresses through 19 discrete checkpoints:

| OMC | Checkpoint | Key Outcomes |
|----------------|--------------------------|------------------------------|
| 1 | Diagnosis (DN) | ASCT eligibility, TFI to L1, OS, Mortality |
| 2 | Line 1 Start (L1S) | Age update, TXR, TXD, OS, Mortality |
| 3 | Line 1 End (L1E) | BCR, ASCT receipt, ASCT BCR, Maintenance, TFI, OS, Mortality |
| 4 | Line 2 Start (L2S) | Age, TXR, TXD, OS, Mortality |
| 5 | Line 2 End (L2E) | BCR, TFI, OS, Mortality |
| ... | Lines 3-9 | Repeated pattern (Start → End) |
| 19 | Line 9 End (L9E) | Final BCR, OS, Terminal mortality |

At each checkpoint, the simulation: 1. **Filters** patients who are alive and have reached this state 2. **Predicts** outcomes using risk equations (ordered logit, parametric survival, multinomial logit) 3. **Updates** patient state matrices with new outcomes 4. **Tracks** cumulative time since diagnosis

------------------------------------------------------------------------

## Data Structures

### Core Patient Matrices (All in Mata Memory Space)

#### **Input Vectors** (Patient Characteristics)

``` stata
vID        : Patient identifiers (n × 1)
vAge       : Age at diagnosis (continuous)
vAge2      : Age squared (for quadratic effects)
vMale      : Sex indicator (0=Female, 1=Male)
vECOG      : Original ECOG category (0, 1, 2+)
vECOG0/1/2 : ECOG dummy variables
vRISS      : Original R-ISS category (1, 2, 3)
vRISS1/2/3 : R-ISS dummy variables
vCMc       : Comorbidity score (0, 1, 2, 3+)
vCM0/1/2/3 : Comorbidity dummy variables
vCKD       : Chronic kidney disease indicator
vAge70/75  : Age threshold indicators
vCons      : Constant vector (all 1s)
vSCT_DN    : ASCT eligibility at diagnosis
vSCT_L1    : ASCT receipt at Line 1
vMNT       : Maintenance therapy receipt
```

#### **Outcome Matrices** (n × 19 columns, one per OMC)

``` stata
mAge       : Age at each checkpoint (continuous time)
mOS        : Overall survival cumulative time
mTNE       : Time to next event (treatment duration or TFI)
mTSD       : Time since diagnosis (cumulative)
mMOR       : Mortality indicator (0=alive, 1=dead)
mOC        : Outcome container [Time, MortalityFlag]
mTXR       : Treatment regimen (1-9 for each line)
mTXD       : Treatment duration (months per line)
mBCR       : Best clinical response (1-6: CR, VGPR, PR, MR, SD, PD)
              Column 10: BCR post-ASCT
mTFI       : Treatment-free intervals (diagnosis + Lines 1-8)
mState     : State tracking for eligibility logic
```

------------------------------------------------------------------------

## Risk Equation Architecture

### 30 Evidence-Based Models

The simulation implements **30 distinct risk equations** derived from MRDR registry data:

| ID | Outcome | Model Type | Key Covariates |
|---------------|---------------|------------------|-------------------------|
| 1-2 | Overall Survival (DN/L1S) | Parametric survival | Age, Sex, ECOG, R-ISS, CM |
| 3 | ASCT Eligibility (DN) | Logistic | Age, ECOG, R-ISS, CM, CKD |
| 4 | Diagnosis to Treatment Interval | Parametric survival | Age, Sex, R-ISS |
| 5-7 | L1 TXD ASCT (3 manual splines) | Parametric survival | Age, Sex, ECOG, TXR |
| 8 | L1 TXD Non-ASCT | Parametric survival | Age, Sex, ECOG, R-ISS, TXR |
| 9 | L1 BCR | Ordered logistic | Age, Sex, ECOG, R-ISS, TXR |
| 10 | ASCT Receipt (L1) | Logistic | Age, ECOG, R-ISS, L1 BCR |
| 11 | ASCT BCR | Ordered logistic | Age, Sex, ECOG, R-ISS, L1 BCR |
| 12 | Maintenance Therapy | Logistic | Age, ECOG, ASCT status, L1 BCR |
| 13-14 | L1 TFI (ASCT vs Non-ASCT) | Parametric survival | Age, Sex, ECOG, R-ISS, L1 BCR |
| 15-30 | L2-L9 TXR, TXD, BCR, TFI, OS | Varied | Line-specific with prior BCR |

### Common Modelling Approaches

#### **1. Ordered Logistic Regression (BCR)**

``` stata
// Multinomial outcome with ordered categories
// Reference: Progressive Disease (PD)
Pr(BCR = CR) = exp(XB_CR) / [1 + Σexp(XB_k)]
Pr(BCR = VGPR) = exp(XB_VGPR) / [1 + Σexp(XB_k)]
...
// Normalisation ensures: Σ Pr(BCR_k) = 1
```

#### **2. Parametric Survival (OS, TXD, TFI)**

``` stata
// Generalised gamma distribution with 3 parameters
log(time) = XB + σ * (κ * z) where z ~ N(0,1)
// Extracted coefficients: XB, σ (sigma), κ (kappa)
// Transform: exp(XB + error_term) with proper error bounds
```

#### **3. Logistic Regression (Binary Outcomes)**

``` stata
// ASCT eligibility, maintenance therapy
Pr(Y = 1) = invlogit(XB) = exp(XB) / [1 + exp(XB)]
// Draw: runiform() < Pr(Y=1)
```

------------------------------------------------------------------------

## Vectorisation Strategy (v2.0 → v2.1 Transformation)

### The Challenge

**v2.0 Loop-Based Approach:**

``` stata
forvalues i = 1/$Obs {
    scalar age = Age[`i']
    scalar male = Male[`i']
    // ... extract all covariates individually
    scalar XB = b1*age + b2*male + ...
    mata: mOS[`i', OMC] = exp(XB + error)
}
```

**Problem:** 10,000 patients × 19 checkpoints × multiple outcomes = \~500,000+ individual loop iterations

### The Solution

**v2.1 Vectorised Approach:**

``` stata
mata {
    // Extract all alive patients at this checkpoint
    idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
    
    // Build patient matrix (n_alive × k_predictors)
    pMat = (vAge[idx], vAge2[idx], vMale[idx], vECOG1[idx], ...)
    
    // Extract coefficient vector (1 × k_predictors)
    coefs = bMatrix[coefCols]
    
    // Matrix multiplication: (n × k) × (k × 1) = (n × 1)
    XB = pMat * coefs'
    
    // Vectorised error generation and transformation
    errors = rnormal(rows(idx), 1, 0, sigma)
    outcomes = exp(XB + bounded(errors, kappa, -maxError, maxError))
    
    // Update outcome matrix in one operation
    mOS[idx, OMC] = outcomes
}
```

**Performance Gain:** - **Before:** 23+ minutes for 10,000 patients - **After:** \<2 seconds for 10,000 patients - **Speedup:** \~700×

------------------------------------------------------------------------

## Key Mata Utility Functions

### `bounded()` - Error Term Constraint

``` stata
real colvector bounded(real colvector values, 
                       real scalar maxValue, 
                       real scalar minValue)
```

**Purpose:** Constrain error terms to prevent implausible predictions - Used in survival models to prevent survival times \<0 or \>age limit - Applies limits only to non-missing values - Returns clamped vector: `min(max(value, minValue), maxValue)`

### `validateDimensions()` - Matrix Multiplication Check

``` stata
void validateDimensions(real matrix pMat, 
                        real rowvector coefs, 
                        string scalar matrixName)
```

**Purpose:** Catch dimension mismatches before matrix multiplication - Critical for debugging coefficient extraction errors - Provides detailed diagnostic output if cols(pMat) ≠ cols(coefs)

### `validateOutcomes()` - Categorical Outcome Validation

``` stata
void validateOutcomes(real colvector outcomes, 
                      real rowvector validValues, 
                      string scalar outcomeName)
```

**Purpose:** Ensure categorical predictions fall within expected ranges - Used for BCR (1-6), TXR regimen codes, binary indicators - Warns if invalid values detected (helps catch logic errors)

------------------------------------------------------------------------

## Treatment Pathway Modelling

### Line 1 Special Handling

**ASCT Pathway:**

```         
Diagnosis → ASCT Eligibility Check
    ↓ (if eligible)
L1 Start → Induction Therapy (TXR, TXD, BCR)
    ↓
L1 End → ASCT Receipt Decision
    ↓ (if receives ASCT)
    → ASCT BCR (separate equation, typically better)
    → Maintenance Therapy Decision
    → Treatment-Free Interval (ASCT-specific coefficients)
```

**Non-ASCT Pathway:**

```         
Diagnosis → ASCT Ineligibility
    ↓
L1 Start → Standard Therapy (TXR, TXD, BCR)
    ↓
L1 End → No ASCT
    → Maintenance Therapy Decision
    → Treatment-Free Interval (non-ASCT coefficients)
```

### Regimen Selection Logic

**Line 1 (3 regimens modelled explicitly):** - **VRd** (Bortezomib-Lenalidomide-Dex): 15% of cohort, superior efficacy - **VCd** (Bortezomib-Cyclophosphamide-Dex): 58% of cohort, standard care - **Other**: 26% of cohort, averaged survival benefit

**Line 2+ (regression to mean approach):** - **DVd** (Daratumumab-Bortezomib-Dex): Modelled explicitly at L2 in some analyses - **Rd** (Lenalidomide-Dex): Common L2 option (16% in MRDR) - **Other**: Averaged coefficients for other regimens at L2+

**Treatment Regimen Codes:**

``` stata
// Encoded as integers in mTXR matrix
1 = VRd
2 = VCd  
3 = Other (L1)
4 = DVd (L2)
5 = Rd (L2)
6-9 = Other regimens (L2+)
```

------------------------------------------------------------------------

## Mortality and Survival Tracking

### Competing Risks Framework

At every checkpoint, patients face **two competing risks:**

1.  **Event-driven mortality**: Death before next treatment/event
2.  **Censoring events**: Progression to next line, treatment discontinuation

``` stata
mata {
    // Predict overall survival time from current age
    vOS_predicted = exp(XB_OS + errors)
    
    // Predict time to next event (TXD or TFI)
    vTNE_predicted = exp(XB_TNE + errors)
    
    // Competing risk: which happens first?
    vDiesFirst = (vOS_predicted :< vTNE_predicted)
    
    // Update mortality matrix
    mMOR[idx, OMC] = vDiesFirst
    
    // Update outcome time container
    mOC[idx, 1] = vDiesFirst :* vOS_predicted + 
                  (1 :- vDiesFirst) :* vTNE_predicted
    mOC[idx, 2] = vDiesFirst  // Mortality flag
}
```

### Terminal Mortality (OMC 19)

``` stata
// All patients still alive at L9E are forced to mortality
idxAlive = selectindex(mMOR[., OMC-1] :== 0)

// Cap OS at age limit (e.g., 100 years)
vExceedsLimit = ((mAge[idxAlive, 1] :+ mOS[idxAlive, OMC]) :> Limit)
mOS[idxAlive, OMC] = vExceedsLimit :* (Limit :- mAge[idxAlive, 1]) :+ 
                      (!vExceedsLimit) :* mOS[idxAlive, OMC]

// Mark all as dead
mMOR[idxAlive, OMC] = J(rows(idxAlive), 1, 1)
mOC[idxAlive, 1] = mOS[idxAlive, OMC]
mOC[idxAlive, 2] = J(rows(idxAlive), 1, 1)
```

------------------------------------------------------------------------

## Best Clinical Response (BCR) Prediction

### Ordered Categories (1-6)

```         
1 = Complete Response (CR)        - Best outcome
2 = Very Good Partial Response (VGPR)
3 = Partial Response (PR)
4 = Minimal Response (MR)
5 = Stable Disease (SD)
6 = Progressive Disease (PD)      - Worst outcome (reference)
```

### Multinomial Logit Implementation

``` stata
mata {
    // Build patient matrix
    pMat = (vAge[idx], vMale[idx], vECOG1[idx], vECOG2[idx], 
            vRISS2[idx], vRISS3[idx], vTXR_VRd[idx], vTXR_Other[idx], vCons[idx])
    
    // Extract coefficients for each BCR category (except reference PD)
    coefCols_CR = (1, 2, 3, 4, 5, 6, 7, 8, 9)      // Columns for CR
    coefCols_VGPR = (10, 11, 12, 13, 14, 15, 16, 17, 18) // Columns for VGPR
    // ... similar for PR, MR, SD
    
    // Calculate XB for each category
    XB_CR = pMat * bBCR[coefCols_CR]'
    XB_VGPR = pMat * bBCR[coefCols_VGPR]'
    // ... etc
    
    // Exponentiate
    exp_CR = exp(XB_CR)
    exp_VGPR = exp(XB_VGPR)
    // ... etc
    exp_PD = J(rows(idx), 1, 1)  // Reference category
    
    // Calculate denominator (normalisation)
    denom = exp_CR + exp_VGPR + exp_PR + exp_MR + exp_SD + exp_PD
    
    // Calculate probabilities
    pr_CR = exp_CR :/ denom
    pr_VGPR = exp_VGPR :/ denom
    // ... etc (probabilities sum to 1.0)
    
    // Cumulative probabilities for drawing
    cum_CR = pr_CR
    cum_VGPR = cum_CR + pr_VGPR
    cum_PR = cum_VGPR + pr_PR
    cum_MR = cum_PR + pr_MR
    cum_SD = cum_MR + pr_SD
    // cum_PD = 1.0 (all remaining probability)
    
    // Draw random uniform and map to category
    u = runiform(rows(idx), 1)
    vBCR_outcome = (u :< cum_CR) :* 1 +
                   ((u :>= cum_CR) :& (u :< cum_VGPR)) :* 2 +
                   ((u :>= cum_VGPR) :& (u :< cum_PR)) :* 3 +
                   ((u :>= cum_PR) :& (u :< cum_MR)) :* 4 +
                   ((u :>= cum_MR) :& (u :< cum_SD)) :* 5 +
                   (u :>= cum_SD) :* 6
    
    // Update BCR matrix
    mBCR[idx, OMC_to_Line_mapping] = vBCR_outcome
}
```

### BCR as Predictor

BCR at Line N becomes a **powerful predictor** for outcomes at Line N+1: - **Better Line 1 BCR** → Better Line 2 BCR, longer TFI, better OS - **CR/VGPR at L1** → 2-3× longer median OS vs PD at L1 - BCR interacts with Line number in OS equations (separate coefficients per line segment)

------------------------------------------------------------------------

## Parametric Survival Models

### Generalised Gamma Distribution

**Why Generalised Gamma?** - Flexible shape: can approximate exponential, Weibull, lognormal, log-logistic - Handles non-monotonic hazards (initially increasing, then decreasing) - Three-parameter control: μ (location), σ (scale), κ (shape)

**Implementation:**

``` stata
mata {
    // Extract coefficients
    XB = pMat * coefs[linearPredictorCols]'      // Linear predictor
    sigma = coefs[sigmaCol]                       // Scale parameter
    kappa = coefs[kappaCol]                       // Shape parameter
    
    // Generate error term
    z = rnormal(rows(idx), 1, 0, 1)              // Standard normal
    error = sigma * (kappa * z)                   // Scaled error
    
    // Bound error to prevent extreme values
    maxError = 3 * abs(sigma)
    minError = -3 * abs(sigma)
    error = bounded(error, maxError, minError)
    
    // Transform to time scale
    log_time = XB + error
    time = exp(log_time)
    
    // Apply clinical constraints
    time = bounded(time, maxTime, minTime)       // E.g., TXD: 1-60 months
    
    // Update outcome matrix
    mTXD[idx, Line] = time
}
```

### Manual Splines for ASCT Patients (L1 TXD)

**Problem:** ASCT patients have bi-phasic treatment duration: 1. **Induction phase**: 3-6 months of chemotherapy 2. **Post-ASCT phase**: Variable duration depending on response

**Solution:** Three separate survival equations

``` stata
// Equation 5: TXD from L1 start to ASCT
// Equation 6: TXD from ASCT to end of induction/consolidation  
// Equation 7: Total L1 TXD (combined model)

// In simulation:
if (vSCT_L1[i] == 1) {
    // Use ASCT-specific coefficients
    TXD_L1 = predict_with_spline(bL1_TXD_ASCT_S1, bL1_TXD_ASCT_S2, bL1_TXD_ASCT_S3)
} else {
    // Use non-ASCT coefficients
    TXD_L1 = predict_parametric(bL1_TXD)
}
```

------------------------------------------------------------------------

## Bootstrap Uncertainty Quantification

### Purpose

- Characterise **parameter uncertainty** in risk equations
- Generate **confidence intervals** for health economic outcomes (QALYs, costs, ICERs)
- Validate **model robustness** (out-of-sample prediction)

### Implementation

**Training Cohort Resampling:**

``` stata
// In estimation repository:
forvalues b = 1/100 {
    // Resample training cohort with replacement
    bsample
    
    // Re-estimate all 30 risk equations
    ologit BCR_L1 age male ECOG RISS TXR
    matrix coef_BCR_L1_B`b' = e(b)
    
    // ... repeat for all 30 equations
    
    // Save coefficients
    mata: mata matsave "coefficients_B`b'.mmat", replace
}
```

**Simulation with Bootstrap Samples:**

``` stata
// In simulation repository:
forvalues b = 1/$max_bs {
    // Load bootstrap coefficient set
    mata: mata matuse "coefficients_B`b'.mmat"
    
    // Run full simulation with this parameter set
    load_patients
    mata_setup
    simulation
    process_data
    
    // Save results
    save "simulated_B`b'.dta", replace
}
```

**Aggregation:**

``` stata
// Combine bootstrap iterations
use "simulated_B1.dta", clear
forvalues b = 2/100 {
    append using "simulated_B`b'.dta"
}

// Calculate 95% confidence intervals
collapse (p50) median_OS = OS (p2.5) lower_OS = OS (p97.5) upper_OS = OS, by(time)
```

------------------------------------------------------------------------

## Validation Strategy

### Out-of-Sample Prediction (70/30 Split)

**Training Cohort (70%):** - Multiple imputation of missing values - Estimate all 30 risk equations - Bootstrap 100 times for uncertainty

**Validation Cohort (30%):** - Independent holdout sample - Diagnostic characteristics provided to simulation - No parameter estimation performed

**Comparison:** - Kaplan-Meier curves: Validation vs Simulated - 95% confidence intervals generated via bootstrap - Monthly p-value tests (H₀: no difference between curves) - Result: 90% of 120 months showed no significant difference (p \> 0.05)

### BCR as Surrogate for OS

**Test:** Do patients with better BCR have better OS?

``` stata
// Stratify by BCR at Line 1
stset OS, failure(death)

sts graph, by(BCR_L1)
// Expected: Monotonically decreasing survival curves from CR → PD

sts test BCR_L1
// Result: Highly significant (p < 0.001)
```

**Finding:** BCR is a valid surrogate for OS - Patients with CR/VGPR have 2-3× longer median survival than PD - Justifies using BCR in cost-effectiveness analysis (intermediate outcome)

------------------------------------------------------------------------

## Command-Line Interface

### Main Entry Point: `main.do`

``` stata
// Syntax:
do "main.do" [analysis] [intervention] [line] [coeffs] [data] ///
             [min_id] [max_id] [bootstrap] [min_bs] [max_bs] [report]
```

### Parameter Specifications

| Arg | Name | Description | Examples |
|---------------|---------------|------------------------|------------------|
| 1 | **analysis** | Analysis configuration folder | `base_model`, `vrd_l1_post`, `dvd_l2_method` |
| 2 | **intervention** | Treatment to simulate | `VRd`, `DVd`, `SoC`, `all` |
| 3 | **line** | Line of therapy (0=all) | `0`, `1`, `2`, `3`, `4` |
| 4 | **coeffs** | Coefficient set | `base_model`, `VRd`, `SoC` |
| 5 | **data** | Patient cohort type | `population`, `predicted` |
| 6 | **min_id** | Minimum patient ID | `1` |
| 7 | **max_id** | Maximum patient ID | `10`, `1000`, `4884` |
| 8 | **bootstrap** | Bootstrap flag | `0` (no), `1` (yes) |
| 9 | **min_bs** | Min bootstrap iteration | `1` |
| 10 | **max_bs** | Max bootstrap iteration | `100`, `500` |
| 11 | **report** | Generate report | `0` (no), `1` (yes) |

### Example Commands

``` stata
// Quick test: 10 patients, no bootstrap
do "main.do" base_model all 0 base_model population 1 10 0 0 0 0

// Full simulation: all 4,884 patients with report
do "main.do" base_model all 0 base_model population 1 4884 0 0 0 1

// Bootstrap analysis: 1000 patients, 100 iterations
do "main.do" base_model VRd 1 base_model predicted 1 1000 1 1 100 1

// Population projection: Line 2 specific
do "main.do" base_model all 2 base_model population 1 4884 0 0 0 1
```

------------------------------------------------------------------------

## Analysis Configurations

### `base_model/`

**Purpose:** Foundation model with all standard regimens - **Regimens:** VRd, VCd, Other (L1); Rd, DVd, Other (L2+) - **Coefficients:** Estimated from full MRDR training cohort - **Use case:** Baseline comparisons, validation

### `vrd_l1_post/`

**Purpose:** Post-market analysis of VRd at Line 1 - **Focus:** VRd-specific outcomes after PBS listing - **Data:** Real-world VRd patients (2014-2023) - **Use case:** Post-market surveillance, effectiveness evaluation

### `dvd_l2_method/`

**Purpose:** DVd transportability methodology at Line 2 - **Focus:** CASTOR trial → MRDR population bridging - **Methods:** Calibrated Transportation with common comparator (Vd) - **Use case:** Pre-market HTA submission, reimbursement decision

------------------------------------------------------------------------

## Health Economic Integration

### Utility Mapping (QALYs)

``` stata
// Map BCR and age to EQ-5D utilities
// Formula: u = baseline - age_decrement + BCR_benefit - treatment_disutility

utility_CR = 0.85 - 0.005*(age-60) + 0.00 - 0.02  // No CR benefit (reference)
utility_VGPR = 0.85 - 0.005*(age-60) - 0.03 - 0.02
utility_PR = 0.85 - 0.005*(age-60) - 0.06 - 0.02
utility_SD = 0.85 - 0.005*(age-60) - 0.12 - 0.02
utility_PD = 0.85 - 0.005*(age-60) - 0.18 - 0.02

// Treatment disutility:
// - On chemotherapy: -0.02
// - Off treatment (TFI): +0.02
```

### Cost Attribution

``` stata
// Per-patient costs by state:
cost_L1_chemo = acquisition_cost + admin_cost + monitoring_cost
cost_L1_ASCT = inpatient_cost + stem_cell_cost + complication_cost
cost_L1_maintenance = acquisition_cost + monitoring_cost
cost_TFI = monitoring_cost  // Minimal

// Aggregate:
total_cost = Σ(state_duration × state_cost) across all lines
```

### Threshold Analysis

``` stata
// Find break-even price where ICER = willingness-to-pay threshold
// E.g., For DVd vs SoC at L4:

ICER(price) = [Cost_DVd(price) - Cost_SoC] / [QALY_DVd - QALY_SoC]

// Solve: ICER(price*) = WTP_threshold
// Implemented via bootstrap confidence intervals
```

------------------------------------------------------------------------

## Common Pitfalls and Debugging

### 1. Dimension Mismatch Errors

**Symptom:**

```         
conformability error
    pMat * coefs'
```

**Cause:** Patient matrix columns ≠ coefficient vector columns

**Solution:**

``` stata
// Always use validateDimensions() before multiplication
mata: validateDimensions(pMat, coefs, "BCR_L1")
```

### 2. Missing Coefficient Extraction

**Symptom:** Some patients get missing outcomes (`.`)

**Cause:** Coefficient column indices don't match matrix structure

**Solution:**

``` stata
// Verify coefficient extraction:
mata: bBCR[coefCols]  // Display extracted coefficients
mata: cols(pMat)      // Should equal cols(bBCR[coefCols])
```

### 3. Probability Sums ≠ 1.0 (BCR)

**Symptom:** Invalid BCR categories or probabilities

**Cause:** Floating point error in cumulative probability calculation

**Solution:**

``` stata
// Force normalisation:
denom = exp_CR + exp_VGPR + ... + exp_PD
pr_CR = exp_CR :/ denom  // Guarantees sum = 1.0
```

### 4. Extreme Survival Predictions

**Symptom:** OS predictions of 500+ years

**Cause:** Unbounded error terms in parametric models

**Solution:**

``` stata
// Always bound error terms:
errors = bounded(rnormal(...), maxError, minError)

// And bound final predictions:
vOS = bounded(exp(XB + errors), maxOS, 0.1)  // E.g., max 30 years
```

### 5. Mortality Not Updating

**Symptom:** Patients "stuck" alive beyond predicted OS

**Cause:** Mortality logic not checking competing risks correctly

**Solution:**

``` stata
// Ensure competing risk check at every OMC:
vDiesFirst = (mOS[idx, OMC] :< mTNE[idx, OMC])
mMOR[idx, OMC] = vDiesFirst
```

------------------------------------------------------------------------

## Performance Optimisation Tips

### 1. Vectorise Everything

**Bad:**

``` stata
forvalues i = 1/$N {
    XB = b1*x1[`i'] + b2*x2[`i'] + ...
    outcome[`i'] = exp(XB)
}
```

**Good:**

``` stata
mata {
    XB = pMat * coefs'          // Single matrix multiplication
    outcomes = exp(XB)           // Vectorised exponentiation
    mOutcome[idx, OMC] = outcomes
}
```

### 2. Pre-filter Indices

**Bad:**

``` stata
mata {
    for (i = 1; i <= rows(mMOR); i++) {
        if (mMOR[i, OMC-1] == 0 & mState[i, 1] <= OMC) {
            // Process patient i
        }
    }
}
```

**Good:**

``` stata
mata {
    idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
    // Now operate on pMat[idx, ], vAge[idx], etc.
}
```

### 3. Avoid Unnecessary Copying

**Bad:**

``` stata
mata {
    tempAge = mAge[., 1]           // Creates copy
    tempAge2 = tempAge :^ 2        // Creates another copy
}
```

**Good:**

``` stata
mata {
    vAge2 = vAge :^ 2              // Direct operation on view
}
```

### 4. Use Mata for All Computations

**Bad:**

``` stata
// Switching between Stata and Mata repeatedly
foreach var in Age Male ECOG {
    mata: v`var' = st_data(., "`var'")
}
```

**Good:**

``` stata
// Load all at once in mata_setup.do
mata {
    vAge = st_data(., "Age_DN")
    vMale = st_data(., "Male")
    vECOG = st_data(., "ECOGcc")
    // ... all variables in one block
}
```

------------------------------------------------------------------------

## Future Development Roadmap

### Planned Enhancements

1.  **Infection Risk Sub-Module**
    - Discrete monthly time steps for infection probability
    - Grade 3-4 infection costs and disutilities
    - Risk factors: Age, neutropenia, regimen intensity
2.  **Progression-Free Survival Integration**
    - Alternative endpoint to BCR
    - Requires PFS data extraction from MRDR
    - Enable comparison: BCR-based vs PFS-based predictions
3.  **R/Python Translation**
    - Broaden accessibility beyond Stata users
    - Leverage modern ML libraries (survival analysis packages)
    - Improved visualisation capabilities
4.  **Smoldering Myeloma Extension**
    - Add pre-diagnosis states (SMM risk stratification)
    - CTX-1 biomarker-guided intervention
    - Compare watch-and-wait vs early treatment strategies
5.  **Performance Profiling**
    - Identify remaining bottlenecks
    - Explore GPU acceleration for bootstrap analyses
    - Target: \<1 second for 10,000 patients

------------------------------------------------------------------------

## Key Insights for Users

### When to Use the Model

**Ideal Use Cases:** - Pre-market HTA submissions (predict RWE from trials) - Post-market surveillance (validate trial predictions) - Cost-effectiveness analysis (QALY estimation) - Budget impact modelling (population projections) - Comparative effectiveness research

**Not Suitable For:** - Individual patient risk prediction (use ISS/R-ISS instead) - Real-time clinical decision support (too computationally intensive) - Non-myeloma cancers (disease-specific model)

### Interpreting Outputs

**Key Variables:** - `OS`: Overall survival from **diagnosis** (not from current line) - `TXD_LN`: Duration of treatment at line N (months) - `TFI_LN`: Treatment-free interval **after** line N (months off treatment) - `BCR_LN`: Best response **during** line N (1=CR to 6=PD)

**Important Caveats:** - Bootstrap CIs reflect **parameter uncertainty**, not **patient heterogeneity** - Model assumes MRDR population characteristics (Australian/NZ) - Extrapolation beyond 10 years increasingly uncertain - Novel regimens not in MRDR require external efficacy assumptions

------------------------------------------------------------------------

## References and Resources

### Key Publications

1.  **Irving et al. (2024)** - "Discrete-event simulation modelling of multiple myeloma using a national clinical registry" (PLOS ONE)
    - DOI: 10.1371/journal.pone.0308812
    - Describes model development, validation, and BCR-OS relationship
2.  **Dahabreh et al. (2019)** - "Extending inferences from a randomized trial to a target population"
    - Framework for transportability and generalisation
3.  **Cole & Stuart (2010)** - "Generalizing evidence from randomized clinical trials to target populations"
    - Foundation for population adjustment methods

### GitHub Repository

- **URL:** `https://github.com/CHE-Monash/EpiMAP-Myeloma`
- **License:** GPL-3.0
- **Version:** 2.1 (October 2025)

### Data Access

- **MRDR:** Not publicly available (patient confidentiality)
- **De-identified data:** Available via application to MRDR Steering Committee
- **Contact:** <https://www.mrdr.net.au/>

------------------------------------------------------------------------

## Conclusion

EpiMAP Myeloma v2.1 represents a **state-of-the-art discrete-event simulation** for multiple myeloma disease modelling. The vectorisation transformation has delivered unprecedented performance whilst maintaining clinical validity. Its 30 evidence-based risk equations, comprehensive treatment pathway modelling, and robust validation make it a powerful tool for health technology assessment, comparative effectiveness research, and health economic evaluation.

The modular architecture enables extension to novel therapies, alternative endpoints, and related disease states, positioning EpiMAP as a flexible platform for future myeloma outcomes research.
