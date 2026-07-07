# Monash Myeloma Model Testing Suite
## Comprehensive Test Specifications for Model Logic Validation

**Version 2.0**  
**Date:** 27 October 2025

---

## Executive Summary

This document outlines a comprehensive testing framework for the Monash Myeloma Model discrete event simulation model. The testing suite is designed to validate model logic across all 50 risk equations, ensure data integrity throughout the simulation, and verify the accuracy of treatment pathway modelling and clinical outcome predictions.

The test suite is organised into seven major categories covering:
1. Input validation
2. Risk equation accuracy
3. Treatment pathway logic
4. Survival calculations
5. Edge case handling
6. Performance validation
7. Output verification

Each test includes specific objectives, expected outcomes, and implementation guidance for Stata do-files.

---

## 1. Input Data Validation Tests

### Overview
These tests ensure that patient input data meets all requirements before simulation begins and that patient characteristics are correctly transformed into model variables.

---

### Test 1.1: Patient Characteristic Range Validation

**Objective:** Verify all patient characteristics fall within valid ranges and follow expected distributions.

**Test Cases:**
- **Age:** 
  - Verify all patients aged 18-100 years
  - Check distribution matches MRDR registry patterns
  - Flag any suspicious outliers (e.g., >100 years or <18 years)
  
- **Sex (Male indicator):**
  - Binary indicator (0/1) only
  - No missing values permitted
  - Check sex distribution approximates 55-60% male
  
- **ECOG Performance Score:**
  - Values must be 0, 1, or 2 only
  - Verify proportions are reasonable (ECOG 0: ~40%, ECOG 1: ~35%, ECOG 2: ~25%)
  - Flag if any invalid values (negative, >2, or missing)
  
- **R-ISS Stage:**
  - Values must be 1, 2, or 3 only
  - Check stage distribution matches expected pattern
  - R-ISS I: ~25-30%, R-ISS II: ~50-55%, R-ISS III: ~15-20%
  
- **Comorbidity Score (CMScore):**
  - Values must be 0-6 (counting: cardiac, pulmonary, diabetes, liver, neuropathy, other malignancy)
  - Verify distribution reasonable (mode typically 0-1)
  - Flag if >50% have CMScore ≥3 (suggests data quality issue)

**Expected Outcome:** 
- 100% of patients pass range checks
- Distributions match MRDR patterns within acceptable tolerance (±10%)
- Clear error messages for any out-of-range values

**Implementation Notes:**
- Create `test_input_validation.do`
- Use `assert` commands for hard constraints
- Generate summary statistics table comparing input data to MRDR benchmarks
- Flag outliers using `kdensity` plots and quantile checks

---

### Test 1.2: Dummy Variable Generation Validation

**Objective:** Validate correct creation of dummy variables for ECOG, R-ISS, and comorbidity categories, ensuring proper vectorisation matches matrix operations.

**Test Cases:**

**ECOG Dummy Variables:**
- Verify `vECOG0`, `vECOG1`, `vECOG2` sum to exactly 1 for each patient
- Check `vECOG0 = (vECOG == 0)`, `vECOG1 = (vECOG == 1)`, `vECOG2 = (vECOG == 2)`
- Cross-validate against `mCore` matrix ECOG columns
- No patients should have all three indicators = 0

**R-ISS Dummy Variables:**
- Verify `vRISS1`, `vRISS2`, `vRISS3` sum to exactly 1 for each patient
- Check `vRISS1 = (vRISS == 1)`, `vRISS2 = (vRISS == 2)`, `vRISS3 = (vRISS == 3)`
- Cross-validate against `mCore` matrix R-ISS columns
- Ensure mutually exclusive (only one =1 per patient)

**Comorbidity Dummy Variables:**
- Verify correct categorisation: `vCMc0 = (CMScore == 0)`, `vCMc1 = (CMScore == 1)`, `vCMc2 = (CMScore == 2)`, `vCMc3 = (CMScore >= 3)`
- Check that exactly one comorbidity dummy = 1 per patient
- Cross-validate against `mCom` matrix

**Age Transformation Variables:**
- `vAge`: Verify matches original age variable
- `vAge2`: Check equals `vAge^2` for all patients
- `vAge70`: Verify equals 1 if `vAge >= 70`, else 0
- `vAge75`: Verify equals 1 if `vAge >= 75`, else 0
- Check against `mCom` matrix columns 1-2

**Expected Outcome:** 
- All dummy variables correctly reflect original patient characteristics
- No summation errors (each set of dummies sums to 1.0)
- Vector operations produce identical results to matrix operations
- Zero tolerance for errors (must be exact match)

**Implementation Notes:**
- Adapt existing `tests/validate_vectors.do` 
- Run on multiple patient cohorts (n=100, n=1000, n=4884)
- Compare scalar vs matrix calculations element-by-element
- Use `assert` statements with `_rc` checks for each validation

---

### Test 1.3: Missing Data Handling

**Objective:** Verify model handles missing data appropriately through multiple imputation or exclusion, following documented procedures.

**Test Cases:**

**Critical Variables (Cannot be Missing):**
- Age, sex, ECOG, R-ISS, CMScore at baseline
- Test dataset with deliberate missing values in each critical variable
- Verify model produces clear error messages identifying which variable is missing
- Confirm simulation does not proceed with missing critical data

**Multiple Imputation Validation:**
- If imputation implemented, verify imputed values fall within valid ranges
- Check that imputation procedure uses appropriate covariates
- Verify imputed datasets have similar distributions to complete-case data
- Test that simulation results are averaged across multiple imputed datasets

**Listwise Deletion:**
- If using complete-case analysis, verify proper exclusion count reported
- Check that patients with any missing critical variables excluded before simulation
- Confirm excluded patients logged with reason for exclusion

**Pattern Detection:**
- Check for systematic patterns in missing data (e.g., older patients missing R-ISS more often)
- Flag if >20% of any critical variable is missing in original data

**Expected Outcome:** 
- Model either properly imputes or appropriately excludes patients with critical missing data
- Clear documentation of handling approach in log files
- No simulation runs with unhandled missing values
- Imputed values (if used) are clinically plausible

**Implementation Notes:**
- Create `test_missing_data.do` with synthetic datasets containing strategic missingness
- Test each critical variable individually and in combination
- Verify error handling produces informative messages
- If using MI, validate against Stata's `mi estimate` output

---

### Test 1.4: Constant Vector Validation

**Objective:** Verify the constant vector (`vCons`) is correctly initialised as all ones for use in risk equations.

**Test Cases:**
- Check `vCons` has same length as number of patients
- Verify every element equals 1.0 exactly
- Confirm `vCons` remains unchanged throughout simulation
- Test that `vCons` properly feeds into coefficient calculations

**Expected Outcome:** 
- `vCons` vector contains only 1.0 values
- No variation across patients
- Consistent across all simulation time points

**Implementation Notes:**
- Quick validation in `mata_setup.do`
- Use `assert min(vCons) == 1 & max(vCons) == 1`

---

## 2. Risk Equation Accuracy Tests

### Overview
These tests validate the 50 risk equations that form the core of the simulation model, ensuring coefficients are correctly applied and predictions match expected values.

---

### Test 2.1: Coefficient Matrix Loading and Integrity

**Objective:** Verify all coefficient matrices load correctly and contain expected values for each analysis scenario.

**Test Cases:**

**Matrix Loading:**
- Verify coefficients load from correct `.mmat` files for each analysis (VRd-Post, DVd-Pre, etc.)
- Check matrix dimensions match expected rows and columns
- Confirm no missing values (`.`) in coefficient matrices
- Test that coefficient labels match variable names in documentation

**Value Validation:**
- Compare loaded coefficients against published values in supplementary materials
- Verify coefficient signs (positive/negative) match clinical expectations
- Check magnitude of coefficients reasonable (no extreme outliers like >50 or <-50)

**Analysis-Specific Coefficients:**
- **VRd analysis:** Verify VRd-specific coefficients load for Line 1 regimen
- **DVd analysis:** Verify DVd-specific coefficients load for Line 2 regimen  
- **Standard of Care:** Verify SoC coefficients reflect pre-VRd era patterns

**Coefficient Consistency:**
- Test that different analysis scenarios load different coefficients appropriately
- Verify Line 1 vs Line 2+ coefficients differ as expected
- Check ASCT vs non-ASCT pathway coefficients where relevant

**Expected Outcome:** 
- All 50 risk equation coefficient sets load without errors
- Coefficients match documented parameter values exactly
- Analysis-specific modifications apply correctly
- No warning messages during matrix operations

**Implementation Notes:**
- Create `test_coefficient_loading.do`
- Use `mata: mata describe` to inspect loaded matrices
- Compare key coefficients to published tables using tolerance of 0.0001
- Test all analysis scenarios: VRd-Post, DVd-Pre, DVd-Post, etc.

---

### Test 2.2: Overall Survival (OS) Risk Equation Validation

**Objective:** Validate the primary overall survival risk equation produces accurate and clinically sensible survival predictions across patient segments.

**Test Cases:**

**Manual Calculation Verification:**
- Create test patients with known characteristics (e.g., Age=60, Male=1, ECOG=0, R-ISS=I, BCR=CR)
- Manually calculate expected survival time using risk equation formula
- Compare model predictions to manual calculations (tolerance <0.01%)
- Test at multiple time segments: Diagnosis, Line 1 Start, Line 1 End, Line 2+

**Prognostic Factor Monotonicity:**
- Better prognostic factors should yield longer survival:
  - Younger age → longer survival
  - Better ECOG (0 > 1 > 2) → longer survival
  - Better R-ISS (I > II > III) → longer survival
  - Better BCR (CR > VGPR > PR > MR > SD > PD) → longer survival
- Test on cohort of 1000 patients to verify trend direction

**Edge Case Testing:**
- Very young patient (age 30) with best prognostic factors → very long survival
- Very old patient (age 90) with worst prognostic factors → short survival
- Test survival time is always positive
- Verify no survival predictions >30 years (unrealistic for MM)

**Segment-Specific Testing:**
- **Segment 0 (Diagnosis/DN):** Baseline survival without treatment
- **Segment 1 (Line 1 No ASCT):** Survival during/after induction without ASCT
- **Segment 2 (Line 1 ASCT):** Survival during/after induction with ASCT
- **Segments 3-7 (Lines 2-6):** Survival at later lines

**BCR Impact Validation:**
- Systematically vary BCR (1-6) holding other factors constant
- Verify survival increases monotonically with better BCR
- Check that BCR effect size matches clinical expectations (CR roughly doubles survival vs PD)

**Expected Outcome:** 
- Model survival predictions match hand-calculated values within 0.01%
- Monotonicity preserved for all prognostic factors
- Edge cases produce clinically plausible results
- Segment-specific coefficients apply correctly based on patient pathway
- BCR-survival relationship robust and clinically valid

**Implementation Notes:**
- Create `test_os_equation.do`
- Use `calcSurvTime()` function directly with test parameters
- Generate Kaplan-Meier style survival curves for visual validation
- Compare to published MRDR survival curves by subgroup
- Test with fixed random seeds for reproducibility

---

### Test 2.3: ASCT Eligibility Prediction Validation

**Objective:** Validate planned ASCT eligibility predictions align with clinical criteria and observed MRDR patterns.

**Test Cases:**

**Age-Based Eligibility:**
- Patients <65 years: High ASCT eligibility (>60% planned ASCT)
- Patients 65-70 years: Moderate eligibility (30-50% planned ASCT)
- Patients 70-75 years: Low eligibility (10-20% planned ASCT)
- Patients ≥75 years: Very low eligibility (<5% planned ASCT)

**ECOG-Based Eligibility:**
- ECOG 0: Highest ASCT probability
- ECOG 1: Moderate reduction in ASCT probability
- ECOG 2: Substantial reduction in ASCT probability
- Within age groups, verify ECOG 0 > ECOG 1 > ECOG 2

**Comorbidity Impact:**
- CMScore 0: Baseline ASCT probability
- CMScore 1-2: Modest reduction in probability
- CMScore 3+: Substantial reduction in probability
- Test that comorbidity coefficient has expected sign (negative)

**Combined Risk Profiles:**
- Best case (Age 55, ECOG 0, CMScore 0, R-ISS I): ~80-90% ASCT planned
- Worst case (Age 75, ECOG 2, CMScore 3+, R-ISS III): <5% ASCT planned
- Average case: Should match MRDR observed rate (~40-50%)

**Population-Level Validation:**
- Aggregate ASCT planning rate across full simulated cohort
- Compare to MRDR observed rate (target: within ±5%)
- Verify distribution by age group matches registry patterns

**Logit Mechanics:**
- Test that predicted probabilities all fall in [0,1] range
- Verify logit transformation applied correctly: p = exp(XB) / (1 + exp(XB))
- Check that extreme XB values don't cause numerical overflow

**Expected Outcome:** 
- ASCT eligibility predictions clinically sensible across all patient profiles
- Age is strongest predictor (as in clinical practice)
- Population ASCT rate matches MRDR within ±5%
- No patients have predicted probability <0 or >1
- Older/sicker patients consistently have lower ASCT probability

**Implementation Notes:**
- Create `test_asct_eligibility.do`
- Generate cross-tabulations of ASCT probability by age/ECOG/R-ISS
- Compare distributions to MRDR using chi-square tests
- Plot predicted probabilities across age continuum

---

### Test 2.4: Best Clinical Response (BCR) Prediction Validation

**Objective:** Validate multinomial logit models for BCR produce appropriate response distributions across lines of therapy and regimens.

**Test Cases:**

**Probability Sum Validation:**
- For every patient at every line: Sum of BCR probabilities (CR, VGPR, PR, MR, SD, PD) must equal 1.0
- Tolerance: <0.0001 deviation from 1.0
- Test across all 1000 test patients at Lines 1, 2, and 3+

**Line 1 BCR by Regimen:**
- **VRd:** Should produce superior response distribution (more CR/VGPR) compared to VCd and Other
- **VCd:** Intermediate response rates
- **Other:** Generally poorer response rates
- Verify VRd has statistically higher ORR (Overall Response Rate = CR+VGPR+PR) than other regimens

**Prognostic Factor Impact:**
- Better baseline characteristics should predict better BCR:
  - Younger age → higher probability of CR/VGPR
  - Better ECOG → higher probability of deep response
  - Better R-ISS → higher probability of CR/VGPR
- Test monotonicity on cohorts stratified by single factor

**ASCT BCR Validation:**
- ASCT BCR should show improvement over induction therapy BCR
- Patients with good induction BCR (CR/VGPR) should maintain or improve with ASCT
- Expected distribution: 60-70% CR/VGPR post-ASCT
- ASCT BCR should be superior to non-ASCT Line 1 end BCR

**Line 2+ BCR Patterns:**
- BCR at Line 2 should depend on Line 1 BCR (prior response predictor)
- Patients with better Line 1 BCR should have better Line 2 BCR probability
- DVd at Line 2 should show superior BCR compared to other Line 2 regimens

**BCR Distribution Benchmarking:**
- Compare aggregate BCR distribution to MRDR observed rates by regimen and line
- Target: Within ±10% for each BCR category
- Generate contingency tables for visual comparison

**Multinomial Logit Mechanics:**
- Reference category correctly set (typically Progressive Disease = reference)
- XB calculations use correct coefficient columns for each BCR outcome
- Exponentiation and normalisation produces valid probability distribution

**Expected Outcome:** 
- All BCR probability vectors sum to 1.0 (no arithmetic errors)
- VRd produces superior response distribution (statistically significant)
- ASCT improves responses appropriately (clinical benefit evident)
- Prognostic factors influence BCR in expected directions
- Population-level BCR distributions match MRDR patterns within ±10%

**Implementation Notes:**
- Create `test_bcr_predictions.do`
- Use ordered logit validation techniques
- Generate cross-tabs of BCR by regimen, line, and patient characteristics
- Perform chi-square tests comparing simulated vs observed BCR distributions
- Visual validation with stacked bar charts by regimen

---

### Test 2.5: Treatment Duration Prediction Validation

**Objective:** Validate parametric survival models for treatment duration produce clinically realistic durations across all lines of therapy.

**Test Cases:**

**Positive Duration Requirement:**
- All predicted treatment durations must be >0 months
- No negative durations permitted
- No missing duration values

**Mean Duration by Regimen and Line:**
- Compare mean simulated duration to MRDR observed mean by regimen and line
- Target: Within ±10% of observed mean
- Test for Line 1 regimens: VRd, VCd, Other
- Test for Line 2 regimens: Rd, DVd, Other

**BCR-Duration Correlation:**
- Better BCR should correlate with longer treatment duration (patients stay on effective therapy)
- CR/VGPR → Longer duration than PR
- PR → Longer than MR/SD/PD
- Calculate Spearman correlation coefficient (expect ρ > 0.3)

**ASCT Pathway Duration:**
- Induction therapy before ASCT typically 3-6 months
- Verify mean induction duration for ASCT patients falls in this range
- ASCT patients should have distinct duration pattern (reflecting pre-transplant protocol)

**Restricted Cubic Splines Validation:**
- For ASCT patients: Three splines used with cutoffs at 3 and 6 months
- Verify smooth duration distribution (no unrealistic spikes at cutoffs)
- Test hazard rate changes are gradual across spline boundaries
- Check that patients surviving past Cutoff 1 use Spline 2 coefficients, past Cutoff 2 use Spline 3

**Duration by Line of Therapy:**
- Line 1 duration typically longest (6-12 months)
- Line 2 duration moderate (4-8 months)  
- Line 3+ duration generally shorter (3-6 months)
- Verify declining trend in mean duration with later lines

**Extreme Value Testing:**
- Flag any durations >36 months (very unusual, suggests model issue)
- Flag any durations <0.5 months (unlikely except for immediate PD)
- Verify <1% of patients have extreme values

**Expected Outcome:** 
- All durations positive and clinically realistic
- Mean durations match registry within ±10%
- BCR-duration relationship positive and significant
- ASCT induction durations fall in expected 3-6 month window
- Spline models produce smooth distributions
- Duration decreases with later lines of therapy

**Implementation Notes:**
- Create `test_treatment_duration.do`
- Compare distributions using kernel density plots
- Use `calcSurvTime()` function with various patient profiles
- Validate spline transitions with boundary case testing
- Generate summary statistics by regimen, line, and BCR

---

### Test 2.6: Treatment-Free Interval (TFI) Validation

**Objective:** Validate TFI predictions between lines of therapy are clinically appropriate and reflect treatment response patterns.

**Test Cases:**

**Non-Negativity Requirement:**
- All TFI values must be ≥0 months
- No negative intervals permitted (would violate time sequence)

**BCR-TFI Relationship:**
- Better BCR on previous line should predict longer TFI
- Test across BCR categories: CR > VGPR > PR > MR > SD/PD
- Calculate correlation coefficient (expect ρ > 0.4)

**ASCT Advantage:**
- ASCT patients should have longer Line 1→Line 2 TFI than non-ASCT patients
- Expected difference: 6-12 months longer TFI for ASCT patients
- Test statistical significance using t-test or Wilcoxon rank-sum

**Maintenance Therapy Effect:**
- Patients receiving maintenance therapy should have extended TFI
- Expected increase: 3-6 months longer TFI with maintenance vs without
- Verify maintenance therapy coefficient has positive sign

**Mean TFI by Line and Pathway:**
- Compare mean TFI to MRDR observed intervals
- **Line 1→2 (No ASCT):** Target ~12-18 months
- **Line 1→2 (ASCT):** Target ~24-36 months
- **Line 2→3:** Target ~8-12 months
- **Line 3→4+:** Target ~6-9 months
- Tolerance: ±15% of observed means

**TFI Distribution Shape:**
- TFI should follow roughly log-normal or Weibull distribution
- Long right tail expected (some patients have very long remissions)
- No unrealistic bimodal distributions

**Parametric Model Validation:**
- Verify selected distribution (Weibull, exponential, Gompertz, etc.) fits data well
- Check AIC/BIC to confirm best-fitting distribution chosen
- Visual inspection of survival curves vs Kaplan-Meier

**Expected Outcome:** 
- All TFI values non-negative
- Better prior response yields significantly longer remissions
- ASCT advantage clearly evident (6-12 month longer TFI)
- Maintenance therapy extends TFI appropriately
- Mean TFI by pathway matches MRDR within ±15%
- Distribution shapes clinically plausible

**Implementation Notes:**
- Create `test_tfi_validation.do`
- Stratify analyses by ASCT status and BCR
- Use `stcox` or `streg` to validate TFI models
- Generate Kaplan-Meier plots by subgroup
- Compare median TFI (more robust than mean for skewed distributions)

---

### Test 2.7: Regimen Selection Validation

**Objective:** Validate multinomial logit models for treatment regimen selection produce appropriate regimen distributions.

**Test Cases:**

**Line 1 Regimen Distribution:**
- Expected proportions (pre-VRd era): VCd ~58%, VRd ~15%, Other ~27%
- Expected proportions (post-VRd era): VRd ~60%, VCd ~25%, Other ~15%
- Verify proportions sum to 100% ±1%

**Line 2 Regimen Distribution:**
- Expected proportions: Rd ~16%, DVd ~11% (post-DVd approval ~40%), Other ~73% (or ~49% post-DVd)
- Verify analysis scenario (DVd-Pre vs DVd-Post) changes distribution appropriately

**BCR Influence on Next Regimen:**
- Better previous BCR should slightly increase probability of more intensive regimens
- Test that CR/VGPR patients have higher probability of triplet therapy at next line
- Correlation expected to be modest (ρ ~ 0.1-0.2)

**Patient Characteristics Impact:**
- Younger patients: More likely to receive intensive regimens
- Better ECOG: Higher probability of complex triplet therapies
- Older/frailer patients: More likely to receive doublet or less intensive regimens

**Regimen Transition Matrices:**
- No patients should receive same regimen twice in sequence (cross-resistance assumed)
- Verify that Line 1 regimen influences Line 2 regimen probabilities
- Expected pattern: VRd→DVd more common than VRd→Rd (in DVd-Post era)

**Expected Outcome:**
- Regimen distributions match MRDR observed patterns within ±10%
- Analysis scenario correctly modifies regimen probabilities
- BCR influences next regimen selection in expected direction
- No same-regimen repeats occur
- Patient characteristics influence regimen choice appropriately

**Implementation Notes:**
- Create `test_regimen_selection.do`
- Generate contingency tables of regimen by line, BCR, age, ECOG
- Use chi-square tests to compare simulated vs observed distributions
- Validate multinomial logit mechanics (probabilities sum to 1.0)

---

### Test 2.8: Maintenance Therapy Receipt Validation

**Objective:** Validate logit model for maintenance therapy receipt produces clinically appropriate proportions.

**Test Cases:**

**Overall Receipt Rate:**
- Expected proportion receiving maintenance: ~40-60% of Line 1 patients
- Compare simulated rate to MRDR observed rate (±10%)

**ASCT Patient Maintenance:**
- Higher maintenance rate among ASCT patients (~70-80%)
- Lower rate among non-ASCT patients (~30-40%)

**BCR Influence:**
- Better BCR to induction therapy increases maintenance probability
- CR/VGPR patients more likely to receive maintenance than PR/MR patients

**Age and Comorbidity Impact:**
- Younger patients: Higher maintenance probability
- Fewer comorbidities: Higher maintenance probability
- Older/frailer patients: Lower maintenance probability (tolerability concerns)

**Expected Outcome:**
- Maintenance receipt rate matches MRDR within ±10%
- ASCT patients have substantially higher maintenance rate
- Better BCR increases maintenance probability
- Patient characteristics influence maintenance in expected directions

**Implementation Notes:**
- Create `test_maintenance_receipt.do`
- Stratify by ASCT status, BCR, age, and comorbidity score
- Use logistic regression diagnostics to validate model

---

## 3. Treatment Pathway Logic Tests

### Overview
These tests verify the simulation correctly implements treatment sequencing, including ASCT pathways, maintenance therapy, and progression through lines of therapy.

---

### Test 3.1: Line of Therapy Progression Sequencing

**Objective:** Verify patients progress through lines of therapy in correct temporal and logical sequence.

**Test Cases:**

**Starting Point:**
- All patients must start at Line 1
- No patients begin simulation at Line 2 or later
- Verify `mState[, 1] == 1` for all patients at simulation start

**Sequential Progression:**
- Patients can only progress Line 1→2→3→4, etc.
- No line skipping (e.g., Line 1 directly to Line 3)
- Verify `mState[, t+1] == mState[, t]` or `mState[, t+1] == mState[, t] + 1`

**Progression Timing:**
- Progression to next line only occurs after: Treatment duration + Treatment-free interval
- Verify: Time at Line N start = Time at Line N-1 start + Duration N-1 + TFI N-1→N
- No temporal overlaps (patient cannot be on two lines simultaneously)

**Death Censoring:**
- Patients who die do not progress to additional lines
- Verify: If `mMOR[i, t] == 1`, then `mState[i, t+1] == mState[i, t]` (state freezes at death)

**Maximum Line Limit:**
- Model supports up to 9 lines of therapy
- Patients reaching Line 9 remain at Line 9 (do not progress to Line 10)
- Very few patients (<5%) should reach Line 6+

**ASCT Patients:**
- ASCT occurs during/after Line 1 induction
- ASCT receipt does not create a separate "line" (remains part of Line 1 treatment)
- Post-ASCT, patients eventually progress to Line 2

**Expected Outcome:**
- 100% of patients start at Line 1
- No line skipping occurs in any patient trajectory
- Progression timing mathematically correct (no negative intervals)
- State properly frozen at death
- <1% of patients progress beyond Line 6
- ASCT fits within Line 1 timeline appropriately

**Implementation Notes:**
- Create `test_line_progression.do`
- Check `mState` matrix for each patient across all time points
- Validate time calculations: `mTSD`, `mTXD`, `mTFI` matrices
- Use visual timelines for sample patients to verify sequence

---

### Test 3.2: ASCT Pathway Logic Validation

**Objective:** Verify ASCT pathway correctly implemented with proper eligibility, timing, and BCR assessment.

**Test Cases:**

**Eligibility Criteria:**
- Only patients with planned ASCT (from eligibility risk equation) are candidates
- Patients with Line 1 BCR = SD or PD are ineligible for ASCT (per clinical protocol)
- Verify: If `vSCT_DN == 1` and `BCR_L1 <= 4`, then eligible for ASCT

**Receipt vs Eligibility:**
- Not all eligible patients receive ASCT (some decline, medical complications)
- Receipt rate among eligible: Expected ~70-80%
- Verify receipt risk equation applied only to eligible patients

**Timing of ASCT:**
- ASCT occurs after induction therapy completion
- Typical timing: 3-6 months post-diagnosis
- Verify: `mTSD[, ASCT_col]` approximately equals `mTSD[, L1_start] + L1_duration`

**ASCT BCR Assessment:**
- ASCT has its own BCR assessment (separate from induction BCR)
- ASCT BCR generally equals or improves upon induction BCR
- Verify ASCT BCR stored in appropriate `mCR` column

**Induction Duration with ASCT:**
- Patients with planned ASCT use three-spline model for induction duration
- Verify correct splines applied based on time: Spline 1 (0-3mo), Spline 2 (3-6mo), Spline 3 (6+mo)

**Post-ASCT Maintenance:**
- Higher proportion of ASCT patients receive maintenance therapy
- Maintenance extends Line 1→2 TFI

**Single ASCT Assumption:**
- Model assumes maximum of one ASCT per patient (per MRDR data: 94% receive ≤1 ASCT)
- No patient receives ASCT at Line 2 or later
- Verify: Sum of ASCT receipts per patient ≤1

**Expected Outcome:**
- ASCT eligibility logic correctly implemented
- Only eligible patients (good BCR, planned ASCT) receive ASCT
- ASCT timing falls within expected 3-6 month window
- ASCT BCR assessment separate and generally superior to induction
- Spline model transitions smooth for ASCT patients
- Each patient receives ≤1 ASCT
- ASCT patients have distinct outcomes (longer TFI, better survival)

**Implementation Notes:**
- Create `test_asct_pathway.do`
- Track ASCT-eligible vs recipient patients
- Validate `mSCT` matrix: planned ASCT and receipt columns
- Check timing: `mTSD` columns for induction end and ASCT
- Compare outcomes: ASCT vs non-ASCT subgroups

---

### Test 3.3: Maintenance Therapy Logic Validation

**Objective:** Verify maintenance therapy correctly extends treatment-free intervals without changing BCR.

**Test Cases:**

**BCR Independence:**
- Maintenance therapy should NOT change BCR
- Verify: BCR at Line 1 end = BCR at Line 2 start (maintenance does not improve response)
- Post-maintenance BCR assessment only upon Line 2 initiation

**TFI Extension:**
- Maintenance therapy should extend Line 1→2 TFI
- Expected increase: 3-6 months longer TFI
- Verify: TFI with maintenance > TFI without maintenance (on average)

**Eligibility and Receipt:**
- Only patients completing Line 1 therapy are candidates
- Not all candidates receive maintenance (patient/physician preference, tolerability)
- Verify receipt aligns with risk equation predictions

**Timing Logic:**
- Maintenance begins after induction therapy completion (or post-ASCT for ASCT patients)
- Maintenance extends duration of Line 1 phase (before Line 2 initiation)
- Verify: Time at Line 2 start = Time at Line 1 end + extended TFI (if maintenance)

**ASCT Interaction:**
- ASCT patients can receive maintenance post-transplant
- Maintenance more common among ASCT patients (~70% vs ~40% non-ASCT)

**Survival Impact:**
- Maintenance extends TFI, which indirectly improves survival (delayed progression)
- Verify: Patients with maintenance have longer overall survival (on average)

**Expected Outcome:**
- Maintenance does not alter BCR (remains constant during maintenance)
- Maintenance extends TFI by 3-6 months on average
- Maintenance receipt rates match clinical patterns (higher in ASCT patients)
- Timing logic preserves temporal sequence
- Survival benefit evident through extended TFI

**Implementation Notes:**
- Create `test_maintenance_logic.do`
- Compare BCR before and during maintenance (should be identical)
- Stratify TFI analysis by maintenance receipt
- Verify maintenance flag in appropriate matrix/variable
- Test interaction with ASCT pathway

---

### Test 3.4: Treatment Regimen Transition Logic

**Objective:** Verify regimen transitions between lines follow clinical patterns and model logic.

**Test Cases:**

**No Same-Regimen Repeats:**
- Patients should not receive the same regimen twice consecutively
- Cross-resistance assumption: If regimen fails, don't reuse immediately
- Verify: `mTXR[i, Line N] ≠ mTXR[i, Line N+1]` for all patients and lines

**BCR Influences Next Regimen:**
- Better previous BCR should slightly increase probability of more intensive next regimen
- Test that regimen selection risk equation includes previous BCR as covariate

**Line-Specific Regimen Options:**
- Line 1: VRd, VCd, Other
- Line 2: Rd, DVd, Other
- Line 3+: Averaged regimen (no specific regimens modelled)
- Verify regimen options match line of therapy

**Analysis Scenario Impact:**
- VRd-Post analysis: Higher proportion receive VRd at Line 1
- DVd-Post analysis: Higher proportion receive DVd at Line 2
- Verify coefficient sets change based on analysis scenario

**Expected Outcome:**
- No same-regimen repeats observed
- BCR→regimen relationship positive (better response → more intensive next regimen)
- Regimen options appropriate for each line
- Analysis scenarios correctly modify regimen distributions

**Implementation Notes:**
- Create `test_regimen_transitions.do`
- Generate transition matrices: Line 1 regimen × Line 2 regimen
- Check diagonal elements (same-regimen repeats) are zero
- Compare regimen distributions by BCR category

---

### Test 3.5: Death Event Logic Validation

**Objective:** Verify death events properly terminate patient trajectories and survival calculations are correct.

**Test Cases:**

**Death Timing:**
- Death can occur at any point: diagnosis, during treatment, during TFI, or at later lines
- Verify death time calculated using OS risk equation
- Death time should be recorded in `mTSD[, Death_col]`

**State Freeze at Death:**
- Once patient dies (`mMOR[i, t] == 1`), patient state freezes
- No further line progressions after death
- No further treatment durations or TFIs calculated
- Verify: `mState[i, t+1] == mState[i, t]` for all t after death

**Competing Events:**
- If death occurs before completing a line of therapy, that line counts as "current line at death"
- If death occurs during TFI, patient does not progress to next line

**Death Flag Consistency:**
- Once `mMOR[i, t] == 1`, it remains 1 for all subsequent t
- Verify: No "resurrection" (death flag does not flip back to 0)

**Survival Time Accuracy:**
- Survival time = Time from diagnosis to death
- For alive patients: Survival time = Time from diagnosis to end of simulation
- Verify: `Survival_time = mTSD[, Death_col]` if died, else `= max(simulation_time)`

**Proportion Dying:**
- Substantial proportion of patients should die during simulation horizon (MM is fatal disease)
- Expected: 40-60% mortality over 10-year simulation
- If <20% die, suggests OS risk equation may be miscalibrated

**Expected Outcome:**
- Death timing calculated correctly using OS risk equation
- Patient trajectories properly terminated at death
- No state changes after death
- Death flags permanent (no reversals)
- Survival times accurately calculated
- Mortality rate clinically realistic (40-60% over 10 years)

**Implementation Notes:**
- Create `test_death_logic.do`
- Check `mMOR` matrix for consistency (once 1, stays 1)
- Verify state freeze: compare `mState` before and after death
- Calculate observed mortality rate and compare to clinical benchmarks
- Validate survival time calculations against hand calculations

---

## 4. Survival Calculation Tests

### Overview
These tests validate the mathematical accuracy of parametric survival functions, time-to-event calculations, and conditional survival probabilities.

---

### Test 4.1: Parametric Survival Function Validation

**Objective:** Verify parametric survival time calculations are mathematically correct for chosen distributions.

**Test Cases:**

**Distribution-Specific Formulas:**
- **Weibull:** Test `S(t) = exp(-λt^γ)` where `λ = exp(XB)` and `γ = exp(aux)`
- **Exponential:** Test `S(t) = exp(-λt)` where `λ = exp(XB)`
- **Gompertz:** Test `S(t) = exp(-λ/γ * (exp(γt) - 1))` where `λ = exp(XB)` and `γ = aux`
- **Log-logistic:** Test `S(t) = 1 / (1 + (λt)^γ)` where `λ = exp(XB)` and `γ = exp(aux)`
- Verify formulas used in `calcSurvTime()` and `calcSurvProb()` match distribution choice

**Inverse CDF Method:**
- Survival time calculation: `t = S^(-1)(U)` where U ~ Uniform(0,1)
- Verify inverse transformation correctly implemented
- Test on known values: If U=0.5, should get median survival time

**Survival Probability Calculation:**
- Given time t and covariates, calculate `S(t) = P(Survival > t)`
- Verify `calcSurvProb()` function produces probabilities in [0,1] range
- Test monotonicity: S(t) should decrease as t increases

**Conditional Survival:**
- For competing events: `S(t2 | survived to t1) = S(t2) / S(t1)`
- Verify conditional probability used for OS calculations after each line
- Test: If patient survived to t1, random draw should be `U * S(t1)` not just `U`

**Numerical Stability:**
- Test extreme XB values (e.g., XB = -10, +10)
- Verify no overflow/underflow in `exp()` calculations
- Confirm survival times finite and positive

**Expected Outcome:**
- Survival functions match parametric distribution formulas exactly
- Inverse CDF correctly calculates survival time from uniform random draw
- Survival probabilities all in [0,1] range and monotonically decreasing
- Conditional survival properly adjusts for prior survival
- No numerical errors for extreme parameter values

**Implementation Notes:**
- Create `test_survival_functions.do`
- Use Mata to test `calcSurvTime()` and `calcSurvProb()` functions
- Compare to Stata's `streg` predictions
- Test with various XB and auxiliary parameter values
- Validate against known survival distributions (e.g., exponential median = ln(2)/λ)

---

### Test 4.2: Random Number Generation Validation

**Objective:** Verify random number generation produces proper uniform and reproducible draws.

**Test Cases:**

**Uniform Distribution:**
- All random draws should be uniformly distributed in [0,1]
- Test using Kolmogorov-Smirnov test: `ksmirnov var = 1`
- Generate 10,000 random draws, verify mean ≈ 0.5, variance ≈ 1/12

**Reproducibility:**
- Setting same random seed should produce identical results
- Verify: `set seed 12345; simulate run 1` ≡ `set seed 12345; simulate run 2`
- Critical for debugging and validation

**Independence:**
- Random draws should be independent (no autocorrelation)
- Test using runs test or Ljung-Box test
- No systematic patterns in sequence of random numbers

**Seed Documentation:**
- If bootstrap enabled, seeds should be documented for each replicate
- Verify seed properly set at simulation start

**Expected Outcome:**
- Random numbers uniformly distributed in [0,1]
- Same seed produces identical simulation results
- Random draws independent (no autocorrelation)
- Seeds properly documented

**Implementation Notes:**
- Create `test_random_numbers.do`
- Generate large sample of random draws
- Use `ksmirnov`, `summarize` to test distribution
- Run simulation twice with same seed, verify identical output

---

### Test 4.3: Time Accumulation Validation

**Objective:** Verify time variables correctly accumulate across events and lines of therapy.

**Test Cases:**

**Time Since Diagnosis (TSD):**
- TSD at any event = Sum of all prior durations and intervals
- Verify: `mTSD[i, event_k] = sum(mTXD[i, 1:k-1]) + sum(mTFI[i, 1:k-1])`
- TSD should be monotonically increasing (never decreases)

**Treatment Duration Accumulation:**
- Each line has a treatment duration
- Verify: Treatment duration properly added to TSD at end of treatment

**Treatment-Free Interval Accumulation:**
- Each TFI adds to TSD before next line starts
- Verify: TFI properly added to TSD at next line start

**Death Time:**
- Death time represents total survival from diagnosis
- Verify: `mTSD[, Death_col]` = survival time calculated via OS risk equation
- Should be consistent with accumulated times from all events

**No Negative Intervals:**
- All duration and interval variables must be ≥0
- Verify: `min(mTXD) >= 0` and `min(mTFI) >= 0`

**Temporal Ordering:**
- Event times should be ordered: Diagnosis < Line 1 start < Line 1 end < Line 2 start < ... < Death
- Verify ordering for every patient

**Expected Outcome:**
- TSD monotonically increases through patient trajectory
- All durations and intervals non-negative
- Event times properly ordered (no temporal violations)
- Death time consistent with event accumulation

**Implementation Notes:**
- Create `test_time_accumulation.do`
- Check `mTSD`, `mTXD`, `mTFI` matrices for each patient
- Calculate expected TSD and compare to model output
- Generate patient timelines for visual verification

---

## 5. Edge Case and Boundary Tests

### Overview
These tests examine model behaviour under extreme or unusual conditions to ensure robustness.

---

### Test 5.1: Extreme Patient Characteristics

**Objective:** Verify model handles extreme patient profiles without errors or implausible results.

**Test Cases:**

**Very Young Patient:**
- Age = 30, best other characteristics (ECOG=0, R-ISS=I, CMScore=0)
- Should have very long predicted survival (>15 years)
- Should have high ASCT probability (>80%)

**Very Old Patient:**
- Age = 95, worst other characteristics (ECOG=2, R-ISS=III, CMScore=6)
- Should have very short predicted survival (<2 years)
- Should have very low ASCT probability (<1%)

**All Poor Prognostic Factors:**
- Age=85, ECOG=2, R-ISS=III, CMScore=5
- Should have poor outcomes across all endpoints
- Model should not crash or produce negative values

**All Good Prognostic Factors:**
- Age=50, ECOG=0, R-ISS=I, CMScore=0
- Should have excellent outcomes across all endpoints
- Survival should be extended but realistic (<30 years)

**Expected Outcome:**
- No errors or crashes for extreme profiles
- Survival predictions clinically plausible (no >30 year predictions)
- Probabilities remain in [0,1] range
- Extreme cases produce extreme (but realistic) outcomes

**Implementation Notes:**
- Create `test_edge_cases.do`
- Manually create patients with extreme characteristics
- Run simulation and check for errors
- Verify outcomes are monotonic with risk profile

---

### Test 5.2: Immediate Disease Progression

**Objective:** Verify model handles patients with immediate progression or death appropriately.

**Test Cases:**

**Immediate Progression (BCR = PD at Line 1):**
- Patient achieves PD at Line 1
- Should be ineligible for ASCT
- Should have very short TFI before Line 2
- Verify model does not assign ASCT to PD patients

**Immediate Death (Death at Diagnosis):**
- Patient dies before starting Line 1 treatment
- Verify: No Line 1 treatment assigned
- Patient should have survival time ≈ 0 (diagnosis to death interval minimal)

**Death During Induction:**
- Patient dies during Line 1 treatment (before completion)
- Verify: Line 1 marked as incomplete
- No Line 2 progression

**Expected Outcome:**
- PD patients correctly excluded from ASCT
- Immediate death properly recorded
- Death during treatment handled correctly (no invalid state transitions)

**Implementation Notes:**
- Create patients with very poor prognosis (force early death)
- Test with manipulated random seeds to trigger early events
- Verify proper handling in state transition logic

---

### Test 5.3: Maximum Line of Therapy

**Objective:** Verify model correctly handles patients reaching maximum line of therapy (Line 9).

**Test Cases:**

**Progression to Line 9:**
- Create very long-surviving patient who progresses through all lines
- Verify patient reaches Line 9 without errors

**Staying at Line 9:**
- Verify patient remains at Line 9 (does not progress to Line 10)
- State should freeze at Line 9 for remaining survival time

**Coefficient Equivalence:**
- Lines 6-9 share same risk equation coefficients
- Verify correct coefficient application at Lines 6, 7, 8, 9

**Expected Outcome:**
- Model supports progression through all 9 lines
- No errors at Line 9
- Patient correctly remains at Line 9 (no Line 10)

**Implementation Notes:**
- Create synthetic patient with excellent prognosis and very long survival
- Force progression through all lines by manipulating TFI
- Verify `mState` matrix maxes out at 9

---

### Test 5.4: Bootstrap and Uncertainty Tests

**Objective:** Verify bootstrap functionality produces appropriate uncertainty estimates.

**Test Cases:**

**Bootstrap Activation:**
- Verify bootstrap runs when `$Bootstrap = 1`
- Verify specified number of bootstrap replicates execute
- Check that bootstrap samples have variation in results

**Coefficient Sampling:**
- If parametric bootstrap implemented, verify coefficients sampled from variance-covariance matrix
- Check that sampled coefficients maintain correlation structure

**Result Aggregation:**
- Verify results properly averaged across bootstrap replicates
- Calculate standard errors and confidence intervals
- Check that CI width is reasonable

**Reproducibility:**
- Setting bootstrap seed should produce reproducible results
- Verify: Same seed → same bootstrap estimates

**Expected Outcome:**
- Bootstrap runs without errors
- Results show appropriate uncertainty (not all replicates identical)
- Confidence intervals have reasonable width
- Bootstrap reproducible with fixed seed

**Implementation Notes:**
- Create `test_bootstrap.do`
- Run small number of bootstrap replicates (5-10) for testing
- Verify result storage and aggregation
- Check seed handling across replicates

---

## 6. Performance and Scalability Tests

### Overview
These tests evaluate computational efficiency and ensure the model can handle large patient cohorts.

---

### Test 6.1: Computational Time Benchmarks

**Objective:** Establish baseline computation times and identify performance bottlenecks.

**Test Cases:**

**Small Cohort (n=100):**
- Measure time to simulate 100 patients
- Should complete in <1 minute
- Benchmark for iterative development/testing

**Medium Cohort (n=1,000):**
- Measure time to simulate 1,000 patients
- Should complete in <5 minutes
- Standard analysis cohort size

**Large Cohort (n=4,884):**
- Measure time to simulate full MRDR-sized cohort
- Should complete in <20 minutes
- Real-world application cohort

**Performance Profiling:**
- Identify which risk equations take longest to calculate
- Determine if bottlenecks in OS, BCR, or duration equations
- Profile Mata vs Stata operations

**Expected Outcome:**
- Linear scaling with cohort size (doubling patients ≈ doubles time)
- Identify any quadratic or exponential scaling issues
- Document baseline performance for future optimization

**Implementation Notes:**
- Create `test_performance.do`
- Use `timer on/off` to profile code sections
- Run on standardised hardware
- Test both matrix and loop-based operations

---

### Test 6.2: Memory Usage Validation

**Objective:** Verify model does not exceed memory limits with large cohorts.

**Test Cases:**

**Matrix Size Estimation:**
- Calculate expected matrix sizes for various cohorts
- Verify matrices fit in available memory
- Check for memory leaks (memory usage should be stable)

**Large Cohort Simulation:**
- Simulate 10,000 patients
- Monitor memory usage throughout simulation
- Verify memory released after simulation completes

**Expected Outcome:**
- Memory usage scales linearly with cohort size
- No memory leaks detected
- Large cohorts (10,000+) successfully simulated

**Implementation Notes:**
- Use `memory` command to check Stata memory settings
- Monitor system memory during execution
- Test on machine with limited RAM to identify issues

---

### Test 6.3: Parallel Processing (if implemented)

**Objective:** Verify parallel processing correctly splits workload and aggregates results.

**Test Cases:**

**Result Consistency:**
- Run same simulation single-core vs multi-core
- Verify identical results (given same random seed)

**Speed-up Factor:**
- Measure computation time with 1, 2, 4, 8 cores
- Verify near-linear speed-up (4 cores ≈ 4× faster)

**Expected Outcome:**
- Parallel processing produces identical results
- Achieves substantial speed-up
- No race conditions or synchronisation issues

**Implementation Notes:**
- If using parallel processing, create `test_parallel.do`
- Test on multi-core machine
- Verify thread-safety of random number generation

---

## 7. Output Verification Tests

### Overview
These tests validate the final simulation output dataset contains correct variables, values, and distributions.

---

### Test 7.1: Output Variable Completeness

**Objective:** Verify all expected output variables are present and populated.

**Test Cases:**

**Patient Identifiers:**
- `patid`: Unique patient identifier (no duplicates)
- `bootstrap`: Bootstrap replicate number (if applicable)

**Baseline Characteristics:**
- `age`, `male`, `ecog`, `riss`, `cmscore`: All baseline variables present
- No missing values in baseline characteristics

**Survival Outcomes:**
- `os_time`: Overall survival time from diagnosis
- `death_flag`: Binary indicator of death (0/1)
- All patients have valid os_time (no missing)

**Treatment Variables:**
- `sct_planned`, `sct_received`: ASCT indicators
- `maintenance_received`: Maintenance therapy indicator
- `regimen_l1`, `regimen_l2`, etc.: Treatment regimen at each line

**Response Variables:**
- `bcr_l1`, `bcr_l2`, etc.: Best clinical response at each line
- Values should be 1-6 (CR to PD) or missing if line not reached

**Duration Variables:**
- `txd_l1`, `txd_l2`, etc.: Treatment duration at each line
- `tfi_l1_l2`, `tfi_l2_l3`, etc.: Treatment-free intervals

**Time Variables:**
- `tsd_l1_start`, `tsd_l1_end`, etc.: Time since diagnosis at each event
- All times should be non-negative

**Expected Outcome:**
- All expected variables present in output dataset
- No unexpected variables (clean output)
- Variable types correct (numeric vs string)
- Variable labels informative

**Implementation Notes:**
- Create `test_output_variables.do`
- Use `describe` and `codebook` to check variables
- Compare to documentation of expected output structure

---

### Test 7.2: Output Value Range Validation

**Objective:** Verify all output values fall within clinically plausible ranges.

**Test Cases:**

**Survival Time:**
- All survival times >0
- <1% survival times >20 years (very long survivals rare)
- Mean survival ~5-7 years for contemporary cohorts

**Treatment Duration:**
- All durations ≥0
- Typical Line 1 duration: 3-12 months
- Flag any durations >36 months

**Treatment-Free Intervals:**
- All TFI ≥0
- Typical Line 1→2 TFI: 6-24 months (higher for ASCT patients)
- Flag any TFI >60 months (very unusual)

**BCR Values:**
- All BCR values in {1,2,3,4,5,6} or missing
- No invalid codes (0, 7, negative values)

**ASCT Indicators:**
- All ASCT variables binary (0/1)
- If ASCT received, must have ASCT planned
- Maximum 1 ASCT per patient

**Expected Outcome:**
- All values within valid ranges
- Extreme values flagged for review (<1% of cohort)
- No impossible values (negative times, invalid codes)

**Implementation Notes:**
- Create `test_output_ranges.do`
- Use `summarize` to check ranges
- Use `assert` statements for hard constraints
- Generate outlier reports

---

### Test 7.3: Output Distribution Validation

**Objective:** Verify output distributions match expected clinical patterns and MRDR benchmarks.

**Test Cases:**

**Survival Distribution:**
- Median OS ~5-7 years
- 5-year survival rate ~50-60%
- 10-year survival rate ~30-40%
- Compare to MRDR Kaplan-Meier curves

**BCR Distribution by Line:**
- Line 1 BCR: CR/VGPR ~50%, PR ~30%, MR/SD/PD ~20%
- Line 2 BCR: Generally poorer than Line 1
- Compare to MRDR observed response rates

**ASCT Rate:**
- Overall ASCT receipt: ~40-50% of eligible patients
- Higher in younger patients (<65 years)
- Compare to MRDR observed ASCT rate

**Line Reached Distribution:**
- ~100% reach Line 1
- ~60-70% reach Line 2
- ~40-50% reach Line 3
- <20% reach Line 4+
- Compare to MRDR progression patterns

**Treatment Regimen Distribution:**
- Line 1: VRd ~60% (post-VRd era), VCd ~25%, Other ~15%
- Line 2: DVd ~40% (post-DVd era), Rd ~15%, Other ~45%
- Compare to MRDR regimen utilisation

**Expected Outcome:**
- All distributions match MRDR patterns within ±15%
- No bimodal or unusual distribution shapes (unless clinically expected)
- Population-level statistics consistent with published MM outcomes

**Implementation Notes:**
- Create `test_output_distributions.do`
- Generate histograms and Kaplan-Meier plots
- Compare to MRDR using chi-square tests and log-rank tests
- Calculate Harrell's C-statistic for survival model discrimination

---

### Test 7.4: Internal Consistency Checks

**Objective:** Verify logical consistency across related output variables.

**Test Cases:**

**ASCT Consistency:**
- If `sct_received == 1`, then `sct_planned == 1`
- If `sct_received == 1`, then `bcr_l1 <= 4` (no SD/PD patients receive ASCT)

**Maintenance Consistency:**
- If `maintenance_received == 1`, patient must have completed Line 1
- Maintenance patients should have longer `tfi_l1_l2` on average

**Line Reached Consistency:**
- If Line 3 variables populated, then Line 2 variables must be populated
- If patient reached Line N, must have reached all prior lines 1 to N-1

**Survival Consistency:**
- If `death_flag == 1`, then `os_time` should equal time to death
- If `death_flag == 0`, then `os_time` should equal censoring time (end of follow-up)

**Time Consistency:**
- `tsd_l2_start` should equal `tsd_l1_end + tfi_l1_l2`
- Verify temporal ordering: all `tsd_*` variables increase monotonically

**Expected Outcome:**
- No logical inconsistencies detected
- All relationships between variables hold
- Zero tolerance for consistency violations

**Implementation Notes:**
- Create `test_output_consistency.do`
- Use `assert` statements for each consistency check
- Generate cross-tabulations to identify violations
- Flag any inconsistent observations for review

---

## 8. Comparative Validation Tests

### Overview
These tests compare simulation outputs to external benchmarks and published results.

---

### Test 8.1: MRDR Validation Cohort Comparison

**Objective:** Validate model predictions against the MRDR validation cohort (30% of registry).

**Test Cases:**

**Overall Survival:**
- Compare predicted vs observed Kaplan-Meier curves
- Use log-rank test for equality
- Calculate integrated Brier score for calibration

**BCR Distribution:**
- Compare predicted vs observed response rates by regimen and line
- Use chi-square tests for each BCR category

**ASCT Rate:**
- Compare predicted vs observed ASCT receipt rate
- Stratify by age, ECOG, and R-ISS

**Treatment Duration:**
- Compare mean predicted vs observed duration by regimen and line
- Use t-tests or Wilcoxon tests

**Line Progression:**
- Compare proportion reaching each line (observed vs predicted)
- Verify similar drop-off pattern as patients progress through lines

**Expected Outcome:**
- Predicted survival curves closely match observed (log-rank p>0.05)
- BCR distributions match within ±10% for each category
- ASCT rate matches within ±5%
- Discrimination metrics (C-statistic) >0.70 for survival models
- Calibration plots show good agreement (slope ≈1, intercept ≈0)

**Implementation Notes:**
- Use MRDR validation cohort (separate from training cohort)
- Generate side-by-side comparisons: predicted vs observed
- Calculate formal validation statistics: C-statistic, integrated Brier score, calibration slope
- Document in validation report

---

### Test 8.2: Published Trial Data Comparison

**Objective:** Compare simulation results to published clinical trial outcomes for key regimens.

**Test Cases:**

**VRd Efficacy (vs SWOG S0777 trial):**
- Median PFS with VRd: Compare simulated TFI to trial result (~43 months)
- ORR with VRd: Compare simulated BCR distribution to trial (ORR ~82%)

**DVd Efficacy (vs CASTOR trial):**
- Median PFS with DVd at relapse: Compare to trial result (~16 months)
- ORR with DVd: Compare to trial (ORR ~83%)

**ASCT Benefit (vs IFM 2009 trial):**
- Compare survival benefit of ASCT vs no ASCT
- Expected PFS advantage: ~5-10 months longer

**Expected Outcome:**
- Simulated outcomes within confidence intervals of published trials
- Direction of effects consistent (e.g., VRd superior to VCd)
- Magnitude of effects clinically plausible (not over/underestimating)

**Implementation Notes:**
- Create `test_trial_comparison.do`
- Restrict simulated cohort to match trial eligibility criteria
- Compare primary endpoints: ORR, PFS, OS
- Acknowledge differences between RCT and real-world populations

---

### Test 8.3: External Economic Model Comparison

**Objective:** Compare model structure and outcomes to previously published MM economic models.

**Test Cases:**

**Survival Predictions:**
- Compare to published partitioned survival models or Markov models
- Verify similar median OS for comparable patient profiles

**QALY Estimates:**
- If utility weights applied, compare QALYs to published economic models
- Verify within expected range for MM population

**Treatment Sequences:**
- Compare proportion of patients receiving various treatment sequences
- Verify alignment with real-world evidence studies

**Expected Outcome:**
- Model produces comparable survival estimates to published models
- Treatment pathway distributions similar to cohort studies
- Any discrepancies explainable by differences in data sources or assumptions

**Implementation Notes:**
- Literature review of published MM economic models
- Document structural differences (DES vs Markov vs PSM)
- Highlight advantages of individual-level DES approach

---

## 9. Reproducibility and Documentation Tests

### Overview
These tests ensure the model is reproducible, well-documented, and accessible to external users.

---

### Test 9.1: Random Seed Reproducibility

**Objective:** Verify setting random seed produces identical results.

**Test Cases:**

**Exact Replication:**
- Set `seed 12345`
- Run simulation on test cohort
- Re-run with same seed
- Verify outputs identical (patient-level and aggregate)

**Different Seed Variation:**
- Run with `seed 12345` and `seed 67890`
- Verify results differ (stochasticity present)
- Verify aggregate results similar (within Monte Carlo error)

**Expected Outcome:**
- Same seed → identical results (bit-for-bit reproducibility)
- Different seeds → different patient-level results but similar aggregates
- Seed properly documented in output/log files

**Implementation Notes:**
- Create `test_reproducibility.do`
- Use `cf` (compare files) command to verify identical datasets
- Document seed setting procedures in user guide

---

### Test 9.2: Documentation Completeness

**Objective:** Verify all code modules are well-documented with clear comments.

**Test Cases:**

**Code Comments:**
- Every do-file has header with purpose, inputs, outputs
- Complex sections have explanatory comments
- All Mata functions documented (purpose, parameters, return values)

**User Guide:**
- Installation instructions clear and complete
- Example commands provided
- All command-line arguments explained

**Technical Specification:**
- All 50 risk equations documented with covariates
- Parametric distributions specified for survival models
- Coefficient matrices described

**Changelog:**
- Version history documented
- Changes from v1.0 to v2.0 clearly explained

**Expected Outcome:**
- External user can understand and run model without author assistance
- All code modules have clear documentation
- Model decisions justified and referenced to clinical literature

**Implementation Notes:**
- Code review for documentation gaps
- User guide tested with naive user
- Technical specification compared to published PLOS ONE paper

---

### Test 9.3: Cross-Platform Compatibility

**Objective:** Verify model runs on different operating systems and Stata versions.

**Test Cases:**

**Operating Systems:**
- Test on Windows, macOS, and Linux
- Verify file paths and directory separators compatible
- Check for OS-specific commands that may fail

**Stata Versions:**
- Test on Stata 15, 16, 17, 18
- Verify Mata code compatible across versions
- Document minimum required Stata version

**Hardware:**
- Test on standard laptop (8GB RAM)
- Test on high-performance workstation (64GB RAM)
- Verify memory requirements documented

**Expected Outcome:**
- Model runs on all major OS platforms
- Minimum Stata version clearly specified (Stata 15+)
- Memory requirements documented (adequate for 10,000 patients)

**Implementation Notes:**
- Coordinate testing across development team members
- Document platform-specific issues and workarounds
- Use relative paths (not absolute) for file references

---

## 10. Implementation Prioritisation

### Recommended Testing Sequence

**Phase 1: Critical Path Tests (Week 1)**
1. Test 1.1: Input validation (ensure data quality)
2. Test 1.2: Dummy variable generation (foundation for all equations)
3. Test 2.1: Coefficient loading (must work before testing equations)
4. Test 2.2: OS risk equation (most critical outcome)
5. Test 7.1-7.2: Output completeness and ranges (basic sanity checks)

**Phase 2: Core Model Logic Tests (Week 2)**
6. Test 2.3: ASCT eligibility
7. Test 2.4: BCR predictions
8. Test 2.5: Treatment duration
9. Test 2.6: TFI validation
10. Test 3.1-3.2: Treatment pathway and ASCT logic

**Phase 3: Comprehensive Validation (Week 3)**
11. Test 4.1-4.3: Survival calculations and time accumulation
12. Test 3.3-3.5: Maintenance, regimen transitions, death logic
13. Test 7.3-7.4: Output distributions and consistency

**Phase 4: Advanced Testing (Week 4)**
14. Test 5.1-5.4: Edge cases and bootstrap
15. Test 8.1-8.3: External validation against MRDR and trials
16. Test 6.1-6.2: Performance benchmarking

**Phase 5: Documentation and Reproducibility (Week 5)**
17. Test 9.1-9.3: Reproducibility, documentation, cross-platform
18. Generate comprehensive validation report
19. Address any outstanding issues from earlier phases

---

## 11. Test Output Requirements

### For Each Test Module

**Test Report Should Include:**

1. **Test Identification**
   - Test number and name
   - Date executed
   - Tester name
   - Stata version and OS

2. **Test Results**
   - Pass/Fail status
   - Number of assertions tested
   - Number of assertions passed/failed
   - Specific failures documented with patient IDs or indices

3. **Summary Statistics**
   - Key outcomes measured (means, SDs, distributions)
   - Comparison to benchmarks or expected values
   - Statistical test results (p-values, confidence intervals)

4. **Visualisations**
   - Plots comparing predicted vs observed (where applicable)
   - Distribution histograms
   - Survival curves

5. **Recommendations**
   - If test failed: specific recommendations for fixes
   - If test passed: any concerns or edge cases noted
   - Suggestions for additional testing

---

## 12. Continuous Testing Strategy

### Ongoing Validation

**After Each Model Modification:**
- Re-run Phase 1 critical path tests
- Re-run tests specific to modified module
- Update test documentation if new tests needed

**Before Each Release:**
- Execute complete test suite (all phases)
- Generate comprehensive validation report
- Obtain sign-off from clinical advisory group

**Quarterly:**
- Re-validate against updated MRDR data
- Update benchmarks as new registry data available
- Revise tests to reflect evolving clinical practice

---

## 13. Appendices

### Appendix A: Test Data Files

**Suggested Test Datasets:**

1. **test_basic_100.dta:** 100 patients with full covariate data
2. **test_extreme_cases.dta:** 50 patients with extreme characteristics
3. **test_missing_data.dta:** 100 patients with strategic missing values
4. **test_large_cohort.dta:** 10,000 patients for scalability testing
5. **mrdr_validation_cohort.dta:** 1,237 MRDR validation patients

### Appendix B: Key Performance Indicators

**Model Performance Targets:**

- **Discrimination:** C-statistic >0.70 for survival models
- **Calibration:** Calibration slope 0.9-1.1, intercept -0.1 to +0.1
- **Accuracy:** Population-level outcomes within ±10% of MRDR
- **Completeness:** 100% of patients have valid outcome predictions
- **Efficiency:** Simulate 1,000 patients in <5 minutes

### Appendix C: Statistical Test Reference

**Suggested Statistical Tests:**

- **Survival comparison:** Log-rank test, Cox proportional hazards
- **Distribution comparison:** Chi-square test, Kolmogorov-Smirnov test
- **Mean comparison:** t-test, Wilcoxon rank-sum test
- **Correlation:** Spearman or Pearson correlation coefficients
- **Calibration:** Integrated Brier score, calibration plots

### Appendix D: Glossary

- **ASCT:** Autologous Stem Cell Transplantation
- **BCR:** Best Clinical Response (6-category IMWG scale)
- **DES:** Discrete Event Simulation
- **ECOG:** Eastern Cooperative Oncology Group performance status
- **LoT:** Line of Therapy
- **MM:** Multiple Myeloma
- **MRDR:** Australia and New Zealand Myeloma and Related Diseases Registry
- **OS:** Overall Survival
- **R-ISS:** Revised International Staging System
- **TFI:** Treatment-Free Interval
- **VRd:** Bortezomib, lenalidomide, dexamethasone regimen

---

## Conclusion

This comprehensive testing suite provides a structured approach to validating the Monash Myeloma Model discrete event simulation model. By systematically testing input validation, risk equations, treatment pathways, survival calculations, edge cases, performance, and outputs, we can ensure the model produces clinically valid, reproducible, and reliable predictions of multiple myeloma outcomes.

The testing framework is designed to be:
- **Comprehensive:** Covering all major model components
- **Modular:** Tests can be run independently or as complete suite
- **Reproducible:** Clear documentation enables replication
- **Iterative:** Supports continuous validation as model evolves

Implementation of this testing suite will provide confidence in model predictions for clinical, policy, and research applications.

---

**Document Version:** 1.0  
**Last Updated:** 27 October 2025  
**Contact:** adam.irving@monash.edu
