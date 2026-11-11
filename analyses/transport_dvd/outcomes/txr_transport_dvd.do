**********
* Monash Myeloma Model - TXR (transport_dvd)
*
* Purpose: Declare every regimen list this analysis models - per-line treatment (TXR_L1..L9), L1
*          maintenance (MNR_L1), and the len-refractory gate (LENREFR_regimens). TXR_L1..L9 are
*          built generically by gen_txr in prep/risk_equations.do; any regimen not listed for a line
*          falls into 0 = 'other'. L2 carries Vd (5, comparator) and DVd (80, intervention)
*          explicitly - the decision-relevant contrast is at second line.
**********

* Regimen codes:
*   4  Bort/Cycl/Dexa    5  Bort/Dexa (Vd)    7  Lena/Dexa      31  Bort/Lena/Dexa
*  49  Carf/Dexa        56  Poma/Dexa        80  Dara/Bort/Dexa (DVd)

global TXR_L1 "4 31"
global TXR_L2 "5 7 80"
global TXR_L3 "7 49"
global TXR_L4 "49 56"
* L5-L9 unset => all 'other'

* Lenalidomide-containing regimen codes, for the refractory gate. The FIT conditions on the true
* drug binary (Lenalidomide == 1); the engine has no drug field, so sim_lenrefr.do counts a line as
* lenalidomide iff the DRAWN regimen is one of these. 'other' (code 0) is treated as non-len, which
* under-counts (docs/refractory.md 5(6)). Same list as the default and car_t analyses, and it covers
* every lenalidomide regimen this analysis models (7 at L2/L3, 31 at L1).
* Saved into the .mmat with the coefficients, so it must be declared at REFIT time, not run time.
global LENREFR_regimens "7 31"

**********
* MAINTENANCE regimens (MNR_L1). gen_mnr builds MNR_L1 from this list; any drug not listed falls to
* 0 = 'other'. Drug codes, per docs/refractory.md 2:
*   0 none/other   1 lenalidomide   2 daratumumab   3 carfilzomib   4 bortezomib   5 thalidomide
*
* LENALIDOMIDE ONLY, as car_t declares and unlike the default analysis's "1 5". Thalidomide
* maintenance ended with the 2020 PBS listing (51.6% of starts in 2019, none from 2021:
* scratch/refractory/mnr_recency.do), and this analysis's cohort enters L2 in 2020-2025, so its L1
* maintenance is a post-listing paradigm. risk_equations.do skips the L1_MND_THAL fit when 5 is
* absent and stores oL1_MNR = 1, so sim_mnr assigns every maintenance patient lenalidomide.
*
* This matters through the POOL BUILD, not the analysis. $line is 2, so maintenance is never
* simulated or costed from L2 onward; but patients/sample_size/cohort_pool.do runs from diagnosis,
* and the len-refractory draw fires only where vMNR == 1. Excluding thalidomide therefore exposes
* every maintenance patient in the pool to that draw rather than skipping the thalidomide share,
* which raises LenRefr prevalence at L2 entry. Check the realised prevalence against the registry's
* 15.5% at L2 after rebuilding the pool.
*
* Lenalidomide is NOT optional: it is the mlogit base outcome and sim_mnr's single-regimen fallback,
* so risk_equations.do exits 2000 on a list that omits it (or an empty list).
global MNR_L1 "1"
