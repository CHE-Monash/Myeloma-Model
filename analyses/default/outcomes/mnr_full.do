**********
* Monash Myeloma Model - MNR maintenance regimens (default analysis)
*
* Purpose: declare the maintenance drug codes to model. gen_mnr in prep/risk_equations.do builds
*          MNR_L1 from this list; any drug not listed falls into 0 = 'other'. Mirrors txr_full.do.
*          This is the CANONICAL maintenance list for the default analysis; the train fit
*          (mnr_train.do) sources it, so the in-sample/out-of-sample validation uses the same
*          regimens.
* Notes:   the list is per-analysis on purpose - the maintenance mix changed abruptly in 2020 and
*          the whole-registry mix describes no year that ever existed. Rationale and the year
*          table: docs/refractory.md 7.4. For the mix by year before choosing a list, see
*          scratch/refractory/mnr_recency.do.
**********

* Maintenance drug codes (MNR_L1 as built by data_extraction.do; see docs/refractory.md 2):
*   0  none/other    1  lenalidomide    2  daratumumab
*   3  carfilzomib   4  bortezomib      5  thalidomide

* ONE list serves both eras, which is why there is no per-analysis switch here:
*   OOS / historical   len + thal   - covers ~85% explicitly. Thalidomide was the MAJORITY
*                                     regimen until 2020 (51.6% of starts in 2019, 0 from 2021).
*   Current paradigm   len only     - thalidomide simply empties out in a modern window, and the
*                                     r(r) == 1 guard in risk_equations.do then assigns everyone
*                                     lenalidomide, rather than needing a different list.
global MNR_L1 "1 5"

* SIMPLE-FIRST: lenalidomide and thalidomide only (docs/refractory.md 4.4). The fits restrict to
* inlist(MNR_L1, 1, 5), so bortezomib, daratumumab and carfilzomib maintenance are EXCLUDED from
* the regimen and duration estimation, and the engine never produces an 'other' maintenance
* regimen - those patients are assigned lenalidomide or thalidomide by the regimen logit. This is
* the accepted simplification: bortezomib is ~5 to 11% of starts but has no PBS maintenance DPMQ to
* price separately anyway (MSAG Clinical Practice Guideline, June 2022), and daratumumab/carfilzomib
* maintenance is a couple of dozen patients, mostly later-line. Adding 'other' back is a later step
* if a modern-paradigm analysis needs bortezomib costed explicitly.
