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
*   OOS / historical   len + thal + other(~bortezomib)   - covers 85% explicitly. Thalidomide was
*                                                           the MAJORITY regimen until 2020 (51.6%
*                                                           of starts in 2019, 0 from 2021).
*   Current paradigm   len + thal(empty) + other(~bort)  - covers 90% explicitly.
* Thalidomide simply empties out in a modern window rather than needing a different list.
global MNR_L1 "1 5"

* Bortezomib falls into 'other' deliberately. It IS current practice (5 to 11% of starts, and 95.4%
* of its episodes are in-gap, so it is genuinely L1 maintenance - it cannot be dismissed as
* later-line the way daratumumab can). But per the MSAG Clinical Practice Guideline (June 2022) it
* is neither TGA-registered nor PBS-reimbursed for maintenance, so there is no maintenance DPMQ to
* price it with separately whatever we do here. Pooling costs little: L1_MND carries a regimen
* slope (7.4), so 'other' gets its own slope, and inside 'other' bortezomib dominates.
*
* Daratumumab and carfilzomib also fall to 'other': 22 and 18 patients, mostly later-line (7.6),
* and they break a 6-level fit (daratumumab coef 7.08, SE 265).
