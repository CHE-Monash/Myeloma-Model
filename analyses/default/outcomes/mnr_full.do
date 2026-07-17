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

* The default analysis fits diagnosis years 1995-2040, and its out-of-sample validation scores
* against observed history - where thalidomide was the MAJORITY maintenance regimen until 2020
* (51.6% of starts in 2019, 0 from 2021). So len + bort + thal are modelled here.
* Daratumumab and carfilzomib are 22 and 18 patients, are mostly later-line (7.6), and break a
* 6-level fit (daratumumab coef 7.08, SE 265), so they fall to 'other'.
global MNR_L1 "1 4 5"

* A current-paradigm analysis would set:  global MNR_L1 "1"
* From 2021 the mix is ~90% lenalidomide, ~7% bortezomib, 0% thalidomide and 0% carfilzomib, so
* 'other' becomes essentially bortezomib. Safe because the regimens inside 'other' are
* interchangeable for the share (bortezomib -0.658 vs thalidomide -0.638, against lenalidomide
* +0.330). Note the PRICES are not interchangeable, so a costed current-paradigm analysis also
* needs the maintenance DPMQ question in docs/refractory.md 7.4 answered first.
