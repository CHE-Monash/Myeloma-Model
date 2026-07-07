**********
* Monash Myeloma Model - TXR OOS
*
* Purpose: Regimen definitions (TXR_L1..L9, built by gen_txr in prep/risk_equations.do) for the OOS
*          70/30 validation. Declared inline so the OOS job is self-contained.
* Notes:   Must be IDENTICAL to analyses/base_model/outcomes/txr_base_model.do -- KEEP IN SYNC. Inline
*          (not sourced) because HPC risk-equation runs ship only analyses/oos/ + prep/ + hpc/ + core/,
*          not analyses/base_model/, so a cross-folder `do` errors r(601).
**********

* Regimen codes:
*   2  Thal/Cycl/Dexa    4  Bort/Cycl/Dexa    7  Lena/Dexa       9  Bort/Thal/Dexa
*  31  Bort/Lena/Dexa   49  Carf/Dexa        56  Poma/Dexa      80  Dara/Bort/Dexa

global TXR_L1 "4 31"
global TXR_L2 "7 80"
global TXR_L3 "49 7"
global TXR_L4 "49 56"
* L5-L9 unset => all 'other'
