**********
* TXR - out-of-sample (70/30) validation
*
* The OOS analysis validates the base-model structure on held-out patients, so its regimen definitions
* must be IDENTICAL to analyses/base_model/outcomes/txr_base_model.do. They are declared inline here
* (not sourced from base_model) so the OOS job is self-contained: HPC risk-equation runs ship only
* analyses/oos/ + prep/ + hpc/ + core/, NOT analyses/base_model/, so a cross-folder `do` errors r(601).
* KEEP IN SYNC with txr_base_model.do. TXR_L1..L9 are built by gen_txr in prep/risk_equations.do.
*
* Regimen codes:
*   2  Thal/Cycl/Dexa    4  Bort/Cycl/Dexa    7  Lena/Dexa       9  Bort/Thal/Dexa
*  31  Bort/Lena/Dexa   49  Carf/Dexa        56  Poma/Dexa      80  Dara/Bort/Dexa
**********

global TXR_L1 "4 31"
global TXR_L2 "7 80"
global TXR_L3 "49 7"
global TXR_L4 "49 56"
* L5-L9 unset => all 'other'
