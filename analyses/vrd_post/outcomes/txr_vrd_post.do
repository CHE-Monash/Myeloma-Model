**********
* Monash Myeloma Model - TXR regimens per line (vrd_post)
*
* Purpose: declare the per-line regimen code lists (MRDR Regimen codes) for the VRd-post scenario, which
*          models a different set than the base model (e.g. Lena/Dexa (7) added at L1/L3/L4). gen_txr in
*          prep/risk_equations.do builds TXR_L1..L9 from these; any regimen not listed falls into 0 = 'other'.
**********

* Regimen codes:
*   2  Thal/Cycl/Dexa    4  Bort/Cycl/Dexa    7  Lena/Dexa       9  Bort/Thal/Dexa
*  31  Bort/Lena/Dexa   49  Carf/Dexa        56  Poma/Dexa      80  Dara/Bort/Dexa

global TXR_L1 "4 7 31"
global TXR_L2 "7 80"
global TXR_L3 "7 49"
global TXR_L4 "7 49"
* L5-L9 unset => all 'other'
