**********
* TXR - VRd Post scenario (regimen definitions)
*
* Declares only the per-line regimen code lists (this scenario models a different set than the base
* model - e.g. Lena/Dexa (7) added at L1/L3/L4). TXR_L1..L9 are built generically by gen_txr in
* prep/risk_equations.do; any regimen not listed for a line falls into 0 = 'other'.
*
* Regimen codes:
*   2  Thal/Cycl/Dexa    4  Bort/Cycl/Dexa    7  Lena/Dexa       9  Bort/Thal/Dexa
*  31  Bort/Lena/Dexa   49  Carf/Dexa        56  Poma/Dexa      80  Dara/Bort/Dexa
**********

global TXR_L1 "4 7 31"
global TXR_L2 "7 80"
global TXR_L3 "7 49"
global TXR_L4 "7 49"
* L5-L9 unset => all 'other'
