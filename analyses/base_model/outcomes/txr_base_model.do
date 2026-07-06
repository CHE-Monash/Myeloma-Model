**********
* TXR - regimens modelled per line (base model)
*
* Declares only the per-line regimen code lists. TXR_L1..L9 are built generically by gen_txr in
* prep/risk_equations.do; any regimen not listed for a line falls into 0 = 'other'.
* To view per-line regimen frequencies before choosing a list, see scratch/regimen_freq.do.
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
