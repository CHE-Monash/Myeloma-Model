**********
* Monash Myeloma Model - TXR regimens per line (default analysis)
*
* Purpose: declare the per-line regimen code lists (MRDR Regimen codes). gen_txr in
*          prep/risk_equations.do builds TXR_L1..L9 from these; any regimen not listed for a line falls
*          into 0 = 'other'. This is the CANONICAL regimen list for the default analysis; the train fit
*          (txr_train.do) sources it, so the in-sample/out-of-sample validation uses the same regimens.
* Notes:   for per-line regimen frequencies before choosing a list, see scratch/regimen_freq.do.
**********

* Regimen codes:
*   2  Thal/Cycl/Dexa    4  Bort/Cycl/Dexa    7  Lena/Dexa       9  Bort/Thal/Dexa
*  31  Bort/Lena/Dexa   49  Carf/Dexa        56  Poma/Dexa      80  Dara/Bort/Dexa

global TXR_L1 "4 7 31"
global TXR_L2 "7 80"
global TXR_L3 "7 49"
global TXR_L4 "7 49 56"
* L5-L9 unset => all 'other'

* Lenalidomide-containing regimen codes: which of the above put the patient on treatment-dose
* lenalidomide, so a line drawn as one of these is eligible to become len-refractory (sim_lenrefr.do).
* Here 7 = Lena/Dexa and 31 = Bort/Lena/Dexa. The 'other' bucket is treated as non-len. This is the
* engine analogue of the fit's Lenalidomide == 1 gate (docs/refractory.md 3.5 / 4).
global LENREFR_regimens "7 31"
