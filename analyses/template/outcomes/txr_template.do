**********
* Monash Myeloma Model - TXR regimens per line (template)
*
* Purpose: declare the per-line regimen code lists (MRDR Regimen codes) this coefficient set models.
*          gen_txr in prep/risk_equations.do builds TXR_L1..L9 from these; any regimen not listed for a
*          line falls into 0 = 'other'. Replace the example lists with your analysis's regimens.
* Usage:   loaded as outcomes/txr_$coeffs.do -- rename to txr_<coeffs>.do to match the $coeffs in
*          simulate.do (here: "template"). Worked example: analyses/transport_dvd/outcomes/txr_transport_dvd.do.
**********

* Regimen code reference (common ones):
*   4  Bort/Cycl/Dexa   5  Bort/Dexa (Vd)   7  Lena/Dexa      31  Bort/Lena/Dexa
*  49  Carf/Dexa       56  Poma/Dexa       80  Dara/Bort/Dexa (DVd)

global TXR_L1 "4 31"
global TXR_L2 "5 7"
global TXR_L3 "7 49"
* L4-L9 unset => all 'other'
