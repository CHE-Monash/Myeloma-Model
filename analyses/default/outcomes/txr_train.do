**********
* Monash Myeloma Model - TXR regimens per line (default analysis, train fit)
*
* Purpose: the 70%-train fit must use the IDENTICAL regimen list as the full fit, so this simply sources
*          the canonical list (txr_full.do) rather than duplicating it. Both files live in the same
*          analysis folder, so this cross-file `do` is HPC-safe (the folder ships together). risk_equations.do
*          loads this via `do "analyses/$analysis/outcomes/txr_$coeffs.do"` when $coeffs == "train".
**********

do "analyses/$analysis/outcomes/txr_full.do"
