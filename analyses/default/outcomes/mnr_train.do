**********
* Monash Myeloma Model - MNR maintenance regimens (default analysis, train fit)
*
* Purpose: the 70%-train fit must use the IDENTICAL maintenance list as the full fit, so this simply
*          sources the canonical list (mnr_full.do) rather than duplicating it. Both files live in the
*          same analysis folder, so this cross-file `do` is HPC-safe (the folder ships together).
*          risk_equations.do loads this via `do "analyses/$analysis/outcomes/mnr_$coeffs.do"` when
*          $coeffs == "train". Mirrors txr_train.do.
**********

do "analyses/$analysis/outcomes/mnr_full.do"
