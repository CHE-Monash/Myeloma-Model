**********
* TXR - out-of-sample (70/30) validation
*
* The OOS analysis validates the base-model structure on held-out patients, so its regimen
* definitions must be IDENTICAL to the base model. Source the base-model spec rather than
* duplicate it (single source of truth, cannot drift). TXR_L* are built by gen_txr in
* prep/risk_equations.do from the $TXR_L* globals this sets.
**********

do "analyses/base_model/outcomes/txr_base_model.do"
