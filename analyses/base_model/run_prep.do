**********
* Monash Myeloma Model - base_model: prep launcher (run_prep.do)
*
* Fits this analysis's inputs. Run from the repository root; uncomment the step
* you need. MRDR data is read via $data_dir (config.do).
* Prep counterpart to the simulation dispatcher (base_model.do).
**********

cap cd "$repo_root"          // ensure repo root (config.do sets $repo_root; load once per session)
capture run "config.do"      // $data_dir (git-ignored; see config.example.do)

* Risk equations -> analyses/base_model/coefficients/coefficients_base_model.mmat
* args: analysis coeffs min_year max_year boot [min_bs max_bs]   (uses outcomes/txr_base_model.do)
* do "prep/risk_equations.do" base_model base_model 1995 2040 0          // deterministic
* do "prep/risk_equations.do" base_model base_model 1995 2040 1 1 500    // bootstrap, iterations 1-500
