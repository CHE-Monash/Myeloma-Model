**********
* Monash Myeloma Model - vrd_post: prep launcher (run_prep.do)
*
* Fits this analysis's inputs (the VRd-excluded coefficient set). Run from the
* repository root; uncomment the step you need. MRDR data is read via $data_dir.
**********

cap cd "$repo_root"          // ensure repo root (config.do sets $repo_root; load once per session)
capture run "config.do"      // $data_dir (git-ignored; see config.example.do)

* Risk equations -> analyses/vrd_post/coefficients/coefficients_vrd_post.mmat
* args: analysis coeffs min_year max_year boot [min_bs max_bs]   (uses outcomes/txr_vrd_post.do)
* do "prep/risk_equations.do" vrd_post vrd_post 1995 2040 0          // deterministic
* do "prep/risk_equations.do" vrd_post vrd_post 1995 2040 1 1 500    // bootstrap, iterations 1-500
