**********
* config.example.do — copy to config.do and set your own paths.
*
* config.do is GIT-IGNORED, so machine-specific paths never reach GitHub.
* The data-prep scripts (e.g. outcomes/calibrated_transport.do) load it with
* `capture run "config.do"` from the repository root.
**********

* Repository root (interactively-run scripts cd here; load config.do once per session):
global repo_root    "/path/to/repo"

* Bootstrap / scratch output:
global scratch_dir  "/path/to/scratch"

* MRDR data cut (folder name + filename suffix):
global data_cut     "251128"

* Restricted MRDR data (request access via the MRDR Steering Committee):
*   epimap_dir   = the EpiMAP project dir on your data drive
*   data_dir     = the dated working dir within it (MRDR Long.dta, MRDR Long MI.dta, ...)
*   mrdr_raw_dir = the raw registry tables (tbl_*.dta)
global epimap_dir   "/path/to/EpiMAP/Myeloma"
global data_dir     "${epimap_dir}/Data/${data_cut}"
global mrdr_raw_dir "/path/to/MRDR/Registry data/MRDR Data/2025/${data_cut}_Data"
