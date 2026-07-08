**********
* Monash Myeloma Model - Config template
*
* Purpose: Template for config.do. Copy to config.do and set your own machine-specific paths.
*          config.do is git-ignored, so those paths never reach GitHub. The data-prep scripts
*          (e.g. outcomes/calibrated_transport.do) load it via `capture run "config.do"` from the
*          repository root.
**********

* Repository root (interactively-run scripts cd here; load config.do once per session):
global repo_path    "/path/to/repo"

* Bootstrap / scratch output:
global scratch_path  "/path/to/scratch"

* MRDR data cut (folder name + filename suffix):
global data_cut     "251128"

* Restricted MRDR data (request access via the MRDR Steering Committee):
*   drive_path   = the EpiMAP project dir on your data drive
*   data_path     = the dated working dir within it (MRDR Long.dta, MRDR Long MI.dta, ...)
*   registry_path = the raw registry tables (tbl_*.dta)
global drive_path   "/path/to/EpiMAP/Myeloma"
global data_path     "${drive_path}/Data/${data_cut}"
global registry_path "/path/to/MRDR/Registry data/MRDR Data/2025/${data_cut}_Data"

* Raw PBS Schedule extract (public; not restricted) consumed by prep/extract_pbs_costs.do.
* Point at the tables_as_csv dir of a dated PBS API CSV download. This is the source for the DEFAULT
* (2026-07-01) schedule only; other price years are selected by passing the schedule date as arg 1
* (e.g. `do prep/extract_pbs_costs.do 2025-07-01`), which resolves the source dir by convention.
global pbs_src      "/path/to/PBS-API-CSV-files/tables_as_csv"
