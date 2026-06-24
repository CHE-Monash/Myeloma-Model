**********
* config.example.do — copy to config.do and set your own paths.
*
* config.do is GIT-IGNORED, so machine-specific paths never reach GitHub.
* The data-prep scripts (e.g. outcomes/calibrated_transport.do) load it with
* `capture run "config.do"` from the repository root.
**********

* Base directory for bootstrap / scratch output:
global scratch_dir  "/path/to/scratch"

* Directory holding the MRDR extract ("MRDR Long MI.dta") and its bootstrap/ subfolder.
* Restricted data — not included in the repo; request access via the MRDR Steering Committee:
global data_dir     "/path/to/MRDR/Data/251128"
