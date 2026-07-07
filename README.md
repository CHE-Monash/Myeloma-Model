# Monash Myeloma Model v3.0

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Stata](https://img.shields.io/badge/Stata-15.0%2B-red.svg)](https://www.stata.com/) [![DOI](https://img.shields.io/badge/DOI-10.1371%2Fjournal.pone.0308812-blue.svg)](https://doi.org/10.1371/journal.pone.0308812)

A comprehensive discrete-event simulation model for multiple myeloma disease outcomes and treatment pathways, developed through collaboration between Monash University's Centre for Health Economics and Transfusion Research Unit.

## What's New in v3.0

- **Calibrated Transport methods**: out-of-trial outcome prediction (e.g. DVd at L2) via the `transport_dvd` analysis
- **Common Random Numbers (CRN)**: aligned RNG across treatment arms for variance-reduced cost-effectiveness comparisons
- **Standardised CSV exports**: machine-readable result surface for downstream and programmatic use
- **Vectorised engine** (incorporated from v2.1): Mata vector/matrix rewrite for dramatically faster large-scale simulation
- **Rebrand**: project renamed from EpiMAP Myeloma to the Monash Myeloma Model

## Model Overview

The model simulates the complete treatment journey of multiple myeloma patients using **30 evidence-based risk equations** derived from the Australia and New Zealand Myeloma and Related Diseases Registry (MRDR).

### Key Features

- **Patient Characteristics**: Age, sex, ECOG performance score, R-ISS staging & co-morbidity score
- **Treatment Pathways**: Comprehensive modelling of up to 9 lines of therapy
- **Clinical Outcomes**: Best Clinical Response (BCR) and Overall Survival (OS)
- **ASCT Modelling**: Separate pathways for transplant-eligible patients
- **Maintenance Therapy**: Post-induction treatment modelling
- **Parametric Survival Models**: Time-to-event analysis for all outcomes
- **High Performance**: Mata implementation for efficient large-scale simulations

## Quick Start

### Prerequisites

- **Stata 15.0 or higher** (valid licence required)
- Windows, macOS, or Linux operating system

### Installation

``` bash
git clone https://github.com/CHE-Monash/Myeloma-Model.git
cd Myeloma-Model
```

### Basic Usage

Each analysis is driven by its own dispatcher in `analyses/<name>/`. Configure a run by editing the globals at the top of the dispatcher, then run it from Stata:

``` stata
cd "path/to/myeloma-model"
do "analyses/base_model/simulate.do"
```

#### Configuration

The dispatcher's configuration block sets the run via globals. Interactive runs need nothing more; `run.do` and the HPC array jobs instead pass a few of these as optional positional arguments (`boot`, `min_bs`, `max_bs`, and — where an analysis has scenarios — `scenario`), which override the globals of the same name:

| Global | Description | Example |
|---|---|---|
| `analysis` | Analysis name (folder under `analyses/`) | `base_model` |
| `int` | Intervention | `all`, `VRd`, `SoC`, `DVd` |
| `line` | Line of therapy assessed (1–9; `0` = all) | `0` |
| `coeffs` | Coefficient set | `base_model` |
| `data` | Patient data | `population`, `predicted` |
| `min_year` / `max_year` | Diagnosis-year range | `1995` / `2040` |
| `min_id` / `max_id` | Patient ID range | `1` / `101212` |
| `boot` | Bootstrap flag (0/1) | `0` |
| `min_bs` / `max_bs` | Bootstrap iteration range | `1` / `100` |
| `cost_year` | Price year for costs (AUD) | `2025` |
| `drate` | Annual discount rate (PBAC = 5%) | `0.05` |
| `report` | Generate PDF report (0/1) | `0` |
| `scenario` | Scenario label | `Base` |

#### Available Analyses

| Dispatcher | Focus |
|---|---|
| `analyses/base_model/simulate.do` | All regimens — current-practice projections |
| `analyses/vrd_post/simulate.do` | VRd at LoT 1, post-market impact |
| `analyses/transport_dvd/simulate.do` | DVd via Calibrated Transport |

Results are written to `analyses/<analysis>/simulated/`.

### Analysis layout convention

Every analysis folder under `analyses/<name>/` follows the same layout, so the same commands work across analyses:

| File / folder | Role |
|---|---|
| `run.do` | **Analysis runbook** — the full ordered pipeline (prep → simulate → validate) as runnable `do` lines, plus the bootstrap HPC plumbing (the canonical record of the `sbatch`/`rsync` sequence, kept in a trailing block comment). Start here to reproduce an analysis. |
| `simulate.do` | **Dispatcher** — one simulation run, configured by the globals block at the top. Also accepts the optional positionals `boot min_bs max_bs [scenario]` (used by `run.do` and the HPC arrays). |
| `outcomes/txr_<coeffs>.do` | Per-line regimen code lists (`$TXR_L1..L9`) for the coefficient set; `prep/risk_equations.do`'s `gen_txr` builds the regimen dummies from them. |
| `outcomes/` | Any analysis-specific outcome overrides (e.g. `sim_bcr_override.do`) and coefficient generators (e.g. `calibrated_transport.do`). |
| `coefficients/` | Fitted coefficients `coefficients_<coeffs>.mmat` (+ `bootstrap/` for the 500 resamples). |
| `patients/` | The simulation cohort `.dta` (git-ignored; restricted data). |
| `simulated/` | Run outputs (per-scenario subfolders where an analysis has scenarios; `bootstrap/` for resamples). |
| `results/` | Analysis-level CSVs/figures + a `results.md` — the canonical read surface for downstream consumers. |
| `README.md` | Analysis-specific notes (research question, methodology, quirks). |

Shared inputs are built once by `prep/` (multiple imputation, population cohorts) and fitted per-analysis by `prep/risk_equations.do <analysis> <coeffs> …`, which loads that analysis's `outcomes/txr_<coeffs>.do`. The bootstrap steps run on the MASSIVE cluster via the **generic** array jobs in `hpc/` (`multiple_imputation.script`, `risk_equations.script`, `simulate.script`); per-analysis arguments ride on each `sbatch --export=` line, recorded in that analysis's `run.do`.

## Result Exports

Each simulation can emit a PDF report (`$report`) and machine-readable **flat CSV** outputs for downstream use (R/Python post-processing, dashboards, manuscript drafting). CSV export is produced for point-estimate runs and skipped under bootstrap.

- **Engine-level** — `core/export_results.do` writes the CSVs every analysis needs (per-patient summary, BCR distribution, mean cost/QALY/LY) into `simulated/<scenario>/`.
- **Analysis-level** — outputs specific to one analysis, plus cross-scenario aggregation, live under `analyses/<name>/results/` with a `results.md` summarising the key figures — the canonical read surface for downstream consumers.

See each analysis README (e.g. `analyses/transport_dvd/README.md`) for specifics.

## Repository Structure

```         
Myeloma-Model/
├── core/                    # Shared simulation engine
│   ├── load_patients.do     # Patient data loading
│   ├── mata_setup.do        # Mata vector/matrix setup
│   ├── mata_functions.do    # Mata utility functions
│   ├── simulation_engine.do # Discrete-event simulation core
│   ├── rng_slots.do         # Common-random-number slot registry
│   ├── run_pipeline.do      # Shared engine pass used by every analysis
│   ├── process_data.do      # Post-simulation processing
│   ├── export_results.do    # Machine-readable CSV exports
│   ├── generate_report.do   # PDF report
│   ├── outcomes/            # Outcome (risk-equation) modules
│   └── tests/               # Engine verification: unit tests + extreme-value
├── analyses/                # Per-analysis dispatchers, data & results
│   ├── base_model/          # All regimens (current practice)
│   ├── vrd_post/            # VRd LoT 1 post-market
│   ├── transport_dvd/       # DVd Calibrated Transport
│   └── oos/                 # Out-of-sample (70/30) validation — mainstay
├── patients/                # Population cohorts (.dta)
├── scratch/                 # Retired validation: OS/TFI debug diagnostics
├── docs/                    # Technical documentation
├── hpc/                     # MASSIVE M3 cluster scripts
├── README.md
├── CHANGELOG.md
└── LICENSE
```

## Performance Improvements in v3.0

The vectorised implementation provides significant performance benefits:

- **Faster Execution**: Matrix operations process all patients simultaneously
- **Better Memory Efficiency**: Reduced overhead from loop-based processing
- **Scalability**: Handles large cohorts (10,000+ patients) more efficiently
- **Maintainability**: Cleaner code structure with modular vector setup

## Validation

The model has been comprehensively validated:

- **Published Validation**: See Irving et al. (2024) in PLOS ONE
- **Out-of-Sample Testing**: 70/30 split validation with 100 bootstrap iterations
- **Survival Curve Accuracy**: No significant difference in 90% of 120 months post-diagnosis
- **Vectorisation Validation**: Comprehensive test suite confirms identical results to original implementation

## Version History

Access previous versions via Git tags:

- **v3.0**: Current version — Calibrated Transport & CRN methods, standardised CSV exports, rebrand to Monash Myeloma Model (incorporates the v2.1 vectorised engine)
- **v2.1**: Vectorised Mata implementation (backward compatible)
- **v2.0**: Reorganised architecture with extended treatment options
- **v1.0**: Initial public release (August 2024) - `git checkout v1.0`

## Citation

### Model Paper

Irving A, Petrie D, Harris A, Fanning L, Wood EM, Moore E, et al. Developing and validating a discrete-event simulation model of multiple myeloma disease outcomes and treatment pathways using a national clinical registry. *PLOS ONE*. 2024;19(8):e0308812. [doi:10.1371/journal.pone.0308812](https://doi.org/10.1371/journal.pone.0308812)

### Software Citation

``` bibtex
@software{monash_myeloma_model_v3_0,
  title = {Monash Myeloma Model (originally published as EpiMAP Myeloma)},
  author = {Irving, Adam and Petrie, Dennis and Harris, Anthony and Fanning, Laura and 
            Wood, Erica M and Moore, Elizabeth and Wellard, Cameron and Waters, Neil and
            Augustson, Bradley and Cook, Gordon and Gay, Francesca and McCaughan, Georgia and
            Mollee, Peter and Spencer, Andrew and McQuilten, Zoe K},
  version = {3.0},
  year = {2026},
  url = {https://github.com/CHE-Monash/Myeloma-Model},
  doi = {10.1371/journal.pone.0308812},
  institution = {Monash University}
}
```

## Research Team

**Health Economists**: Adam Irving, Dennis Petrie, Anthony Harris, Laura Fanning

**Clinical Experts**: Zoe K McQuilten, Erica M Wood, Bradley Augustson, Gordon Cook, Francesca Gay, Georgia McCaughan, Peter Mollee, Andrew Spencer

**Registry Team**: Elizabeth Moore, Cameron Wellard, Neil Waters

**Consumer Representative**: Andrew Marks

## Data Access

### MRDR Patient Data

For access to genuine patient data from the Australia and New Zealand Myeloma and Related Diseases Registry: - Submit applications to the MRDR Steering Committee - Visit: [mrdr.net.au](https://www.mrdr.net.au/)

## Support

- **Model Questions**: [adam.irving\@monash.edu](mailto:adam.irving@monash.edu)
- **Technical Issues**: [Create an issue](https://github.com/CHE-Monash/Myeloma-Model/issues)
- **Collaboration Enquiries**: Contact the research team

## Licence

This project is licenced under the GNU General Public Licence v3.0 - see the [LICENSE](LICENSE) file for details.

## Related Resources

- [Australia and New Zealand Myeloma and Related Diseases Registry](https://www.mrdr.net.au/)
- [Monash Centre for Health Economics](https://www.monash.edu/business/che)
- [Monash Transfusion Research Unit](https://www.monash.edu/medicine/sphpm/units/transfusion-research)

## Acknowledgements

The Monash Myeloma Model project (originally EpiMAP Myeloma) was supported by Medical Research Future Fund GNT1200706 and GNT2017480 and by National Health and Medical Research Council GNT1189490, GNT2024876 and GNT2036025. We thank patients, clinicians, and research staff at participating centres for their invaluable contributions to the MRDR.

------------------------------------------------------------------------

**Important**: This model is designed for research purposes. Clinical decisions should always involve qualified healthcare professionals and consider individual patient circumstances.
