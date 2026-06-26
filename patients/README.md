# Population Datasets

## Overview

This directory contains 10 different population realisations for multiple myeloma patients diagnosed in Australia between 1995 and 2040. These are base datasets that can be used across different analyses without modification. Each dataset contains 101,212 patients.

## Files

- `population_1995_2040_1.dta` to `population_1995_2040_10.dta` - Patient population datasets
- `metadata/population_characteristics.dta` - Summary statistics for all populations
- `metadata/generation_log.txt` - Documentation of how populations were generated

## Population Characteristics

Each population dataset contains 101,212 patients with:

| Variable | Description | Range/Categories |
|----------|-------------|------------------|
| **ID** | Patient identifier | 1 to N |
| **Age** | Age at diagnosis | 18-100 years |
| **Male** | Sex indicator | 0=Female, 1=Male |
| **ECOGcc** | ECOG performance status | 0, 1, 2 |
| **RISS** | Revised International Staging | 1, 2, 3 |
| **CMc** | Comorbidity Score | 0, 1, 2 |
| **State** | Initial disease state | Usually 1 |
| **DateDN** | Date of diagnosis | Date format |

## Usage in Analyses

### Option 1: Use specific population
```stata
use "patients/population_1995_2040_3.dta", clear
```

### Option 2: Programmatic selection
```stata
local pop_number = 1
use "patients/population_1995_2040_`pop_number'.dta", clear
```

### Option 3: Parameter-driven selection
```stata
// In analysis file, use global parameter
global PopulationNumber = 5
use "patients/population_1995_2040_${PopulationNumber}.dta", clear
```

## Population Differences

Each population represents a different random realisation of the same underlying demographic and clinical distributions. They can be used for:

1. **Sensitivity Analysis** - Test model robustness across different populations
2. **Uncertainty Quantification** - Multiple runs with different base populations  
3. **Scenario Testing** - Different starting populations for interventions
4. **Validation** - Cross-validation using different population samples

## Integration with Analyses

Analyses can reference these populations in several ways:

1. **Base Model**: Uses populations as-is for standard projections
2. **Intervention Studies**: Applies intervention logic to population data
3. **Comparative Studies**: Runs same population through different treatment scenarios
4. **Sensitivity Analysis**: Tests intervention across multiple populations

## Selection Recommendations

- **Base Model**: Use Population 1 as default
- **Sensitivity Analysis**: Use Populations 1-5 for primary sensitivity
- **Extensive Testing**: Use all 10 populations for comprehensive analysis
- **Quick Testing**: Use Population 1 with reduced sample size

## Maintenance

- Populations should be regenerated when underlying demographic assumptions change
- Version control should track when populations were last updated
- Summary statistics should be regenerated after any population updates
