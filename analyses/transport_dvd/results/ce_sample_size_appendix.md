# Monte Carlo precision and the simulated cohort size

*Supplementary methods + figure captions. Numbers from `ce_precision_sigma.csv`
(per-scenario σ_pp, deterministic two-arm runs) and `ce_sample_size.csv`
(500-iteration PSA per scenario, N = 50,000). Reusable per analysis: re-run
`ce_precision.do` then `ce_sample_size.do` and substitute the values.*

## Methods text

Reported uncertainty in the incremental QALY should reflect parameter uncertainty,
not the first-order (stochastic) noise of the patient-level simulation. We therefore
chose the simulated cohort size so that the Monte Carlo contribution to the reported
uncertainty was negligible. We did **not** adopt a fixed acceptance threshold for the
Monte Carlo error — none is established in the methodological guidance (NICE DSU TSD 15
recommends judging the simulated sample size by whether the residual error is acceptably
small for the decision, rather than against a fixed ratio) — and instead report the
contribution Monte Carlo error makes to the reported intervals directly.

The Monte Carlo standard deviation of the incremental QALY for a cohort of *N* patients
is σ_pp / √N, where σ_pp is the per-patient standard deviation of the *paired*
incremental QALY and is independent of *N*. Both arms were simulated under common random
numbers, which align each patient's stochastic draws across arms and cancel the shared
lifetime-trajectory variance from the increment (ISPOR-SMDM discrete-event-simulation
good-practice recommendation; Karnon et al. 2012); this reduced σ_pp roughly fourfold
(e.g. 2.1 to 0.51 for the calibrated-transport scenario), and hence the patients required
for a given precision by roughly an order of magnitude.

For each scenario the parameter standard deviation was estimated from its 500-iteration
PSA and corrected for residual Monte Carlo content via the variance decomposition
SD²_param = Var_PSA − σ_pp²/N (O'Hagan, Stevenson & Madan 2007). At N = 50,000 the Monte
Carlo standard deviation of the incremental QALY is at most 0.0025 across the three
scenarios. Expressed as a fraction of each scenario's parameter SD this is 7.2%
(traditional, A), 4.6% (calibrated transport, B) and 2.3% (observed, C); equivalently,
Monte Carlo error inflates the reported 95% intervals by at most **0.26%** (A; 0.10% and
0.03% for B and C) and accounts for under 0.6% of the total variance of each estimate.
Because the incremental-QALY intervals for the calibrated-transport and observed scenarios
already cross zero, the cost-effectiveness conclusions are governed entirely by parameter
uncertainty, and Monte Carlo noise of this magnitude cannot alter them. We therefore
simulated N = 50,000 patients per arm (a single windowed cohort shared across the three
scenarios). Supplementary Figures X–Y document the precision.

## Per-scenario figures

| Scenario | Approach | σ_pp(ΔQALY) | Parameter SD | MC SD @ 50k | MC as % of param SD | 95% interval inflation |
|---|---|---:|---:|---:|---:|---:|
| A_trial | Traditional trial-based | 0.551 | 0.034 | 0.0025 | 7.2% | 0.26% |
| B_transport | Calibrated Transport | 0.514 | 0.050 | 0.0023 | 4.6% | 0.10% |
| C_mrdr | Observed | 0.181 | 0.035 | 0.0008 | 2.3% | 0.03% |

ΔQALY (95% PSA interval): A 0.180 (0.111 to 0.245); B 0.071 (−0.034 to 0.166);
C 0.016 (−0.054 to 0.086).

## Figure captions

**Supplementary Figure X. Convergence of the incremental outcomes with simulated cohort
size, by scenario (DVd vs Vd).** For each scenario, incremental cost, incremental QALY and
the cost per QALY gained plotted against the number of simulated patients on a log scale,
each shown as the full-sample estimate with a ±1.96 × σ_pp/√N Monte Carlo band that narrows
as the cohort grows (the jackknife-equivalent band for the ratio). The dotted line marks the
cohort size used (N = 50,000). After NICE DSU TSD 15 (Davis et al. 2014), Figures 5–7.

**Supplementary Figure Y. Monte Carlo error relative to parameter uncertainty, by scenario.**
The Monte Carlo standard deviation of the incremental QALY expressed as a fraction of each
scenario's parameter (PSA) standard deviation, against simulated cohort size *N* (log scale).
The grey dashed line marks 5% of the parameter SD as a visual reference (not an acceptance
threshold); the dotted vertical line marks N = 50,000. Both arms were simulated under common
random numbers.

## Notes

- σ_pp differs by scenario because the L2 best-response model differs; the traditional
  scenario (A) has the largest σ_pp and the smallest parameter SD, so it is the most
  demanding and sets the binding cohort size.
- The parameter SD is N-independent (the Monte Carlo correction at N = 50,000 is < 1% of
  the variance), so the reported intervals would be essentially unchanged at a larger N;
  a larger cohort would only reduce the already-negligible Monte Carlo contribution.
