**********
* DEPRECATED - split into two scripts (June 2026)
*
* This monolith rebuilt the pool every time it ran a size sweep. It is now:
*   build_cohort_pool.do  - builds the reusable L2-entry pool (expensive, once)
*   ce_convergence.do     - draws sizes from the pool and runs the sweep (cheap to repeat)
*
* Run build_cohort_pool.do once, then ce_convergence.do as often as you like.
**********

di as error "generate_ce_cohort.do is deprecated. Use build_cohort_pool.do then ce_convergence.do."
exit 198
