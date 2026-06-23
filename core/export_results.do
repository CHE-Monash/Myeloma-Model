**********
* Monash Myeloma Model - Export Results (CSV)
*
* Purpose: Engine-level export of machine-readable CSV results for downstream
*          access (post-processing, dashboards, manuscript drafting). This is
*          the standard hand-off format for results going forward.
*
*          Model-wide: produces the CSVs common to EVERY analysis from a single
*          run. Analysis-specific CSVs and cross-scenario aggregation live under
*          analyses/<name>/ instead (no single run can see other scenarios).
*
* Behaviour:
*   - Runs on the in-memory processed dataset, immediately after process_data,
*     once per arm. Does NOT modify the data in memory (file write +
*     preserve/restore), so the dispatcher's subsequent save is unaffected.
*   - Point-estimate runs only; bootstrap is aggregated downstream (Tier 3),
*     so this exits early when $boot == 1.
*
* Output: $simulated_path/$scenario/
*   bcr_<stub>.csv       BCR distribution at the assessed line ($line)
*   econ_<stub>.csv      Mean cost / QALY / LY (discounted and undiscounted)
*   patients_<stub>.csv  Per-patient key outcomes (flat)
*   where <stub> = ${int}_${line}_${data}_${min_id}_${max_id}
*
* Note: reads the variable names produced by core/process_data.do
*       (cost_total_d, qaly_total_d, BCR_L#, ...), NOT the older names used in
*       core/generate_report.do (cTotald, qTotald). Keep this aligned with
*       process_data.do if those variables are renamed.
**********

capture program drop export_results
program export_results

	* Point-estimate runs only; bootstrap aggregation happens downstream
	if ("$boot" == "1") exit

	local L = $line
	local dir "$simulated_path/$scenario"
	capture mkdir "`dir'"
	local stub "${int}_${line}_${data}_${min_id}_${max_id}"

	di as text "Exporting CSV results -> `dir' (`stub')"

	**********
	* BCR distribution at the assessed line
	**********
	tempname fb
	file open `fb' using "`dir'/bcr_`stub'.csv", write replace
	file write `fb' "scenario,intervention,line,bcr_code,bcr_label,n,pct" _n

	qui count if !missing(BCR_L`L')
	local denom = r(N)

	local bcr_labs `" "CR" "VGPR" "PR" "MR" "SD" "PD" "'
	forval b = 1/6 {
		qui count if BCR_L`L' == `b'
		local nb = r(N)
		local pct = cond(`denom' > 0, 100 * `nb' / `denom', .)
		local lab : word `b' of `bcr_labs'
		local row "$scenario,$int,`L',`b',`lab',`nb',`=trim(string(`pct', "%9.2f"))'"
		file write `fb' `"`row'"' _n
	}
	file close `fb'

	**********
	* Economic summary (means)
	**********
	tempname fe
	file open `fe' using "`dir'/econ_`stub'.csv", write replace
	file write `fe' "scenario,intervention,line,n,mean_cost_disc,mean_cost_undisc,mean_qaly_disc,mean_qaly_undisc,mean_ly" _n

	qui count
	local n = r(N)
	qui summarize cost_total_d, meanonly
	local mcd = r(mean)
	qui summarize cost_total, meanonly
	local mcu = r(mean)
	qui summarize qaly_total_d, meanonly
	local mqd = r(mean)
	qui summarize qaly_total, meanonly
	local mqu = r(mean)
	qui summarize OC_TIME_L, meanonly
	local mly = r(mean) / 12

	local erow "$scenario,$int,`L',`n'"
	local erow "`erow',`=trim(string(`mcd', "%15.2f"))',`=trim(string(`mcu', "%15.2f"))'"
	local erow "`erow',`=trim(string(`mqd', "%15.4f"))',`=trim(string(`mqu', "%15.4f"))'"
	local erow "`erow',`=trim(string(`mly', "%15.4f"))'"
	file write `fe' `"`erow'"' _n
	file close `fe'

	**********
	* Per-patient key outcomes (flat)
	**********
	preserve
	keep ID Male ECOGcc RISS BCR_L`L' TXR_L`L' cost_total cost_total_d qaly_total qaly_total_d OC_TIME OC_TIME_L
	export delimited using "`dir'/patients_`stub'.csv", replace
	restore

end
