**********
* Monash Myeloma Model - Generate Report
*
* Purpose: Build the per-run PDF report (putpdf) of patient, survival, treatment and
*          economic outcomes from the processed dataset. Called by the dispatcher.
* Notes:   Output -> `report_dir'/${int}_${line}_${data}.pdf.
**********


**********
* Paths — simulated output lives under an optional scenario subfolder (empty
* $scenario => simulated_path itself), matching simulate.do / export_results.do.
local sim_out = cond("$scenario" == "", "$simulated_path", "$simulated_path/$scenario")
capture mkdir "`sim_out'"
capture mkdir "`sim_out'/report"
local report_dir "`sim_out'/report"

**********
* TWO-ARM MODE ($report_twoarm == "1"): a landscape comparison report that puts the
* intervention ($int1) and comparator ($int0) outcomes side-by-side. Reloads each arm's
* saved .dta itself, so it does not depend on what is in memory. The single-arm report
* (the else branch below) is unchanged.
**********
if ("$report_twoarm" == "1") {

    local arm1 "$int1"
    local arm0 "$int0"
    local L0   = ${line}                 // decision line: ignore everything before it
    local drate_pct = string(${drate} * 100, "%3.1f")
    * BCR shown at the decision line and the next two lines (e.g. L2 -> L2, L3, L4)
    local rl1 = `L0'
    local rl2 = `L0' + 1
    local rl3 = `L0' + 2

    * ---- Per-arm summaries -> locals suffixed _1 (intervention) / _0 (comparator) ----
    foreach a in 1 0 {
        local arm "`arm`a''"
        capture confirm file "`sim_out'/`arm'_${line}_${data}.dta"
        if _rc {
            di as error "Two-arm report: arm file not found -> `sim_out'/`arm'_${line}_${data}.dta"
            exit 601
        }
        qui use "`sim_out'/`arm'_${line}_${data}.dta", clear
        qui count
        local n_`a' = r(N)

        * Overall survival FROM the decision line: OC_TIME_L = OC_TIME - time-to-L`L0' (months).
        * This is the decision-relevant clock (matches the os_l`L0' benchmark and the discounted QALY/
        * cost basis). Denominator = patients with OC_TIME_L observed (i.e. who reached the line).
        capture confirm variable OC_TIME_L
        if _rc {
            gen double OC_TIME_L = OC_TIME
            if `L0' > 1 {
                capture replace OC_TIME_L = OC_TIME - TSD_L`L0'S
            }
        }
        qui count if !missing(OC_TIME_L)
        local nosl_`a' = r(N)
        qui sum OC_TIME_L, detail
        local osmean_`a' = r(mean)/12
        local osmed_`a'  = r(p50)/12
        local osp25_`a'  = r(p25)/12
        local osp75_`a'  = r(p75)/12
        foreach y in 1 2 3 5 10 {
            qui count if OC_TIME_L/12 >= `y' & !missing(OC_TIME_L)
            local surv`y'_`a' = cond(`nosl_`a'' > 0, 100 * r(N) / `nosl_`a'', .)
        }
        * From-diagnosis mean OS kept as an extra reference row (OC_TIME = months from diagnosis)
        qui sum OC_TIME
        local osmeandx_`a' = r(mean)/12

        * QALYs (discounted): total + components (decision-line treatment + everything after)
        qui sum qaly_total_d
        local qaly_`a' = r(mean)
        local qtxd_`a' = .
        capture confirm variable qaly_txd_L`L0'_d
        if !_rc {
            qui sum qaly_txd_L`L0'_d
            local qtxd_`a' = r(mean)
        }
        local qpost_`a' = .
        capture confirm variable qaly_post_L`L0'_d
        if !_rc {
            qui sum qaly_post_L`L0'_d
            local qpost_`a' = r(mean)
        }

        * Costs (discounted): total + treatment / non-treatment segments + 5-year (undisc.)
        * Cost means are taken over patients with a non-missing TOTAL, so the segment parts add
        * up to the total exactly and match the headline / ICER (which use cost_total_d).
        qui sum cost_total_d
        local cost_`a' = r(mean)
        qui sum cost_tx_d if !missing(cost_total_d)
        local ctx_`a' = r(mean)
        qui sum cost_nt_d if !missing(cost_total_d)
        local cnt_`a' = r(mean)
        foreach seg in hosp mbs emer {
            local cnt_`seg'_`a' = .
            capture confirm variable cost_nt_`seg'_d
            if !_rc {
                qui sum cost_nt_`seg'_d if !missing(cost_total_d)
                local cnt_`seg'_`a' = r(mean)
            }
        }
        local c5yr_`a' = .
        capture confirm variable cost_5yr
        if !_rc {
            qui sum cost_5yr
            local c5yr_`a' = r(mean)
        }

        * Discounted treatment cost by line: txc = mean per patient TREATED at the line;
        * txcpop = mean per ARM patient over the complete-cost population (r(sum)/N_complete, = 0 for
        * uncosted lines). Using the same denominator as the segment/headline makes it sum to the
        * treatment total exactly.
        qui count if !missing(cost_total_d)
        local ncomp_`a' = r(N)
        forvalues L = `L0'/9 {
            local txc`L'_`a' = .
            local txcpop`L'_`a' = 0
            capture confirm variable cost_tx_L`L'_d
            if !_rc {
                qui sum cost_tx_L`L'_d
                local txc`L'_`a' = r(mean)
                qui sum cost_tx_L`L'_d if !missing(cost_total_d)
                local txcpop`L'_`a' = cond(`ncomp_`a'' > 0, r(sum) / `ncomp_`a'', 0)
            }
        }

        * Mean treatment duration (TXD, months among patients treated at the line)
        forvalues L = `L0'/9 {
            local txd`L'_`a' = .
            capture confirm variable TXD_L`L'
            if !_rc {
                qui sum TXD_L`L' if TXD_L`L' > 0 & !missing(TXD_L`L')
                local txd`L'_`a' = r(mean)
            }
        }

        * Mean treatment-free interval (TFI, months to next line among those who progress)
        forvalues L = `L0'/8 {
            local tfi`L'_`a' = .
            capture confirm variable TFI_L`L'
            if !_rc {
                qui sum TFI_L`L' if !missing(TFI_L`L')
                local tfi`L'_`a' = r(mean)
            }
        }

        * BCR distribution at each of the three lines, as % of patients REACHING that line
        foreach ln in `rl1' `rl2' `rl3' {
            capture confirm variable BCR_L`ln'
            if !_rc {
                qui count if BCR_L`ln' >= 1 & BCR_L`ln' < .
                local nresp`ln'_`a' = r(N)
                forvalues b = 1/6 {
                    qui count if BCR_L`ln' == `b'
                    local bcr`ln'_`b'_`a' = cond(`nresp`ln'_`a'' > 0, 100 * r(N) / `nresp`ln'_`a'', .)
                }
            }
        }

        * Pathways: % reaching each line and max-line-reached distribution (% of arm)
        capture drop LOT_MAX
        gen LOT_MAX = 0
        forvalues l = 1/9 {
            qui replace LOT_MAX = `l' if !missing(TXR_L`l')
            qui count if !missing(TXR_L`l')
            local reach`l'_`a' = 100 * r(N) / `n_`a''
        }
        forvalues l = 1/9 {
            qui count if LOT_MAX == `l'
            local maxl`l'_`a' = 100 * r(N) / `n_`a''
        }
    }

    * ---- Incrementals (intervention - comparator) ----
    local d_os   = `osmean_1' - `osmean_0'
    local d_osm  = `osmed_1'  - `osmed_0'
    local d_osdx = `osmeandx_1' - `osmeandx_0'
    local d_qaly = `qaly_1'   - `qaly_0'
    local d_cost = `cost_1'   - `cost_0'
    local icer   = `d_cost' / `d_qaly'

    * ---- Overlaid OS Kaplan-Meier (build the figure before opening the PDF) ----
    set graphics off
    capture mkdir "$simulated_path/report"
    capture mkdir "$simulated_path/report/figures"
    foreach a in 1 0 {
        qui use "`sim_out'/`arm`a''_${line}_${data}.dta", clear
        capture confirm variable OC_TIME_L
        if _rc {
            gen double OC_TIME_L = OC_TIME
            if `L0' > 1 {
                capture replace OC_TIME_L = OC_TIME - TSD_L`L0'S
            }
        }
        keep OC_TIME_L OC_MORT
        gen arm = `a'
        tempfile _os`a'
        qui save `_os`a''
    }
    use `_os1', clear
    append using `_os0'
    stset OC_TIME_L if OC_TIME_L < 240 & !missing(OC_TIME_L), failure(OC_MORT)
    sts graph, by(arm) ///
        xtitle("Months from Line `L0'") ytitle("Survival probability") title("") ///
        ylabel(0(0.2)1, angle(0) format(%3.1f)) xlabel(0(24)240) ///
        legend(order(1 "`arm0'" 2 "`arm1'") rows(1) pos(6)) ///
        graphregion(color(white)) name(os_cmp, replace)
    graph export "$simulated_path/report/figures/os_compare.png", replace width(1600)

    * ================= BUILD PDF (landscape) =================
    capture putpdf clear
    putpdf begin, landscape

    putpdf paragraph, halign(center)
    putpdf text ("Monash Myeloma Model v3.0"), bold font(,18)
    putpdf paragraph, halign(center)
    putpdf text ("Two-Arm Comparison Report"), bold font(,15)
    putpdf paragraph, halign(center)
    putpdf text ("`arm1' (intervention) vs `arm0' (comparator) - Line `L0' onwards"), font(,12)

    * Run settings
    putpdf paragraph
    local n_fmt = trim(string(`n_1', "%12.0fc"))
    putpdf table s = (6,2), border(all)
    putpdf table s(1,1) = ("Setting"), bold
    putpdf table s(1,2) = ("Value"), bold
    putpdf table s(2,1) = ("Analysis")
    putpdf table s(2,2) = ("$analysis")
    putpdf table s(3,1) = ("Line / Data")
    putpdf table s(3,2) = ("${line} / ${data}")
    putpdf table s(4,1) = ("Scenario")
    putpdf table s(4,2) = ("$scenario")
    putpdf table s(5,1) = ("N per arm")
    putpdf table s(5,2) = ("`n_fmt'")
    putpdf table s(6,1) = ("Report date")
    putpdf table s(6,2) = ("`c(current_date)'")

    * ---- Headline outcomes: Outcome | Intervention | Comparator | Incremental ----
    putpdf paragraph
    putpdf text ("Outcomes per patient"), bold font(,13)
    putpdf table h = (6,4), border(all)
    putpdf table h(1,1) = ("Outcome"), bold
    putpdf table h(1,2) = ("`arm1'"), bold
    putpdf table h(1,3) = ("`arm0'"), bold
    putpdf table h(1,4) = ("Incremental"), bold
    local v1 = string(`osmean_1', "%4.2f")
    local v0 = string(`osmean_0', "%4.2f")
    local vd = string(`d_os', "%4.2f")
    putpdf table h(2,1) = ("Mean OS from L`L0' (years, undisc.)")
    putpdf table h(2,2) = ("`v1'")
    putpdf table h(2,3) = ("`v0'")
    putpdf table h(2,4) = ("`vd'")
    local v1 = string(`osmed_1', "%4.2f")
    local v0 = string(`osmed_0', "%4.2f")
    local vd = string(`d_osm', "%4.2f")
    putpdf table h(3,1) = ("Median OS from L`L0' (years, undisc.)")
    putpdf table h(3,2) = ("`v1'")
    putpdf table h(3,3) = ("`v0'")
    putpdf table h(3,4) = ("`vd'")
    local v1 = string(`qaly_1', "%5.2f")
    local v0 = string(`qaly_0', "%5.2f")
    local vd = string(`d_qaly', "%5.2f")
    putpdf table h(4,1) = ("Mean QALYs (disc.)")
    putpdf table h(4,2) = ("`v1'")
    putpdf table h(4,3) = ("`v0'")
    putpdf table h(4,4) = ("`vd'")
    local v1 = trim(string(`cost_1', "%12.0fc"))
    local v0 = trim(string(`cost_0', "%12.0fc"))
    local vd = trim(string(`d_cost', "%12.0fc"))
    putpdf table h(5,1) = ("Mean total cost (disc., AUD)")
    putpdf table h(5,2) = ("$`v1'")
    putpdf table h(5,3) = ("$`v0'")
    putpdf table h(5,4) = ("$`vd'")
    local vi = trim(string(`icer', "%12.0fc"))
    putpdf table h(6,1) = ("ICER (disc., \$/QALY)"), bold
    putpdf table h(6,4) = ("$`vi'"), bold
    putpdf paragraph
    putpdf text ("OS is undiscounted (descriptive); QALYs, costs and the ICER are discounted at `drate_pct'%."), font(,9) italic

    * ---- Treatment pathways ----
    putpdf pagebreak
    putpdf paragraph
    putpdf text ("Treatment pathways (Line `L0' onwards)"), bold font(,14)
    local nprows = 9 - `L0' + 2
    putpdf table pw = (`nprows',5), border(all)
    putpdf table pw(1,1) = ("Line"), bold
    putpdf table pw(1,2) = ("Reached % (`arm1')"), bold
    putpdf table pw(1,3) = ("Reached % (`arm0')"), bold
    putpdf table pw(1,4) = ("Max line % (`arm1')"), bold
    putpdf table pw(1,5) = ("Max line % (`arm0')"), bold
    local prow 1
    forvalues l = `L0'/9 {
        local prow = `prow' + 1
        local r1 = string(`reach`l'_1', "%4.1f")
        local r0 = string(`reach`l'_0', "%4.1f")
        local m1 = string(`maxl`l'_1', "%4.1f")
        local m0 = string(`maxl`l'_0', "%4.1f")
        putpdf table pw(`prow',1) = ("L`l'")
        putpdf table pw(`prow',2) = ("`r1'%")
        putpdf table pw(`prow',3) = ("`r0'%")
        putpdf table pw(`prow',4) = ("`m1'%")
        putpdf table pw(`prow',5) = ("`m0'%")
    }
    putpdf paragraph
    putpdf text ("Reached % = of the arm; Max line % = of the arm whose last treated line is that line."), font(,9) italic

    * ---- Overall survival ----
    putpdf pagebreak
    putpdf paragraph
    putpdf text ("Overall survival (from Line `L0')"), bold font(,14)
    putpdf table os = (9,3), border(all)
    putpdf table os(1,1) = ("Statistic"), bold
    putpdf table os(1,2) = ("`arm1'"), bold
    putpdf table os(1,3) = ("`arm0'"), bold
    local v1 = string(`osmean_1', "%4.2f")
    local v0 = string(`osmean_0', "%4.2f")
    putpdf table os(2,1) = ("Mean OS from L`L0' (years)")
    putpdf table os(2,2) = ("`v1'")
    putpdf table os(2,3) = ("`v0'")
    local v1 = string(`osmed_1', "%4.2f") + " [" + string(`osp25_1', "%4.2f") + "-" + string(`osp75_1', "%4.2f") + "]"
    local v0 = string(`osmed_0', "%4.2f") + " [" + string(`osp25_0', "%4.2f") + "-" + string(`osp75_0', "%4.2f") + "]"
    putpdf table os(3,1) = ("Median OS from L`L0' [IQR], years")
    putpdf table os(3,2) = ("`v1'")
    putpdf table os(3,3) = ("`v0'")
    local v1 = string(`osmeandx_1', "%4.2f")
    local v0 = string(`osmeandx_0', "%4.2f")
    putpdf table os(4,1) = ("Mean OS from diagnosis (years)")
    putpdf table os(4,2) = ("`v1'")
    putpdf table os(4,3) = ("`v0'")
    local orow 4
    foreach y in 1 2 3 5 10 {
        local orow = `orow' + 1
        local v1 = string(`surv`y'_1', "%4.1f")
        local v0 = string(`surv`y'_0', "%4.1f")
        putpdf table os(`orow',1) = ("`y'-year survival (from L`L0')")
        putpdf table os(`orow',2) = ("`v1'%")
        putpdf table os(`orow',3) = ("`v0'%")
    }
    putpdf paragraph
    putpdf text ("OS is measured from Line `L0' (the decision point), matching the discounted QALY/cost basis; mean OS from diagnosis is shown as a reference row. From-L`L0' survival denominator = patients reaching L`L0'."), font(,9) italic
    putpdf paragraph
    putpdf image "$simulated_path/report/figures/os_compare.png", width(7)

    * ---- Best clinical response at the decision line, L+1 and L+2 ----
    putpdf pagebreak
    putpdf paragraph
    putpdf text ("Best clinical response by line (% of patients reaching each line)"), bold font(,14)
    putpdf table bc = (7,7), border(all)
    putpdf table bc(1,1) = ("Response"), bold
    putpdf table bc(1,2) = ("L`rl1' `arm1'"), bold
    putpdf table bc(1,3) = ("L`rl1' `arm0'"), bold
    putpdf table bc(1,4) = ("L`rl2' `arm1'"), bold
    putpdf table bc(1,5) = ("L`rl2' `arm0'"), bold
    putpdf table bc(1,6) = ("L`rl3' `arm1'"), bold
    putpdf table bc(1,7) = ("L`rl3' `arm0'"), bold
    local rlab1 "CR"
    local rlab2 "VGPR"
    local rlab3 "PR"
    local rlab4 "MR"
    local rlab5 "SD"
    local rlab6 "PD"
    forvalues b = 1/6 {
        local row = `b' + 1
        putpdf table bc(`row',1) = ("`rlab`b''")
        local col 1
        foreach ln in `rl1' `rl2' `rl3' {
            local col1 = `col' + 1
            local col2 = `col' + 2
            local col = `col' + 2
            local p1 = string(`bcr`ln'_`b'_1', "%4.1f")
            local p0 = string(`bcr`ln'_`b'_0', "%4.1f")
            putpdf table bc(`row',`col1') = ("`p1'%")
            putpdf table bc(`row',`col2') = ("`p0'%")
        }
    }

    * ---- QALYs (discounted) ----
    putpdf pagebreak
    putpdf paragraph
    putpdf text ("Quality-adjusted life years (discounted at `drate_pct'%)"), bold font(,14)
    putpdf table q = (4,4), border(all)
    putpdf table q(1,1) = ("QALY component"), bold
    putpdf table q(1,2) = ("`arm1'"), bold
    putpdf table q(1,3) = ("`arm0'"), bold
    putpdf table q(1,4) = ("Incremental"), bold
    local v1 = string(`qtxd_1', "%5.3f")
    local v0 = string(`qtxd_0', "%5.3f")
    local vd = string(`qtxd_1' - `qtxd_0', "%5.3f")
    putpdf table q(2,1) = ("Line `L0' treatment")
    putpdf table q(2,2) = ("`v1'")
    putpdf table q(2,3) = ("`v0'")
    putpdf table q(2,4) = ("`vd'")
    local v1 = string(`qpost_1', "%5.3f")
    local v0 = string(`qpost_0', "%5.3f")
    local vd = string(`qpost_1' - `qpost_0', "%5.3f")
    putpdf table q(3,1) = ("Post-Line `L0' (all later lines)")
    putpdf table q(3,2) = ("`v1'")
    putpdf table q(3,3) = ("`v0'")
    putpdf table q(3,4) = ("`vd'")
    local v1 = string(`qaly_1', "%5.2f")
    local v0 = string(`qaly_0', "%5.2f")
    local vd = string(`d_qaly', "%5.2f")
    putpdf table q(4,1) = ("Total"), bold
    putpdf table q(4,2) = ("`v1'"), bold
    putpdf table q(4,3) = ("`v0'"), bold
    putpdf table q(4,4) = ("`vd'"), bold
    putpdf paragraph
    putpdf text ("The model does not decompose QALYs by line beyond Line `L0'; later lines are pooled into Post-Line `L0'."), font(,9) italic

    * ---- Costs (discounted) ----
    putpdf pagebreak
    putpdf paragraph
    putpdf text ("Costs (discounted at `drate_pct'%)"), bold font(,14)

    * Cost by segment
    putpdf paragraph
    putpdf text ("Cost by segment (mean, AUD)"), bold font(,13)
    putpdf table cs = (7,4), border(all)
    putpdf table cs(1,1) = ("Segment"), bold
    putpdf table cs(1,2) = ("`arm1'"), bold
    putpdf table cs(1,3) = ("`arm0'"), bold
    putpdf table cs(1,4) = ("Incremental"), bold
    local v1 = trim(string(`ctx_1', "%12.0fc"))
    local v0 = trim(string(`ctx_0', "%12.0fc"))
    local vd = trim(string(`ctx_1' - `ctx_0', "%12.0fc"))
    putpdf table cs(2,1) = ("Treatment (regimens)")
    putpdf table cs(2,2) = ("$`v1'")
    putpdf table cs(2,3) = ("$`v0'")
    putpdf table cs(2,4) = ("$`vd'")
    local v1 = trim(string(`cnt_1', "%12.0fc"))
    local v0 = trim(string(`cnt_0', "%12.0fc"))
    local vd = trim(string(`cnt_1' - `cnt_0', "%12.0fc"))
    putpdf table cs(3,1) = ("Non-treatment (total)")
    putpdf table cs(3,2) = ("$`v1'")
    putpdf table cs(3,3) = ("$`v0'")
    putpdf table cs(3,4) = ("$`vd'")
    local segrow 3
    local slab_hosp "  Hospital"
    local slab_mbs  "  MBS"
    local slab_emer "  Emergency"
    foreach seg in hosp mbs emer {
        local segrow = `segrow' + 1
        local v1 = trim(string(`cnt_`seg'_1', "%12.0fc"))
        local v0 = trim(string(`cnt_`seg'_0', "%12.0fc"))
        local vd = trim(string(`cnt_`seg'_1' - `cnt_`seg'_0', "%12.0fc"))
        putpdf table cs(`segrow',1) = ("`slab_`seg''")
        putpdf table cs(`segrow',2) = ("$`v1'")
        putpdf table cs(`segrow',3) = ("$`v0'")
        putpdf table cs(`segrow',4) = ("$`vd'")
    }
    local v1 = trim(string(`cost_1', "%12.0fc"))
    local v0 = trim(string(`cost_0', "%12.0fc"))
    local vd = trim(string(`d_cost', "%12.0fc"))
    putpdf table cs(7,1) = ("Total"), bold
    putpdf table cs(7,2) = ("$`v1'"), bold
    putpdf table cs(7,3) = ("$`v0'"), bold
    putpdf table cs(7,4) = ("$`vd'"), bold

    * Treatment cost by line - per patient treated at the line
    putpdf paragraph
    putpdf text ("Treatment cost by line - per patient treated at the line (mean, AUD)"), bold font(,13)
    local nclrows = 9 - `L0' + 2
    putpdf table cl = (`nclrows',4), border(all)
    putpdf table cl(1,1) = ("Line"), bold
    putpdf table cl(1,2) = ("`arm1'"), bold
    putpdf table cl(1,3) = ("`arm0'"), bold
    putpdf table cl(1,4) = ("Incremental"), bold
    local crow 1
    forvalues L = `L0'/9 {
        local crow = `crow' + 1
        local v1 = cond(missing(`txc`L'_1'), "n/a", "$" + trim(string(`txc`L'_1', "%12.0fc")))
        local v0 = cond(missing(`txc`L'_0'), "n/a", "$" + trim(string(`txc`L'_0', "%12.0fc")))
        local vd = cond(missing(`txc`L'_1') | missing(`txc`L'_0'), "n/a", "$" + trim(string(`txc`L'_1' - `txc`L'_0', "%12.0fc")))
        putpdf table cl(`crow',1) = ("L`L'")
        putpdf table cl(`crow',2) = ("`v1'")
        putpdf table cl(`crow',3) = ("`v0'")
        putpdf table cl(`crow',4) = ("`vd'")
    }
    putpdf paragraph
    putpdf text ("n/a = line not costed (e.g. a forced/unmodelled regimen at the decision line has no drug cost)."), font(,9) italic

    * Treatment cost by line - per ARM patient (population; = per-treated cost x proportion reaching)
    putpdf paragraph
    putpdf text ("Treatment cost by line - per arm patient (population mean, AUD; sums to the treatment total)"), bold font(,13)
    putpdf table pl = (`nclrows',4), border(all)
    putpdf table pl(1,1) = ("Line"), bold
    putpdf table pl(1,2) = ("`arm1'"), bold
    putpdf table pl(1,3) = ("`arm0'"), bold
    putpdf table pl(1,4) = ("Incremental"), bold
    local plrow 1
    forvalues L = `L0'/9 {
        local plrow = `plrow' + 1
        local v1 = trim(string(`txcpop`L'_1', "%12.0fc"))
        local v0 = trim(string(`txcpop`L'_0', "%12.0fc"))
        local vd = trim(string(`txcpop`L'_1' - `txcpop`L'_0', "%12.0fc"))
        putpdf table pl(`plrow',1) = ("L`L'")
        putpdf table pl(`plrow',2) = ("$`v1'")
        putpdf table pl(`plrow',3) = ("$`v0'")
        putpdf table pl(`plrow',4) = ("$`vd'")
    }
    putpdf paragraph
    putpdf text ("Population cost = per-treated cost x the fraction of the arm reaching that line, so it reflects both dynamics at once (e.g. an intervention that treats fewer patients downstream but each for longer)."), font(,9) italic

    * Diagnosis-to-5-year cost (undiscounted) - only if the analysis populates it
    if !missing(`c5yr_1') | !missing(`c5yr_0') {
        putpdf paragraph
        putpdf text ("Diagnosis-to-5-year cost (undiscounted, mean AUD)"), bold font(,13)
        putpdf table c5 = (2,3), border(all)
        putpdf table c5(1,1) = ("`arm1'"), bold
        putpdf table c5(1,2) = ("`arm0'"), bold
        putpdf table c5(1,3) = ("Incremental"), bold
        local v1 = trim(string(`c5yr_1', "%12.0fc"))
        local v0 = trim(string(`c5yr_0', "%12.0fc"))
        local vd = trim(string(`c5yr_1' - `c5yr_0', "%12.0fc"))
        putpdf table c5(2,1) = ("$`v1'")
        putpdf table c5(2,2) = ("$`v0'")
        putpdf table c5(2,3) = ("$`vd'")
    }

    * ---- Treatment duration & treatment-free interval by line (explains the per-line cost gap) ----
    putpdf pagebreak
    putpdf paragraph
    putpdf text ("Treatment duration and treatment-free interval by line"), bold font(,14)

    putpdf paragraph
    putpdf text ("Mean treatment duration (TXD, months among patients treated at the line)"), bold font(,13)
    local ntdrows = 9 - `L0' + 2
    putpdf table td = (`ntdrows',4), border(all)
    putpdf table td(1,1) = ("Line"), bold
    putpdf table td(1,2) = ("`arm1'"), bold
    putpdf table td(1,3) = ("`arm0'"), bold
    putpdf table td(1,4) = ("Incremental"), bold
    local trow 1
    forvalues L = `L0'/9 {
        local trow = `trow' + 1
        local v1 = cond(missing(`txd`L'_1'), "n/a", string(`txd`L'_1', "%5.1f"))
        local v0 = cond(missing(`txd`L'_0'), "n/a", string(`txd`L'_0', "%5.1f"))
        local vd = cond(missing(`txd`L'_1') | missing(`txd`L'_0'), "n/a", string(`txd`L'_1' - `txd`L'_0', "%5.1f"))
        putpdf table td(`trow',1) = ("L`L'")
        putpdf table td(`trow',2) = ("`v1'")
        putpdf table td(`trow',3) = ("`v0'")
        putpdf table td(`trow',4) = ("`vd'")
    }

    putpdf paragraph
    putpdf text ("Mean treatment-free interval (TFI, months to next line among those who progress)"), bold font(,13)
    local ntfrows = 8 - `L0' + 2
    putpdf table tf = (`ntfrows',4), border(all)
    putpdf table tf(1,1) = ("Interval"), bold
    putpdf table tf(1,2) = ("`arm1'"), bold
    putpdf table tf(1,3) = ("`arm0'"), bold
    putpdf table tf(1,4) = ("Incremental"), bold
    local frow 1
    forvalues L = `L0'/8 {
        local frow = `frow' + 1
        local nx = `L' + 1
        local v1 = cond(missing(`tfi`L'_1'), "n/a", string(`tfi`L'_1', "%5.1f"))
        local v0 = cond(missing(`tfi`L'_0'), "n/a", string(`tfi`L'_0', "%5.1f"))
        local vd = cond(missing(`tfi`L'_1') | missing(`tfi`L'_0'), "n/a", string(`tfi`L'_1' - `tfi`L'_0', "%5.1f"))
        putpdf table tf(`frow',1) = ("L`L' to L`nx'")
        putpdf table tf(`frow',2) = ("`v1'")
        putpdf table tf(`frow',3) = ("`v0'")
        putpdf table tf(`frow',4) = ("`vd'")
    }
    putpdf paragraph
    putpdf text ("Longer TXD at a line reflects better response carried forward (higher BCR -> longer time on treatment), which raises per-patient treatment cost at that line even when fewer patients reach it."), font(,9) italic

    set graphics on
    putpdf save "`report_dir'/compare_${line}_${data}.pdf", replace
    di as text "Two-arm comparison report -> `report_dir'/compare_${line}_${data}.pdf"
}
else {

* Load Data
qui use "`sim_out'/${int}_${line}_${data}.dta", clear


**********
* Start PDF
**********

capture putpdf clear
set graphics off
putpdf begin

// Title page
putpdf paragraph, halign(center)
putpdf text ("Monash Myeloma Model v3.0"), bold font(,18)

putpdf paragraph
putpdf text ("Simulation Report"), bold font(,16)

putpdf table settings = (8, 2), border(all)
putpdf table settings(1,1) = ("Setting"), bold
putpdf table settings(1,2) = ("Value"), bold
putpdf table settings(2,1) = ("Analysis")
putpdf table settings(2,2) = ("$analysis")
putpdf table settings(3,1) = ("Intervention")
putpdf table settings(3,2) = ("$int")
putpdf table settings(4,1) = ("Line")
putpdf table settings(4,2) = ("$line")
putpdf table settings(5,1) = ("Data")
putpdf table settings(5,2) = ("$data")
putpdf table settings(6,1) = ("Patient IDs")
putpdf table settings(6,2) = ("$min_id to $max_id")
putpdf table settings(7,1) = ("Scenario")
putpdf table settings(7,2) = ("$scenario")
putpdf table settings(8,1) = ("Report Date")
putpdf table settings(8,2) = ("`c(current_date)'")

**********
* Patients
**********

putpdf paragraph
putpdf text ("Patients"), bold font(,16)

// Sample size
quietly count
local total_n = string(r(N), "%9.0fc")

// Sex
quietly count if Male == 1
local male_n = string(r(N), "%9.0fc")
local male_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("Sample Size: `total_n'"), linebreak
putpdf text ("Male: `male_n' (`male_pct'%)")

// Age
putpdf paragraph
putpdf text ("Age at Diagnosis"), bold

quietly summarize Age_DN, detail
local mean_age = string(r(mean), "%4.1f")
local sd_age = string(r(sd), "%4.1f")
local median_age = string(r(p50), "%4.1f")
local p25_age = string(r(p25), "%4.1f")
local p75_age = string(r(p75), "%4.1f")

putpdf table age_sum = (3, 2), border(all)
putpdf table age_sum(1,1) = ("Statistic"), bold
putpdf table age_sum(1,2) = ("Years"), bold
putpdf table age_sum(2,1) = ("Mean (SD)")
putpdf table age_sum(2,2) = ("`mean_age' (`sd_age')")
putpdf table age_sum(3,1) = ("Median [IQR]")
putpdf table age_sum(3,2) = ("`median_age' [`p25_age' - `p75_age']")

putpdf paragraph

// ECOG
quietly count if ECOGcc == 0
local ecog0_n = string(r(N), "%9.0fc")
local ecog0_pct = string(100*r(N)/_N, "%4.1f")
quietly count if ECOGcc == 1
local ecog1_n = string(r(N), "%9.0fc")
local ecog1_pct = string(100*r(N)/_N, "%4.1f")
quietly count if ECOGcc == 2
local ecog2_n = string(r(N), "%9.0fc")
local ecog2_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("ECOG"), bold

putpdf table ecog_tbl = (4, 3), border(all)
putpdf table ecog_tbl(1,1) = ("ECOG Score"), bold
putpdf table ecog_tbl(1,2) = ("N"), bold
putpdf table ecog_tbl(1,3) = ("%"), bold
putpdf table ecog_tbl(2,1) = ("0")
putpdf table ecog_tbl(2,2) = ("`ecog0_n'")
putpdf table ecog_tbl(2,3) = ("`ecog0_pct'")
putpdf table ecog_tbl(3,1) = ("1")
putpdf table ecog_tbl(3,2) = ("`ecog1_n'")
putpdf table ecog_tbl(3,3) = ("`ecog1_pct'")
putpdf table ecog_tbl(4,1) = ("2+")
putpdf table ecog_tbl(4,2) = ("`ecog2_n'")
putpdf table ecog_tbl(4,3) = ("`ecog2_pct'")


// R-ISS
quietly count if RISS == 1
local riss1_n = string(r(N), "%9.0fc")
local riss1_pct = string(100*r(N)/_N, "%4.1f")
quietly count if RISS == 2
local riss2_n = string(r(N), "%9.0fc")
local riss2_pct = string(100*r(N)/_N, "%4.1f")
quietly count if RISS == 3
local riss3_n = string(r(N), "%9.0fc")
local riss3_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("R-ISS"), bold

putpdf table riss_tbl = (4, 3), border(all)
putpdf table riss_tbl(1,1) = ("R-ISS Stage"), bold
putpdf table riss_tbl(1,2) = ("N"), bold
putpdf table riss_tbl(1,3) = ("%"), bold
putpdf table riss_tbl(2,1) = ("I")
putpdf table riss_tbl(2,2) = ("`riss1_n'")
putpdf table riss_tbl(2,3) = ("`riss1_pct'")
putpdf table riss_tbl(3,1) = ("II")
putpdf table riss_tbl(3,2) = ("`riss2_n'")
putpdf table riss_tbl(3,3) = ("`riss2_pct'")
putpdf table riss_tbl(4,1) = ("III")
putpdf table riss_tbl(4,2) = ("`riss3_n'")
putpdf table riss_tbl(4,3) = ("`riss3_pct'")

putpdf pagebreak

**********
* Treatment
**********

putpdf paragraph
putpdf text ("Treatments"), bold font(,16)

// Line 1 Regimen × BCR Cross-tabulation
quietly count if TXR_L1 < .
local l1_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL1_TXR)
mata: st_numscalar("n_regimens", cols(oL1_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
    // Assign name based on code
    if `code' == 0  local reg_`i'_name "Other"
    if `code' == 2  local reg_`i'_name "TCd"
    if `code' == 4  local reg_`i'_name "VCd"
    if `code' == 7  local reg_`i'_name "Rd"
    if `code' == 9  local reg_`i'_name "VTd"
    if `code' == 31 local reg_`i'_name "VRd"
    if `code' == 49 local reg_`i'_name "Kd"
    if `code' == 56 local reg_`i'_name "Pd"
    if `code' == 80 local reg_`i'_name "DVd"
    
    // Count N for this regimen
    quietly count if TXR_L1 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l1_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD)
    forval b = 1/6 {
        quietly count if TXR_L1 == `code' & BCR_L1 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 8 rows (header + N + 6 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l1 = (8, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l1(1,1) = ("Line 1"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l1(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l1(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l1(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR" "VGPR" "PR" "MR" "SD" "PD" "'
forval b = 1/6 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l1(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l1(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

// ASCT BCR Table with Patient Characteristics
quietly count if TXR_L1 < .
local l1_total = r(N)

quietly count if SCT_L1 == 1
local sct_n = r(N)
local sct_n_fmt = string(r(N), "%9.0fc")
local sct_pct = string(100*r(N)/`l1_total', "%4.1f")

quietly count if SCT_L1 == 0
local nosct_n = r(N)
local nosct_n_fmt = string(r(N), "%9.0fc")
local nosct_pct = string(100*r(N)/`l1_total', "%4.1f")

// BCR within ASCT patients (BCR_SCT: 1=CR, 2=VGPR, 3=PR, 4=MR)
forval b = 1/4 {
    quietly count if SCT_L1 == 1 & BCR_SCT == `b'
    if `sct_n' > 0 {
        local sct_bcr`b' = string(100*r(N)/`sct_n', "%4.1f")
    }
    else {
        local sct_bcr`b' = "—"
    }
}

// Age breakdown for ASCT patients
quietly count if SCT_L1 == 1 & Age_L1S < 65
local sct_age1 = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 1 & Age_L1S >= 65 & Age_L1S < 70
local sct_age2 = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 1 & Age_L1S >= 70 & Age_L1S < 75
local sct_age3 = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 1 & Age_L1S >= 75
local sct_age4 = string(100*r(N)/`sct_n', "%4.1f")

// Age breakdown for No ASCT patients
quietly count if SCT_L1 == 0 & Age_L1S < 65
local nosct_age1 = string(100*r(N)/`nosct_n', "%4.1f")
quietly count if SCT_L1 == 0 & Age_L1S >= 65 & Age_L1S < 70
local nosct_age2 = string(100*r(N)/`nosct_n', "%4.1f")
quietly count if SCT_L1 == 0 & Age_L1S >= 70 & Age_L1S < 75
local nosct_age3 = string(100*r(N)/`nosct_n', "%4.1f")
quietly count if SCT_L1 == 0 & Age_L1S >= 75
local nosct_age4 = string(100*r(N)/`nosct_n', "%4.1f")

// BCR L1 CR/VGPR for ASCT vs No ASCT
quietly count if SCT_L1 == 1 & (BCR_L1 == 1 | BCR_L1 == 2)
local sct_crvgpr = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 0 & (BCR_L1 == 1 | BCR_L1 == 2)
local nosct_crvgpr = string(100*r(N)/`nosct_n', "%4.1f")

// Create table: 6 rows × 5 columns
putpdf table txr_bcr_sct = (6, 5), border(all)

// Header row
putpdf table txr_bcr_sct(1,1) = ("ASCT"), bold
putpdf table txr_bcr_sct(1,2) = (""), bold
putpdf table txr_bcr_sct(1,3) = (""), bold
putpdf table txr_bcr_sct(1,4) = ("ASCT"), bold
putpdf table txr_bcr_sct(1,5) = ("No ASCT"), bold

// Row 2: N / Age < 65
putpdf table txr_bcr_sct(2,1) = ("N"), bold
putpdf table txr_bcr_sct(2,2) = ("`sct_n_fmt' (`sct_pct'%)")
putpdf table txr_bcr_sct(2,3) = ("Age < 65")
putpdf table txr_bcr_sct(2,4) = ("`sct_age1'%")
putpdf table txr_bcr_sct(2,5) = ("`nosct_age1'%")

// Row 3: CR / Age 65-69
putpdf table txr_bcr_sct(3,1) = ("CR")
putpdf table txr_bcr_sct(3,2) = ("`sct_bcr1'%")
putpdf table txr_bcr_sct(3,3) = ("Age >= 65 & < 70")
putpdf table txr_bcr_sct(3,4) = ("`sct_age2'%")
putpdf table txr_bcr_sct(3,5) = ("`nosct_age2'%")

// Row 4: VGPR / Age 70-74
putpdf table txr_bcr_sct(4,1) = ("VGPR")
putpdf table txr_bcr_sct(4,2) = ("`sct_bcr2'%")
putpdf table txr_bcr_sct(4,3) = ("Age >= 70 & < 75")
putpdf table txr_bcr_sct(4,4) = ("`sct_age3'%")
putpdf table txr_bcr_sct(4,5) = ("`nosct_age3'%")

// Row 5: PR / Age 75+
putpdf table txr_bcr_sct(5,1) = ("PR")
putpdf table txr_bcr_sct(5,2) = ("`sct_bcr3'%")
putpdf table txr_bcr_sct(5,3) = ("Age >= 75")
putpdf table txr_bcr_sct(5,4) = ("`sct_age4'%")
putpdf table txr_bcr_sct(5,5) = ("`nosct_age4'%")

// Row 6: MR / BCR L1 CR/VGPR
putpdf table txr_bcr_sct(6,1) = ("MR")
putpdf table txr_bcr_sct(6,2) = ("`sct_bcr4'%")
putpdf table txr_bcr_sct(6,3) = ("BCR L1 CR/VGPR")
putpdf table txr_bcr_sct(6,4) = ("`sct_crvgpr'%")
putpdf table txr_bcr_sct(6,5) = ("`nosct_crvgpr'%")

// Line 2 Regimen × BCR Cross-tabulation
quietly count if TXR_L2 < .
local l2_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL2_TXR)
mata: st_numscalar("n_regimens", cols(oL2_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
    // Assign name based on code
    if `code' == 0  local reg_`i'_name "Other"
    if `code' == 2  local reg_`i'_name "TCd"
    if `code' == 4  local reg_`i'_name "VCd"
    if `code' == 7  local reg_`i'_name "Rd"
    if `code' == 9  local reg_`i'_name "VTd"
    if `code' == 31 local reg_`i'_name "VRd"
    if `code' == 49 local reg_`i'_name "Kd"
    if `code' == 56 local reg_`i'_name "Pd"
    if `code' == 80 local reg_`i'_name "DVd"
    
    // Count N for this regimen
    quietly count if TXR_L2 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l2_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD)
    forval b = 1/6 {
        quietly count if TXR_L2 == `code' & BCR_L2 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 8 rows (header + N + 6 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l2 = (8, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l2(1,1) = ("Line 2"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l2(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l2(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l2(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR" "VGPR" "PR" "MR" "SD" "PD" "'
forval b = 1/6 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l2(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l2(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

// Line 3 Regimen × BCR Cross-tabulation
quietly count if TXR_L3 < .
local l3_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL3_TXR)
mata: st_numscalar("n_regimens", cols(oL3_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
    // Assign name based on code
    if `code' == 0  local reg_`i'_name "Other"
    if `code' == 2  local reg_`i'_name "TCd"
    if `code' == 4  local reg_`i'_name "VCd"
    if `code' == 7  local reg_`i'_name "Rd"
    if `code' == 9  local reg_`i'_name "VTd"
    if `code' == 31 local reg_`i'_name "VRd"
    if `code' == 49 local reg_`i'_name "Kd"
    if `code' == 56 local reg_`i'_name "Pd"
    if `code' == 80 local reg_`i'_name "DVd"
    
    // Count N for this regimen
    quietly count if TXR_L3 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l3_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD)
    forval b = 1/6 {
        quietly count if TXR_L3 == `code' & BCR_L3 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 8 rows (header + N + 6 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l3 = (8, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l3(1,1) = ("Line 3"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l3(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l3(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l3(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR" "VGPR" "PR" "MR" "SD" "PD" "'
forval b = 1/6 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l3(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l3(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

// Line 4 Regimen × BCR Cross-tabulation
quietly count if TXR_L4 < .
local l4_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL4_TXR)
mata: st_numscalar("n_regimens", cols(oL4_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
    // Assign name based on code
    if `code' == 0  local reg_`i'_name "Other"
    if `code' == 2  local reg_`i'_name "TCd"
    if `code' == 4  local reg_`i'_name "VCd"
    if `code' == 7  local reg_`i'_name "Rd"
    if `code' == 9  local reg_`i'_name "VTd"
    if `code' == 31 local reg_`i'_name "VRd"
    if `code' == 49 local reg_`i'_name "Kd"
    if `code' == 56 local reg_`i'_name "Pd"
    if `code' == 80 local reg_`i'_name "DVd"
    
    // Count N for this regimen
    quietly count if TXR_L4 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l4_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD)
    forval b = 1/6 {
        quietly count if TXR_L4 == `code' & BCR_L4 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 8 rows (header + N + 6 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l4 = (8, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l4(1,1) = ("Line 4"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l4(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l4(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l4(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR" "VGPR" "PR" "MR" "SD" "PD" "'
forval b = 1/6 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l4(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l4(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

**********
* Overall Survival
**********

putpdf paragraph
putpdf text ("Overall Survival Results"), bold font(,16)

// Summary statistics
quietly summarize OC_TIME, detail
local mean = string(r(mean)/12, "%6.2f")
local sd = string(r(sd)/12, "%5.2f")
local median = string(r(p50)/12, "%6.2f")
local p25 = string(r(p25)/12, "%6.2f")
local p75 = string(r(p75)/12, "%6.2f")

putpdf paragraph
putpdf text ("Summary Statistics"), bold

putpdf table os_sum = (3, 2), border(all)
putpdf table os_sum(1,1) = ("Statistic"), bold
putpdf table os_sum(1,2) = ("Years"), bold
putpdf table os_sum(2,1) = ("Mean (SD)")
putpdf table os_sum(2,2) = ("`mean' (`sd')")
putpdf table os_sum(3,1) = ("Median [IQR]")
putpdf table os_sum(3,2) = ("`median' [`p25'-`p75']")

// Survival at key time points
putpdf paragraph
putpdf text ("Survival at Key Time Points"), bold

putpdf table surv_time = (6, 2), border(all)
putpdf table surv_time(1,1) = ("Time Point"), bold
putpdf table surv_time(1,2) = ("Survival % (95% CI)"), bold

local row = 2
foreach year in 1 2 3 5 10 {
    quietly count if OC_TIME/12 >= `year'
    local pct = (r(N) / _N) * 100    
    local pct_str = string(`pct', "%5.1f")
    
    putpdf table surv_time(`row',1) = ("`year'-year")
    putpdf table surv_time(`row',2) = ("`pct_str'%")
    local row = `row' + 1
}

// Generate and insert KM curves
capture mkdir "$simulated_path/report/figures"

// Overall KM
preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
sts graph, ///
    xtitle("Months") ytitle("Probability") title("") ///
    ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
    xlabel(0(24)240) ci risktable legend(off) ///
    graphregion(color(white)) name(os, replace)
graph export "$simulated_path/report/figures/os.png", replace width(1200)

putpdf paragraph
putpdf text ("Overall Survival"), bold linebreak(2)
putpdf image "$simulated_path/report/figures/os.png", width(7)

// By ASCT
gen asct = SCT_L1
label define asct_lbl 0 "No ASCT" 1 "ASCT"
label values asct asct_lbl
    
sts graph, by(asct) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "No ASCT") label(2 "ASCT")) ///
	name(os_asct, replace)
graph export "$simulated_path/report/figures/os_asct.png", replace width(1200)
    
putpdf paragraph
putpdf text ("Overall Survival by ASCT Status"), bold
putpdf image "$simulated_path/report/figures/os_asct.png", width(6)
restore

// By BCR L1 / ASCT
putpdf pagebreak

preserve
gen bcr_group = .
replace bcr_group = 1 if SCT_L1 == 0 & BCR_L1 == 1 | SCT_L1 == 1 & BCR_SCT == 1
replace bcr_group = 2 if SCT_L1 == 0 & BCR_L1 == 2 | SCT_L1 == 1 & BCR_SCT == 2
replace bcr_group = 3 if SCT_L1 == 0 & BCR_L1 == 3 | SCT_L1 == 1 & BCR_SCT == 3
replace bcr_group = 4 if SCT_L1 == 0 & BCR_L1 == 4 | SCT_L1 == 1 & BCR_SCT == 4
replace bcr_group = 5 if BCR_L1 == 5
replace bcr_group = 6 if BCR_L1 == 6
    
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
sts graph, by(bcr_group) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "CR") label(2 "VGPR") label(3 "PR") label(4 "MR") label(5 "SD") label(6 "PD") rows(2)) ///
	name(os_bcr, replace)
graph export "$simulated_path/report/figures/os_bcr.png", replace width(1200)
   
putpdf paragraph
putpdf text ("Overall Survival by BCR LoT 1 / ASCT"), bold
putpdf image "$simulated_path/report/figures/os_bcr.png", width(6)
restore

// By Age
preserve
gen age_group = .
replace age_group = 1 if Age_DN < 65
replace age_group = 2 if Age_DN >= 65 & Age_DN < 75
replace age_group = 3 if Age_DN >= 75 & Age_DN < .

stset OC_TIME if OC_TIME < 240, failure(OC_MORT)

sts graph, by(age_group) ///
    xtitle("Months") ytitle("Probability") title("") ///
    ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
    graphregion(color(white)) ///
    legend(label(1 "<65") label(2 "65-74") label(3 "≥75") rows(1)) ///
    name(os_age, replace)
graph export "$simulated_path/report/figures/os_age.png", replace width(1200)

putpdf paragraph
putpdf text ("Overall Survival by Age Group"), bold
putpdf image "$simulated_path/report/figures/os_age.png", width(6)
restore

// By R-ISS
putpdf pagebreak

preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
sts graph, by(RISS) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "Stage I") label(2 "Stage II") label(3 "Stage III") rows(1)) ///
	name(os_riss, replace)
graph export "$simulated_path/report/figures/os_riss.png", replace width(1200)

putpdf paragraph
putpdf text ("Overall Survival by R-ISS Stage"), bold
putpdf image "$simulated_path/report/figures/os_riss.png", width(6)
restore

// By ECOG
preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
sts graph, by(ECOGcc) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "ECOG 0") label(2 "ECOG 1") label(3 "ECOG 2+") rows(1)) ///
	name(os_ecog, replace)
graph export "$simulated_path/report/figures/os_ecog.png", replace width(1200)
    
putpdf paragraph
putpdf text ("Overall Survival by ECOG Status"), bold
putpdf image "$simulated_path/report/figures/os_ecog.png", width(6)
restore

**********
* Lines of Therapy
**********

putpdf pagebreak

putpdf paragraph
putpdf text ("Lines of Therapy Distribution"), bold font(,16)

// Count patients receiving each line using TXR_L* variables
// If TXR_L# is not missing, patient received that line

// Calculate N receiving each line
local n_total = _N
forvalues l = 1/9 {
	quietly count if !missing(TXR_L`l')
	local n_l`l' = r(N)
	local pct_l`l' = (`n_l`l'' / `n_total') * 100
}

// Determine max line reached per patient
gen LOT_MAX = 0
forvalues l = 1/9 {
	replace LOT_MAX = `l' if !missing(TXR_L`l')
}

// Count patients by their maximum line reached
forvalues l = 1/9 {
	quietly count if LOT_MAX == `l'
	local n_max_l`l' = r(N)
	local pct_max_l`l' = (`n_max_l`l'' / `n_total') * 100
}

// Display to log
di _col(5) "Line" _col(15) "Received" _col(30) "%" _col(40) "Max Line" _col(55) "%"
di "{hline 60}"
forvalues l = 1/9 {
	di _col(5) "L`l'" _col(15) %9.0fc `n_l`l'' _col(30) %5.1f `pct_l`l'' _col(40) %9.0fc `n_max_l`l'' _col(55) %5.1f `pct_max_l`l''
}

// Create PDF table
putpdf paragraph

putpdf table lot_tbl = (11, 5), border(all)
putpdf table lot_tbl(1,1) = ("Line"), bold
putpdf table lot_tbl(1,2) = ("N Received"), bold
putpdf table lot_tbl(1,3) = ("% Received"), bold
putpdf table lot_tbl(1,4) = ("N Max Line"), bold
putpdf table lot_tbl(1,5) = ("% Max Line"), bold

forvalues l = 1/9 {
	local row = `l' + 1
	putpdf table lot_tbl(`row',1) = ("L`l'")
	putpdf table lot_tbl(`row',2) = ("`=string(`n_l`l'', "%9.0fc")'")
	putpdf table lot_tbl(`row',3) = ("`=string(`pct_l`l'', "%5.1f")'%")
	putpdf table lot_tbl(`row',4) = ("`=string(`n_max_l`l'', "%9.0fc")'")
	putpdf table lot_tbl(`row',5) = ("`=string(`pct_max_l`l'', "%5.1f")'%")
}

putpdf table lot_tbl(11,1) = ("Total"), bold
putpdf table lot_tbl(11,2) = ("`=string(`n_total', "%9.0fc")'"), bold
putpdf table lot_tbl(11,3) = ("—")
putpdf table lot_tbl(11,4) = ("`=string(`n_total', "%9.0fc")'"), bold
putpdf table lot_tbl(11,5) = ("100.0%"), bold

// Summary statistics
quietly summarize LOT_MAX, detail
local mean_lot = string(r(mean), "%4.2f")
local median_lot = string(r(p50), "%4.0f")

putpdf paragraph
putpdf text ("Mean lines received: `mean_lot'; Median: `median_lot'")

// Mortality by maximum line reached

putpdf paragraph
putpdf text ("Mortality by Maximum Line Reached"), bold

putpdf table mort_lot = (11, 4), border(all)
putpdf table mort_lot(1,1) = ("Line"), bold
putpdf table mort_lot(1,2) = ("Deaths"), bold
putpdf table mort_lot(1,3) = ("Deaths during TXD"), bold
putpdf table mort_lot(1,4) = ("Deaths during TFI"), bold

// Row for DN (patients who died before L1)
quietly count if MOR_DN == 1
local n_died_dn = r(N)
putpdf table mort_lot(2,1) = ("DN")
putpdf table mort_lot(2,2) = ("`=string(`n_died_dn', "%9.0fc")'")
putpdf table mort_lot(2,3) = ("—")
putpdf table mort_lot(2,4) = ("`=string(`n_died_dn', "%9.0fc")'")

// L1 to L9
forvalues l = 1/9 {
	local row = `l' + 2
	
	// N with max line = L`l'
	quietly count if LOT_MAX == `l'
	local n_lot = r(N)
	
	// Deaths on treatment (MOR_L#S)
	quietly count if LOT_MAX == `l' & MOR_L`l'S == 1
	local n_died_s = r(N)
	
	// Deaths at exit/TFI (MOR_L#E)
	quietly count if LOT_MAX == `l' & MOR_L`l'E == 1
	local n_died_e = r(N)
	
	putpdf table mort_lot(`row',1) = ("L`l'")
	putpdf table mort_lot(`row',2) = ("`=string(`n_lot', "%9.0fc")'")
	putpdf table mort_lot(`row',3) = ("`=string(`n_died_s', "%9.0fc")'")
	putpdf table mort_lot(`row',4) = ("`=string(`n_died_e', "%9.0fc")'")
}

**********
* Treatment Duration
**********

putpdf pagebreak
putpdf paragraph
putpdf text ("Treatment Duration (TXD)"), bold font(,16)

// Summary statistics for TXD by line
putpdf paragraph
putpdf text ("TXD Summary Statistics (months)"), bold

// Create summary table for TXD
local nrows = 10
putpdf table txd_sum = (`nrows', 6), border(all)
putpdf table txd_sum(1,1) = ("Line"), bold
putpdf table txd_sum(1,2) = ("N"), bold
putpdf table txd_sum(1,3) = ("Mean"), bold
putpdf table txd_sum(1,4) = ("SD"), bold
putpdf table txd_sum(1,5) = ("Median"), bold
putpdf table txd_sum(1,6) = ("IQR"), bold

forvalues l = 1/9 {
	local row = `l' + 1
	quietly summarize TXD_L`l' if TXD_L`l' > 0, detail
	local n = r(N)
	local mean = string(r(mean), "%5.1f")
	local sd = string(r(sd), "%5.1f")
	local median = string(r(p50), "%5.1f")
	local iqr = "[" + string(r(p25), "%4.1f") + "-" + string(r(p75), "%4.1f") + "]"
		
	putpdf table txd_sum(`row',1) = ("L`l'")
	putpdf table txd_sum(`row',2) = ("`=string(`n', "%9.0fc")'")
	putpdf table txd_sum(`row',3) = ("`mean'")
	putpdf table txd_sum(`row',4) = ("`sd'")
	putpdf table txd_sum(`row',5) = ("`median'")
	putpdf table txd_sum(`row',6) = ("`iqr'")
		
	di "TXD L`l': N=" %9.0fc `n' ", Mean=" %5.1f r(mean) ", Median=" %5.1f r(p50)
}

// Generate KM curves for TXD (Lines 1-9)
preserve

gen patient_id = _n
tempfile base_data
save `base_data'

// Stack all lines into long format
clear
local first = 1
forvalues l = 1/9 {
	use `base_data', clear
	keep patient_id TXD_L`l'
	rename TXD_L`l' txd_time
	gen line = `l'
	keep if !missing(txd_time) & txd_time > 0
	gen txd_event = 1
			
	if `first' {
		tempfile txd_stacked
		save `txd_stacked'
		local first = 0
	}
	else {
		append using `txd_stacked'
	save `txd_stacked', replace
	}
}

label define line_lbl 1 "L1" 2 "L2" 3 "L3" 4 "L4" 5 "L5" 6 "L6" 7 "L7" 8 "L8" 9 "L9"
label values line line_lbl

stset txd_time, failure(txd_event)

sts graph, by(line) ///
	ytitle("Proportion on treatment") xtitle("Months") title("") ///
	ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
	legend(order(1 "L1" 2 "L2" 3 "L3" 4 "L4" 5 "L5" ///
           6 "L6" 7 "L7" 8 "L8" 9 "L9") ///
			rows(2) size(small) pos(6)) ///
	scheme(s2color) ///
	graphregion(color(white)) ///
	name(txd, replace)

graph export "$simulated_path/report/figures/txd.png", replace width(1600)

restore

putpdf paragraph
putpdf text ("Treatment Duration by Line of Therapy"), bold
putpdf image "$simulated_path/report/figures/txd.png", width(8)


**********
* Treatment-free Intervals
**********

putpdf pagebreak
putpdf paragraph
putpdf text ("Treatment-Free Interval (TFI)"), bold font(,16)

// Summary statistics for TFI
putpdf paragraph
putpdf text ("TFI Summary Statistics (months)"), bold

// Create summary table for TFI (DN→L1 through L8→L9)
local nrows = 10  // Header + DN + L1-L8
putpdf table tfi_sum = (`nrows', 7), border(all)
putpdf table tfi_sum(1,1) = ("Interval"), bold
putpdf table tfi_sum(1,2) = ("N"), bold
putpdf table tfi_sum(1,3) = ("Mean"), bold
putpdf table tfi_sum(1,4) = ("SD"), bold
putpdf table tfi_sum(1,5) = ("Median"), bold
putpdf table tfi_sum(1,6) = ("IQR"), bold
putpdf table tfi_sum(1,7) = ("Zero TFI %"), bold

// Row 2: DN→L1
capture confirm variable TFI_DN
if !_rc {
	quietly summarize TFI_DN if !missing(TFI_DN), detail
	local n = r(N)
	local mean = string(r(mean), "%5.1f")
	local sd = string(r(sd), "%5.1f")
	local median = string(r(p50), "%5.1f")
	local iqr = "[" + string(r(p25), "%4.1f") + "-" + string(r(p75), "%4.1f") + "]"
	
	quietly count if TFI_DN == 0 & !missing(TFI_DN)
	local n_zero = r(N)
	local pct_zero = string((`n_zero' / `n') * 100, "%4.1f")
	
	putpdf table tfi_sum(2,1) = ("DN→L1")
	putpdf table tfi_sum(2,2) = ("`=string(`n', "%9.0fc")'")
	putpdf table tfi_sum(2,3) = ("`mean'")
	putpdf table tfi_sum(2,4) = ("`sd'")
	putpdf table tfi_sum(2,5) = ("`median'")
	putpdf table tfi_sum(2,6) = ("`iqr'")
	putpdf table tfi_sum(2,7) = ("`pct_zero'%")
	
	di "TFI DN→L1: N=" %9.0fc `n' ", Mean=" `mean' ", Median=" `median' ", Zero TFI=" `pct_zero' "%"
}

// Rows 3-10: L1→L2 through L8→L9
forvalues l = 1/8 {
	local row = `l' + 2
	local next_line = `l' + 1
	capture confirm variable TFI_L`l'
	if !_rc {
		quietly summarize TFI_L`l' if !missing(TFI_L`l'), detail
		local n = r(N)
		local mean = string(r(mean), "%5.1f")
		local sd = string(r(sd), "%5.1f")
		local median = string(r(p50), "%5.1f")
		local iqr = "[" + string(r(p25), "%4.1f") + "-" + string(r(p75), "%4.1f") + "]"
		
		quietly count if TFI_L`l' == 0 & !missing(TFI_L`l')
		local n_zero = r(N)
		local pct_zero = string((`n_zero' / `n') * 100, "%4.1f")
		
		putpdf table tfi_sum(`row',1) = ("L`l'→L`next_line'")
		putpdf table tfi_sum(`row',2) = ("`=string(`n', "%9.0fc")'")
		putpdf table tfi_sum(`row',3) = ("`mean'")
		putpdf table tfi_sum(`row',4) = ("`sd'")
		putpdf table tfi_sum(`row',5) = ("`median'")
		putpdf table tfi_sum(`row',6) = ("`iqr'")
		putpdf table tfi_sum(`row',7) = ("`pct_zero'%")
		
		di "TFI L`l'→L`next_line': N=" %9.0fc `n' ", Mean=" `mean' ", Median=" `median' ", Zero TFI=" `pct_zero' "%"
	}
}
	
// Generate KM curves for TFI (DN→L1 through L8→L9)
preserve

gen patient_id = _n
tempfile base_data
save `base_data'

// Stack all TFI intervals into long format
clear
local first = 1

// TFI DN
use `base_data', clear
keep patient_id TFI_DN
rename TFI_DN tfi_time
gen interval = 0
keep if !missing(tfi_time)
gen tfi_event = 1
			
tempfile tfi_stacked
save `tfi_stacked'
local first = 0

// TFI L1 to L8
forvalues l = 1/8 {
use `base_data', clear
	keep patient_id TFI_L`l'
	rename TFI_L`l' tfi_time
	gen interval = `l'
	keep if !missing(tfi_time)
	gen tfi_event = 1
				
	if `first' {
		tempfile tfi_stacked
		save `tfi_stacked'
		local first = 0
	}
	else {
	append using `tfi_stacked'
	save `tfi_stacked', replace
	}
}

label define int_lbl 0 "DN" 1 "L1" 2 "L2" 3 "L3" 4 "L4" ///
						 5 "L5" 6 "L6" 7 "L7" 8 "L8"
label values interval int_lbl

stset tfi_time, failure(tfi_event)

sts graph, by(interval) ///
	ytitle("Proportion not yet starting next line") xtitle("Months") title("") ///
	ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
	legend(order(1 "DN" 2 "L1" 3 "L2" 4 "L3" 5 "L4" ///
           6 "L5" 7 "L6" 8 "L7" 9 "L8") ///
			rows(2) size(small) pos(6)) ///
	scheme(s2color) ///
	graphregion(color(white)) ///
	name(tfi, replace)

graph export "$simulated_path/report/figures/tfi.png", replace width(1600)

restore
	
putpdf paragraph
putpdf text ("Treatment-Free Interval by LoT"), bold
putpdf image "$simulated_path/report/figures/tfi.png", width(8)

**********
* Economic Outcomes
**********

putpdf pagebreak
putpdf paragraph
putpdf text ("Economic Outcomes"), bold font(,16)

// Get discount rate for display
local drate_pct = string($drate * 100, "%3.1f")

**********
* Costs
**********

putpdf paragraph
putpdf text ("Costs (Discounted at `drate_pct'%)"), bold font(,14)

// Total Costs Summary
quietly summarize cost_total_d, detail
local n_cost = string(r(N), "%9.0fc")
local mean_cost = string(r(mean), "%12.0fc")
local sd_cost = string(r(sd), "%12.0fc")
local median_cost = string(r(p50), "%12.0fc")
local p25_cost = string(r(p25), "%12.0fc")
local p75_cost = string(r(p75), "%12.0fc")

putpdf paragraph
putpdf text ("Total Costs Summary"), bold

putpdf table cost_sum = (3, 2), border(all)
putpdf table cost_sum(1,1) = ("Statistic"), bold
putpdf table cost_sum(1,2) = ("AUD"), bold
putpdf table cost_sum(2,1) = ("Mean (SD)")
putpdf table cost_sum(2,2) = ("$`mean_cost' ($`sd_cost')")
putpdf table cost_sum(3,1) = ("Median [IQR]")
putpdf table cost_sum(3,2) = ("$`median_cost' [$`p25_cost' - $`p75_cost']")

// Cost Components (all means over all patients, so the parts are additive)
putpdf paragraph
putpdf text ("Cost Components (Mean over all patients)"), bold

// ASCT is missing for non-recipients; set to 0 so the population mean is well-defined and additive.
capture drop _asct_d0
gen double _asct_d0 = cond(missing(cost_tx_asct_d), 0, cost_tx_asct_d)
quietly summarize cost_tx_d
scalar m_tx = r(mean)
quietly summarize cost_nt_d
scalar m_nt = r(mean)
quietly summarize _asct_d0
scalar m_asct = r(mean)
scalar m_pbs = m_tx - m_asct           // treatment less ASCT = PBS drugs (regimens + maintenance)
scalar m_total = m_tx + m_nt
drop _asct_d0

putpdf table cost_comp = (6, 2), border(all)
putpdf table cost_comp(1,1) = ("Component"), bold
putpdf table cost_comp(1,2) = ("Mean (AUD)"), bold
putpdf table cost_comp(2,1) = ("Treatment costs (total)")
putpdf table cost_comp(2,2) = ("$" + string(m_tx, "%12.0fc"))
putpdf table cost_comp(3,1) = ("  PBS drugs (regimens + maintenance)")
putpdf table cost_comp(3,2) = ("$" + string(m_pbs, "%12.0fc"))
putpdf table cost_comp(4,1) = ("  ASCT (AR-DRG)")
putpdf table cost_comp(4,2) = ("$" + string(m_asct, "%12.0fc"))
putpdf table cost_comp(5,1) = ("Non-treatment costs (total)")
putpdf table cost_comp(5,2) = ("$" + string(m_nt, "%12.0fc"))
putpdf table cost_comp(6,1) = ("Total"), bold
putpdf table cost_comp(6,2) = ("$" + string(m_total, "%12.0fc")), bold

// Per-recipient context: ASCT and maintenance are received by only a subset of patients.
quietly summarize cost_tx_asct_d if SCT_L1 == 1
local mean_asct_r = cond(r(N) > 0, string(r(mean), "%12.0fc"), "N/A")
local n_asct = string(r(N), "%9.0fc")
quietly summarize cost_tx_mnt_d if MNT == 1
local mean_mnt_r = cond(r(N) > 0, string(r(mean), "%12.0fc"), "N/A")
local n_mnt = string(r(N), "%9.0fc")
putpdf paragraph
putpdf text ("Mean among recipients (discounted): ASCT $`mean_asct_r' (n=`n_asct'); maintenance $`mean_mnt_r' (n=`n_mnt')."), font(,9) italic

// Diagnosis-to-5-year Cost (undiscounted) - table identical to the Total Costs Summary
quietly count if cost_5yr < .
if r(N) > 0 {
	quietly summarize cost_5yr, detail
	local mean_5yr = string(r(mean), "%12.0fc")
	local sd_5yr = string(r(sd), "%12.0fc")
	local median_5yr = string(r(p50), "%12.0fc")
	local p25_5yr = string(r(p25), "%12.0fc")
	local p75_5yr = string(r(p75), "%12.0fc")

	putpdf paragraph
	putpdf text ("Diagnosis-to-5-year Cost (undiscounted)"), bold

	putpdf table cost_5yr_sum = (3, 2), border(all)
	putpdf table cost_5yr_sum(1,1) = ("Statistic"), bold
	putpdf table cost_5yr_sum(1,2) = ("AUD"), bold
	putpdf table cost_5yr_sum(2,1) = ("Mean (SD)")
	putpdf table cost_5yr_sum(2,2) = ("$`mean_5yr' ($`sd_5yr')")
	putpdf table cost_5yr_sum(3,1) = ("Median [IQR]")
	putpdf table cost_5yr_sum(3,2) = ("$`median_5yr' [$`p25_5yr' - $`p75_5yr']")
}

// Treatment Costs by Line of Therapy
putpdf paragraph
putpdf text ("Treatment Costs by Line of Therapy"), bold

// Count rows needed (only include lines with patients)
local n_lines = 0
forval l = 1/9 {
	quietly count if cost_tx_L`l'_d != . & cost_tx_L`l'_d > 0
	if r(N) > 0 local n_lines = `l'
}

// Create table with header + lines
local n_rows = `n_lines' + 1
putpdf table cost_line = (`n_rows', 3), border(all)
putpdf table cost_line(1,1) = ("Line"), bold
putpdf table cost_line(1,2) = ("N Treated"), bold
putpdf table cost_line(1,3) = ("Mean Cost (AUD)"), bold

local row = 2
forval l = 1/`n_lines' {
	quietly count if cost_tx_L`l'_d != . & cost_tx_L`l'_d > 0
	local n_l = string(r(N), "%9.0fc")
	quietly summarize cost_tx_L`l'_d if cost_tx_L`l'_d > 0
	if r(N) > 0 {
		local mean_l = string(r(mean), "%12.0fc")
	}
	else {
		local mean_l = "0"
	}
	
	putpdf table cost_line(`row',1) = ("Line `l'")
	putpdf table cost_line(`row',2) = ("`n_l'")
	putpdf table cost_line(`row',3) = ("$`mean_l'")
	local row = `row' + 1
}

**********
* QALYs
**********

putpdf pagebreak
putpdf paragraph
putpdf text ("Quality-Adjusted Life Years (Discounted at `drate_pct'%)"), bold font(,14)

// Total QALYs Summary
quietly summarize qaly_total_d, detail
local n_qaly = string(r(N), "%9.0fc")
local mean_qaly = string(r(mean), "%5.2f")
local sd_qaly = string(r(sd), "%5.2f")
local median_qaly = string(r(p50), "%5.2f")
local p25_qaly = string(r(p25), "%5.2f")
local p75_qaly = string(r(p75), "%5.2f")

putpdf paragraph
putpdf text ("Total QALYs Summary"), bold

putpdf table qaly_sum = (3, 2), border(all)
putpdf table qaly_sum(1,1) = ("Statistic"), bold
putpdf table qaly_sum(1,2) = ("QALYs"), bold
putpdf table qaly_sum(2,1) = ("Mean (SD)")
putpdf table qaly_sum(2,2) = ("`mean_qaly' (`sd_qaly')")
putpdf table qaly_sum(3,1) = ("Median [IQR]")
putpdf table qaly_sum(3,2) = ("`median_qaly' [`p25_qaly' - `p75_qaly']")

// QALY Components by Health State
putpdf paragraph
putpdf text ("QALYs by Health State (Mean)"), bold

quietly summarize qaly_tfi_DN_d
local q_tfi_dn = string(r(mean), "%5.3f")
quietly summarize qaly_txd_L1_d
local q_txd_l1 = string(r(mean), "%5.3f")
quietly summarize qaly_tfi_L1_d
local q_tfi_l1 = string(r(mean), "%5.3f")
quietly summarize qaly_txd_L2_d
local q_txd_l2 = string(r(mean), "%5.3f")
quietly summarize qaly_post_L2_d
local q_post = string(r(mean), "%5.3f")

putpdf table qaly_comp = (6, 3), border(all)
putpdf table qaly_comp(1,1) = ("Health State"), bold
putpdf table qaly_comp(1,2) = ("Utility Weight"), bold
putpdf table qaly_comp(1,3) = ("Mean QALYs"), bold
putpdf table qaly_comp(2,1) = ("TFI Pre-L1")
putpdf table qaly_comp(2,2) = ("0.72")
putpdf table qaly_comp(2,3) = ("`q_tfi_dn'")
putpdf table qaly_comp(3,1) = ("L1 Treatment")
putpdf table qaly_comp(3,2) = ("0.63")
putpdf table qaly_comp(3,3) = ("`q_txd_l1'")
putpdf table qaly_comp(4,1) = ("TFI Post-L1")
putpdf table qaly_comp(4,2) = ("0.72")
putpdf table qaly_comp(4,3) = ("`q_tfi_l1'")
putpdf table qaly_comp(5,1) = ("L2 Treatment")
putpdf table qaly_comp(5,2) = ("0.67")
putpdf table qaly_comp(5,3) = ("`q_txd_l2'")
putpdf table qaly_comp(6,1) = ("Post-L2")
putpdf table qaly_comp(6,2) = ("0.63")
putpdf table qaly_comp(6,3) = ("`q_post'")

**********
* Undiscounted vs Discounted Comparison
**********

putpdf paragraph
putpdf text ("Discounted vs Undiscounted Comparison"), bold font(,14)

// Get undiscounted values
quietly summarize cost_total
local mean_cost_undisc = string(r(mean), "%12.0fc")
quietly summarize cost_total_d
local mean_cost_disc = string(r(mean), "%12.0fc")

quietly summarize qaly_total
local mean_qaly_undisc = string(r(mean), "%5.2f")
quietly summarize qaly_total_d
local mean_qaly_disc = string(r(mean), "%5.2f")

putpdf table disc_comp = (3, 3), border(all)
putpdf table disc_comp(1,1) = ("Outcome"), bold
putpdf table disc_comp(1,2) = ("Undiscounted"), bold
putpdf table disc_comp(1,3) = ("Discounted (`drate_pct'%)"), bold
putpdf table disc_comp(2,1) = ("Mean Total Cost (AUD)")
putpdf table disc_comp(2,2) = ("$`mean_cost_undisc'")
putpdf table disc_comp(2,3) = ("$`mean_cost_disc'")
putpdf table disc_comp(3,1) = ("Mean QALYs")
putpdf table disc_comp(3,2) = ("`mean_qaly_undisc'")
putpdf table disc_comp(3,3) = ("`mean_qaly_disc'")

**********
* Costs and QALYs by Subgroup
**********

putpdf paragraph
putpdf text ("Economic Outcomes by Subgroup"), bold font(,14)

// By ASCT Status
putpdf paragraph
putpdf text ("By ASCT Status"), bold

quietly summarize cost_total_d if SCT_L1 == 0
local cost_noasct = string(r(mean), "%12.0fc")
quietly summarize qaly_total_d if SCT_L1 == 0
local qaly_noasct = string(r(mean), "%5.2f")
quietly count if SCT_L1 == 0
local n_noasct = string(r(N), "%9.0fc")

quietly summarize cost_total_d if SCT_L1 == 1
local cost_asct = string(r(mean), "%12.0fc")
quietly summarize qaly_total_d if SCT_L1 == 1
local qaly_asct = string(r(mean), "%5.2f")
quietly count if SCT_L1 == 1
local n_asct = string(r(N), "%9.0fc")

putpdf table asct_econ = (3, 4), border(all)
putpdf table asct_econ(1,1) = ("ASCT Status"), bold
putpdf table asct_econ(1,2) = ("N"), bold
putpdf table asct_econ(1,3) = ("Mean Cost (AUD)"), bold
putpdf table asct_econ(1,4) = ("Mean QALYs"), bold
putpdf table asct_econ(2,1) = ("No ASCT")
putpdf table asct_econ(2,2) = ("`n_noasct'")
putpdf table asct_econ(2,3) = ("$`cost_noasct'")
putpdf table asct_econ(2,4) = ("`qaly_noasct'")
putpdf table asct_econ(3,1) = ("ASCT")
putpdf table asct_econ(3,2) = ("`n_asct'")
putpdf table asct_econ(3,3) = ("$`cost_asct'")
putpdf table asct_econ(3,4) = ("`qaly_asct'")

// By Age Group
putpdf paragraph
putpdf text ("By Age Group"), bold

forval a = 1/3 {
	if `a' == 1 {
		local age_cond "Age_DN < 65"
		local age_lab "<65"
	}
	if `a' == 2 {
		local age_cond "Age_DN >= 65 & Age_DN < 75"
		local age_lab "65-74"
	}
	if `a' == 3 {
		local age_cond "Age_DN >= 75 & Age_DN < ."
		local age_lab "≥75"
	}
	
	quietly count if `age_cond'
	local n_age`a' = string(r(N), "%9.0fc")
	quietly summarize cost_total_d if `age_cond'
	local cost_age`a' = string(r(mean), "%12.0fc")
	quietly summarize qaly_total_d if `age_cond'
	local qaly_age`a' = string(r(mean), "%5.2f")
}

putpdf table age_econ = (4, 4), border(all)
putpdf table age_econ(1,1) = ("Age Group"), bold
putpdf table age_econ(1,2) = ("N"), bold
putpdf table age_econ(1,3) = ("Mean Cost (AUD)"), bold
putpdf table age_econ(1,4) = ("Mean QALYs"), bold
putpdf table age_econ(2,1) = ("<65")
putpdf table age_econ(2,2) = ("`n_age1'")
putpdf table age_econ(2,3) = ("$`cost_age1'")
putpdf table age_econ(2,4) = ("`qaly_age1'")
putpdf table age_econ(3,1) = ("65-74")
putpdf table age_econ(3,2) = ("`n_age2'")
putpdf table age_econ(3,3) = ("$`cost_age2'")
putpdf table age_econ(3,4) = ("`qaly_age2'")
putpdf table age_econ(4,1) = ("≥75")
putpdf table age_econ(4,2) = ("`n_age3'")
putpdf table age_econ(4,3) = ("$`cost_age3'")
putpdf table age_econ(4,4) = ("`qaly_age3'")

// By R-ISS Stage
putpdf paragraph
putpdf text ("By R-ISS Stage"), bold

forval r = 1/3 {
	quietly count if RISS == `r'
	local n_riss`r' = string(r(N), "%9.0fc")
	quietly summarize cost_total_d if RISS == `r'
	local cost_riss`r' = string(r(mean), "%12.0fc")
	quietly summarize qaly_total_d if RISS == `r'
	local qaly_riss`r' = string(r(mean), "%5.2f")
}

putpdf table riss_econ = (4, 4), border(all)
putpdf table riss_econ(1,1) = ("R-ISS Stage"), bold
putpdf table riss_econ(1,2) = ("N"), bold
putpdf table riss_econ(1,3) = ("Mean Cost (AUD)"), bold
putpdf table riss_econ(1,4) = ("Mean QALYs"), bold
putpdf table riss_econ(2,1) = ("Stage I")
putpdf table riss_econ(2,2) = ("`n_riss1'")
putpdf table riss_econ(2,3) = ("$`cost_riss1'")
putpdf table riss_econ(2,4) = ("`qaly_riss1'")
putpdf table riss_econ(3,1) = ("Stage II")
putpdf table riss_econ(3,2) = ("`n_riss2'")
putpdf table riss_econ(3,3) = ("$`cost_riss2'")
putpdf table riss_econ(3,4) = ("`qaly_riss2'")
putpdf table riss_econ(4,1) = ("Stage III")
putpdf table riss_econ(4,2) = ("`n_riss3'")
putpdf table riss_econ(4,3) = ("$`cost_riss3'")
putpdf table riss_econ(4,4) = ("`qaly_riss3'")

**********
* Save PDF
**********

set graphics on
local output_file "`report_dir'/${int}_${line}_${data}.pdf"
putpdf save "`output_file'", replace
n di as result _n "Report saved at `output_file'."

}
* end of single-arm ($report_twoarm != "1") report branch
