**********
* Monash Myeloma Model - Sim BCR Override (template)
*
* Purpose: replace the best-clinical-response (BCR) draw at line $line for this analysis -- e.g. inject a
*          trial-calibrated response distribution for a drug the registry model can't fit. Inert until you
*          fill the dispatch branches. Delete this file if the analysis uses the standard BCR model.
* Usage:   auto-run by core/outcomes/sim_bcr.do when the current Line equals $line (no registration).
*          Set $line in simulate.do; place this file in analyses/<analysis>/outcomes/.
* Notes:   Mata state (per current Line/OMC): mBCR response matrix (write mBCR[idx, Line]); mMOR
*          mortality-by-OMC; mState pathway state; vXB default linear predictor; mPat default design;
*          rn_bcr(Line) the BCR CRN slot (reuse so arms stay aligned). Worked example:
*          analyses/transport_dvd/outcomes/sim_bcr_override.do.
**********

* ---- Example A: assign BCR from a fixed 6-category distribution, ranked by prognosis (vXB) ----
*   Pass the name of a 1x6 probability matrix (loaded via mata matuse below).
capture program drop bcr_from_distribution
program define bcr_from_distribution
    args bcr_matrix
    mata {
        // Alive, non-prevalent patients at this stage
        idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
        if (rows(idx) > 0) {
            nPatients = rows(idx)
            mBCR_probs = *findexternal(st_local("bcr_matrix"))     // 1x6 category proportions

            // Rank by linear predictor (1 = worst prognosis) -> percentile -> category
            ranks       = invorder(order(vXB, -1))
            percentiles = ranks :/ nPatients
            cumProbs    = runningsum(mBCR_probs)
            vOC = 1 :* (percentiles :<= cumProbs[1]) +
                  2 :* (percentiles :>  cumProbs[1] :& percentiles :<= cumProbs[2]) +
                  3 :* (percentiles :>  cumProbs[2] :& percentiles :<= cumProbs[3]) +
                  4 :* (percentiles :>  cumProbs[3] :& percentiles :<= cumProbs[4]) +
                  5 :* (percentiles :>  cumProbs[4] :& percentiles :<= cumProbs[5]) +
                  6 :* (percentiles :>  cumProbs[5])
            mBCR[idx, Line] = vOC                                   // <- the override
        }
    }
end

* ---- Example B: assign BCR from a fitted ordered logit (indicator design) ----
*   Pass the name of a coefficient row vector = (betas | cutpoints).
capture program drop bcr_from_ologit
program define bcr_from_ologit
    args coef_matrix
    mata {
        idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
        if (rows(idx) > 0) {
            nPatients = rows(idx)
            // Build your design matrix here (indicator example: intercept + one arm flag)
            armFlag = (st_global("int") == "<YOUR_INT>")           // <- edit
            mPatO   = (J(nPatients,1,1), J(nPatients,1,armFlag))

            vCoef_full  = *findexternal(st_local("coef_matrix"))
            vCoef       = vCoef_full[1, 1..cols(mPatO)]'
            vXBo        = mPatO * vCoef
            nCutPoints  = 5                                          // 6 categories -> 5 cutpoints
            cutPoints   = vCoef_full[1, (cols(vCoef_full)-nCutPoints+1)..cols(vCoef_full)]

            cumProbs = calcOrdLogitProbs(vXBo, cutPoints)
            vRN      = rnDraw(idx, rn_bcr(Line))                     // CRN: aligned across arms
            vOC      = assignOrdOutcome(vRN, cumProbs, (1,2,3,4,5,6))
            mBCR[idx, Line] = vOC
        }
    }
end

* ---- Dispatch: pick a program + coefficient file per $scenario / $int / $boot ----
* Fill these branches in. Deterministic files live under $outcomes_path/<scenario>/; bootstrap files
* under $outcomes_path/<scenario>/bootstrap/ suffixed _B$b (see transport_dvd for the full pattern).
if ($boot == 0) {
    if ("$scenario" == "<YOUR_SCENARIO>") {
        // mata: mata matuse "$outcomes_path/<scenario>/<your_probs>.mmat", replace
        // bcr_from_distribution <your_probs>
    }
}
else if ($boot == 1) {
    if ("$scenario" == "<YOUR_SCENARIO>") {
        // mata: mata matuse "$outcomes_path/<scenario>/bootstrap/<your_probs>_B${b}.mmat", replace
        // bcr_from_distribution <your_probs>
    }
}
