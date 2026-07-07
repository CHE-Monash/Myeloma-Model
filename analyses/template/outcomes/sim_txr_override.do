**********
* Monash Myeloma Model - Sim TXR Override (template)
*
* Purpose: force the treatment regimen (TXR) at line $line -- e.g. assign the intervention drug to the
*          whole cohort at the decision line instead of letting the regimen model draw it. Inert until
*          you fill the $int branches. Delete this file if the analysis uses the standard TXR model.
* Usage:   auto-run by core/outcomes/sim_txr.do when the current Line equals $line (no registration).
*          Set $line in simulate.do; place this file in analyses/<analysis>/outcomes/.
* Notes:   Mata state: mTXR regimen matrix (write mTXR[idx, Line]); mMOR, mState, OMC, Line. Regimen
*          codes are MRDR Regimen values (see outcomes/txr_template.do). Worked example:
*          analyses/transport_dvd/outcomes/sim_txr_override.do.
**********

capture program drop txr_override
program define txr_override
    args tx_code
    mata {
        // Alive, eligible patients at this stage
        idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
        if (rows(idx) > 0) {
            mTXR[idx, Line] = J(rows(idx), 1, strtoreal(st_local("tx_code")))   // <- forced regimen
        }
    }
end

* ---- Dispatch: choose the regimen code per intervention arm ----
* if ("$int" == "<YOUR_INTERVENTION>") txr_override <regimen_code>
* if ("$int" == "<YOUR_COMPARATOR>")   txr_override <regimen_code>
