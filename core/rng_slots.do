**********
* EpiMAP Myeloma - RNG Slot Registry (Common Random Numbers)
*
* Purpose: Single source of truth for the column layout of the pre-drawn
*          random-number matrix mRN (Obs x K) used for common random numbers
*          (CRN). Every stochastic event reads its own FIXED column for its own
*          FIXED patient (row), so the same uniform feeds the same event for the
*          same patient in every arm of a comparison. See CRN_implementation_plan.md.
*
* What this file provides:
*   - rn_K()                         total number of columns (K) to allocate
*   - rn_<event>(point) accessors    the column index for each (event, point)
*   - rnDraw(idx, slot)              the shared uniforms mRN[idx, slot] for those
*                                    patients at that event. CRN is unconditional -
*                                    there is NO runtime toggle and NO runiform
*                                    fallback (plan Q5, 15 Jun 2026).
*
* What it deliberately does NOT do (done in mata_setup.do):
*   - build mRN                      set seed (base + b); mRN = runiform(Obs, rn_K())
*
* Loading: run once before the simulation pass, alongside the other Mata
*          definitions, i.e. add  run "core/rng_slots.do"  to run_pipeline (after
*          core/mata_functions.do). Functions compiled in a mata{} block persist
*          for the session, exactly like core/mata_functions.do.
*
* Design rules (from the plan, section 4):
*   - One slot per (event, pathway point). Modules that draw at several lines/OMC
*     (OS, BCR, TXR, TXD, TFI) get a DISTINCT column per line/OMC.
*   - Sequential draws for the SAME patient need distinct slots (e.g. the L1 TXD
*     splines S1->S2->S3). Disjoint patient subsets at one point (ASCT vs no-ASCT)
*     are given SEPARATE slots here for readability, though they could share.
*   - An override that REPLACES a core event must reuse that event's slot
*     (e.g. the L2 BCR override calls rn_bcr($line)); an override that introduces
*     a genuinely NEW stochastic event takes a column from the reserved block
*     (rn_override(i)).
*
* Column map (1-indexed; K = 74). Engine call sites are the live draws in
*   core/outcomes/ inventoried in section 5 of the plan.
*
*   cols  event              keyed by   accessor             core call site(s)
*   ----  -----------------  ---------  -------------------  --------------------------------
*    1-19 OS                 OMC 1..19  rn_os(omc)           sim_os.do:108 (DN), :105 (OMC>=2, * vPR)
*   20-28 BCR                Line 1..9  rn_bcr(line)         sim_bcr.do:96  (+ overrides reuse)
*      29 BCR after ASCT     -          rn_bcr_asct()        sim_bcr_asct.do:47
*   30-38 TXR (regimen)      Line 1..9  rn_txr(line)         sim_txr.do:96
*   39-43 TXD L1 (5 draws)   sub 1..5   rn_txd_l1(sub)       sim_txd_l1.do:52,73,100,147,174
*   44-51 TXD L2+            Line 2..9  rn_txd(line)         sim_txd.do:81
*      52 TFI at DN          -          rn_tfi_dn()          sim_tfi_dn.do:32
*   53-54 TFI at L1 (2)      branch 1,2 rn_tfi_l1(branch)    sim_tfi_l1.do:44 (ASCT), :80 (noASCT)
*   55-62 TFI L2+            Line 2..9  rn_tfi(line)         sim_tfi.do:69
*      63 ASCT at DN         -          rn_asct_dn()         sim_asct_dn.do:36
*      64 ASCT at L1         -          rn_asct_l1()         sim_asct_l1.do:44
*   65-66 Maintenance (2)    branch 1,2 rn_mnt(branch)       sim_mnt.do:69 (ASCT), :119 (noASCT)
*   67-74 Reserved override  i 1..8     rn_override(i)       analysis overrides introducing new draws
*
*   TXD_L1 sub-index:  1 = ASCT spline 1, 2 = ASCT spline 2 (cond.),
*                      3 = ASCT spline 3 (cond.), 4 = no-ASCT, 5 = continuous therapy
*   TFI_L1 branch:     1 = ASCT, 2 = no-ASCT
*   MNT branch:        1 = ASCT, 2 = no-ASCT
*
* Author: Adam Irving
* Date: June 2026
**********

mata:

// ---- Tiny assert helper (defined first; loud failure beats a silent mis-slot) ----
void rn_assert(real scalar cond, string scalar msg) {
	if (cond == 0) _error(msg)
}

// ---- Block bases (0-indexed; a block of width w occupies base+1 .. base+w) ----
real scalar rn_base_os()       return(0)    // 1..19
real scalar rn_base_bcr()      return(19)   // 20..28
real scalar rn_base_bcrasct()  return(28)   // 29
real scalar rn_base_txr()      return(29)   // 30..38
real scalar rn_base_txdl1()    return(38)   // 39..43
real scalar rn_base_txd()      return(43)   // 44..51  (line 2..9)
real scalar rn_base_tfidn()    return(51)   // 52
real scalar rn_base_tfil1()    return(52)   // 53..54
real scalar rn_base_tfi()      return(54)   // 55..62  (line 2..9)
real scalar rn_base_asctdn()   return(62)   // 63
real scalar rn_base_asctl1()   return(63)   // 64
real scalar rn_base_mnt()      return(64)   // 65..66
real scalar rn_base_override() return(66)   // 67..74

// ---- Total columns to allocate ----
real scalar rn_K() return(74)

// ---- Accessors: each returns the absolute column index for (event, point) ----

// Overall survival, one column per pathway point (OMC = 1=DN, 2=L1S, ... 19=L9E)
real scalar rn_os(real scalar omc) {
	rn_assert(omc >= 1 & omc <= 19, "rn_os: OMC out of range 1..19")
	return(rn_base_os() + omc)
}

// Best clinical response, one column per line (1..9). Overrides reuse this slot.
real scalar rn_bcr(real scalar line) {
	rn_assert(line >= 1 & line <= 9, "rn_bcr: line out of range 1..9")
	return(rn_base_bcr() + line)
}

// Best clinical response after ASCT (L1E only)
real scalar rn_bcr_asct() return(rn_base_bcrasct() + 1)

// Treatment regimen, one column per line (1..9)
real scalar rn_txr(real scalar line) {
	rn_assert(line >= 1 & line <= 9, "rn_txr: line out of range 1..9")
	return(rn_base_txr() + line)
}

// Treatment duration at L1 (5 sequential/disjoint sub-draws; see header for codes)
real scalar rn_txd_l1(real scalar sub) {
	rn_assert(sub >= 1 & sub <= 5, "rn_txd_l1: sub-index out of range 1..5")
	return(rn_base_txdl1() + sub)
}

// Treatment duration at L2+ (line 2..9)
real scalar rn_txd(real scalar line) {
	rn_assert(line >= 2 & line <= 9, "rn_txd: line out of range 2..9")
	return(rn_base_txd() + (line - 1))
}

// Treatment-free interval at diagnosis
real scalar rn_tfi_dn() return(rn_base_tfidn() + 1)

// Treatment-free interval at L1E (branch 1=ASCT, 2=no-ASCT)
real scalar rn_tfi_l1(real scalar branch) {
	rn_assert(branch == 1 | branch == 2, "rn_tfi_l1: branch must be 1 or 2")
	return(rn_base_tfil1() + branch)
}

// Treatment-free interval at L2+ (line 2..9)
real scalar rn_tfi(real scalar line) {
	rn_assert(line >= 2 & line <= 9, "rn_tfi: line out of range 2..9")
	return(rn_base_tfi() + (line - 1))
}

// ASCT decision at diagnosis / at L1
real scalar rn_asct_dn() return(rn_base_asctdn() + 1)
real scalar rn_asct_l1() return(rn_base_asctl1() + 1)

// Maintenance therapy at L1E (branch 1=ASCT, 2=no-ASCT)
real scalar rn_mnt(real scalar branch) {
	rn_assert(branch == 1 | branch == 2, "rn_mnt: branch must be 1 or 2")
	return(rn_base_mnt() + branch)
}

// Reserved columns for overrides that introduce a NEW stochastic event (i = 1..8)
real scalar rn_override(real scalar i) {
	rn_assert(i >= 1 & i <= 8, "rn_override: index out of range 1..8")
	return(rn_base_override() + i)
}

// ---- Draw accessor: shared uniforms for patients `idx` at event column `slot` ----
//   CRN is the model's behaviour: no runtime toggle, no runiform fallback. To
//   switch generators in future (e.g. a counter-based hash) change only this.
//   mRN is the external matrix built once per cohort in core/mata_setup.do.
real colvector rnDraw(real colvector idx, real scalar slot) {
	external real matrix mRN
	return(mRN[idx, slot])
}

end
