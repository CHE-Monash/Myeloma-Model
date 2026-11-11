"""Poster recolour of Figure 1 (predicted vs observed BCR).

Same data, same panels and same whisker rule as the manuscript figure
(analyses/transport_dvd/draw_figure1.R). Only the palette and the typeface change:
a single-hue sequential ramp, light to dark across Trial-based -> Calibrated
Transport -> Observed, anchored on the poster's green so the figure sits with the
section headings. Sequential, not categorical, because the three approaches are
ordered along the response gradient.
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import font_manager

# Resolve relative to this file, so the script runs from anywhere. It previously carried an
# absolute path into a session sandbox (/mnt/user-data/uploads/...), which exists on no machine.
_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(_HERE, 'results')
OUT = os.path.join(RES, 'figure1_bcr_predicted_vs_observed_poster.png')

# Carlito is metric-compatible with Calibri, which the poster body now uses.
FONT = 'Carlito' if any('Carlito' in f.name for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams.update({'font.family': FONT, 'axes.linewidth': 1.0})

LIGHT, MID, DARK = '#C7E08F', '#7AB800', '#3A5A1B'   # validated: monotone L, CVD dE 15.2
EDGE = '#2A2A2A'                                      # dark outline: the contrast relief
INK, MUTED = '#1A1A1A', '#4A4A4A'

CATS = ['CR', 'VGPR', 'PR', 'MR', 'SD', 'PD']

# CASTOR Table 2 counts (Palumbo 2016), as hard-coded in the R script and in
# outcomes/calibrated_transport.do. The trial reported these, so bootstrap noise
# has no business moving the bar.
CASTOR = {'vd':  dict(CR=21, VGPR=47, PR=80, MR=20, SD=47, PD=16),
          'dvd': dict(CR=46, VGPR=96, PR=57, MR=10, SD=24, PD=5)}

def reported(arm):
    c = CASTOR[arm]; n = sum(c.values())
    return [100.0 * c[k] / n for k in CATS]

def observed(fname):
    d = {r['bcr']: float(r['pct']) for r in csv.DictReader(open(os.path.join(RES, fname))) if r['bcr'] != 'N'}
    return [d[k] for k in CATS]

def boot(scenario, arm):
    rows = [r for r in csv.DictReader(open(os.path.join(RES, 'bcr_distributions.csv')))
            if r['scenario'] == scenario and r['intervention'] == arm]
    d = {r['bcr']: (float(r['mean']), float(r['lo95']), float(r['hi95'])) for r in rows}
    return ([d[k][0] for k in CATS], [d[k][1] for k in CATS], [d[k][2] for k in CATS])

vd_trial, vd_obs = reported('vd'), observed('observed_vd_bcr.csv')
dvd_trial = reported('dvd')
_, dvd_trial_lo, dvd_trial_hi = boot('A_trial', 'dvd')
dvd_ct, dvd_ct_lo, dvd_ct_hi = boot('B_transport', 'dvd')
dvd_obs = observed('observed_dvd_bcr.csv')

fig, (axL, axR) = plt.subplots(1, 2, figsize=(15.6, 6.0), sharey=True,
                               gridspec_kw=dict(wspace=0.06))
fig.patch.set_facecolor('white')

BARW, STEP = 0.24, 0.255
x = list(range(len(CATS)))

def bars(ax, offset, vals, color, lo=None, hi=None, label=None):
    pos = [i + offset for i in x]
    ax.bar(pos, vals, BARW, color=color, edgecolor=EDGE, linewidth=1.1,
           label=label, zorder=3)
    if lo is not None:
        ax.errorbar(pos, vals, yerr=[[v - l for v, l in zip(vals, lo)],
                                     [h - v for v, h in zip(vals, hi)]],
                    fmt='none', ecolor=INK, elinewidth=1.6, capsize=0, zorder=4)

# Left: the comparator. Two bars, no whiskers - this is the evidence the
# calibration is estimated from, not a prediction being judged.
bars(axL, -STEP / 2, vd_trial, LIGHT)
bars(axL, +STEP / 2, vd_obs, DARK)

# Right: the novel treatment. Whiskers on the two predictions only.
bars(axR, -STEP, dvd_trial, LIGHT, dvd_trial_lo, dvd_trial_hi, 'Trial-based')
bars(axR, 0.0, dvd_ct, MID, dvd_ct_lo, dvd_ct_hi, 'Calibrated Transport')
bars(axR, +STEP, dvd_obs, DARK, label='Observed')

for ax, title in ((axL, 'Vd (common comparator)'), (axR, 'DVd (novel treatment)')):
    ax.set_title(title, fontsize=21, color=INK, pad=14)
    ax.set_xticks(x); ax.set_xticklabels(CATS, fontsize=19, color=INK)
    ax.set_ylim(0, 50); ax.set_xlim(-0.55, len(CATS) - 0.45)
    ax.yaxis.grid(True, color='#DDDDDD', linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.spines['left'].set_color('#9A9A9A'); ax.spines['bottom'].set_color('#9A9A9A')
    ax.tick_params(axis='both', length=0, colors=INK, labelsize=19)

axL.set_ylabel('Patients (%)', fontsize=20, color=INK, labelpad=10)
axR.legend(frameon=False, fontsize=19, loc='upper right', labelcolor=INK,
           handlelength=1.5, handleheight=1.1, borderaxespad=0.4)
fig.text(0.5, 0.015, 'Best clinical response', ha='center', fontsize=20, color=INK)

fig.subplots_adjust(left=0.055, right=0.99, top=0.91, bottom=0.13)
fig.savefig(OUT, dpi=260, facecolor='white')
print('wrote', OUT, os.path.getsize(OUT))
