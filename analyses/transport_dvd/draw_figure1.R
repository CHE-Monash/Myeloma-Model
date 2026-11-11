##########
# Monash Myeloma Model - Figure 1 (predicted vs observed BCR)
#
# Purpose: Two-panel grouped bar chart of second-line best clinical response.
#          LEFT  - Vd, the common comparator: trial (CASTOR) vs observed (registry).
#                  This is the efficacy-effectiveness gap the calibration is estimated FROM.
#          RIGHT - DVd, the novel treatment: Trial-based / Calibrated Transport / Observed.
#                  This is the gap the method has to PREDICT.
#          Colour encodes the approach and is consistent across panels, so the
#          navy-to-green drop means the same thing in both. Bootstrap 95% percentile
#          whiskers on the two DVd predictions only: the observed marginals are fixed
#          benchmarks (as in Table 2 and in the MAE definition) and the Vd bars are the
#          evidence the calibration uses, not predictions being judged for accuracy.
#
# Inputs (all written by the Stata pipeline; bar heights/CIs are the bootstrap mean and
# 95% percentile interval from bootstrap_summary.do):
#   results/bcr_distributions.csv  - scenario x arm x BCR. Supplies the Calibrated
#                                    Transport bar and the DVd trial whisker. The two
#                                    TRIAL bar heights are the published CASTOR counts
#                                    (hard-coded below), not the bootstrap mean.
#   results/observed_dvd_bcr.csv   - registry DVd held-out benchmark, n=399
#   results/observed_vd_bcr.csv    - registry Vd calibration cohort,  n=93
#     (written by outcomes/calibrated_transport.do on the deterministic run)
#
# Outputs: results/figure1_bcr_predicted_vs_observed.{png,pdf}
#
# Usage:   Rscript analyses/transport_dvd/draw_figure1.R      # from the repository root
#
# NOTE: do NOT take the Vd series from the "vd" rows of bcr_distributions.csv for
#       scenarios B/C. Those are the simulated Vd arm of the CE cohort: they imply a
#       mean cut-point shift of about -0.02 against the fitted beta_MRDR of 0.66, and
#       put real-world Vd ABOVE trial Vd, which contradicts the fitted model. The
#       observed Vd marginal is observed_vd_bcr.csv.
#
# Requires: ggplot2 (>= 3.4). No other packages.
##########

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)      # unit(); base-shipped, so no extra install
})

## ---- 0. settings ----------------------------------------------------------
# Times New Roman, to match the manuscript. On a machine without it (most Linux
# boxes) change to "Liberation Serif" or "Nimbus Roman", both Times-metric.
FONT   <- "Times New Roman"

# Greyscale, not colour. Health Economics prefers line figures supplied in black and
# white so they stay legible if a reader prints them, charges for colour in print, and
# otherwise converts colour figures to B&W. The previous navy/maroon/green palette
# converted to greys 55, 67 and 82 of 255, i.e. near-indistinguishable in print.
# Shading runs light -> dark across trial-based -> Calibrated Transport -> observed,
# which also tracks the response gradient. Whiskers fall only on the two LIGHTEST
# bars, so plain black reads clearly against them.
LIGHT  <- "#DEDEDE"   # trial-based
MID    <- "#9A9A9A"   # Calibrated Transport
DARK   <- "#3F3F3F"   # observed
EDGE   <- "#1A1A1A"   # bar outline, so the lightest fill still reads
BARW   <- 0.23        # one bar width for BOTH panels
STEP   <- 0.235       # centre-to-centre spacing (bars just touch)

res <- "analyses/transport_dvd/results"                 # repo root, as run.do assumes
if (!dir.exists(res) && dir.exists("results")) res <- "results"   # or this script's folder
stopifnot(dir.exists(res))

CATS <- c("CR", "VGPR", "PR", "MR", "SD", "PD")
APPR <- c("Trial-based", "Calibrated Transport", "Observed")
PANL <- c("Vd (common comparator)", "DVd (novel treatment)")

## ---- 1. read ---------------------------------------------------------------
bcr <- read.csv(file.path(res, "bcr_distributions.csv"), stringsAsFactors = FALSE)
odv <- read.csv(file.path(res, "observed_dvd_bcr.csv"),  stringsAsFactors = FALSE)
ovd <- read.csv(file.path(res, "observed_vd_bcr.csv"),   stringsAsFactors = FALSE)
odv <- odv[odv$bcr != "N", ]
ovd <- ovd[ovd$bcr != "N", ]

pick <- function(scn, arm) {
  d <- bcr[bcr$scenario == scn & bcr$intervention == arm, ]
  d[match(CATS, d$bcr), c("bcr", "mean", "lo95", "hi95")]
}
obs_row <- function(d) data.frame(bcr = CATS,
                                  mean = d$pct[match(CATS, d$bcr)],
                                  lo95 = NA_real_, hi95 = NA_real_)

# ---- trial arms: plot the PUBLISHED distribution, not the bootstrap mean ----
# CASTOR Table 2 counts (Palumbo et al. 2016), also hard-coded in
# outcomes/calibrated_transport.do (see "CASTOR Vd arm" / "CASTOR DVd arm").
# The trial reported these numbers, so resampling noise has no business moving the
# bar: the bootstrap mean sits up to 0.3 pp away, enough to shift DVd CR from 19% to
# 20% once rounded. The DVd whisker is still the bootstrap interval, because that is
# the uncertainty the analysis actually propagated.
CASTOR <- list(
  vd  = c(CR = 21, VGPR = 47, PR = 80, MR = 20, SD = 47, PD = 16),   # n = 231
  dvd = c(CR = 46, VGPR = 96, PR = 57, MR = 10, SD = 24, PD =  5)    # n = 238
)
rep_pct <- function(arm) {
  cts <- CASTOR[[arm]][CATS]
  data.frame(bcr = CATS, mean = 100 * cts / sum(cts),
             lo95 = NA_real_, hi95 = NA_real_, row.names = NULL)
}

# guard: the bootstrap mean is an unbiased estimate of the reported value, so a gap of
# more than a point means the hard-coded counts have drifted from the Stata side.
.chk <- merge(rep_pct("dvd"), pick("A_trial", "dvd"), by = "bcr", suffixes = c(".rep", ".boot"))
if (any(abs(.chk$mean.rep - .chk$mean.boot) > 1, na.rm = TRUE)) {
  warning("CASTOR counts here differ from the bootstrap mean by >1 pp; ",
          "check them against outcomes/calibrated_transport.do")
}

# The Vd panel carries no intervals at all: it is the evidence the calibration is
# estimated from, not a prediction being judged for accuracy.
vd_trial  <- rep_pct("vd")
dvd_trial <- rep_pct("dvd")
.boot     <- pick("A_trial", "dvd")
dvd_trial$lo95 <- .boot$lo95[match(dvd_trial$bcr, .boot$bcr)]
dvd_trial$hi95 <- .boot$hi95[match(dvd_trial$bcr, .boot$bcr)]

dat <- rbind(
  cbind(vd_trial,                    panel = PANL[1], approach = APPR[1]),
  cbind(obs_row(ovd),                panel = PANL[1], approach = APPR[3]),
  cbind(dvd_trial,                   panel = PANL[2], approach = APPR[1]),
  cbind(pick("B_transport", "dvd"),  panel = PANL[2], approach = APPR[2]),
  cbind(obs_row(odv),                panel = PANL[2], approach = APPR[3])
)
dat$panel    <- factor(dat$panel,    levels = PANL)
dat$approach <- factor(dat$approach, levels = APPR)
dat$cat      <- match(dat$bcr, CATS)

## ---- 2. explicit x positions ----------------------------------------------
# Bar width is identical in both panels, so the Vd panel simply carries more air
# between response categories. Positions are computed rather than dodged, because
# position_dodge(preserve = "single") left-aligns the 2-bar panel instead of centring it.
dat$xpos <- NA_real_
off <- c(-STEP/2, STEP/2)                       # Vd: 2 bars, centred
dat$xpos[dat$panel == PANL[1]] <-
  dat$cat[dat$panel == PANL[1]] +
  off[match(as.character(dat$approach[dat$panel == PANL[1]]), APPR[c(1, 3)])]
off3 <- c(-STEP, 0, STEP)                       # DVd: 3 bars, centred
dat$xpos[dat$panel == PANL[2]] <-
  dat$cat[dat$panel == PANL[2]] +
  off3[as.integer(dat$approach[dat$panel == PANL[2]])]

## ---- 3. draw ---------------------------------------------------------------
p <- ggplot(dat, aes(x = xpos, y = mean, fill = approach)) +
  geom_col(width = BARW, colour = EDGE, linewidth = 0.35) +
  geom_linerange(aes(ymin = lo95, ymax = hi95),
                 na.rm = TRUE, colour = "black", linewidth = 0.55) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = setNames(c(LIGHT, MID, DARK), APPR), drop = FALSE) +
  scale_x_continuous(breaks = 1:6, labels = CATS,
                     limits = c(0.5, 6.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Best clinical response", y = "Patients (%)", fill = NULL) +
  theme_minimal(base_size = 12.5, base_family = FONT) +
  theme(
    # Times is a light face, so everything is bold and near-black: theme_minimal's
    # default grey30 text reads as washed out at figure sizes.
    text = element_text(face = "bold", colour = "#111111"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.45),
    axis.line.x  = element_line(colour = "#333333", linewidth = 0.5),
    axis.line.y  = element_line(colour = "#333333", linewidth = 0.5),
    axis.ticks   = element_blank(),
    axis.text    = element_text(face = "bold", colour = "#111111"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
    strip.text   = element_text(face = "bold", size = 13, margin = margin(b = 8)),
    panel.spacing = unit(1.4, "lines"),
    legend.background = element_blank(),
    legend.key.size   = unit(0.95, "lines"),
    legend.text       = element_text(face = "bold", size = 11),
    plot.margin       = margin(6, 10, 4, 6)
  )

# Legend inside the top-right of the DVd panel. ggplot2 3.5 renamed this argument;
# nudge the two coordinates if it sits awkwardly on your version.
if (utils::packageVersion("ggplot2") >= "3.5.0") {
  p <- p + theme(legend.position = "inside",
                 legend.position.inside = c(0.995, 0.99),
                 legend.justification = c(1, 1))
} else {
  p <- p + theme(legend.position = c(0.995, 0.99),
                 legend.justification = c(1, 1))
}

## ---- 4. save ---------------------------------------------------------------
png_out <- file.path(res, "figure1_bcr_predicted_vs_observed.png")
pdf_out <- file.path(res, "figure1_bcr_predicted_vs_observed.pdf")
ggsave(png_out, p, width = 11.6, height = 4.7, dpi = 300, bg = "white")
# cairo_pdf embeds system fonts properly; fall back to the base pdf device (which uses
# the built-in Times) if this build of R has no cairo support.
tryCatch(
  ggsave(pdf_out, p, width = 11.6, height = 4.7, device = cairo_pdf, bg = "white"),
  error = function(e) {
    message("cairo_pdf unavailable (", conditionMessage(e), "); using base pdf device.")
    ggsave(pdf_out, p, width = 11.6, height = 4.7, device = grDevices::pdf, family = "Times")
  }
)
cat("Figure 1 written to:\n  ", png_out, "\n  ", pdf_out, "\n", sep = "")
