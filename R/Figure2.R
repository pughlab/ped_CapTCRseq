#!/usr/bin/env Rscript

# Title: Figure 2

# Inputs:
#   - meta_div_goodsamples.rds (from OneDrive datapath as per notebooks)
#   - Color and theme helpers from `R/functions/color_schemes.R` and `R/functions/ggplot2_theme.R`
#   - Misc helpers for p-value formatting from `R/functions/Misc_functions.R`
# Outputs:
#   - Plots saved to `plots/` (e.g., Fig2_Naive_X01.pdf)
# Usage:
#   Rscript R/Figure2.R

# --- Setup --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggdist)
  library(ggsignif)
})
temp <- "~/git/ped_CapTCRseq"
# Load helper palettes and themes
source(file.path(temp, "R/functions", "color_schemes.R"))
source(file.path(temp, "R/functions", "ggplot2_theme.R"))
source(file.path(temp, "R/functions", "Misc_functions.R"))

# Project-rooted paths
project_root <- temp
plots_dir <- file.path(project_root, "plots")
input_rds <- file.path(project_root, "data", "metadata.rds")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)


# --- Reusable plotting helpers -------------------------------------------------

baseraincloud_plot.fx <- function(mydf, yvar, xvar, fillvar, colpal,
                                 scaledots = 0.5, justdots = 1.2, scaleslab = 0.5,
                                 justslab = -0.2, adjustslab = 0.5, binwidth. = ggplot2::unit(0.01, "npc")) {
  p0 <- ggplot(data = mydf, aes(x = eval(parse(text = xvar)), y = eval(parse(text = yvar)), fill = eval(parse(text = fillvar)))) +
    ggdist::geom_dots(side = "left", scale = scaledots, justification = justdots, color = "transparent", overlaps = "nudge", binwidth = binwidth.) +
    ggdist::stat_slab(scale = scaleslab, adjust = adjustslab, justification = justslab) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
    scale_fill_manual(values = colpal) +
    myaxis +
    myplot +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, 110), breaks = c(0, 25, 50, 75, 100)) +
    labs(x = "", y = gsub("`", "", yvar))
  return(p0)
}


generate_raincloud_with_ks <- function(
  df,
  y_col,
  x_col = "cancergroup",
  fill_col = x_col,
  colpal = group_col,
  scaledots = 0.5,
  justdots = 1.2,
  scaleslab = 0.5,
  justslab = -0.2,
  adjustslab = 0.5,
  binwidth. = 2,
  pairs = NULL,
  y_positions = NULL,
  title_label = NULL,
  textsize = 5,
  tip_length = 0.01
) {
  # Build strings for baseraincloud_plot.fx (handles special characters via backticks)
  needs_backticks <- grepl("[^A-Za-z0-9_]", y_col)
  yvar_str <- if (needs_backticks) paste0("`", y_col, "`") else y_col

  # Base plot using the notebook-style helper
  p_base <- baseraincloud_plot.fx(
    mydf = df, yvar = yvar_str, xvar = x_col, fillvar = fill_col, colpal = colpal,
    scaledots = scaledots, justdots = justdots, scaleslab = scaleslab,
    justslab = justslab, adjustslab = adjustslab, binwidth. = binwidth.
  )

  # Summary by group
  message(paste0("Summary for ", y_col, " by ", x_col, ":"))
  print(tapply(df[[y_col]], df[[x_col]], summary, na.rm = TRUE))

  # Determine group level order and default pairs if not provided
  x_levels <- levels(factor(df[[x_col]]))
  if (is.null(pairs)) {
    pairs <- combn(x_levels, 2, simplify = FALSE)
  }

  # KS tests
  message(paste0("KS test results for ", y_col, ":"))
  pvals <- lapply(pairs, function(pr) {
    g1 <- pr[1]; g2 <- pr[2]
    x1 <- df[[y_col]][df[[x_col]] == g1]
    x2 <- df[[y_col]][df[[x_col]] == g2]
    suppressWarnings(stats::ks.test(x1, x2)$p.value)
  })

  # Print results with headers; format only if p < 0.05
  for (i in seq_along(pairs)) {
    pr <- pairs[[i]]
    pv <- pvals[[i]]
    pv_disp <- if (is.finite(pv) && pv < 0.05) round_and_format(pv) else sprintf("%.3f", pv)
    message(sprintf("  %s vs %s: p = %s", pr[1], pr[2], pv_disp))
  }

  # Add annotations for significant comparisons only (p < 0.05)
  sig_idx <- which(vapply(pvals, function(pv) is.finite(pv) && pv < 0.05, logical(1)))
  p_final <- p_base
  if (length(sig_idx) > 0) {
    # Compute x positions
    xmin <- sapply(sig_idx, function(i) match(pairs[[i]][1], x_levels))
    xmax <- sapply(sig_idx, function(i) match(pairs[[i]][2], x_levels))
    annotations <- paste0("p = ", vapply(pvals[sig_idx], round_and_format, character(1)))

    # Y positions: stack slightly above current data max
    y_max <- suppressWarnings(max(df[[y_col]], na.rm = TRUE))
    if (!is.finite(y_max)) y_max <- 1
    if (is.null(y_positions)) {
      # Stagger if multiple sigs
      base <- y_max * 0.92
      step <- (y_max * 0.06) / max(1, (length(sig_idx) - 1))
      y_positions <- base + step * seq_along(sig_idx)
    } else {
      y_positions <- y_positions[seq_along(sig_idx)]
    }

    p_final <- p_final +
      ggsignif::geom_signif(
        y_position = y_positions,
        xmin = xmin,
        xmax = xmax,
        annotation = annotations,
        tip_length = tip_length,
        textsize = textsize
      )
  }

  # Title
  if (is.null(title_label)) title_label <- y_col
  p_final <- p_final + ggplot2::ggtitle(as.expression(bquote(underline(.(title_label)))))

  return(p_final)
}


# --- Load inputs ---------------------------------------------------------------

metadata <- readr::read_rds(input_rds)

# --- Filters --------------------------------------------------

# Exclude missing CD3
metadata <- metadata[!is.na(metadata$CD3), ]
# Exclude T-cell malignancies
metadata <- metadata[metadata$cancergroup != "T-cell malignancies", ]
# PBMC only
pbmc <- metadata[metadata$sampletype == "PBMC", ]
# Remove zero CD3
pbmc <- pbmc[pbmc$CD3 > 0, ]
# Baseline X01
pbmc_01 <- pbmc[pbmc$cycle == "X01", ]


# --- Figure 2A ---------------------------------

p0 <- generate_raincloud_with_ks(
  df = pbmc_01,
  y_col = "PD1%",
  x_col = "cancergroup",
  fill_col = "cancergroup",
  colpal = group_col,
  scaledots = 0.5,
  justdots = 1.2,
  scaleslab = 0.5,
  justslab = -0.2,
  adjustslab = 1,
  binwidth. = NA,
  y_positions = c(80, 90, 100),
  pairs = list(c("Leukemia", "Lymphoma"),
  c("Solid tumors", "Lymphoma"), 
  c("Solid tumors", "Leukemia"))
)

# Save
ggsave(file.path(plots_dir, "Fig2A.pdf"), p0, width = 6, height = 6)
message("Figure 2A saved to ", file.path(plots_dir, "Fig2A.pdf"))

# --- Figure 2B ---------------------------------
p1 <- generate_raincloud_with_ks(
  df = pbmc_01,
  y_col = "LAG3%",
  x_col = "cancergroup",
  fill_col = "cancergroup",
  colpal = group_col,
  scaledots = 0.5,
  justdots = 1.2,
  scaleslab = 0.5,
  justslab = -0.2,
  adjustslab = 0.5,
  binwidth. = NA,
  y_positions = c(85, 100),
  pairs = list(c("Leukemia", "Lymphoma"),
  c("Solid tumors", "Lymphoma"))
)

# Save
ggsave(file.path(plots_dir, "Fig2B.pdf"), p1, width = 6, height = 6)
message("Figure 2B saved to ", file.path(plots_dir, "Fig2B.pdf"))

# --- Figure 2C ---------------------------------
p2 <- generate_raincloud_with_ks(
  df = pbmc_01,
  y_col = "TIM3%",
  x_col = "cancergroup",
  fill_col = "cancergroup",
  colpal = group_col,
  scaledots = 0.5,
  justdots = 1.2,
  scaleslab = 0.5,
  justslab = -0.2,
  adjustslab = 0.5,
  binwidth. = NA,
  y_positions = c(78, 88,105),
  pairs = list(c("Leukemia", "Lymphoma"),
  c("Solid tumors", "Lymphoma"), 
  c("Solid tumors", "Leukemia"))
)

# Save
ggsave(file.path(plots_dir, "Fig2C.pdf"), p2, width = 6, height = 6)
message("Figure 2C saved to ", file.path(plots_dir, "Fig2C.pdf"))

# --- Figure 2D ---------------------------------
# --- Figure 2E ---------------------------------

# --- Figure 2F ---------------------------------
p3 <- generate_raincloud_with_ks(
  df = pbmc_01,
  y_col = "CM%",
  x_col = "cancergroup",
  fill_col = "cancergroup",
  colpal = group_col,
  scaledots = 0.5,
  justdots = 1.2,
  scaleslab = 0.5,
  justslab = -0.2,
  adjustslab = 0.5,
  binwidth. = NA,
  y_positions = c(75, 85, 100),
  pairs = list(c("Leukemia", "Lymphoma"),
  c("Solid tumors", "Lymphoma"), 
  c("Solid tumors", "Leukemia"))
)

# Save
ggsave(file.path(plots_dir, "Fig2F.pdf"), p3, width = 6, height = 6)
message("Figure 2F saved to ", file.path(plots_dir, "Fig2F.pdf"))

# --- Figure 2G ---------------------------------
p4 <- generate_raincloud_with_ks(
  df = pbmc_01,
  y_col = "Naïve%",
  x_col = "cancergroup",
  fill_col = "cancergroup",
  colpal = group_col,
  scaledots = 0.5,
  justdots = 1.2,
  scaleslab = 0.5,
  justslab = -0.2,
  adjustslab = 0.5,
  binwidth. = NA,
  y_positions = 100,
  pairs = list(c("Solid tumors", "Leukemia"))
)

# Save
ggsave(file.path(plots_dir, "Fig2G.pdf"), p4, width = 6, height = 6)
message("Figure 2G saved to ", file.path(plots_dir, "Fig2G.pdf"))

# --- Figure 2H ---------------------------------
p5 <- generate_raincloud_with_ks(
  df = pbmc_01,
  y_col = "TE%",
  x_col = "cancergroup",
  fill_col = "cancergroup",
  colpal = group_col,
  scaledots = 0.5,
  justdots = 1.2,
  scaleslab = 0.5,
  justslab = -0.2,
  adjustslab = 0.5,
  binwidth. = NA,
  y_positions = 105,
  pairs = list(c("Solid tumors", "Leukemia"))
)

# Save
ggsave(file.path(plots_dir, "Fig2H.pdf"), p5, width = 6, height = 6)
message("Figure 2H saved to ", file.path(plots_dir, "Fig2H.pdf"))

# --- Figure 2I ---------------------------------