#!/usr/bin/env Rscript

# Title: Figure 2

# Inputs:
#   - metadata.rds 
#   - Color and theme helpers from `R/functions/color_schemes.R` and `R/functions/ggplot2_theme.R`
#   - Misc helpers for p-value formatting from `R/functions/Misc_functions.R`
#   - Plotting helpers from `R/functions/Plotting_functions.R`
# Outputs:
#   - Plots saved to `plots/` (e.g., Fig2A.pdf)
#   - Tables saved to `tables/` (e.g., PD1percent_cohort_median_summary.xlsx)
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

# Load helper palettes and themes
source(file.path(getwd(), "R/functions", "color_schemes.R"))
source(file.path(getwd(), "R/functions", "ggplot2_theme.R"))
source(file.path(getwd(), "R/functions", "Misc_functions.R"))
source(file.path(getwd(), "R/functions", "Plotting_functions.R"))

# Project-rooted paths
project_root <- getwd()
plots_dir <- file.path(project_root, "plots")
tables_dir <- file.path(project_root, "tables")
input_rds <- file.path(project_root, "data", "metadata.rds")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)


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

build_marker_summary_plot <- function(
  df,
  marker,
  disease_col = "Disease_type",
  cancergroup_col = "cancergroup",
  group_col,
  median_fun = median.cohorts.fx,
  sort_fun = sort.df.fx,
  splot_fun = Splot.fx,
  ylim_range = c(0, 100),
  y_label = NULL,
  title_underlined = TRUE,
  dummy_groups = c("EMPTY1", "EMPTY2", "EMPTY3"),
  dummy_cancergroups = c("Lymphoma", "Solid tumors", "T-cell malignancies"),
  export_path = NULL,
  export_engine = c("auto", "writexl", "openxlsx"),
  export_overwrite = TRUE
) {
  export_engine <- match.arg(export_engine)

  # Defensive: ensure needed columns exist and are not factors (for safe assignment)
  if (!cancergroup_col %in% names(df)) {
    df[[cancergroup_col]] <- NA
  }
  if (is.factor(df[[disease_col]])) df[[disease_col]] <- as.character(df[[disease_col]])
  if (is.factor(df[[cancergroup_col]])) df[[cancergroup_col]] <- as.character(df[[cancergroup_col]])

  # Create 3 NA dummy rows preserving column classes
  if (nrow(df) == 0) stop("Input df has 0 rows; cannot infer column classes.")
  dummy_rows <- df[rep(1, length(dummy_groups)), , drop = FALSE]
  dummy_rows[] <- NA

  # Fill dummy rows
  for (i in seq_along(dummy_groups)) {
    dummy_rows[i, disease_col] <- dummy_groups[i]
    dummy_rows[i, marker] <- -1
    dummy_rows[i, cancergroup_col] <- dummy_cancergroups[i]
  }

  # Prepend dummy rows
  df_augmented <- rbind(dummy_rows, df)

  # Compute medians per cohort
  mymed <- median_fun(df_augmented, marker, disease_col)

  # Attach cancergroup to medians by matching on group vs disease_col
  mymed[[cancergroup_col]] <- df_augmented[[cancergroup_col]][
    match(mymed$group, df_augmented[[disease_col]])
  ]

  # Order by cancergroup then median
  mymed <- mymed[order(mymed[[cancergroup_col]], mymed$median), , drop = FALSE]

  # Exclude dummy rows for export/return
  mymed_no_dummy <- mymed[!(mymed$group %in% dummy_groups), , drop = FALSE]

  # Sort full data using provided helper
  sorted_df_list <- sort_fun(df_augmented, mymed, marker, disease_col)

  # Colors for dummy entries (white for empties)
  dummy_colors <- rep("black", nrow(mymed))
  dummy_colors[mymed$group %in% dummy_groups] <- "white"

  # Quote marker for non-syntactic names when passing as string to plotting helper
  is_standard_name <- grepl("^[A-Za-z\\.][A-Za-z0-9_\\.]*$", marker)
  marker_for_plot <- if (is_standard_name) marker else paste0("`", marker, "`")

  # Build plot
  sp <- splot_fun(sorted_df_list, marker_for_plot, cancergroup_col, group_col, "", rmEMPTY = dummy_colors)

  # Y label and title
  final_ylabel <- if (is.null(y_label)) marker else y_label
  if (isTRUE(title_underlined)) {
    plot_title <- bquote(underline(.(marker)))
  } else {
    plot_title <- marker
  }
  sp_limited <- sp + ylim(ylim_range) + labs(y = final_ylabel) + ggtitle(plot_title)

  # Optional export to Excel
  if (!is.null(export_path)) {
    if (export_engine == "auto") {
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(mymed_no_dummy, path = export_path)
      } else if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(mymed_no_dummy, file = export_path, overwrite = export_overwrite)
      } else {
        stop("No Excel writer available. Install either 'writexl' or 'openxlsx'.")
      }
    } else if (export_engine == "writexl") {
      if (!requireNamespace("writexl", quietly = TRUE)) {
        stop("Package 'writexl' not installed.")
      }
      writexl::write_xlsx(mymed_no_dummy, path = export_path)
    } else if (export_engine == "openxlsx") {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("Package 'openxlsx' not installed.")
      }
      openxlsx::write.xlsx(mymed_no_dummy, file = export_path, overwrite = export_overwrite)
      message("Exported to ", export_path)
    }
  }

  # Return everything useful
  list(
    data_augmented = df_augmented,
    medians = mymed,
    medians_no_dummy = mymed_no_dummy,
    colors_dummy = dummy_colors,
    sorted = sorted_df_list,
    plot = sp,
    plot_final = sp_limited + theme(axis.text = element_text(size = 13), 
                                    axis.title = element_text(size = 13), 
                                    plot.title = element_text(size = 13), 
                                    plot.margin = unit(c(0, 0, 0, 0), "cm")),
    export_path = export_path
  )
}
# --- Load inputs ---------------------------------------------------------------

metadata <- readr::read_rds(input_rds)

# --- Filters --------------------------------------------------
## For raincloud plots
# Exclude missing CD3
metadata <- metadata[!is.na(metadata$CD3), ]
# Exclude T-cell malignancies
metadata_raincloud <- metadata[metadata$cancergroup != "T-cell malignancies", ]
# PBMC only
pbmc <- metadata_raincloud[metadata_raincloud$sampletype == "PBMC", ]
# Remove zero CD3
pbmc <- pbmc[pbmc$CD3 > 0, ]
# Baseline X01
pbmc_01 <- pbmc[pbmc$cycle == "X01", ]

## For S plots - include T-cell malignancies
# PBMC only
pbmc <- metadata[metadata$sampletype == "PBMC", ]
# Remove zero CD3
pbmc <- pbmc[pbmc$CD3 > 0, ]
# Baseline X01
pbmc_01_splot <- pbmc[pbmc$cycle == "X01", ]

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
p3 <- build_marker_summary_plot(
  pbmc_01_splot,
  marker = "PD1%",
  disease_col = "Disease_type",
  cancergroup_col = "cancergroup",
  group_col = group_col,
  median_fun = median.cohorts.fx,
  sort_fun = sort.df.fx,
  splot_fun = Splot.fx,
  ylim_range = c(0, 100),
  y_label = "PD1%",
  title_underlined = TRUE,
  dummy_groups = c("EMPTY1", "EMPTY2", "EMPTY3"),
  dummy_cancergroups = c("Lymphoma", "Solid tumors", "T-cell malignancies"),
  export_path = file.path(tables_dir, "PD1percent_cohort_median_summary.xlsx"),
  export_engine = "auto",
  export_overwrite = TRUE
) 

ggsave(file.path(plots_dir, "Fig2D.pdf"), p3$plot_final, width = 6, height = 4)
message("Figure 2D saved to ", file.path(plots_dir, "Fig2D.pdf"))

# --- Figure 2E ---------------------------------

p4 <- build_marker_summary_plot(
  pbmc_01_splot,
  marker = "LAG3%",
  disease_col = "Disease_type",
  cancergroup_col = "cancergroup",
  group_col = group_col,
  median_fun = median.cohorts.fx,
  sort_fun = sort.df.fx,
  splot_fun = Splot.fx,
  ylim_range = c(0, 100),
  y_label = "LAG3%",
  title_underlined = TRUE,
  dummy_groups = c("EMPTY1", "EMPTY2", "EMPTY3"),
  dummy_cancergroups = c("Lymphoma", "Solid tumors", "T-cell malignancies"),
  export_path = file.path(tables_dir, "LAG3percent_cohort_median_summary.xlsx"),
  export_engine = "auto",
  export_overwrite = TRUE
) 
ggsave(file.path(plots_dir, "Fig2E.pdf"), p4$plot_final, width = 6, height = 4)
message("Figure 2E saved to ", file.path(plots_dir, "Fig2E.pdf"))

# --- Figure 2F ---------------------------------
p5 <- generate_raincloud_with_ks(
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
ggsave(file.path(plots_dir, "Fig2F.pdf"), p5, width = 6, height = 6)
message("Figure 2F saved to ", file.path(plots_dir, "Fig2F.pdf"))

# --- Figure 2G ---------------------------------
p6 <- generate_raincloud_with_ks(
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
ggsave(file.path(plots_dir, "Fig2G.pdf"), p6, width = 6, height = 6)
message("Figure 2G saved to ", file.path(plots_dir, "Fig2G.pdf"))

# --- Figure 2H ---------------------------------
p7 <- generate_raincloud_with_ks(
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
ggsave(file.path(plots_dir, "Fig2H.pdf"), p7, width = 6, height = 6)
message("Figure 2H saved to ", file.path(plots_dir, "Fig2H.pdf"))

# --- Figure 2I ---------------------------------

p8 <- build_marker_summary_plot(
  pbmc_01_splot,
  marker = "CM%",
  disease_col = "Disease_type",
  cancergroup_col = "cancergroup",
  group_col = group_col,
  median_fun = median.cohorts.fx,
  sort_fun = sort.df.fx,
  splot_fun = Splot.fx,
  ylim_range = c(0, 100),
  y_label = "CM%",
  title_underlined = TRUE,
  dummy_groups = c("EMPTY1", "EMPTY2", "EMPTY3"),
  dummy_cancergroups = c("Lymphoma", "Solid tumors", "T-cell malignancies"),
  export_path = file.path(tables_dir, "CMpercent_cohort_median_summary.xlsx"),
  export_engine = "auto",
  export_overwrite = TRUE
) 

ggsave(file.path(plots_dir, "Fig2I.pdf"), p8$plot_final, width = 6, height = 4)
message("Figure 2I saved to ", file.path(plots_dir, "Fig2I.pdf"))