#!/usr/bin/env Rscript

# Title: Figure 3
#
# Inputs:
#   - `data/metadata.rds`
#   - `data/soc_good_with_baseregimen.rds`
#   - Color and theme helpers from `R/functions/color_schemes.R` and `R/functions/ggplot2_theme.R`
#   - Misc helpers for p-value formatting and delta calculation from `R/functions/Misc_functions.R`
#   - Plotting helpers from `R/functions/Plotting_functions.R`
# Outputs:
#   - Plots saved to `plots/` (Fig3A.pdf ... Fig3E.pdf)
# Usage:
#   Rscript R/Figure3.R

# --- Setup --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

# Load helper palettes and themes
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "color_schemes.R"))
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "ggplot2_theme.R"))
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "Misc_functions.R"))
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "Plotting_functions.R"))

# Project-rooted paths
project_root <- "/Users/anabbi/git/ped_CapTCRseq/"
plots_dir <- file.path(project_root, "plots")
input_rds <- file.path(project_root, "data", "metadata.rds")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- Reusable helpers specific to Fig3 ----------------------------------------

#' Segment-annotated spider plot where line segments are colored by next-cycle regimen
#' @param df_diff Long or wide df containing Patient, cycle, Diff (delta from baseline), exposed_regimen
#' @param x_var Cycle column name (default "cycle")
#' @param color_var Column for coloring segments (default "exposed_regimen")
#' @param colpal Named palette mapping regimen strings to colors
#' @param facet_formula Optional facet formula string (e.g., "~ marker + Disease_type")
#' @param y_lab Y-axis label
#' @return ggplot object
segment_spider_plot <- function(df_diff,
                                x_var = "cycle",
                                color_var = "exposed_regimen",
                                colpal = NULL,
                                facet_formula = NULL,
                                y_lab = NULL) {
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    warning("Package 'ggh4x' not installed; facet_wrap2 will be unavailable.")
  }

  # Prepare segment data frame to color by next segment regimen
  # Group by Patient and marker (if marker exists) to handle multi-marker plots correctly
  group_vars <- "Patient"
  arrange_vars <- c("Patient", x_var)
  if ("marker" %in% names(df_diff)) {
    group_vars <- c("Patient", "marker")
    arrange_vars <- c("Patient", x_var)  # Keep Patient and cycle for arrangement
  }
  
  segments_data <- df_diff %>%
    arrange(across(all_of(arrange_vars))) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      x_start = .data[[x_var]],
      y_start = .data$Diff,
      x_end   = dplyr::lead(.data[[x_var]]),
      y_end   = dplyr::lead(.data$Diff),
      segment_regimen = dplyr::lead(.data[[color_var]])
    ) %>%
    filter(!is.na(x_end)) %>%
    ungroup()

  p <- ggplot(df_diff, aes(x = .data[[x_var]], y = .data$Diff)) +
    geom_segment(
      data = segments_data,
      aes(x = .data$x_start, y = .data$y_start, xend = .data$x_end, yend = .data$y_end, color = .data$segment_regimen),
      linewidth = 1
    ) +
    geom_point(cex = 2) +
    scale_color_manual(values = colpal) +
    myplot + myaxis +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")

  if (!is.null(facet_formula)) {
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::facet_wrap2(as.formula(facet_formula), axes = "all", remove_labels = "all") +
        theme(
          axis.title.x = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(size = 13),
          strip.background.x = element_blank(),
          strip.placement = "outside",
          axis.line = element_line(colour = "black")
        )
    } else {
      p <- p + facet_wrap(as.formula(facet_formula), scales = "free")
    }
  }

  if (!is.null(y_lab)) {
    p <- p + ylab(y_lab)
  }

  p
}

# Default regimen palette
regimen_palette <- c(
  "Vincristine" = "#E41A1C",
  "Methotrexate" = "#377EB8",
  "Cyclophosphamide" = "#4DAF4A",
  "Vincristine/Methotrexate" = "#984EA3",
  "Vincristine/Cyclophosphamide" = "#fa9734",
  "Vincristine/Methotrexate/Cyclophosphamide" = "#f781bf",
  "No therapy" = "grey",
  "Others" = "#6c6c6c"
)

# Override delta_basespiderplot.fx to use segments for coloring by next cycle regimen
delta_basespiderplot.fx <- function(df_diff, var1, clrby, colpal) {
  # Create segments data for coloring lines by regimen
  segments_data <- df_diff %>%
    arrange(Patient, eval(parse(text = var1))) %>%
    group_by(Patient) %>%
    mutate(
      x_start = eval(parse(text = var1)),
      y_start = Difference,
      x_end = lead(eval(parse(text = var1))),
      y_end = lead(Difference),
      segment_regimen = lead(exposed_regimen)
    ) %>%
    filter(!is.na(x_end)) %>%
    ungroup()
  
  p0 <- ggplot(df_diff, aes(x = eval(parse(text = var1)), y = Difference)) +
    geom_segment(
      data = segments_data,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = segment_regimen),
      linewidth = 1
    ) +
    geom_point(cex = 2) +
    scale_color_manual(values = colpal) +
    myplot +
    myaxis +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  return(p0)
}

# --- Load inputs ---------------------------------------------------------------

metadata <- readr::read_rds(input_rds)

# Load regimen data and merge with metadata
datapath <- file.path(project_root, "data/")
soc_good <- readr::read_rds(paste0(datapath, "soc_good_with_baseregimen.rds"))
soc_good_long <- soc_good %>%
  pivot_longer(
    cols = starts_with("base_regimen_"),
    names_to = "cycle",
    values_to = "regimen",
    names_prefix = "base_regimen_"
  ) %>%
  select(Patient, cycle, regimen)

head(soc_good_long)
# Join soc_good_long to metadata
# First, create a mapping between cycle formats (X05 -> cycle5)
metadata$cycle_mapped <- paste0("cycle", gsub("X0?", "", metadata$cycle))

# Join the datasets
metadata <- metadata %>%
  left_join(soc_good_long, by = c("Patient" = "Patient", "cycle_mapped" = "cycle"))

# Clean up the temporary column
metadata$cycle_mapped <- NULL

# Create exposed_regimen column, for each patient, if cycle is X01, then it is baseline, otherwise it is the previous cycle's regimen
## Use the soc_good dataframe to get the regimen for the previous cycle
## First, create a mapping between cycle formats (X05 -> cycle5)
## be explicit, eg get regimen of X04 for exposed_regimen of X05
# Create exposed_regimen column
metadata$exposed_regimen <- NA

# For each sample, determine the exposed regimen
for (i in 1:nrow(metadata)) {
  current_patient <- metadata$Patient[i]
  current_cycle <- metadata$cycle[i]
  
  # If it's X01 (first cycle), exposed_regimen is NA (baseline)
  if (current_cycle == "X01") {
    metadata$exposed_regimen[i] <- "baseline"
  } else {
    # Extract cycle number and get previous cycle
    cycle_num <- as.numeric(gsub("X0?", "", current_cycle))
    prev_cycle_num <- cycle_num - 1
    
    # Format previous cycle (e.g., 4 -> "cycle4")
    prev_cycle_col <- paste0("cycle", prev_cycle_num)
    
    # Get the base regimen for the previous cycle from soc_good
    prev_regimen <- as.character(soc_good[soc_good$Patient == current_patient, paste0("base_regimen_", prev_cycle_col), drop = T])
    
    if (length(prev_regimen) > 0 && !is.na(prev_regimen)) {
      metadata$exposed_regimen[i] <- prev_regimen
    }
  }
}

# Derive cfDNA/PBMC subsets similar to notebook
# Exclude T-cell malignancies
metadata <- metadata[metadata$cancergroup != "T-cell malignancies", ]

# PBMC only for cellular markers; ensure CD3 > 0 and non-missing
pbmc <- metadata[metadata$sampletype == "PBMC", ]
pbmc <- pbmc[!is.na(pbmc$CD3) & pbmc$CD3 > 0, ]

# cfDNA for diversity
cfdna <- metadata[metadata$sampletype == "cfDNA", ]

# --- Prepare deltas ------------------------------------------------------------

# Fig3B-C: cfDNA diversity deltas per disease panels
# Use Misc_functions::calculate_delta.fx which computes per-patient delta from X01 baseline
if (!"log10shann_scaled" %in% names(cfdna) && "log10shann" %in% names(cfdna)) {
  # Scale per overall to mimic notebook
  cfdna$log10shann_scaled <- scale(cfdna$log10shann)
}

diff_log10shann_cfdna <- tryCatch({
  calculate_delta.fx(cfdna, "cycle", "log10shann_scaled")
}, error = function(e) {
  message("calculate_delta.fx on cfDNA failed: ", conditionMessage(e))
  cfdna
})

diff_log10shann_cfdna <- diff_log10shann_cfdna[!is.na(diff_log10shann_cfdna$Difference), ]
# Add grp column for faceting
diff_log10shann_cfdna$grp <- "cfTCR Diversity"

# Fig3A/D/E: PBMC deltas for cell subsets and immune checkpoints
pbmc1 <- pbmc

# Define markers per notebook (using backticked names to match column names)
mycells <- c("`Naïve%`", "`SCM%`", "`CM%`", "`EM%`", "`TE%`", "`PD1%`", "`LAG3%`", "`TIM3%`")
cells_markers <- c("Naïve%", "SCM%", "CM%", "EM%", "TE%")
ic_markers <- c("PD1%", "LAG3%", "TIM3%")
all_markers <- mycells

# Create scaled columns and per-sample differences vs baseline for each marker
for (i in 1:length(mycells)) {
  pbmc1$grp <- gsub("`", "", mycells[i])
  mycol <- paste0("Diff_", gsub("%", "", pbmc1$grp[1]))
  myvar <- paste0("scaled_", gsub("%", "", pbmc1$grp[1]))
  pbmc1[[myvar]] <- NA
  pbmc1[[myvar]] <- scale(pbmc1[[pbmc1$grp[1]]])
  mydiff <- calculate_delta.fx(pbmc1, "cycle", myvar)
  pbmc1[[mycol]] <- mydiff$Difference[match(pbmc1$sample_id, mydiff$sample_id)]
}

# Long format for plotting
pbmc1_long <- pbmc1 %>%
  select(
    sample_id, Patient, cycle, cancergroup, Disease_type, Age, regimen, exposed_regimen,
    Relapse, Diff_Naïve, Diff_SCM, Diff_CM, Diff_EM, Diff_TE,
    Diff_PD1, Diff_LAG3, Diff_TIM3
  ) %>%
  tidyr::pivot_longer(cols = starts_with("Diff_"), names_to = "marker", values_to = "Diff")

pbmc1_long <- pbmc1_long[!is.na(pbmc1_long$Diff), ]
pbmc1_long$marker <- gsub("Diff_", "", pbmc1_long$marker)
pbmc1_long$marker <- paste0(pbmc1_long$marker, "%")

pbmc1_long_cells <- pbmc1_long[pbmc1_long$marker %in% cells_markers, ]
pbmc1_long_ic <- pbmc1_long[pbmc1_long$marker %in% ic_markers, ]

# Convert marker to character for filtering
pbmc1_long_cells$marker <- as.character(pbmc1_long_cells$marker)
pbmc1_long_ic$marker <- as.character(pbmc1_long_ic$marker)

# --- Figures ------------------------------------------------------------------

# Fig3A: ALL/HR ALL naive+CM deltas, spider segments colored by next regimen
fig3a_df <- pbmc1_long_cells[
  pbmc1_long_cells$Disease_type %in% c("ALL", "HR ALL") & pbmc1_long_cells$marker %in% c("Naïve%", "CM%"),
]
# Factorize marker for correct ordering
fig3a_df$marker <- factor(fig3a_df$marker, levels = c("Naïve%", "CM%"))

p_fig3a <- segment_spider_plot(
  fig3a_df,
  x_var = "cycle",
  color_var = "exposed_regimen",
  colpal = regimen_palette,
  facet_formula = "~ marker + Disease_type",
  y_lab = "D T-cell subsets"
)

ggsave(file.path(plots_dir, "Fig3A.pdf"), cowplot::plot_grid(p_fig3a, labels = "A"), width = 4, height = 6)
message("Figure 3A saved to ", file.path(plots_dir, "Fig3A.pdf"))

# Fig3B: Leukemia panel (ALL/HR ALL/AML) cfDNA diversity
fig3b_df <- diff_log10shann_cfdna[diff_log10shann_cfdna$Disease_type %in% c("ALL", "HR ALL", "AML"), ]

p_fig3b <- delta_basespiderplot.fx(fig3b_df, "cycle", "exposed_regimen", regimen_palette) +
  theme(legend.position = "none")
if (requireNamespace("ggh4x", quietly = TRUE)) {
  p_fig3b <- p_fig3b + ggh4x::facet_wrap2(~ grp + Disease_type, axes = "all", remove_labels = "all") +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.title.x = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(size = 13),
      strip.background.x = element_blank(),
      strip.placement = "outside"
    ) +
    ylab("D cfTCR diversity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_y_continuous(limits = c(-2, 4), breaks = c(-2, 0, 2, 4))
}

ggsave(file.path(plots_dir, "Fig3B.pdf"), cowplot::plot_grid(p_fig3b, labels = "B"), width = 6, height = 3)
message("Figure 3B saved to ", file.path(plots_dir, "Fig3B.pdf"))

# Fig3C: Lymphoma panel (BL, HD) cfDNA diversity
fig3c_df <- diff_log10shann_cfdna[diff_log10shann_cfdna$Disease_type %in% c("BL", "HD"), ]

p_fig3c <- delta_basespiderplot.fx(fig3c_df, "cycle", "exposed_regimen", regimen_palette) +
  theme(legend.position = "none")
if (requireNamespace("ggh4x", quietly = TRUE)) {
  p_fig3c <- p_fig3c + ggh4x::facet_wrap2(~ grp + Disease_type, axes = "all", remove_labels = "all") +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.title.x = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(size = 13),
      strip.background.x = element_blank(),
      strip.placement = "outside"
    ) +
    ylab("D cfTCR diversity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")
}

ggsave(file.path(plots_dir, "Fig3C.pdf"), cowplot::plot_grid(p_fig3c, labels = "C"), width = 4, height = 3)
message("Figure 3C saved to ", file.path(plots_dir, "Fig3C.pdf"))

# Fig3D: Solid tumors (ERMS, EWS, HB, NB, RMS) CM% deltas
fig3d_df <- pbmc1_long_cells[
  pbmc1_long_cells$Disease_type %in% c("ERMS", "EWS", "HB", "NB", "RMS") &
    pbmc1_long_cells$marker == "CM%",
]

p_fig3d <- segment_spider_plot(
  fig3d_df,
  x_var = "cycle",
  color_var = "exposed_regimen",
  colpal = regimen_palette,
  facet_formula = "~ marker + Disease_type",
  y_lab = "D T-cell subsets"
)

ggsave(file.path(plots_dir, "Fig3D.pdf"), cowplot::plot_grid(p_fig3d, labels = "D"), width = 12, height = 3)
message("Figure 3D saved to ", file.path(plots_dir, "Fig3D.pdf"))

# Fig3E: Solid subset (ERMS, OS, HB) PD1% & LAG3% deltas
fig3e_df <- pbmc1_long_ic[
  pbmc1_long_ic$Disease_type %in% c("ERMS", "OS", "HB") &
    pbmc1_long_ic$marker %in% c("PD1%", "LAG3%"),
]

p_fig3e <- segment_spider_plot(
  fig3e_df,
  x_var = "cycle",
  color_var = "exposed_regimen",
  colpal = regimen_palette,
  facet_formula = "~ marker + Disease_type",
  y_lab = "D immune checkpoints"
)

ggsave(file.path(plots_dir, "Fig3E.pdf"), cowplot::plot_grid(p_fig3e, labels = "E"), width = 6, height = 6)
message("Figure 3E saved to ", file.path(plots_dir, "Fig3E.pdf"))
