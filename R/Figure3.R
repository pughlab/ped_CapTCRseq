#!/usr/bin/env Rscript

# Title: Figure 3
#
# Inputs:
#   - `data/metadata.rds`
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
source(file.path(getwd(), "R/functions", "color_schemes.R"))
source(file.path(getwd(), "R/functions", "ggplot2_theme.R"))
source(file.path(getwd(), "R/functions", "Misc_functions.R"))
source(file.path(getwd(), "R/functions", "Plotting_functions.R"))

# Project-rooted paths
project_root <- getwd()
plots_dir <- file.path(project_root, "plots")
input_rds <- file.path(project_root, "data", "metadata.rds")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- Reusable helpers specific to Fig3 ----------------------------------------

#' Build regimen exposure per sample using previous cycle
#' @param df Soc/metadata-like df with columns Patient, cycle (X01..), and regimen
#' @return df with exposed_regimen column ("baseline" for X01 else previous cycle regimen)
compute_exposed_regimen <- function(df) {
  df$exposed_regimen <- NA_character_
  for (i in seq_len(nrow(df))) {
    current_patient <- df$Patient[i]
    current_cycle <- df$cycle[i]
    if (is.na(current_patient) || is.na(current_cycle)) next
    if (current_cycle == "X01") {
      df$exposed_regimen[i] <- "baseline"
    } else {
      cycle_num <- suppressWarnings(as.numeric(gsub("X0?", "", current_cycle)))
      if (!is.finite(cycle_num)) next
      prev_cycle <- paste0("X0", cycle_num - 1)
      prev_row <- which(df$Patient == current_patient & df$cycle == prev_cycle)
      if (length(prev_row) > 0) {
        prev_regimen <- as.character(df$regimen[prev_row[1]])
        if (!is.na(prev_regimen) && nzchar(prev_regimen)) {
          df$exposed_regimen[i] <- prev_regimen
        }
      }
    }
  }
  df
}

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
  segments_data <- df_diff %>%
    arrange(.data$Patient, .data[[x_var]]) %>%
    group_by(.data$Patient) %>%
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

# Default regimen palette from notebook
regimen_palette <- c(
  "Vincristine" = "#E41A1C",
  "Methotrexate" = "#377EB8",
  "Cyclophosphamide" = "#4DAF4A",
  "Vincristine/Methotrexate" = "#984EA3",
  "Vincristine/Cyclophosphamide" = "#fa9734",
  "Vincristine/Methotrexate/Cyclophosphamide" = "#f781bf",
  "No therapy" = "grey",
  "Others" = "black",
  "baseline" = "black"
)

# --- Load inputs ---------------------------------------------------------------

metadata <- readr::read_rds(input_rds)

# Derive cfDNA/PBMC subsets similar to notebook
# Exclude T-cell malignancies
metadata <- metadata[metadata$cancergroup != "T-cell malignancies", ]

# PBMC only for cellular markers; ensure CD3 > 0 and non-missing
pbmc <- metadata[metadata$sampletype == "PBMC", ]
pbmc <- pbmc[!is.na(pbmc$CD3) & pbmc$CD3 > 0, ]

# cfDNA for diversity
cfdna <- metadata[metadata$sampletype == "cfDNA", ]

# Ensure regimen exists; if absent, create from available columns if present
if (!"regimen" %in% names(metadata)) {
  metadata$regimen <- NA
}

# Compute exposed_regimen
metadata <- compute_exposed_regimen(metadata)
pbmc <- metadata[metadata$sampletype == "PBMC", ]
pbmc <- pbmc[!is.na(pbmc$CD3) & pbmc$CD3 > 0, ]
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

# Fig3A/D/E: PBMC deltas for cell subsets and immune checkpoints
pbmc1 <- pbmc

# Define markers per notebook
cells_markers <- c("Naïve%", "SCM%", "CM%", "EM%", "TE%")
ic_markers <- c("PD1%", "LAG3%", "TIM3%")
all_markers <- c(cells_markers, ic_markers)

# Create scaled columns and per-sample differences vs baseline for each marker
for (mk in all_markers) {
  # scaled column name and diff column name
  mk_clean <- gsub("%", "", mk)
  scaled_col <- paste0("scaled_", mk_clean)
  diff_col <- paste0("Diff_", mk_clean)

  if (!scaled_col %in% names(pbmc1)) pbmc1[[scaled_col]] <- NA_real_
  # Some markers have non-syntactic names; use get with backticks via [[ ]]
  if (mk %in% names(pbmc1)) {
    pbmc1[[scaled_col]] <- as.numeric(scale(pbmc1[[mk]]))
  } else if (paste0("`", mk, "`") %in% names(pbmc1)) {
    pbmc1[[scaled_col]] <- as.numeric(scale(pbmc1[[paste0("`", mk, "`")]]))
  }
  mydiff <- tryCatch({
    calculate_delta.fx(pbmc1, "cycle", scaled_col)
  }, error = function(e) {
    message("calculate_delta.fx failed for ", mk, ": ", conditionMessage(e))
    NULL
  })
  if (!is.null(mydiff)) {
    pbmc1[[diff_col]] <- mydiff$Difference[match(pbmc1$sample_id, mydiff$sample_id)]
  }
}

# Long format for plotting
pbmc1_long <- pbmc1 %>%
  select(
    sample_id, Patient, cycle, cancergroup, Disease_type, Age, regimen, exposed_regimen,
    Relapse, starts_with("Diff_")
  ) %>%
  tidyr::pivot_longer(cols = starts_with("Diff_"), names_to = "marker", values_to = "Diff") %>%
  filter(!is.na(Diff))

pbmc1_long$marker <- gsub("Diff_", "", pbmc1_long$marker)
pbmc1_long$marker <- paste0(pbmc1_long$marker, "%")

pbmc1_long_cells <- pbmc1_long[pbmc1_long$marker %in% cells_markers, ]
pbmc1_long_ic <- pbmc1_long[pbmc1_long$marker %in% ic_markers, ]

# --- Figures ------------------------------------------------------------------

# Fig3A: ALL/HR ALL naive+CM deltas, spider segments colored by next regimen
fig3a_df <- pbmc1_long_cells[
  pbmc1_long_cells$Disease_type %in% c("ALL", "HR ALL") & pbmc1_long_cells$marker %in% c("Naïve%", "CM%"),
]

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
fig3b_df$grp <- "cfTCR Diversity"

p_fig3b <- delta_basespiderplot.fx(fig3b_df, "cycle", "exposed_regimen", regimen_palette) +
  theme(legend.position = "none")
if (requireNamespace("ggh4x", quietly = TRUE)) {
  p_fig3b <- p_fig3b + ggh4x::facet_wrap2(~ grp + Disease_type, axes = "all", remove_labels = "all") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x = element_blank(), legend.position = "none")
}

ggsave(file.path(plots_dir, "Fig3B.pdf"), cowplot::plot_grid(p_fig3b, labels = "B"), width = 6, height = 3)
message("Figure 3B saved to ", file.path(plots_dir, "Fig3B.pdf"))

# Fig3C: Lymphoma panel (BL, HD) cfDNA diversity
fig3c_df <- diff_log10shann_cfdna[diff_log10shann_cfdna$Disease_type %in% c("BL", "HD"), ]
fig3c_df$grp <- "cfTCR Diversity"

p_fig3c <- delta_basespiderplot.fx(fig3c_df, "cycle", "exposed_regimen", regimen_palette) +
  theme(legend.position = "none")
if (requireNamespace("ggh4x", quietly = TRUE)) {
  p_fig3c <- p_fig3c + ggh4x::facet_wrap2(~ grp + Disease_type, axes = "all", remove_labels = "all") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x = element_blank(), legend.position = "none")
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
