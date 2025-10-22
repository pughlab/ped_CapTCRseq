#!/usr/bin/env Rscript

# Title: Figure 4
#
# Inputs:
#   - data/metadata.rds (master metadata table for PBMC and cfDNA)
#   - data/chemo_cycles_cleaned.rds (chemotherapy cycles to derive base regimens)
#   - Helper palettes and themes from `R/functions/color_schemes.R` and `R/functions/ggplot2_theme.R`
#   - Misc helpers from `R/functions/Misc_functions.R` (optional)
# Outputs:
#   - Plots saved to `plots/` (e.g., Fig4A.pdf, Fig4B.pdf, ...)
# Usage:
#   Rscript R/Figure4.R

# --- Setup --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

# Load helper palettes and themes
source(file.path(getwd(), "R/functions", "color_schemes.R"))
source(file.path(getwd(), "R/functions", "ggplot2_theme.R"))
# Optional helpers used in the project
misc_path <- file.path(getwd(), "R/functions", "Misc_functions.R")
if (file.exists(misc_path)) source(misc_path)
plotfx_path <- file.path(getwd(), "R/functions", "Plotting_functions.R")
if (file.exists(plotfx_path)) source(plotfx_path)

# Project-rooted paths
project_root <- getwd()
plots_dir <- file.path(project_root, "plots")
input_metadata <- file.path(project_root, "data", "metadata.rds")
input_chemo <- file.path(project_root, "data", "chemo_cycles_cleaned.rds")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- Reusable helpers ----------------------------------------------------------

#' Derive base_regimen per cycle from chemo data (Cycle N Drugs)
#' @param soc data.frame chemo cycles with columns `Patient ID`, `Cancergroup`, `Disease`, and "Cycle N Drugs"
#' @return original df with base_regimen_cycle1..5 columns added
compute_base_regimens <- function(soc) {
  normalize_drug_names <- function(x) {
    x <- gsub("dactinomycin", "Dactinomycin", x, ignore.case = TRUE)
    x <- gsub("Cisplatinum", "Cisplatin", x, ignore.case = TRUE)
    x <- gsub("Cyatarabine", "Cytarabine", x, ignore.case = TRUE)
    x <- sub("\\*.*$", "", x)
    x <- sub("held during radiation", "", x, ignore.case = TRUE)
    trimws(x)
  }

  drug_cols <- grep("^Cycle [1-5] Drugs$", names(soc), value = TRUE)
  # Create base_regimen_cycleN columns
  for (cycle in 1:5) {
    col_cycle <- paste0("Cycle ", cycle, " Drugs")
    if (!col_cycle %in% drug_cols) {
      soc[[paste0("base_regimen_cycle", cycle)]] <- NA_character_
      next
    }
    cycle_data <- soc
    cycle_data[[col_cycle]] <- normalize_drug_names(cycle_data[[col_cycle]])

    vin_patients <- unique(cycle_data$`Patient ID`[grepl("(^|,)\\s*Vincristine(,|$)", cycle_data[[col_cycle]], ignore.case = TRUE)])
    mtx_patients <- unique(cycle_data$`Patient ID`[grepl("(^|,)\\s*Methotrexate(,|$)", cycle_data[[col_cycle]], ignore.case = TRUE)])
    cyclo_patients <- unique(cycle_data$`Patient ID`[grepl("(^|,)\\s*Cyclophosphamide(,|$)", cycle_data[[col_cycle]], ignore.case = TRUE)])
    norx_patients <- unique(cycle_data$`Patient ID`[grepl("(^|,)\\s*N/A(,|$)", cycle_data[[col_cycle]], ignore.case = TRUE)])

    col_name <- paste0("base_regimen_cycle", cycle)
    soc[[col_name]] <- sapply(soc$`Patient ID`, function(pat) {
      regimens <- character(0)
      if (pat %in% vin_patients) regimens <- c(regimens, "Vincristine")
      if (pat %in% mtx_patients) regimens <- c(regimens, "Methotrexate")
      if (pat %in% cyclo_patients) regimens <- c(regimens, "Cyclophosphamide")
      if (pat %in% norx_patients) regimens <- c(regimens, "No therapy")
      if (length(regimens) == 0) return(NA_character_)
      paste(regimens, collapse = "/")
    })
    soc[[col_name]][is.na(soc[[col_name]])] <- "Others"
  }
  soc
}

#' Build long table of base_regimen by cycle for joining to metadata
#' @param soc df from compute_base_regimens
#' @return long df with Patient (CHP_<id>), cycle (cycle1..), regimen
build_regimen_long <- function(soc) {
  soc %>%
    mutate(Patient = paste0("CHP_", `Patient ID`)) %>%
    pivot_longer(
      cols = starts_with("base_regimen_cycle"),
      names_to = "cycle",
      names_prefix = "base_regimen_",
      values_to = "regimen"
    ) %>%
    mutate(cycle = gsub("cycle", "cycle", cycle)) %>%
    select(Patient, cycle, regimen)
}

#' Compute "exposed_regimen" as baseline for X01 else previous cycle base_regimen
#' @param df metadata joined with regimen_long
#' @return df with exposed_regimen column
compute_exposed_regimen <- function(df, soc) {
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
      prev_cycle <- paste0("cycle", cycle_num - 1)
      prev_reg <- soc$regimen[soc$Patient == current_patient & soc$cycle == prev_cycle]
      if (length(prev_reg) > 0 && !is.na(prev_reg[1])) df$exposed_regimen[i] <- prev_reg[1]
    }
  }
  df
}

#' Prepare PCA input matrix from PBMC features per notebook
#' @param df PBMC-only df
#' @return list(mat_narm, df_narm)
prepare_pca_inputs <- function(df) {
  cols_flow <- c('ATC', 'Naïve%','SCM%','CM%','EM%','TE%','PD1%','LAG3%','TIM3%')
  cols_others <- c("observed_Shannon", "cfShannon")
  mymat <- df[, intersect(c(cols_flow, cols_others), names(df)), drop = FALSE]
  rownames(mymat) <- df$sample_id
  # remove rows all NA or zero
  mymat <- mymat[rowSums(is.na(mymat)) != ncol(mymat), , drop = FALSE]
  mymat <- mymat[rowSums(mymat, na.rm = TRUE) != 0, , drop = FALSE]
  # keep complete cases
  mymat_narm <- mymat[rowSums(is.na(mymat)) == 0, , drop = FALSE]
  df_narm <- df[df$sample_id %in% rownames(mymat_narm), , drop = FALSE]
  # rename for plotting
  colnames(mymat_narm)[colnames(mymat_narm) == "observed_Shannon"] <- "TCR Diversity"
  colnames(mymat_narm)[colnames(mymat_narm) == "cfShannon"] <- "cfTCR Diversity"
  list(mat = mymat_narm, df = df_narm)
}

#' Build biplot colored by exposed_regimen
#' @param pca prcomp object
#' @param df_narm dataframe aligned with pca$x rows
#' @param palette named vector for regimen colors (baseline included)
#' @return ggplot
build_biplot_by_regimen <- function(pca, df_narm, palette) {
  # PC scores
  df_narm$Dim1 <- pca$x[, 1][match(df_narm$sample_id, rownames(pca$x))]
  df_narm$Dim2 <- pca$x[, 2][match(df_narm$sample_id, rownames(pca$x))]

  p <- ggplot(df_narm, aes(x = Dim1, y = Dim2)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(colour = exposed_regimen), size = 3, alpha = 0.6, shape = 16) +
    myplot + myaxis +
    scale_color_manual(values = palette) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.position = "right"
    ) +
    labs(x = "Inverted Immune Checkpoint axis (PC1)", y = "Diversity axis (PC2)")
  p
}

#' Build ridgeline plots for PC1 and PC2 by regimen
#' @return list(ridge_x, ridge_y)
build_ridge_plots <- function(df_narm, palette) {
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Package 'ggridges' is required for ridge plots.")
  }
  ridge_x <- ggplot(df_narm, aes(x = Dim1, y = exposed_regimen, fill = exposed_regimen, point_color = exposed_regimen)) +
    ggridges::geom_density_ridges(
      alpha = 0.4,
      scale = 1,
      panel_scaling = TRUE,
      quantile_lines = FALSE,
      jittered_points = TRUE, point_shape = "|", point_size = 2.5,
      position = ggridges::position_points_jitter(height = 0)
    ) +
    geom_point(shape = "|", size = 2.5) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(x = "", y = "") +
    ggridges::theme_ridges(font_size = 11, grid = TRUE) +
    theme(legend.position = "none", axis.text.y = element_blank())

  ridge_y <- ggplot(df_narm, aes(x = Dim2, y = exposed_regimen, fill = exposed_regimen, point_color = exposed_regimen)) +
    ggridges::geom_density_ridges(
      alpha = 0.4,
      scale = 1,
      panel_scaling = TRUE,
      quantile_lines = FALSE,
      jittered_points = TRUE, point_shape = "-", point_size = 4,
      position = ggridges::position_points_jitter(height = 0)
    ) +
    geom_point(shape = "-", size = 4) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(x = "", y = "") +
    ggridges::theme_ridges(font_size = 11, grid = TRUE) +
    theme(legend.position = "none", axis.text.x = element_blank()) +
    coord_flip()

  list(ridge_x = ridge_x, ridge_y = ridge_y)
}

#' Combine ridge plots with biplot similar to notebook cowplot layout
#' @return cowplot object
combine_ridge_and_biplot <- function(ridge_x, biplot, ridge_y) {
  cowplot::plot_grid(
    ridge_x + theme(axis.title.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")), NULL, NULL,
    NULL, NULL, NULL,
    biplot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
    NULL,
    ridge_y + theme(axis.title.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")),
    ncol = 3, nrow = 3,
    align = "hv",
    rel_widths = c(1.5, -0.2, 1),
    rel_heights = c(1, -0.2, 1.5)
  )
}

# Default regimen palette including baseline
regimen_palette <- c(
  "baseline" = "black",
  "Vincristine" = "#E41A1C",
  "Methotrexate" = "#377EB8",
  "Cyclophosphamide" = "#4DAF4A",
  "Vincristine/Methotrexate" = "#984EA3",
  "Vincristine/Cyclophosphamide" = "#fa9734",
  "Vincristine/Methotrexate/Cyclophosphamide" = "#f781bf",
  "No therapy" = "grey",
  "Others" = "#6c6c6c"
)

# --- Load inputs and prepare data ---------------------------------------------

metadata <- readr::read_rds(input_metadata)
chemo <- readr::read_rds(input_chemo)

# Create base regimen by cycle from chemo and derive long mapping
chemo <- compute_base_regimens(chemo)
regimen_long <- build_regimen_long(chemo)

# PBMC/cfDNA subsets and Shannon mapping
pbmc <- metadata[metadata$sampletype == "PBMC", ]
cfdna <- metadata[metadata$sampletype == "cfDNA", ]
# Map cfDNA Shannon to PBMC by sample_id when available
if ("observed_Shannon" %in% names(cfdna)) {
  pbmc$cfShannon <- cfdna$observed_Shannon[match(pbmc$sample_id, cfdna$sample_id)]
}

# Join regimen by cycle_mapped (X05 -> cycle5)
pbmc$cycle_mapped <- paste0("cycle", gsub("X0?", "", pbmc$cycle))
pbmc <- pbmc %>% left_join(regimen_long, by = c("Patient" = "Patient", "cycle_mapped" = "cycle"))
pbmc$cycle_mapped <- NULL
names(pbmc)[names(pbmc) == "regimen"] <- "regimen"

# Compute exposed_regimen from previous cycle regimen
pbmc <- compute_exposed_regimen(pbmc, regimen_long)

# Ensure exposed_regimen factor order baseline first then named palette order
pbmc$exposed_regimen <- factor(pbmc$exposed_regimen, levels = names(regimen_palette))

# Prepare PCA inputs
pca_in <- prepare_pca_inputs(pbmc)
mymat_narm <- pca_in$mat
df_narm <- pca_in$df

# PCA
set.seed(123)
if (ncol(mymat_narm) < 2 || nrow(mymat_narm) < 2) {
  stop("Insufficient non-missing features for PCA in Figure 4.")
}
pca <- prcomp(mymat_narm, center = TRUE, scale. = TRUE)

# Build plots mirroring notebook (biplot by regimen and ridge plots)
p_biplot <- build_biplot_by_regimen(pca, df_narm, regimen_palette)
ridges <- build_ridge_plots(df_narm, regimen_palette)
combined <- combine_ridge_and_biplot(ridges$ridge_x, p_biplot, ridges$ridge_y)

# --- Save outputs --------------------------------------------------------------

ggsave(file.path(plots_dir, "Fig4A_biplot_regimen.pdf"), p_biplot, width = 7, height = 7)
# Combined layout similar to pca_regimen_ridgeplot_cowplot.pdf in notebook
ggsave(file.path(plots_dir, "Fig4B_pca_ridge_combined.pdf"), combined, width = 10, height = 8)

ggsave(file.path(plots_dir, "Fig4C_ridge_PC1.pdf"), ridges$ridge_x, width = 8, height = 4)
# ridge_y is rotated by coord_flip(), so invert dimensions
ggsave(file.path(plots_dir, "Fig4D_ridge_PC2.pdf"), ridges$ridge_y, width = 4, height = 8)

message("Figure 4 outputs saved to ", plots_dir)
