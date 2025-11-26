#!/usr/bin/env Rscript

# Title: Figure 4
#
# Inputs:
#   - `data/metadata.rds`
#   - `data/soc_good_with_baseregimen.rds`
#   - Color and theme helpers from `R/functions/color_schemes.R` and `R/functions/ggplot2_theme.R`
#   - Misc helpers for p-value formatting and delta calculation from `R/functions/Misc_functions.R`
# Outputs:
#   - Plots saved to `plots/` 
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
  library(ggbeeswarm)
  library(ggridges)
  library(factoextra)
  library(ggh4x)
  library(grid)
})

# Set locale
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

# Load helper palettes and themes
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "color_schemes.R"))
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "ggplot2_theme.R"))
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "Misc_functions.R"))
source(file.path("/Users/anabbi/git/ped_CapTCRseq/", "R/functions", "Plotting_functions.R"))

# Project-rooted paths
project_root <- "/Users/anabbi/git/ped_CapTCRseq/"
# project_root <- getwd()
plots_dir <- file.path(project_root, "plots")
input_rds <- file.path(project_root, "data", "metadata.rds")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- Reusable helpers specific to Fig4 ----------------------------------------

# Default regimen palette
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

# Split by sample type
metadata_pbmc <- metadata[metadata$sampletype == "PBMC", ]
metadata_cfdna <- metadata[metadata$sampletype == "cfDNA", ]

# Map cfDNA Shannon to PBMC
metadata_pbmc$cfShannon <- metadata_cfdna$observed_Shannon[match(metadata_pbmc$sample_id, metadata_cfdna$sample_id)]

# --- PCA Analysis --------------------------------------------------------------

# Define feature columns
cols_flow <- c('ATC', 'Naïve%','SCM%','CM%','EM%','TE%','PD1%','LAG3%','TIM3%')
cols_others <- c("observed_Shannon", "cfShannon")

# Create feature matrix
mymat <- metadata_pbmc[, c(cols_flow, cols_others)]
rownames(mymat) <- metadata_pbmc[["sample_id"]]

# Remove rows with all NAs or only zeros
mymat <- mymat[rowSums(is.na(mymat)) != ncol(mymat), ]
mymat <- mymat[rowSums(mymat, na.rm = TRUE) != 0, ]

# Keep only complete cases
mymat_narm <- mymat[rowSums(is.na(mymat)) == 0, ]

# Get corresponding metadata
df_narm <- metadata_pbmc[metadata_pbmc$sample_id %in% rownames(mymat_narm), ]

# Rename columns for plotting
colnames(mymat_narm)[colnames(mymat_narm) == "observed_Shannon"] <- "TCR Diversity"
colnames(mymat_narm)[colnames(mymat_narm) == "cfShannon"] <- "cfTCR Diversity"

# load PCA from RDS
pca <- readr::read_rds(file.path(datapath, "pca.rds"))


df_narm$Dim1 <- pca$x[,1][match( df_narm$sample_id, rownames(pca$x))]
df_narm$Dim2 <- pca$x[,2][match( df_narm$sample_id, rownames(pca$x))]

# Factor exposed_regimen with baseline first
df_narm$exposed_regimen <- factor(df_narm$exposed_regimen, levels = names(regimen_palette))

# Factor exposed_regimen with baseline first and Others last
df_narm$exposed_regimen <- factor(df_narm$exposed_regimen, 
                                  levels = names(regimen_palette))
table(df_narm$exposed_regimen)


# Factor Disease_type for faceting
df_narm$Disease_type <- factor(df_narm$Disease_type, 
                              levels = c("ALL", "AML", "CML", "HR ALL",
                                        "BL", "DLBCL", "HD", "PMBCL",
                                        "ARMS", "ERMS", "EWS", "HB", "NB", "OS", "WILMS", " ",
                                        "ALCL", "T-ALL"))

summary(pca)$importance[2, 1:5] * 100
# --- Figures -------------------------------------------------------------------

# Variable plot (arrows showing feature contributions)
p_var <- fviz_pca_var(pca, repel = FALSE, geom = "arrow", axes.linetype = NA) +
  geom_text_repel(aes(label = rownames(pca$rotation)),
                  point.padding = unit(0.5, "lines"),
                  min.segment.length = unit(0, "lines"),
                  nudge_x = 0.01, nudge_y = 0.01,
                  direction = "y") +
  myaxis + myplot

ggsave(file.path(plots_dir, "Variable_plot.pdf"), p_var, width = 7, height = 7)
message("Variable plot saved to ", file.path(plots_dir, "Variable_plot.pdf"))

# Biplot colored by regimen
p_biplot <- fviz_pca_biplot(pca,
    geom.var = c("arrow", "text"),
    geom.ind = "point",
    alpha.ind = 0, alpha.var = 1,
    col.var = "black",
    repel = TRUE,
    labelsize = 6,
    arrowsize = 0.5,
    addEllipses = FALSE, 
    mean.point = FALSE,
    title = "") +
    geom_point(aes(colour = df_narm$exposed_regimen), size = 4, alpha = 0.5, shape = 16) +
    scale_color_manual(values = regimen_palette) + myplot + myaxis +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20), 
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20)) +
    scale_x_continuous(limits=c(-6,3)) + scale_y_continuous(limits=c(-3,4))
# increase size of the legends
p_biplot <- p_biplot + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10))
p_biplot <- p_biplot + labs(x = "Inverted Immune Checkpoint axis\n(PC1 = 31.9%)", y = "Diversity axis\n(PC2 = 20%)")
ggsave(file.path(plots_dir, "Biplot.pdf"), cowplot::plot_grid(p_biplot + 
theme(legend.position = "none"), labels = "A"), width = 7, height = 7)
message("Biplot saved to ", file.path(plots_dir, "Biplot.pdf"))

# Legend for biplot
pdf(file.path(plots_dir, "Fig4A_legend.pdf"), width = 5, height = 5)
grid.draw(get_legend(p_biplot + theme(legend.background = element_blank(), legend.position = "right") + labs(colour = "Regimen")) )
dev.off()


# Fig4A : Ridge plots for PC1 and PC2
ridge_x <- ggplot(df_narm, aes(x = Dim1, y = exposed_regimen, 
                               fill = exposed_regimen, point_color = exposed_regimen)) +
  geom_density_ridges(
    alpha = 0.4,
    scale = 1,
    panel_scaling = TRUE,
    quantile_lines = FALSE, 
    jittered_points = TRUE, 
    point_shape = "|", 
    point_size = 3, 
    position = position_points_jitter(height = 0)
  ) +
  geom_point(shape = "|", size = 3) +
  scale_fill_manual(values = regimen_palette) +
  scale_color_manual(values = regimen_palette) +
  labs(x = "", y = "") +
  theme_ridges(font_size = 13, grid = TRUE) +
  theme(legend.position = "none", axis.text.y = element_blank()) + 
  scale_x_continuous(limits = c(-6, 3))

ridge_y <- ggplot(df_narm, aes(x = Dim2, y = exposed_regimen, 
                               fill = exposed_regimen, point_color = exposed_regimen)) +
  geom_density_ridges(
    alpha = 0.4,
    scale = 1,
    panel_scaling = TRUE,
    quantile_lines = FALSE, 
    jittered_points = TRUE, 
    point_shape = "-", 
    point_size = 5, 
    position = position_points_jitter(height = 0)
  ) +
  geom_point(shape = "-", size = 5) +
  scale_fill_manual(values = regimen_palette) +
  scale_color_manual(values = regimen_palette) +
  labs(x = "", y = "") +
  theme_ridges(font_size = 13, grid = TRUE) +
  theme(legend.position = "none", axis.text.x = element_blank()) + 
  scale_x_continuous(limits = c(-3, 4)) + 
  coord_flip()

# Fig4A: Combined plot with ridges and biplot
combined_plot <- plot_grid(
  ridge_x + theme(axis.title.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")),
  NULL, 
  NULL, 
  NULL, 
  NULL, 
  NULL, 
  p_biplot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
  NULL,
  ridge_y + theme(axis.title.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")),
  ncol = 3, nrow = 3,
  align = "hv",
  rel_widths = c(1.5, -0.2, 1),
  rel_heights = c(1, -0.2, 1.5),
  labels = "A"
)

ggsave(file.path(plots_dir, "Fig4A.pdf"), combined_plot, width = 10, height = 8)
message("Figure 4A saved to ", file.path(plots_dir, "Fig4A.pdf"))



# Fig4B: Faceted plot by disease type
design <- "
 ABCD
 EFGH
 IJKL
 MNO#
 PQ##
"

baseplot <- ggplot(df_narm, aes(x = Dim1, y = Dim2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(colour = exposed_regimen), size = 4, alpha = 0.5, shape = 16) +
  myplot + myaxis +
  scale_color_manual(values = regimen_palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20), 
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.position = "right",
        strip.text.x = element_text(size = 13),
        strip.background.x = element_blank(),
        strip.placement = "outside",
        axis.line = element_line()) +
  labs(x = "Inverted Immune Checkpoint axis\n(PC1 = 31.9%)", 
       y = "Diversity axis\n(PC2 = 20%)") + 
  ggh4x::facet_manual(~ Disease_type, design, scales = "free") +
  scale_x_continuous(limits = c(-6, 3)) + 
  scale_y_continuous(limits = c(-3, 4))

ggsave(file.path(plots_dir, "Fig4B.pdf"), cowplot::plot_grid(baseplot + theme(legend.position = "none"), labels = "B"), width = 12, height = 12)
message("Figure 4B saved to ", file.path(plots_dir, "Fig4B.pdf"))



# --- Statistical tests (optional output to console) ---------------------------

# # Filter out regimens with <= 2 cases for testing
# regimen_counts <- table(df_narm$exposed_regimen)
# regimens_to_keep <- names(regimen_counts[regimen_counts > 2])
# df_filtered <- df_narm[df_narm$exposed_regimen %in% regimens_to_keep, ]

# # T-tests for PC1 (Dim1)
# baseline_vs_other_dim1 <- lapply(
#   setdiff(unique(df_filtered$exposed_regimen), "baseline"),
#   function(reg) {
#     group1 <- df_filtered$Dim1[df_filtered$exposed_regimen == "baseline"]
#     group2 <- df_filtered$Dim1[df_filtered$exposed_regimen == reg]
#     if (length(group2) > 2) {
#       t_res <- t.test(group2, group1)
#       data.frame(
#         regimen = as.character(reg),
#         p.value = t_res$p.value,
#         mean.diff = mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)
#       )
#     } else {
#       data.frame(regimen = as.character(reg), p.value = NA, mean.diff = NA)
#     }
#   }
# )
# baseline_vs_other_dim1_df <- do.call(rbind, baseline_vs_other_dim1)

# # T-tests for PC2 (Dim2)
# baseline_vs_other_dim2 <- lapply(
#   setdiff(unique(df_filtered$exposed_regimen), "baseline"),
#   function(reg) {
#     group1 <- df_filtered$Dim2[df_filtered$exposed_regimen == "baseline"]
#     group2 <- df_filtered$Dim2[df_filtered$exposed_regimen == reg]
#     if (length(group2) > 2) {
#       t_res <- t.test(group2, group1)
#       data.frame(
#         regimen = as.character(reg),
#         p.value = t_res$p.value,
#         mean.diff = mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)
#       )
#     } else {
#       data.frame(regimen = as.character(reg), p.value = NA, mean.diff = NA)
#     }
#   }
# )
# baseline_vs_other_dim2_df <- do.call(rbind, baseline_vs_other_dim2)

# # Linear models adjusting for Disease_type
# if (requireNamespace("broom", quietly = TRUE)) {
#   lm_dim1 <- lm(Dim1 ~ exposed_regimen + Disease_type, data = df_filtered)
#   coef_tbl_dim1 <- broom::tidy(lm_dim1, conf.int = TRUE)
  
#   lm_dim2 <- lm(Dim2 ~ exposed_regimen + Disease_type, data = df_filtered)
#   coef_tbl_dim2 <- broom::tidy(lm_dim2, conf.int = TRUE)
# }
