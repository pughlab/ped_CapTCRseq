#!/usr/bin/env Rscript

# Title: Figure 1

# Inputs:
#   - `chemo_cycles_cleaned.rds`
#   - Color and theme helpers from `R/color_schemes.R` and `R/ggplot2_theme.R`
# Outputs:
#   - Plots saved to `plots/`
#   - Tables saved to `tables/`
# Usage:
#   Rscript R/Figure1.R

# --- Paths and setup ---------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(gtsummary)
  library(gt)
  library(swimplot)
})

# Load helper palettes and themes
source(file.path(getwd(), "R/functions", "color_schemes.R"))
source(file.path(getwd(), "R/functions", "ggplot2_theme.R"))

# Project-rooted paths
project_root <- getwd()
input_rds <- file.path(project_root, "data", "chemo_cycles_cleaned.rds")
plots_dir <- file.path(project_root, "plots")
tables_dir <- file.path(project_root, "tables")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)


# --- Load inputs -------------------------------------------------------------

# Read data
soc <- readr::read_rds(input_rds)

# Clean drug name inconsistencies similar to notebook
normalize_drug_names <- function(x) {
  x <- gsub("dactinomycin", "Dactinomycin", x, ignore.case = TRUE)
  x <- gsub("Cisplatinum", "Cisplatin", x, ignore.case = TRUE)
  x <- gsub("Cyatarabine", "Cytarabine", x, ignore.case = TRUE)
  x <- sub("\\*.*$", "", x)
  x <- sub("held during radiation", "", x, ignore.case = TRUE)
  x <- trimws(x)
  x
}

# Long-format for bubble plot across all cycles
drug_cols <- c(
  "Cycle 1 Drugs", "Cycle 2 Drugs", "Cycle 3 Drugs", "Cycle 4 Drugs", "Cycle 5 Drugs"
)

soc_long <- soc |>
  tidyr::pivot_longer(
    cols = matches("Cycle \\d+ Drugs"),
    names_to = "cycle",
    values_to = "drugs",
    names_pattern = "Cycle (\\d+) Drugs"
  ) %>%
  # Keep all other columns
  select(everything()) %>%
  # Clean up cycle numbers
  mutate(cycle = as.numeric(cycle)) |>
  filter(!is.na(drugs), drugs != "N/A") |>
  mutate(drugs = normalize_drug_names(drugs)) |>
  tidyr::separate_rows(drugs, sep = ", ") |>
  mutate(drugs = trimws(drugs)) 

soc_long$drug_type <- NA
soc_long$drug_type[grepl("Vincristine|Daunorubicin|Asparaginase|Mercaptopurine|Cyclophosphamide|Cytarabine|Methotrexate|Doxorubicin|Thioguanine|Cisplatin|Bleomyocin|Etoposide|5-Fluorouracil|Mitoxantrone|Ifosphamide|Topotecan|Dactinomycin", soc_long$drugs, ignore.case=TRUE)] <- "Chemotherapy"
soc_long$drug_type[grepl("Dexamethasone|Prednisone", soc_long$drugs, ignore.case=TRUE)] <- "Steroid"
soc_long$drug_type[grepl("Imatinib", soc_long$drugs, ignore.case=TRUE)] <- "Small molecule"

  tab <- soc_long %>%
    group_by(Disease, Cancergroup, cycle, drugs) %>% 
    summarise(Count = n(), .groups = "drop")  
  # Add drug_type column by matching drugs from data
  drug_type_mapping <- soc_long %>%
    select(drugs, drug_type) %>%
    distinct()  
  tab <- tab %>%
    left_join(drug_type_mapping, by = c("drugs" = "drugs"))
  tab <- tab %>%
    arrange(Cancergroup, Disease, Count) %>%
    mutate(Disease = factor(Disease, levels = unique(Disease)))
  # Create drug_type_summary for faceting order
  drug_type_summary <- tab %>%
    group_by(drug_type) %>%
    summarise(total_count = sum(Count)) %>%
    arrange(desc(total_count))
  # Factorize drugs based on total count within each drug type
  tab <- tab %>%
    group_by(drugs, drug_type) %>%
    mutate(total_count = sum(Count)) %>%
    arrange(total_count) %>%
    #mutate(drugs = factor(drugs, levels = unique(drugs))) %>%
    ungroup()
    tab$drugs <- factor(tab$drugs, levels = unique(tab$drugs))
  tab$cycle <- paste0("Cycle ", tab$cycle)
  tab$cycle <- factor(tab$cycle, levels = c("Cycle 1", "Cycle 2", "Cycle 3", "Cycle 4", "Cycle 5"))

# --- Fig1A ------------------------------------
# Bubble plot of drug usage by disease
p <- ggplot(tab, aes(x = Disease, y = drugs, size = Count, color = Cancergroup)) + 
    geom_point(alpha = 0.6, shape = 16) + 
    scale_size(range = c(1,6)) +
    scale_color_manual(values = group_col) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.title = element_blank(),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")) +
    coord_cartesian(clip = "off") +
    labs(title = paste("Drug Usage by Disease"))
  
  # Add facets
  p <- p + facet_grid(vars(drug_type), vars(cycle), scales = "free_y", space = "free", switch = "y") +
    theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside")

ggsave(file.path(plots_dir, "Fig1A.pdf"), p, width = 15, height = 6)
message("Figure 1A saved to ", file.path(plots_dir, "Fig1A.pdf"))

# --- Drug usage summary table ------------------------------------

# Table note: tab[order(tab$total_count, decreasing = TRUE), ] # convert to gtsummary
# Since total_count is not explicitly computed in the notebook state, compute by drug
tab$perc <- tab$total_count/sum(tab$Count)
writexl::write_xlsx(
  tab[order(tab$perc, decreasing = TRUE), ],
  path = file.path(tables_dir, "drug_usage_summary.xlsx")
)


# --- Fig1B ------------------------------------
# Derive base_regimen per cycle 
for (cycle in 1:5) {
  cycle_data <- soc_long[soc_long$cycle == cycle, ]
  vin_patients <- unique(cycle_data$`Patient ID`[cycle_data$drugs == "Vincristine"])
  mtx_patients <- unique(cycle_data$`Patient ID`[cycle_data$drugs == "Methotrexate"])
  cyclo_patients <- unique(cycle_data$`Patient ID`[cycle_data$drugs == "Cyclophosphamide"])

  col_name <- paste0("base_regimen_cycle", cycle)
  soc[[col_name]] <- sapply(soc$`Patient ID`, function(pat) {
    regimens <- c()
    if (pat %in% vin_patients) regimens <- c(regimens, "Vincristine")
    if (pat %in% mtx_patients) regimens <- c(regimens, "Methotrexate")
    if (pat %in% cyclo_patients) regimens <- c(regimens, "Cyclophosphamide")
    if (length(regimens) == 0) return(NA)
    paste(regimens, collapse = "/")
  })
  soc[[col_name]][is.na(soc[[col_name]])] <- "Others"
}


drugs_df <- soc |>
  transmute(
    Patient = paste0("CHP_", `Patient ID`),
    Cancergroup = `Cancergroup`,
    Disease = `Disease`,
    `Cycle 1 Duration (days)` = suppressWarnings(as.numeric(`Cycle 1 Duration (days)`)),
    `Cycle 2 Duration (days)` = suppressWarnings(as.numeric(`Cycle 2 Duration (days)`)),
    `Cycle 3 Duration (days)` = suppressWarnings(as.numeric(`Cycle 3 Duration (days)`)),
    `Cycle 4 Duration (days)` = suppressWarnings(as.numeric(`Cycle 4 Duration (days)`)),
    `Cycle 5 Duration (days)` = suppressWarnings(as.numeric(`Cycle 5 Duration (days)`)),
    base_regimen_cycle1 = `base_regimen_cycle1`,
    base_regimen_cycle2 = `base_regimen_cycle2`,
    base_regimen_cycle3 = `base_regimen_cycle3`,
    base_regimen_cycle4 = `base_regimen_cycle4`,
    base_regimen_cycle5 = `base_regimen_cycle5`
  ) |>
  mutate(
    trtend = rowSums(across(contains("Duration (days)")), na.rm = TRUE) + 50
  )


# --- Regimen summary table ------------------------------------

all_regimens <- c(
  as.character(drugs_df$base_regimen_cycle1),
  as.character(drugs_df$base_regimen_cycle2),
  as.character(drugs_df$base_regimen_cycle3),
  as.character(drugs_df$base_regimen_cycle4),
  as.character(drugs_df$base_regimen_cycle5)
)
# Create a table of counts
regimen_counts <- table(all_regimens, useNA = "ifany")
# Calculate total count (excluding NA)
total_count <- sum(regimen_counts[!is.na(names(regimen_counts))])
# Calculate percentages
percent <- (as.numeric(regimen_counts) / total_count) * 100
# Make dataframe
regimen_summary <- data.frame(
  Regimen = names(regimen_counts),
  Count = as.numeric(regimen_counts),
  Percentage = round(percent, 2)
)

# Add the total as a new column (same value for all rows)
regimen_summary$Total = total_count

# Save regimen_summary as an Excel file
writexl::write_xlsx(regimen_summary, path = file.path(tables_dir, "regimen_summary.xlsx"))

# Create a vector of the "Cycle X Duration (days)" column names
cycle_duration_cols <- grep("Cycle [1-5] Duration \\(days\\)", colnames(drugs_df), value = TRUE)
# Summarize each Cycle X Duration column
cycle_summaries <- lapply(cycle_duration_cols, function(col) {
  message(paste0("\nSummary for ", col, ":\n"))
  print(summary(drugs_df[[col]]))
})
# Summarize all Cycle X Duration columns together (stacked)
all_durations <- unlist(drugs_df[cycle_duration_cols])
message("\nSummary for all Cycle Duration columns combined:\n")
print(summary(all_durations))

# --- Swimmer plot ------------------------------------

# Order patients by cancer group and trtend then by cycle1 regimen
drugs_df <- drugs_df |>
  arrange(Cancergroup, trtend, base_regimen_cycle1) |>
  mutate(Patient = factor(Patient, levels = unique(Patient))) |>
  as.data.frame()

# Compute cycle start and end markers with 10-day gaps as in notebook
drugs_df$cycle2start <- drugs_df$`Cycle 1 Duration (days)`
drugs_df$cycle2start <- drugs_df$cycle2start + 10
drugs_df$cycle2end <- rowSums(drugs_df[, c("cycle2start", "Cycle 2 Duration (days)")], na.rm = TRUE)

drugs_df$cycle3start <- rowSums(drugs_df[, c("cycle2start", "Cycle 2 Duration (days)")], na.rm = TRUE)
drugs_df$cycle3start <- drugs_df$cycle3start + 10
drugs_df$cycle3end <- rowSums(drugs_df[, c("cycle3start", "Cycle 3 Duration (days)")], na.rm = TRUE)

drugs_df$cycle4start <- rowSums(drugs_df[, c("cycle3start", "Cycle 3 Duration (days)")], na.rm = TRUE)
drugs_df$cycle4start <- drugs_df$cycle4start + 10
drugs_df$cycle4end <- rowSums(drugs_df[, c("cycle4start", "Cycle 4 Duration (days)")], na.rm = TRUE)

drugs_df$cycle5start <- rowSums(drugs_df[, c("cycle4start", "Cycle 4 Duration (days)")], na.rm = TRUE)
drugs_df$cycle5start <- drugs_df$cycle5start + 10
drugs_df$cycle5end <- rowSums(drugs_df[, c("cycle5start", "Cycle 5 Duration (days)")], na.rm = TRUE)

# Convert NAs for absent cycles and add existence flags
drugs_df$cycle2start[is.na(drugs_df$`Cycle 2 Duration (days)`)] <- NA
drugs_df$cycle2end[is.na(drugs_df$`Cycle 2 Duration (days)`)] <- NA
drugs_df$cycle3start[is.na(drugs_df$`Cycle 3 Duration (days)`)] <- NA
drugs_df$cycle3end[is.na(drugs_df$`Cycle 3 Duration (days)`)] <- NA
drugs_df$cycle4start[is.na(drugs_df$`Cycle 4 Duration (days)`)] <- NA
drugs_df$cycle4end[is.na(drugs_df$`Cycle 4 Duration (days)`)] <- NA
drugs_df$cycle5start[is.na(drugs_df$`Cycle 5 Duration (days)`)] <- NA
drugs_df$cycle5end[is.na(drugs_df$`Cycle 5 Duration (days)`)] <- NA

drugs_df$cycle1exist <- ifelse(!is.na(drugs_df$`Cycle 1 Duration (days)`), "Yes", NA)
drugs_df$cycle2exist <- ifelse(!is.na(drugs_df$`Cycle 2 Duration (days)`), "Yes", NA)
drugs_df$cycle3exist <- ifelse(!is.na(drugs_df$`Cycle 3 Duration (days)`), "Yes", NA)
drugs_df$cycle4exist <- ifelse(!is.na(drugs_df$`Cycle 4 Duration (days)`), "Yes", NA)
drugs_df$cycle5exist <- ifelse(!is.na(drugs_df$`Cycle 5 Duration (days)`), "Yes", NA)

# Base swimmer
swim <- swimmer_plot(
  df = drugs_df, id = "Patient", end = "trtend", id_order = drugs_df$Patient,
  col = "#cdcaca", fill = "#cdcaca", alpha = 0.5, width = 0.6
)

# Add cycle segments and points similar to notebook
myswim <- swim +
  swimmer_points(df_points = drugs_df, id = "Patient", time = 0, name_col = "cycle1exist", size = 1) +
  swimmer_lines(df_lines = drugs_df, id = "Patient", start = "0", end = "`Cycle 1 Duration (days)`", name_col = "base_regimen_cycle1", size = 0.5) +
  swimmer_points(df_points = drugs_df, id = "Patient", time = "cycle2start", name_col = "cycle2exist", size = 1) +
  swimmer_lines(df_lines = drugs_df, id = "Patient", start = "cycle2start", end = "cycle2end", name_col = "base_regimen_cycle2", size = 0.5) +
  swimmer_points(df_points = drugs_df, id = "Patient", time = "cycle3start", name_col = "cycle3exist", size = 1) +
  swimmer_lines(df_lines = drugs_df, id = "Patient", start = "cycle3start", end = "cycle3end", name_col = "base_regimen_cycle3", size = 0.5) +
  swimmer_points(df_points = drugs_df, id = "Patient", time = "cycle4start", name_col = "cycle4exist", size = 1) +
  swimmer_lines(df_lines = drugs_df, id = "Patient", start = "cycle4start", end = "cycle4end", name_col = "base_regimen_cycle4", size = 0.5) +
  swimmer_points(df_points = drugs_df, id = "Patient", time = "cycle5start", name_col = "cycle5exist", size = 1) +
  swimmer_lines(df_lines = drugs_df, id = "Patient", start = "cycle5start", end = "cycle5end", name_col = "base_regimen_cycle5", size = 0.5)

myp <- myswim + 
scale_color_manual(name = "Regimen", 
values = c("Vincristine" = "#E41A1C", 
           "Methotrexate" = "#377EB8",
           "Cyclophosphamide" = "#4DAF4A",
           "Vincristine/Methotrexate" = "#984EA3",
           "Vincristine/Cyclophosphamide" = "#fa9734",
           "Vincristine/Methotrexate/Cyclophosphamide" = "#f781bf",
           "Others" = "black",
           "Yes" = "black",
           na.value=NA),
breaks = c("Vincristine", "Methotrexate", "Cyclophosphamide", 
           "Vincristine/Methotrexate", "Vincristine/Cyclophosphamide",
           "Vincristine/Methotrexate/Cyclophosphamide",
           "Others")) +
guides(color = guide_legend(override.aes = list(shape = NA))) +
labs(y = "Days after Cycle 1 therapy start", x = "") 


# Save swimmer plot and version labeled by disease on x
ggsave(file.path(plots_dir, "Fig1B.pdf"), 
myp+ scale_x_discrete(labels = drugs_df$Disease[match(drugs_df$Patient, drugs_df$Patient)]) , width = 12, height = 12)

message("Figure 1B saved to ", file.path(plots_dir, "Fig1B.pdf"))