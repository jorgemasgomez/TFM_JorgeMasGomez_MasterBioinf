############################################################
# Clean Template for Merging BLUEs/BLUPs, Heritability Outputs
# and Generating Publication-Quality Plots
#
# Master's Thesis - Bioinformatics
#
# Purpose:
#   This script merges variance components, BLUEs and BLUPs generated
#   from different trait datasets, formats BLUEs for downstream analyses,
#   and generates summary plots for variance components, BLUP correlations
#   and BLUP-based PCA.
#
# Main outputs:
#   - ALL_variance_components_h2.txt
#   - ALL_BLUEs_traits.txt
#   - ALL_BLUPs_traits.txt
#   - BLUEs_final_formatted.txt
#   - variance_components_stacked_selected_traits.png / .pdf
#   - BLUP_correlations_pairwise_selected_traits.txt
#   - BLUP_correlation_heatmap_selected_traits.png / .pdf
#   - BLUP_PCA_scores.txt
#   - BLUP_PCA_variance_explained.txt
#   - BLUP_PCA.png / .pdf
############################################################


############################################################
# 0. Packages
############################################################

load_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      paste0("Package '", package_name, "' is required but is not installed."),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages(
    library(package_name, character.only = TRUE)
  )
}

load_package("dplyr")
load_package("tidyr")
load_package("stringr")
load_package("ggplot2")
load_package("Hmisc")
load_package("reshape2")
load_package("scales")


############################################################
# 1. User configuration
############################################################

input_dir <- "."
output_dir <- "merged_BLUEs_BLUPs_outputs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Files used to filter the final BLUE table.
ids_keep_file <- "IDs_CHIP.txt"

# Output file names.
variance_output_file <- file.path(output_dir, "ALL_variance_components_h2.txt")
blues_output_file <- file.path(output_dir, "ALL_BLUEs_traits.txt")
blups_output_file <- file.path(output_dir, "ALL_BLUPs_traits.txt")
blues_wide_output_file <- file.path(output_dir, "BLUEs_final_formatted.txt")


############################################################
# 2. Trait labels
############################################################

# Labels used for the variance component plot.
trait_labels_variance <- c(
  "Kernel_Area" = "Kernel Area",
  "Kernel_Circularity" = "Kernel Circularity",
  "Kernel_EF_PC1" = "Kernel EF PC1",
  "Kernel_EF_PC3" = "Kernel EF PC3",
  "Kernel_EF_PC4" = "Kernel EF PC4",
  "Kernel_EF_PC5" = "Kernel EF PC5",
  "Kernel_EF_PC6" = "Kernel EF PC6",
  "Kernel_EF_PC7" = "Kernel EF PC7",
  "Kernel_EF_PC8" = "Kernel EF PC8",
  "Kernel_EF_PC9" = "Kernel EF PC9",
  "Kernel_EF_PC10" = "Kernel EF PC10",
  "Kernel_Ellipse_Ratio" = "Kernel Ellipse Ratio",
  "Kernel_Length" = "Kernel Length",
  "Kernel_Perimeter" = "Kernel Perimeter",
  "Kernel_Ratio_Length_Width" = "Kernel Ratio Length-Width",
  "Kernel_Ratio_Length_Width25" = "Kernel Ratio Length-Width25",
  "Kernel_Ratio_Length_Width75" = "Kernel Ratio Length-Width75",
  "Kernel_Ratio_Width25_Width_75" = "Kernel Ratio Width25-Width 75",
  "Kernel_Ratio_Width_Width_25" = "Kernel Ratio Width-Width 25",
  "Kernel_Ratio_Width_Width_75" = "Kernel Ratio Width-Width 75",
  "Kernel_Shoulder_symmetry" = "Kernel Shoulder symmetry",
  "Kernel_Symmetry_h" = "Kernel Symmetry h",
  "Kernel_Symmetry_v" = "Kernel Symmetry v",
  "Kernel_Mean_weight_per_grain" = "Kernel Weight",
  "Kernel_Width" = "Kernel Width",
  "Kernel_Width_25" = "Kernel Width 25",
  "Kernel_Width_75" = "Kernel Width 75",

  "Ratio_Area" = "KS Ratio Area",
  "Ratio_Circularity" = "KS Ratio Circularity",
  "Ratio_Ellipse_Ratio" = "KS Ratio Ellipse Ratio",
  "Ratio_Length" = "KS Ratio Length",
  "Ratio_Ratio_Length_Width" = "KS Ratio Ratio Length-Width",
  "Ratio_Ratio_Length_Width25" = "KS Ratio Ratio Length-Width25",
  "Ratio_Ratio_Length_Width75" = "KS Ratio Ratio Length-Width75",
  "Ratio_Ratio_Width25_Width_75" = "KS Ratio Ratio Width25-Width 75",
  "Ratio_Ratio_Width_Width_25" = "KS Ratio Ratio Width-Width 25",
  "Ratio_Ratio_Width_Width_75" = "KS Ratio Ratio Width-Width 75",
  "Ratio_Symmetry_h" = "KS Ratio Symmetry h",
  "Ratio_Symmetry_v" = "KS Ratio Symmetry v",
  "Ratio_Width" = "KS Ratio Width",
  "Ratio_Width_25" = "KS Ratio Width 25",
  "Ratio_Width_75" = "KS Ratio Width 75",
  "Ratio_Mean_weight_per_grain" = "Ratio KS Weight",

  "Pred_Fats" = "Pred Fats",
  "Pred_Fiber" = "Pred Fiber",
  "Pred_Protein" = "Pred Protein",
  "Pred_Sucrose" = "Pred Sucrose",

  "Shell_Area" = "Shell Area",
  "Shell_Circularity" = "Shell Circularity",
  "Shell_EF_PC1" = "Shell EF PC1",
  "Shell_EF_PC2" = "Shell EF PC2",
  "Shell_EF_PC3" = "Shell EF PC3",
  "Shell_EF_PC4" = "Shell EF PC4",
  "Shell_EF_PC5" = "Shell EF PC5",
  "Shell_EF_PC6" = "Shell EF PC6",
  "Shell_EF_PC7" = "Shell EF PC7",
  "Shell_EF_PC8" = "Shell EF PC8",
  "Shell_EF_PC9" = "Shell EF PC9",
  "Shell_EF_PC10" = "Shell EF PC10",
  "Shell_Ellipse_Ratio" = "Shell Ellipse Ratio",
  "Shell_Length" = "Shell Length",
  "Shell_Perimeter" = "Shell Perimeter",
  "Shell_Ratio_Length_Width" = "Shell Ratio Length-Width",
  "Shell_Ratio_Length_Width25" = "Shell Ratio Length-Width25",
  "Shell_Ratio_Length_Width75" = "Shell Ratio Length-Width75",
  "Shell_Ratio_Width25_Width_75" = "Shell Ratio Width25-Width 75",
  "Shell_Ratio_Width_Width_25" = "Shell Ratio Width-Width 25",
  "Shell_Ratio_Width_Width_75" = "Shell Ratio Width-Width 75",
  "Shell_Shoulder_symmetry" = "Shell Shoulder symmetry",
  "Shell_Symmetry_h" = "Shell Symmetry h",
  "Shell_Symmetry_v" = "Shell Symmetry v",
  "Shell_Mean_weight_per_grain" = "Shell Weight",
  "Shell_Width" = "Shell Width",
  "Shell_Width_25" = "Shell Width 25",
  "Shell_Width_75" = "Shell Width 75",

  "Thickness_pred" = "Thickness pred"
)

# Reduced trait list used for the BLUP correlation heatmap.
trait_labels_correlation <- c(
  "Kernel_Mean_weight_per_grain" = "Kernel Weight",
  "Kernel_Width" = "Kernel Width",
  "Shell_Width_25" = "Shell Width 25",
  "Shell_Symmetry_v" = "Shell Symmetry v",
  "Kernel_Width_25" = "Kernel Width 25",
  "Kernel_Ratio_Length_Width25" = "Kernel Ratio Length-Width25",
  "Shell_Area" = "Shell Area",
  "Shell_Shoulder_symmetry" = "Shell Shoulder symmetry",
  "Shell_Width" = "Shell Width",
  "Shell_EF_PC4" = "Shell EF PC4",
  "Shell_Ratio_Width25_Width_75" = "Shell Ratio Width25-Width 75",
  "Kernel_Ratio_Length_Width75" = "Kernel Ratio Length-Width75",
  "Ratio_Circularity" = "KS Ratio Circularity",
  "Ratio_Width_25" = "KS Ratio Width 25",
  "Shell_EF_PC3" = "Shell EF PC3",
  "Shell_Ratio_Width_Width_25" = "Shell Ratio Width-Width 25",
  "Shell_Length" = "Shell Length",
  "Ratio_Symmetry_v" = "KS Ratio Symmetry v",
  "Shell_Symmetry_h" = "Shell Symmetry h",
  "Shell_EF_PC1" = "Shell EF PC1",
  "Shell_Ellipse_Ratio" = "Shell Ellipse Ratio",
  "Kernel_Shoulder_symmetry" = "Kernel Shoulder symmetry",
  "Shell_Ratio_Length_Width" = "Shell Ratio Length-Width",
  "Ratio_Length" = "KS Ratio Length",
  "Shell_EF_PC9" = "Shell EF PC9",
  "Kernel_EF_PC1" = "Kernel EF PC1",
  "Shell_Perimeter" = "Shell Perimeter",
  "Shell_Ratio_Width_Width_75" = "Shell Ratio Width-Width 75",
  "Kernel_Perimeter" = "Kernel Perimeter",
  "Shell_Width_75" = "Shell Width 75",
  "Shell_Ratio_Length_Width25" = "Shell Ratio Length-Width25",
  "Shell_Ratio_Length_Width75" = "Shell Ratio Length-Width75",
  "Kernel_Area" = "Kernel Area",
  "Shell_Circularity" = "Shell Circularity",
  "Ratio_Area" = "KS Ratio Area",
  "Pred_Fiber" = "Pred Fiber",
  "Kernel_Ratio_Length_Width" = "Kernel Ratio Length-Width",
  "Ratio_Ratio_Length_Width" = "KS Ratio Ratio Length-Width",
  "Ratio_Width" = "KS Ratio Width",
  "Kernel_Width_75" = "Kernel Width 75",
  "Kernel_Circularity" = "Kernel Circularity",
  "Ratio_Ellipse_Ratio" = "KS Ratio Ellipse Ratio",
  "Ratio_Ratio_Width_Width_75" = "KS Ratio Ratio Width-Width 75",
  "Ratio_Ratio_Length_Width75" = "KS Ratio Ratio Length-Width75",
  "Ratio_Width_75" = "KS Ratio Width 75",
  "Ratio_Ratio_Width25_Width_75" = "KS Ratio Ratio Width25-Width 75",
  "Pred_Fats" = "Pred Fats",
  "Kernel_Length" = "Kernel Length",
  "Shell_EF_PC5" = "Shell EF PC5",
  "Pred_Protein" = "Pred Protein",
  "Pred_Sucrose" = "Pred Sucrose",
  "Kernel_Ellipse_Ratio" = "Kernel Ellipse Ratio",
  "Kernel_Symmetry_v" = "Kernel Symmetry v",
  "Thickness_pred" = "Thickness pred",
  "Kernel_EF_PC5" = "Kernel EF PC5",
  "Ratio_Mean_weight_per_grain" = "Ratio KS Weight",
  "Kernel_Ratio_Width25_Width_75" = "Kernel Ratio Width25-Width 75",
  "Kernel_EF_PC4" = "Kernel EF PC4",
  "Ratio_Ratio_Length_Width25" = "KS Ratio Ratio Length-Width25",
  "Kernel_Ratio_Width_Width_75" = "Kernel Ratio Width-Width 75",
  "Ratio_Ratio_Width_Width_25" = "KS Ratio Ratio Width-Width 25",
  "Ratio_Symmetry_h" = "KS Ratio Symmetry h",
  "Shell_Mean_weight_per_grain" = "Shell Weight",
  "Kernel_Symmetry_h" = "Kernel Symmetry h",
  "Kernel_EF_PC3" = "Kernel EF PC3",
  "Kernel_Ratio_Width_Width_25" = "Kernel Ratio Width-Width 25"
)


############################################################
# 3. Helper functions
############################################################

detect_trait_prefix <- function(file_path) {
  file_lower <- tolower(basename(file_path))

  if (str_detect(file_lower, "kernel")) {
    return("Kernel_")
  }

  if (str_detect(file_lower, "shell")) {
    return("Shell_")
  }

  return("")
}

merge_trait_result_files <- function(
    input_dir,
    file_pattern,
    output_file,
    add_prefix = TRUE
) {
  files <- list.files(
    path = input_dir,
    pattern = file_pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop(
      paste0("No files found matching pattern: ", file_pattern),
      call. = FALSE
    )
  }

  cat("\nFiles detected for pattern:", file_pattern, "\n")
  print(files)

  df_list <- lapply(files, function(file_path) {
    df <- read.table(
      file_path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (!"Trait" %in% colnames(df)) {
      stop(
        paste0("Column 'Trait' was not found in file: ", file_path),
        call. = FALSE
      )
    }

    if (add_prefix) {
      trait_prefix <- detect_trait_prefix(file_path)
      df$Trait <- paste0(trait_prefix, df$Trait)
    }

    df$Source_file <- basename(file_path)

    return(df)
  })

  merged_df <- bind_rows(df_list)

  write.table(
    merged_df,
    output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  cat("\nMerged table saved to:", output_file, "\n")
  cat("Dimensions:", paste(dim(merged_df), collapse = " x "), "\n")

  return(merged_df)
}

filter_and_rename_traits <- function(df, trait_labels) {
  df$Trait <- trimws(as.character(df$Trait))

  missing_traits <- setdiff(names(trait_labels), unique(df$Trait))

  cat("\nTraits requested but not found:\n")
  print(missing_traits)

  df_filtered <- df %>%
    filter(Trait %in% names(trait_labels)) %>%
    mutate(
      Trait_oldname = Trait,
      Trait = unname(trait_labels[Trait])
    )

  cat("\nNumber of traits kept:\n")
  print(length(unique(df_filtered$Trait)))

  cat("\nTraits kept:\n")
  print(unique(df_filtered$Trait))

  return(df_filtered)
}

read_ids_to_keep <- function(ids_keep_file) {
  ids_keep <- read.table(
    ids_keep_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ids_keep <- unique(ids_keep[[1]])

  cat("\nNumber of IDs to keep:", length(ids_keep), "\n")

  return(ids_keep)
}


############################################################
# 4. Merge variance components, BLUEs and BLUPs
############################################################

merge_all_outputs <- function() {
  variance_components <- merge_trait_result_files(
    input_dir = input_dir,
    file_pattern = "Variance_components_h2.txt$",
    output_file = variance_output_file,
    add_prefix = TRUE
  )

  blues_all <- merge_trait_result_files(
    input_dir = input_dir,
    file_pattern = "_BLUEs_traits.txt$",
    output_file = blues_output_file,
    add_prefix = TRUE
  )

  blups_all <- merge_trait_result_files(
    input_dir = input_dir,
    file_pattern = "_BLUPs_traits.txt$",
    output_file = blups_output_file,
    add_prefix = TRUE
  )

  return(
    list(
      variance_components = variance_components,
      blues = blues_all,
      blups = blups_all
    )
  )
}


############################################################
# 5. Format BLUEs for downstream analyses
############################################################

format_blues_for_downstream_analysis <- function(
    blues_all,
    ids_keep_file,
    output_file
) {
  if (!file.exists(ids_keep_file)) {
    warning(
      paste0(
        "ID filter file was not found: ", ids_keep_file,
        ". The BLUE table will not be filtered by ID."
      ),
      call. = FALSE
    )

    blues_filtered <- blues_all
  } else {
    ids_keep <- read_ids_to_keep(ids_keep_file)

    blues_filtered <- blues_all %>%
      filter(ID %in% ids_keep)
  }

  blues_wide <- blues_filtered %>%
    select(ID, Trait, predicted.value) %>%
    pivot_wider(
      names_from = Trait,
      values_from = predicted.value
    )

  colnames(blues_wide)[colnames(blues_wide) == "ID"] <- "<Phenotype>"

  write.table(
    blues_wide,
    output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  cat("\nFormatted BLUE table saved to:", output_file, "\n")
  cat("Dimensions:", paste(dim(blues_wide), collapse = " x "), "\n")

  return(blues_wide)
}


############################################################
# 6. Variance component plot
############################################################

plot_selected_variance_components <- function(
    variance_components,
    trait_labels,
    output_prefix = "variance_components_stacked_selected_traits",
    include_residual = FALSE,
    width = 8,
    height = 16
) {
  df_filtered <- filter_and_rename_traits(
    df = variance_components,
    trait_labels = trait_labels
  )

  trait_order <- df_filtered %>%
    group_by(Trait) %>%
    summarise(
      h2 = mean(h2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(h2) %>%
    pull(Trait)

  df_filtered$Trait <- factor(df_filtered$Trait, levels = trait_order)

  if (include_residual) {
    df_long <- df_filtered %>%
      mutate(
        Total = Genotype + Environment + Residual,
        Genotype = Genotype / Total,
        Environment = Environment / Total,
        Residual = Residual / Total
      ) %>%
      select(Trait, Genotype, Environment, Residual) %>%
      pivot_longer(
        cols = c(Genotype, Environment, Residual),
        names_to = "Component",
        values_to = "Proportion"
      ) %>%
      mutate(
        Component = factor(
          Component,
          levels = c("Residual", "Environment", "Genotype")
        )
      )

    fill_values <- c(
      "Residual" = "#B0B0B0",
      "Environment" = "#4C956C",
      "Genotype" = "#4E79A7"
    )
  } else {
    df_long <- df_filtered %>%
      mutate(
        Total = Genotype + Environment + Residual,
        Genotype = Genotype / Total,
        Environment = Environment / Total
      ) %>%
      select(Trait, Genotype, Environment) %>%
      pivot_longer(
        cols = c(Genotype, Environment),
        names_to = "Component",
        values_to = "Proportion"
      ) %>%
      mutate(
        Component = factor(
          Component,
          levels = c("Environment", "Genotype")
        )
      )

    fill_values <- c(
      "Environment" = "#4C956C",
      "Genotype" = "#4E79A7"
    )
  }

  p <- ggplot(
    df_long,
    aes(
      x = Trait,
      y = Proportion,
      fill = Component
    )
  ) +
    geom_bar(
      stat = "identity",
      width = 1,
      color = "grey85",
      linewidth = 0.15
    ) +
    geom_text(
      aes(
        label = ifelse(
          Proportion >= 0.03,
          scales::percent(Proportion, accuracy = 1),
          ""
        )
      ),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "white",
      fontface = "bold"
    ) +
    coord_flip() +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(values = fill_values) +
    labs(
      x = "Trait",
      y = "Variance explained (%)",
      fill = "Variance component"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(
        size = 14,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 14,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 11,
        margin = margin(t = 10)
      ),
      axis.title.y = element_text(
        size = 11,
        margin = margin(r = 10)
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.25
      ),
      legend.position = "top",
      legend.title = element_text(
        size = 14,
        face = "bold"
      ),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )

  print(p)

  png_file <- file.path(output_dir, paste0(output_prefix, ".png"))
  pdf_file <- file.path(output_dir, paste0(output_prefix, ".pdf"))

  ggsave(
    filename = png_file,
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    units = "in",
    bg = "white"
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  cat("\nVariance component plot saved to:\n")
  cat(png_file, "\n")
  cat(pdf_file, "\n")

  return(p)
}


############################################################
# 7. BLUP correlation heatmap
############################################################

plot_blup_correlation_heatmap <- function(
    blups_all,
    trait_labels,
    output_prefix = "BLUP_correlation_heatmap_selected_traits",
    width = 18,
    height = 18
) {
  blups_all$ID <- as.character(blups_all$ID)
  blups_all$Trait <- trimws(as.character(blups_all$Trait))

  df_filtered <- filter_and_rename_traits(
    df = blups_all,
    trait_labels = trait_labels
  )

  df_unique <- df_filtered %>%
    group_by(ID, Trait) %>%
    summarise(
      predicted.value = mean(predicted.value, na.rm = TRUE),
      .groups = "drop"
    )

  wide <- df_unique %>%
    pivot_wider(
      names_from = Trait,
      values_from = predicted.value
    )

  mat <- wide %>%
    select(-ID)

  mat <- data.frame(
    lapply(mat, function(x) as.numeric(as.character(x)))
  )

  mat <- mat[, colSums(!is.na(mat)) > 5, drop = FALSE]
  mat <- mat[, apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

  cat("\nTraits used in the correlation matrix:\n")
  print(colnames(mat))

  res <- Hmisc::rcorr(as.matrix(mat), type = "pearson")

  cor_mat <- res$r
  p_mat <- res$P

  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "average")

  trait_order <- colnames(cor_mat)[hc$order]

  cor_mat <- cor_mat[trait_order, trait_order]
  p_mat <- p_mat[trait_order, trait_order]

  cor_melt <- reshape2::melt(cor_mat, na.rm = TRUE)
  p_melt <- reshape2::melt(p_mat, na.rm = TRUE)

  cor_df <- merge(
    cor_melt,
    p_melt,
    by = c("Var1", "Var2")
  )

  colnames(cor_df) <- c("Trait1", "Trait2", "cor", "p_value")

  cor_df$Trait1 <- factor(cor_df$Trait1, levels = trait_order)
  cor_df$Trait2 <- factor(cor_df$Trait2, levels = rev(trait_order))

  cor_pairs <- cor_df %>%
    filter(as.character(Trait1) != as.character(Trait2)) %>%
    mutate(
      pair = ifelse(
        as.character(Trait1) < as.character(Trait2),
        paste(Trait1, Trait2),
        paste(Trait2, Trait1)
      )
    ) %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(Trait1, Trait2, cor, p_value) %>%
    arrange(desc(abs(cor)))

  pairwise_file <- file.path(
    output_dir,
    "BLUP_correlations_pairwise_selected_traits.txt"
  )

  write.table(
    cor_pairs,
    pairwise_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p <- ggplot(
    cor_df,
    aes(
      Trait1,
      Trait2,
      fill = cor
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.1
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson r"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 14,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 14,
        color = "black"
      ),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(
        size = 18,
        face = "bold"
      ),
      legend.text = element_text(size = 14),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    coord_fixed()

  print(p)

  png_file <- file.path(output_dir, paste0(output_prefix, ".png"))
  pdf_file <- file.path(output_dir, paste0(output_prefix, ".pdf"))

  ggsave(
    filename = png_file,
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    width = width,
    height = height,
    bg = "white"
  )

  cat("\nCorrelation outputs saved to:\n")
  cat(pairwise_file, "\n")
  cat(png_file, "\n")
  cat(pdf_file, "\n")

  return(
    list(
      plot = p,
      pairwise_correlations = cor_pairs,
      correlation_matrix = cor_mat,
      p_value_matrix = p_mat
    )
  )
}


############################################################
# 8. BLUP PCA
############################################################

run_blup_pca <- function(
    blups_all,
    output_prefix = "BLUP_PCA",
    width = 8,
    height = 6
) {
  blups_all$ID <- as.character(blups_all$ID)
  blups_all$Trait <- trimws(as.character(blups_all$Trait))

  wide <- blups_all %>%
    group_by(ID, Trait) %>%
    summarise(
      predicted.value = mean(predicted.value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Trait,
      values_from = predicted.value
    )

  ids <- wide$ID

  mat <- wide %>%
    select(-ID)

  mat <- data.frame(
    lapply(mat, function(x) as.numeric(as.character(x)))
  )

  # Remove columns with all missing values.
  keep_cols <- colSums(is.na(mat)) < nrow(mat)
  mat <- mat[, keep_cols, drop = FALSE]

  # Remove columns with zero variance.
  keep_cols <- apply(mat, 2, sd, na.rm = TRUE) != 0
  mat <- mat[, keep_cols, drop = FALSE]

  # Remove rows with missing or infinite values, while keeping IDs aligned.
  mat[is.infinite(as.matrix(mat))] <- NA
  complete_rows <- complete.cases(mat)

  mat_clean <- mat[complete_rows, , drop = FALSE]
  ids_clean <- ids[complete_rows]

  cat("\nNumber of individuals used in PCA:", nrow(mat_clean), "\n")
  cat("Number of traits used in PCA:", ncol(mat_clean), "\n")

  pca <- prcomp(
    mat_clean,
    center = TRUE,
    scale. = TRUE
  )

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

  pca_df <- data.frame(pca$x)
  pca_df$ID <- ids_clean

  scores_file <- file.path(output_dir, paste0(output_prefix, "_scores.txt"))
  variance_file <- file.path(output_dir, paste0(output_prefix, "_variance_explained.txt"))

  write.table(
    pca_df,
    scores_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(
      PC = paste0("PC", seq_along(var_exp)),
      Variance_explained = var_exp,
      Variance_explained_percent = var_exp * 100
    ),
    variance_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p <- ggplot(
    pca_df,
    aes(
      PC1,
      PC2
    )
  ) +
    geom_point(
      size = 3,
      alpha = 0.85
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    xlab(paste0("PC1 (", round(var_exp[1] * 100, 1), "%)")) +
    ylab(paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"))

  print(p)

  png_file <- file.path(output_dir, paste0(output_prefix, ".png"))
  pdf_file <- file.path(output_dir, paste0(output_prefix, ".pdf"))

  ggsave(
    filename = png_file,
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    width = width,
    height = height,
    bg = "white"
  )

  cat("\nPCA outputs saved to:\n")
  cat(scores_file, "\n")
  cat(variance_file, "\n")
  cat(png_file, "\n")
  cat(pdf_file, "\n")

  return(
    list(
      pca = pca,
      scores = pca_df,
      variance_explained = var_exp,
      plot = p
    )
  )
}


############################################################
# 9. Main execution
############################################################

# Merge output files.
merged_outputs <- merge_all_outputs()

variance_components <- merged_outputs$variance_components
blues_all <- merged_outputs$blues
blups_all <- merged_outputs$blups

# Create a filtered and wide BLUE table for downstream analyses.
blues_wide <- format_blues_for_downstream_analysis(
  blues_all = blues_all,
  ids_keep_file = ids_keep_file,
  output_file = blues_wide_output_file
)

# Plot selected variance components.
p_variance <- plot_selected_variance_components(
  variance_components = variance_components,
  trait_labels = trait_labels_variance,
  output_prefix = "variance_components_stacked_selected_traits",
  include_residual = FALSE,
  width = 8,
  height = 16
)

# BLUP correlation heatmap.
correlation_results <- plot_blup_correlation_heatmap(
  blups_all = blups_all,
  trait_labels = trait_labels_correlation,
  output_prefix = "BLUP_correlation_heatmap_selected_traits",
  width = 18,
  height = 18
)

# BLUP PCA.
pca_results <- run_blup_pca(
  blups_all = blups_all,
  output_prefix = "BLUP_PCA",
  width = 8,
  height = 6
)

############################################################
# End of script
############################################################
