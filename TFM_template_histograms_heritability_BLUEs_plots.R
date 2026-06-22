############################################################
# Template for Heritability, BLUEs, BLUPs and Trait Plots
# Master's Thesis - Bioinformatics
#
# Project:
# Genome-wide association study for morphological and nutritional
# traits in almond using high-throughput genotyping and phenotyping
#
# Author: Dr. Jorge Mas Gómez
#
# Description:
# This script provides a reusable template to:
#   1. Load and prepare phenotypic data.
#   2. Generate descriptive histograms by environment or year.
#   3. Detect outliers using the MAD modified Z-score method.
#   4. Estimate BLUEs and BLUPs using ASReml.
#   5. Estimate variance components and broad-sense heritability.
#   6. Generate publication-quality plots for traits, BLUEs and variance components.
#
# Notes:
#   - Edit the "User configuration" section before running the script.
#   - ASReml is required only for BLUEs, BLUPs, variance components and heritability.
#   - The script assumes tab-delimited input files.
############################################################


############################################################
# 0. Packages
############################################################

load_package <- function(package_name, required = TRUE) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    message <- paste0(
      "Package '", package_name, "' is not installed."
    )

    if (required) {
      stop(message, call. = FALSE)
    } else {
      warning(message, call. = FALSE)
      return(FALSE)
    }
  }

  suppressPackageStartupMessages(
    library(package_name, character.only = TRUE)
  )

  return(TRUE)
}

# Core packages.
load_package("dplyr")
load_package("tidyr")
load_package("ggplot2")
load_package("grid")

# Optional packages used by specific sections.
has_patchwork <- load_package("patchwork", required = FALSE)
has_ggh4x <- load_package("ggh4x", required = FALSE)
has_asreml <- load_package("asreml", required = FALSE)


############################################################
# 1. User configuration
############################################################

# Input files.
phenotype_file <- "input_calidad.txt"
families_file <- "families.txt"

# Output directory.
output_dir <- "TFM_heritability_BLUEs_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Column names in the phenotype table.
id_col <- "ID_CHIP"
location_col <- "Location"
year_col <- "Year"
environment_col <- "Environment"

# Column names in the family/group table.
family_id_col <- "Sample ID"
family_group_col <- "Group"

# Traits to analyse.
traits <- c(
  "Pred_Fiber",
  "Pred_Fats",
  "Pred_Sucrose",
  "Pred_Protein"
)

# Trait groups for plots.
trait_groups <- list(
  Quality_traits = c(
    "Pred_Fats",
    "Pred_Fiber",
    "Pred_Protein",
    "Pred_Sucrose"
  )
)

# Order used for family/group plots.
family_order <- c(
  "Antoñeta_x_Penta",
  "Antoñeta_x_Tardona",
  "Antoñeta_x_Marcona",
  "Florida_x_Marcona",
  "Marcona_x_S4017",
  "R1000_x_Desmayo",
  "Selections",
  "Germplasm"
)

# ASReml memory settings.
asreml_workspace <- 128e6
asreml_pworkspace <- 128e6

# Set to FALSE if you only want descriptive statistics and plots.
run_asreml_analysis <- TRUE


############################################################
# 2. Helper functions
############################################################

clean_label <- function(x) {
  x <- gsub("_x_", " × ", x)
  x <- gsub("_", " ", x)
  return(x)
}

make_named_palette <- function(levels_vector, base_palette) {
  n <- length(levels_vector)

  if (n <= length(base_palette)) {
    cols <- base_palette[seq_len(n)]
  } else {
    cols <- grDevices::colorRampPalette(base_palette)(n)
  }

  setNames(cols, levels_vector)
}

get_varcomp_component <- function(varcomp_table, component_name) {
  if (component_name %in% rownames(varcomp_table)) {
    return(varcomp_table[component_name, "component"])
  }

  warning(
    paste0("Variance component '", component_name, "' was not found."),
    call. = FALSE
  )

  return(NA_real_)
}


############################################################
# 3. Data loading and preparation
############################################################

read_phenotype_data <- function(
    phenotype_file,
    id_col,
    location_col = "Location",
    year_col = "Year",
    environment_col = "Environment"
) {
  df <- read.csv(
    phenotype_file,
    sep = "\t",
    na.strings = c("", "NA"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(id_col, location_col, year_col, environment_col)
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required columns are missing from the phenotype file: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  df$ID <- as.factor(df[[id_col]])
  df$Location <- as.factor(df[[location_col]])
  df$Year <- as.factor(df[[year_col]])
  df$Environment <- as.factor(df[[environment_col]])

  return(df)
}

read_family_data <- function(
    families_file,
    family_id_col = "Sample ID",
    family_group_col = "Group",
    family_order = NULL
) {
  families <- read.table(
    families_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(family_id_col, family_group_col)
  missing_cols <- setdiff(required_cols, colnames(families))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required columns are missing from the family file: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  families <- families %>%
    rename(
      ID = all_of(family_id_col),
      Family = all_of(family_group_col)
    ) %>%
    mutate(
      ID = as.character(ID),
      Family = as.character(Family),
      Family_label = clean_label(Family)
    )

  if (!is.null(family_order)) {
    families$Family_label <- factor(
      families$Family_label,
      levels = clean_label(family_order)
    )
  }

  return(families)
}


############################################################
# 4. Descriptive histograms
############################################################

plot_trait_histograms <- function(
    df,
    traits,
    group_col,
    output_name,
    plot_title = NULL,
    bins = 30,
    width = 12,
    height = 8
) {
  traits_available <- traits[traits %in% colnames(df)]
  traits_missing <- setdiff(traits, colnames(df))

  cat("\nPlot:", output_name, "\n")
  cat("Traits included:\n")
  print(traits_available)

  if (length(traits_missing) > 0) {
    cat("Traits not found and skipped:\n")
    print(traits_missing)
  }

  if (length(traits_available) == 0) {
    stop("None of the requested traits were found in the phenotype table.")
  }

  if (!group_col %in% colnames(df)) {
    stop(paste0("Grouping column '", group_col, "' was not found."))
  }

  df_long <- df %>%
    select(all_of(c(group_col, traits_available))) %>%
    pivot_longer(
      cols = all_of(traits_available),
      names_to = "Trait",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value)) %>%
    mutate(
      Trait_label = factor(
        clean_label(Trait),
        levels = clean_label(traits_available)
      ),
      Group_label = clean_label(as.character(.data[[group_col]]))
    )

  if (is.factor(df[[group_col]])) {
    df_long$Group_label <- factor(
      df_long$Group_label,
      levels = clean_label(levels(df[[group_col]]))
    )
  }

  mean_df <- df_long %>%
    group_by(Trait_label, Group_label) %>%
    summarise(
      Mean_value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(df_long, aes(x = Value)) +
    geom_histogram(
      bins = bins,
      fill = "#4E79A7",
      color = "grey92",
      linewidth = 0.15
    ) +
    geom_vline(
      data = mean_df,
      aes(xintercept = Mean_value),
      color = "#08306B",
      linewidth = 0.5,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = mean_df,
      aes(
        x = Mean_value,
        y = Inf,
        label = sprintf("%.2f", Mean_value)
      ),
      color = "#08306B",
      angle = 90,
      vjust = 1.3,
      hjust = 1.1,
      size = 2.5,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    facet_grid(
      Group_label ~ Trait_label,
      scales = "free"
    ) +
    labs(
      x = "Phenotypic value",
      y = "Frequency",
      title = plot_title
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.background = element_rect(
        fill = "#4E79A7",
        color = "#4E79A7",
        linewidth = 0.25
      ),
      strip.text.x = element_text(
        size = 8,
        face = "bold",
        angle = 90,
        color = "white",
        margin = ggplot2::margin(t = 4, r = 2, b = 4, l = 2)
      ),
      strip.text.y = element_text(
        size = 8,
        face = "bold",
        color = "white",
        margin = ggplot2::margin(t = 2, r = 4, b = 2, l = 4)
      ),
      axis.text.x = element_text(
        size = 6,
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 6,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 10,
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 10,
        margin = ggplot2::margin(r = 8)
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.spacing.x = grid::unit(0.35, "lines"),
      panel.spacing.y = grid::unit(0.35, "lines"),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      ),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  print(p)

  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".png")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  return(p)
}


############################################################
# 5. BLUEs, BLUPs, variance components and heritability
############################################################

run_asreml_trait_models <- function(
    df,
    traits,
    id_col = "ID_CHIP",
    environment_col = "Environment",
    output_dir = ".",
    workspace = 128e6,
    pworkspace = 128e6
) {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop(
      "ASReml is required for BLUE, BLUP, variance component and heritability estimation.",
      call. = FALSE
    )
  }

  df$ID <- as.factor(df[[id_col]])
  df$Environment <- as.factor(df[[environment_col]])

  blue_list <- list()
  blup_list <- list()
  varcomp_list <- data.frame()
  cv_list <- data.frame()

  for (trait in traits) {
    if (!trait %in% colnames(df)) {
      warning(
        paste0("Trait '", trait, "' was not found and will be skipped."),
        call. = FALSE
      )
      next
    }

    df_sub <- df[!is.na(df[[trait]]), ]

    if (nrow(df_sub) == 0) {
      warning(
        paste0("Trait '", trait, "' has no non-missing values and will be skipped."),
        call. = FALSE
      )
      next
    }

    cat("\nFitting models for trait:", trait, "\n")

    ########################################################
    # Model 1: BLUEs
    ########################################################

    model_blue <- asreml::asreml(
      fixed = as.formula(paste(trait, "~ Environment + ID")),
      data = df_sub,
      workspace = workspace,
      pworkspace = pworkspace
    )

    blue <- predict(model_blue, classify = "ID")$pvals
    blue$Trait <- trait
    blue_list[[trait]] <- blue

    ########################################################
    # Residual coefficient of variation
    ########################################################

    vc_blue <- summary(model_blue)$varcomp

    var_R <- get_varcomp_component(vc_blue, "units!R")
    mean_trait <- mean(df_sub[[trait]], na.rm = TRUE)
    cv_error <- (sqrt(var_R) / mean_trait) * 100

    cv_list <- rbind(
      cv_list,
      data.frame(
        Trait = trait,
        Residual_Variance = var_R,
        Mean = mean_trait,
        CV_percent = cv_error
      )
    )

    ########################################################
    # Model 2: Variance components and broad-sense heritability
    ########################################################

    model_var <- asreml::asreml(
      fixed = as.formula(paste(trait, "~ 1")),
      random = ~ ID + Environment,
      data = df_sub,
      workspace = workspace,
      pworkspace = pworkspace
    )

    vc_var <- summary(model_var)$varcomp

    var_ID <- get_varcomp_component(vc_var, "ID")
    var_Environment <- get_varcomp_component(vc_var, "Environment")
    var_Residual <- get_varcomp_component(vc_var, "units!R")

    # Broad-sense heritability approximation.
    # This formula partitions the total phenotypic variance into genotype,
    # environment and residual components.
    h2 <- var_ID / (var_ID + var_Environment + var_Residual)

    varcomp_list <- rbind(
      varcomp_list,
      data.frame(
        Trait = trait,
        Genotype = var_ID,
        Environment = var_Environment,
        Residual = var_Residual,
        h2 = h2
      )
    )

    ########################################################
    # Model 3: BLUPs
    ########################################################

    model_blup <- asreml::asreml(
      fixed = as.formula(paste(trait, "~ Environment")),
      random = ~ ID,
      data = df_sub,
      workspace = workspace,
      pworkspace = pworkspace
    )

    blup <- predict(model_blup, classify = "ID")$pvals
    blup$Trait <- trait
    blup_list[[trait]] <- blup
  }

  ##########################################################
  # Export results
  ##########################################################

  blues_all <- do.call(rbind, blue_list)
  blups_all <- do.call(rbind, blup_list)

  if (!is.null(blues_all) && nrow(blues_all) > 0) {
    blues_all$CV_BLUE <- (
      blues_all$std.error / abs(blues_all$predicted.value)
    ) * 100

    write.table(
      blues_all,
      file = file.path(output_dir, "BLUEs_traits.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  if (!is.null(blups_all) && nrow(blups_all) > 0) {
    write.table(
      blups_all,
      file = file.path(output_dir, "BLUPs_traits.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  write.table(
    varcomp_list,
    file = file.path(output_dir, "Variance_components_h2.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  write.table(
    cv_list,
    file = file.path(output_dir, "Residual_CV_by_trait.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  return(
    list(
      BLUEs = blues_all,
      BLUPs = blups_all,
      Variance_components = varcomp_list,
      Residual_CV = cv_list
    )
  )
}



############################################################
# 6. BLUE histograms by family/group
############################################################

prepare_blues_family_table <- function(
    blues_all,
    families,
    family_order = NULL
) {
  blues_all <- blues_all %>%
    mutate(ID = as.character(ID))

  blues_family <- blues_all %>%
    left_join(
      families %>% select(ID, Family, Family_label),
      by = "ID"
    )

  missing_family <- blues_family %>%
    filter(is.na(Family)) %>%
    distinct(ID)

  cat(
    "Number of IDs without family/group information:",
    nrow(missing_family),
    "\n"
  )

  if (nrow(missing_family) > 0) {
    print(missing_family)
  }

  blues_family <- blues_family %>%
    filter(!is.na(Family))

  if (!is.null(family_order)) {
    blues_family$Family_label <- factor(
      blues_family$Family_label,
      levels = clean_label(family_order)
    )
  }

  return(blues_family)
}

plot_blue_histograms_by_group <- function(
    blues_family,
    trait_group,
    output_name,
    plot_title = NULL,
    bins = 30,
    width = 12,
    height = 9
) {
  traits_available <- trait_group[trait_group %in% unique(blues_family$Trait)]
  traits_missing <- setdiff(trait_group, unique(blues_family$Trait))

  cat("\nPlot:", output_name, "\n")
  cat("Traits included:\n")
  print(traits_available)

  if (length(traits_missing) > 0) {
    cat("Traits not found in BLUEs and skipped:\n")
    print(traits_missing)
  }

  if (length(traits_available) == 0) {
    stop("No traits from trait_group were found in blues_family.")
  }

  plot_df <- blues_family %>%
    filter(Trait %in% traits_available) %>%
    mutate(
      Trait_label = factor(
        clean_label(Trait),
        levels = clean_label(traits_available)
      )
    )

  mean_df <- plot_df %>%
    group_by(Trait_label, Family_label) %>%
    summarise(
      Mean_BLUE = mean(predicted.value, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(plot_df, aes(x = predicted.value)) +
    geom_histogram(
      bins = bins,
      fill = "#4E79A7",
      color = "grey92",
      linewidth = 0.15
    ) +
    geom_vline(
      data = mean_df,
      aes(xintercept = Mean_BLUE),
      color = "#08306B",
      linewidth = 0.5,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = mean_df,
      aes(
        x = Mean_BLUE,
        y = Inf,
        label = sprintf("%.2f", Mean_BLUE)
      ),
      color = "#08306B",
      angle = 90,
      vjust = 1.3,
      hjust = 1.1,
      size = 2.5,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    facet_grid(
      Family_label ~ Trait_label,
      scales = "free"
    ) +
    labs(
      x = "BLUE",
      y = "Frequency",
      title = plot_title
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.background = element_rect(
        fill = "#4E79A7",
        color = "#4E79A7",
        linewidth = 0.25
      ),
      strip.text.x = element_text(
        size = 8,
        face = "bold",
        angle = 90,
        color = "white",
        margin = ggplot2::margin(t = 4, r = 2, b = 4, l = 2)
      ),
      strip.text.y = element_text(
        size = 8,
        face = "bold",
        color = "white",
        margin = ggplot2::margin(t = 2, r = 4, b = 2, l = 4)
      ),
      axis.text.x = element_text(
        size = 6,
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 6,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 10,
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 10,
        margin = ggplot2::margin(r = 8)
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.spacing.x = grid::unit(0.35, "lines"),
      panel.spacing.y = grid::unit(0.35, "lines"),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      ),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  print(p)

  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".png")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  return(p)
}


############################################################
# 7. Violin plots: environment and BLUEs by group
############################################################

plot_violin_environment_and_blues_group <- function(
    df,
    blues_family,
    trait_group,
    output_name,
    plot_title = NULL,
    width = 14,
    height = 10
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for this combined plot.")
  }

  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Package 'ggh4x' is required for this combined plot.")
  }

  traits_available_env <- trait_group[trait_group %in% colnames(df)]
  traits_available_blue <- trait_group[trait_group %in% unique(blues_family$Trait)]

  traits_available <- intersect(traits_available_env, traits_available_blue)
  traits_missing <- setdiff(trait_group, traits_available)

  cat("\nPlot:", output_name, "\n")
  cat("Traits included:\n")
  print(traits_available)

  if (length(traits_missing) > 0) {
    cat("Traits skipped because they were not found in both tables:\n")
    print(traits_missing)
  }

  if (length(traits_available) == 0) {
    stop("No traits from trait_group were found in both df and blues_family.")
  }

  env_long <- df %>%
    select(Environment, all_of(traits_available)) %>%
    pivot_longer(
      cols = all_of(traits_available),
      names_to = "Trait",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value)) %>%
    mutate(
      Trait_label = factor(
        clean_label(Trait),
        levels = clean_label(traits_available)
      ),
      Environment_label = clean_label(as.character(Environment))
    )

  env_long$Environment_label <- factor(
    env_long$Environment_label,
    levels = clean_label(levels(df$Environment))
  )

  mean_env <- env_long %>%
    group_by(Trait_label, Environment_label) %>%
    summarise(
      Mean_value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )

  blue_long <- blues_family %>%
    filter(Trait %in% traits_available) %>%
    mutate(
      Trait_label = factor(
        clean_label(Trait),
        levels = clean_label(traits_available)
      )
    ) %>%
    filter(!is.na(predicted.value))

  mean_blue <- blue_long %>%
    group_by(Trait_label, Family_label) %>%
    summarise(
      Mean_value = mean(predicted.value, na.rm = TRUE),
      .groups = "drop"
    )

  paper_palette_env <- c(
    "#3B4252",
    "#5E81AC",
    "#81A1C1",
    "#88C0D0",
    "#8FBCBB",
    "#A3BE8C",
    "#B48EAD",
    "#D8A657"
  )

  paper_palette_family <- c(
    "#3B4252",
    "#6B8E73",
    "#5E81AC",
    "#A3BE8C",
    "#B48EAD",
    "#D08770",
    "#88C0D0",
    "#7A8C99"
  )

  env_levels <- levels(env_long$Environment_label)
  group_levels <- levels(blue_long$Family_label)
  group_levels <- group_levels[!is.na(group_levels)]

  env_colors <- make_named_palette(
    levels_vector = env_levels,
    base_palette = paper_palette_env
  )

  group_colors <- make_named_palette(
    levels_vector = group_levels,
    base_palette = paper_palette_family
  )

  theme_paper <- theme_classic(base_size = 10) +
    theme(
      strip.background = element_rect(
        fill = "#44546A",
        color = "#44546A",
        linewidth = 0.25
      ),
      strip.text = element_text(
        size = 8,
        face = "bold",
        color = "white"
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(
        size = 6,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 10,
        margin = margin(r = 8)
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks.length = unit(1.5, "mm"),
      panel.spacing = unit(0.55, "lines"),
      legend.position = "bottom",
      legend.title = element_text(
        size = 9,
        face = "bold"
      ),
      legend.text = element_text(
        size = 8
      ),
      legend.key.size = unit(0.35, "cm"),
      legend.box = "vertical",
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0
      ),
      plot.margin = margin(6, 8, 6, 8)
    )

  p_env <- ggplot(
    env_long,
    aes(
      x = Environment_label,
      y = Value,
      fill = Environment_label
    )
  ) +
    geom_violin(
      color = "grey30",
      linewidth = 0.3,
      trim = FALSE,
      scale = "width",
      alpha = 0.9
    ) +
    geom_boxplot(
      width = 0.14,
      fill = "white",
      color = "grey15",
      linewidth = 0.3,
      outlier.shape = NA
    ) +
    geom_point(
      data = mean_env,
      aes(
        x = Environment_label,
        y = Mean_value
      ),
      color = "#1F2A44",
      size = 1.5,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = mean_env,
      aes(
        x = Environment_label,
        y = Mean_value,
        label = sprintf("%.2f", Mean_value)
      ),
      fill = "white",
      color = "#1F2A44",
      size = 2.25,
      fontface = "bold",
      label.size = 0.12,
      label.padding = unit(0.08, "lines"),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = env_colors,
      name = "Environment"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    ) +
    ggh4x::facet_wrap2(
      ~ Trait_label,
      scales = "free_y",
      axes = "all",
      nrow = 1
    ) +
    labs(
      x = NULL,
      y = "Phenotypic value",
      title = "Environment"
    ) +
    theme_paper

  p_blue <- ggplot(
    blue_long,
    aes(
      x = Family_label,
      y = predicted.value,
      fill = Family_label
    )
  ) +
    geom_violin(
      color = "grey30",
      linewidth = 0.3,
      trim = FALSE,
      scale = "width",
      alpha = 0.9
    ) +
    geom_boxplot(
      width = 0.14,
      fill = "white",
      color = "grey15",
      linewidth = 0.3,
      outlier.shape = NA
    ) +
    geom_point(
      data = mean_blue,
      aes(
        x = Family_label,
        y = Mean_value
      ),
      color = "#1F2A44",
      size = 1.5,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = mean_blue,
      aes(
        x = Family_label,
        y = Mean_value,
        label = sprintf("%.2f", Mean_value)
      ),
      fill = "white",
      color = "#1F2A44",
      size = 2.25,
      fontface = "bold",
      label.size = 0.12,
      label.padding = unit(0.08, "lines"),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = group_colors,
      name = "Group"
    ) +
    guides(
      fill = guide_legend(
        nrow = 2,
        byrow = TRUE
      )
    ) +
    ggh4x::facet_wrap2(
      ~ Trait_label,
      scales = "free_y",
      axes = "all",
      nrow = 1
    ) +
    labs(
      x = NULL,
      y = "BLUE",
      title = "BLUEs by group"
    ) +
    theme_paper

  p_final <- p_env / p_blue +
    patchwork::plot_layout(
      heights = c(1, 1.15),
      guides = "keep"
    ) +
    patchwork::plot_annotation(
      title = plot_title,
      theme = theme(
        plot.title = element_text(
          size = 13,
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 8)
        )
      )
    )

  print(p_final)

  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".png")),
    plot = p_final,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".pdf")),
    plot = p_final,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  return(p_final)
}


############################################################
# 8. Main execution
############################################################

# Read phenotype data.
df <- read_phenotype_data(
  phenotype_file = phenotype_file,
  id_col = id_col,
  location_col = location_col,
  year_col = year_col,
  environment_col = environment_col
)

# Inspect data structure.
str(df)

############################################################
# 8.1. Descriptive plots
############################################################

p_hist_environment <- plot_trait_histograms(
  df = df,
  traits = trait_groups$Quality_traits,
  group_col = "Environment",
  output_name = "Fig_histograms_quality_traits_by_environment",
  plot_title = "Quality traits by environment",
  width = 12,
  height = 9
)

p_hist_year <- plot_trait_histograms(
  df = df,
  traits = trait_groups$Quality_traits,
  group_col = "Year",
  output_name = "Fig_histograms_quality_traits_by_year",
  plot_title = "Quality traits by year",
  width = 12,
  height = 9
)


############################################################
# 8.2. BLUEs, BLUPs, variance components and heritability
############################################################

if (run_asreml_analysis) {
  if (!has_asreml) {
    stop(
      "run_asreml_analysis is TRUE, but ASReml is not available. ",
      "Install/load ASReml or set run_asreml_analysis <- FALSE.",
      call. = FALSE
    )
  }

  asreml_results <- run_asreml_trait_models(
    df = df,
    traits = traits,
    id_col = id_col,
    environment_col = "Environment",
    output_dir = output_dir,
    workspace = asreml_workspace,
    pworkspace = asreml_pworkspace
  )

  p_varcomp <- plot_variance_components(
    varcomp_list = asreml_results$Variance_components,
    output_name = "Fig_variance_components_h2",
    plot_title = "Variance components and broad-sense heritability",
    width = 10,
    height = 6
  )
}

############################################################
# 8.3. BLUE plots by family/group
############################################################

if (run_asreml_analysis && file.exists(families_file)) {
  families <- read_family_data(
    families_file = families_file,
    family_id_col = family_id_col,
    family_group_col = family_group_col,
    family_order = family_order
  )

  blues_family <- prepare_blues_family_table(
    blues_all = asreml_results$BLUEs,
    families = families,
    family_order = family_order
  )

  p_blue_hist_group <- plot_blue_histograms_by_group(
    blues_family = blues_family,
    trait_group = trait_groups$Quality_traits,
    output_name = "Fig_BLUE_histograms_by_group_quality_traits",
    plot_title = "Quality traits",
    width = 12,
    height = 9
  )

  if (has_patchwork && has_ggh4x) {
    p_violin <- plot_violin_environment_and_blues_group(
      df = df,
      blues_family = blues_family,
      trait_group = trait_groups$Quality_traits,
      output_name = "Fig_violin_environment_BLUEs_group_quality_traits",
      plot_title = "Quality traits",
      width = 14,
      height = 10
    )
  }
} else {
  cat(
    "Family/group BLUE plots were skipped. ",
    "Check run_asreml_analysis and families_file.\n"
  )
}


