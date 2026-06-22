############################################################
# Script 04 - Allelic Effect Boxplots for GWAS Markers
#
# Purpose:
#   Generate violin + boxplot figures showing allelic/genotypic
#   effects for selected marker-trait combinations.
#
# Expected inputs:
#   - myG.txt                 : HapMap genotype file
#   - myY.txt                 : phenotype or BLUE table
#   - marcadores_selected.txt : table with marker-trait combinations
#
# Required columns in marcadores_selected.txt:
#   - ID
#   - Trait
#
# Main outputs:
#   - one PNG and one PDF per marker-trait combination
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
load_package("readr")
load_package("ggplot2")


############################################################
# 1. User configuration
############################################################

genotype_file <- "myG.txt"
phenotype_file <- "myY.txt"
marker_trait_file <- "marcadores_selected.txt"

output_dir <- "allelic_effect_boxplots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Marker and trait columns in marker_trait_file.
marker_col <- "ID"
trait_col <- "Trait"

# Missing genotype codes to exclude.
missing_genotypes <- c("NN", "N", "--", "", NA)

# Genotype color palette.
genotype_palette <- c(
  "AA" = "#a6cee3",
  "CC" = "#1f78b4",
  "GG" = "#b2df8a",
  "TT" = "#33a02c",
  "AC" = "#fb9a99",
  "AG" = "#e31a1c",
  "AT" = "#fdbf6f",
  "CG" = "#ff7f00",
  "CT" = "#cab2d6",
  "GT" = "#6a3d9a"
)


############################################################
# 2. Helper functions
############################################################

read_input_tables <- function(
    genotype_file,
    phenotype_file,
    marker_trait_file
) {
  genotype <- readr::read_delim(
    genotype_file,
    delim = "\t",
    show_col_types = FALSE
  )

  phenotype <- readr::read_delim(
    phenotype_file,
    delim = "\t",
    show_col_types = FALSE
  )

  marker_traits <- readr::read_delim(
    marker_trait_file,
    delim = "\t",
    show_col_types = FALSE
  )

  colnames(phenotype)[1] <- "ID"

  required_marker_cols <- c(marker_col, trait_col)
  missing_cols <- setdiff(required_marker_cols, colnames(marker_traits))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The marker-trait file is missing the following columns: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  return(
    list(
      genotype = genotype,
      phenotype = phenotype,
      marker_traits = marker_traits
    )
  )
}

extract_marker_genotypes <- function(genotype, marker_id) {
  marker_row <- genotype %>%
    filter(`rs#` == marker_id)

  if (nrow(marker_row) == 0) {
    warning(
      paste0("Marker not found in genotype file: ", marker_id),
      call. = FALSE
    )
    return(NULL)
  }

  marker_info <- marker_row %>%
    select(chrom, pos)

  genotype_table <- marker_row %>%
    select(-chrom, -pos, -strand, -assembly, -center, -protLSID,
           -assayLSID, -panelLSID, -QCcode, -alleles, `-rs#`) %>%
    t() %>%
    as.data.frame()

  # Fallback for HapMap files where some metadata columns may not exist.
  if (ncol(genotype_table) == 0) {
    genotype_table <- marker_row %>%
      select(-matches("^(rs#|alleles|chrom|pos|strand|assembly|center|protLSID|assayLSID|panelLSID|QCcode)$")) %>%
      t() %>%
      as.data.frame()
  }

  colnames(genotype_table) <- "Genotype"
  genotype_table$ID <- rownames(genotype_table)

  return(
    list(
      marker_info = marker_info,
      genotype_table = genotype_table
    )
  )
}

prepare_allelic_effect_data <- function(
    genotype,
    phenotype,
    marker_id,
    trait_name
) {
  if (!trait_name %in% colnames(phenotype)) {
    warning(
      paste0("Trait not found in phenotype file: ", trait_name),
      call. = FALSE
    )
    return(NULL)
  }

  marker_data <- extract_marker_genotypes(
    genotype = genotype,
    marker_id = marker_id
  )

  if (is.null(marker_data)) {
    return(NULL)
  }

  plot_data <- marker_data$genotype_table %>%
    left_join(phenotype, by = "ID") %>%
    select(ID, Genotype, all_of(trait_name)) %>%
    filter(
      !is.na(Genotype),
      !is.na(.data[[trait_name]]),
      !Genotype %in% missing_genotypes
    ) %>%
    mutate(
      Genotype = as.character(Genotype)
    )

  if (nrow(plot_data) == 0) {
    warning(
      paste0(
        "No valid observations for marker ", marker_id,
        " and trait ", trait_name
      ),
      call. = FALSE
    )
    return(NULL)
  }

  return(
    list(
      data = plot_data,
      marker_info = marker_data$marker_info
    )
  )
}

plot_allelic_effect <- function(
    genotype,
    phenotype,
    marker_id,
    trait_name,
    output_dir = "allelic_effect_boxplots",
    genotype_palette = genotype_palette
) {
  prepared <- prepare_allelic_effect_data(
    genotype = genotype,
    phenotype = phenotype,
    marker_id = marker_id,
    trait_name = trait_name
  )

  if (is.null(prepared)) {
    return(NULL)
  }

  plot_data <- prepared$data
  marker_info <- prepared$marker_info

  genotype_counts <- plot_data %>%
    group_by(Genotype) %>%
    summarise(
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(label = paste("n =", n))

  present_genotypes <- sort(unique(plot_data$Genotype))
  available_colors <- genotype_palette[
    names(genotype_palette) %in% present_genotypes
  ]

  missing_colors <- setdiff(present_genotypes, names(available_colors))
  if (length(missing_colors) > 0) {
    extra_colors <- grDevices::rainbow(length(missing_colors))
    names(extra_colors) <- missing_colors
    available_colors <- c(available_colors, extra_colors)
  }

  y_min <- min(plot_data[[trait_name]], na.rm = TRUE)
  y_max <- max(plot_data[[trait_name]], na.rm = TRUE)
  y_range <- y_max - y_min

  if (y_range == 0) {
    y_label_position <- y_min
  } else {
    y_label_position <- y_min - y_range * 0.08
  }

  chromosome_label <- marker_info$chrom[1]
  position_label <- marker_info$pos[1]

  p <- ggplot(
    plot_data,
    aes(
      x = Genotype,
      y = .data[[trait_name]],
      fill = Genotype
    )
  ) +
    geom_violin(
      alpha = 0.45,
      color = "grey30",
      trim = TRUE
    ) +
    geom_boxplot(
      width = 0.15,
      color = "black",
      alpha = 0.75,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.10,
      alpha = 0.50,
      size = 1,
      color = "black"
    ) +
    geom_text(
      data = genotype_counts,
      aes(
        x = Genotype,
        y = y_label_position,
        label = label
      ),
      size = 3.5,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = available_colors) +
    labs(
      title = paste(marker_id, "-", trait_name),
      subtitle = paste("Chr:", chromosome_label, "| Pos:", position_label),
      x = "Genotype",
      y = paste(trait_name, "(BLUE)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )

  safe_marker <- gsub("[^A-Za-z0-9_.-]", "_", marker_id)
  safe_trait <- gsub("[^A-Za-z0-9_.-]", "_", trait_name)

  png_file <- file.path(output_dir, paste0(safe_marker, "_", safe_trait, ".png"))
  pdf_file <- file.path(output_dir, paste0(safe_marker, "_", safe_trait, ".pdf"))

  ggsave(
    filename = png_file,
    plot = p,
    width = 6,
    height = 5,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    width = 6,
    height = 5,
    bg = "white"
  )

  cat("Saved:", png_file, "\n")

  return(p)
}


############################################################
# 3. Main execution
############################################################

input_tables <- read_input_tables(
  genotype_file = genotype_file,
  phenotype_file = phenotype_file,
  marker_trait_file = marker_trait_file
)

genotype <- input_tables$genotype
phenotype <- input_tables$phenotype
marker_traits <- input_tables$marker_traits

for (i in seq_len(nrow(marker_traits))) {
  marker_id <- marker_traits[[marker_col]][i]
  trait_name <- marker_traits[[trait_col]][i]

  plot_allelic_effect(
    genotype = genotype,
    phenotype = phenotype,
    marker_id = marker_id,
    trait_name = trait_name,
    output_dir = output_dir,
    genotype_palette = genotype_palette
  )
}

cat("\nAllelic effect boxplots completed.\n")
cat("Output directory:", output_dir, "\n")

############################################################
# End of script
############################################################
