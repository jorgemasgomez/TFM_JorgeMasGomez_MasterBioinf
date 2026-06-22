############################################################
# Script 03 - Multivariate GWAS Circos Plot
#
# Purpose:
#   Generate a publication-quality circular GWAS plot for almond,
#   including:
#     - 8 Prunus dulcis chromosomes
#     - one Manhattan-style track per trait group
#     - SNP density track
#     - hotspot links between genomic regions
#     - compact legend with trait colors and R2 size scale
#
# Expected input:
#   - resultado_filtrado_por_grupos.txt
#
# Main outputs:
#   - Circos_GWAS_Final_Report.png
#   - Circos_GWAS_Final_Report_paper.png
#   - filtered_GWAS_for_circos.csv
#   - hotspots_by_group.txt
#   - hotspot_thresholds_by_group.txt
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

load_package("circlize")
load_package("dplyr")
load_package("grid")
load_package("png")


############################################################
# 1. User configuration
############################################################

input_file <- "resultado_filtrado_por_grupos.txt"
output_dir <- "GWAS_circos_outputs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Almond chromosome sizes.
chr_data <- data.frame(
  chr = c("Pd01", "Pd02", "Pd03", "Pd04", "Pd05", "Pd06", "Pd07", "Pd08"),
  start = rep(0, 8),
  end = c(
    50452767,
    30972344,
    30498384,
    28068711,
    22241271,
    32371708,
    24859559,
    27108964
  )
)

# Significance and hotspot settings.
lod_threshold <- 3.0
hotspot_window_size <- 250000
min_prop_traits_hotspot <- 0.20
min_abs_traits_hotspot <- 2
max_hotspots_per_group <- 10

# Filter requiring the same association to be detected by at least this many methods.
min_methods_per_trait <- 2


############################################################
# 2. Helper functions
############################################################

prepare_gwas_data <- function(
    input_file,
    min_methods_per_trait = 2
) {
  raw_data <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    sep = "\t",
    check.names = FALSE
  )

  required_cols <- c(
    "Trait.name",
    "Group",
    "Chromosome",
    "Marker.position..bp.",
    "LOD.score",
    "r2....",
    "Method"
  )

  missing_cols <- setdiff(required_cols, colnames(raw_data))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required columns are missing: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  gwas_prepared <- raw_data %>%
    rename(
      Trait = Trait.name,
      Group = Group,
      Chromosome_raw = Chromosome,
      Position = Marker.position..bp.,
      LOD = LOD.score,
      R2 = r2....,
      Method = Method
    ) %>%
    mutate(
      Chromosome = sprintf("Pd%02d", as.numeric(Chromosome_raw)),
      LOD = as.numeric(LOD),
      R2 = as.numeric(R2),
      Position = as.numeric(Position)
    ) %>%
    filter(
      !is.na(Chromosome),
      !is.na(Position),
      !is.na(LOD),
      !is.na(R2)
    )

  gwas_filtered <- gwas_prepared %>%
    group_by(Trait, Chromosome, Position) %>%
    mutate(Num_Methods = n_distinct(Method)) %>%
    ungroup() %>%
    filter(Num_Methods >= min_methods_per_trait)

  if (nrow(gwas_filtered) == 0) {
    stop("The dataset is empty after filtering by method count.", call. = FALSE)
  }

  return(gwas_filtered)
}

assign_trait_colors_and_sizes <- function(gwas_data) {
  high_contrast_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
    "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#000080"
  )

  groups <- sort(unique(gwas_data$Group))
  color_map <- c()

  for (group_name in groups) {
    group_traits <- sort(unique(gwas_data$Trait[gwas_data$Group == group_name]))
    assigned_colors <- rep(high_contrast_colors, length.out = length(group_traits))
    names(assigned_colors) <- group_traits
    color_map <- c(color_map, assigned_colors)
  }

  gwas_data$Color <- color_map[gwas_data$Trait]

  min_r2 <- min(gwas_data$R2, na.rm = TRUE)
  max_r2 <- max(gwas_data$R2, na.rm = TRUE)

  if (max_r2 == min_r2) {
    gwas_data$Cex_Size <- 0.8
  } else {
    gwas_data <- gwas_data %>%
      mutate(
        Cex_Size = 0.3 + 1.1 * ((R2 - min_r2) / (max_r2 - min_r2))
      )
  }

  return(
    list(
      data = gwas_data,
      color_map = color_map,
      groups = groups,
      min_r2 = min_r2,
      max_r2 = max_r2
    )
  )
}

detect_hotspots <- function(
    gwas_data,
    chr_data,
    lod_threshold,
    hotspot_window_size,
    min_prop_traits_hotspot,
    min_abs_traits_hotspot,
    max_hotspots_per_group
) {
  significant_snps <- gwas_data %>%
    filter(LOD >= lod_threshold)

  traits_by_group <- gwas_data %>%
    group_by(Group) %>%
    summarise(
      n_traits_group = n_distinct(Trait),
      min_traits_hotspot = pmax(
        min_abs_traits_hotspot,
        ceiling(n_traits_group * min_prop_traits_hotspot)
      ),
      .groups = "drop"
    )

  if (nrow(significant_snps) == 0) {
    return(
      list(
        significant_snps = significant_snps,
        traits_by_group = traits_by_group,
        hotspots_by_group = data.frame()
      )
    )
  }

  hotspots_by_group <- significant_snps %>%
    left_join(
      chr_data %>% rename(Chromosome = chr, chr_end = end),
      by = "Chromosome"
    ) %>%
    mutate(
      Hotspot_window = floor(Position / hotspot_window_size),
      Hotspot_start = Hotspot_window * hotspot_window_size,
      Hotspot_end = pmin(Hotspot_start + hotspot_window_size, chr_end),
      Hotspot_mid = (Hotspot_start + Hotspot_end) / 2
    ) %>%
    group_by(
      Group,
      Chromosome,
      Hotspot_window,
      Hotspot_start,
      Hotspot_end,
      Hotspot_mid
    ) %>%
    summarise(
      n_assoc = n(),
      n_traits = n_distinct(Trait),
      traits_in_hotspot = paste(sort(unique(Trait)), collapse = "; "),
      max_LOD = max(LOD, na.rm = TRUE),
      mean_LOD = mean(LOD, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(traits_by_group, by = "Group") %>%
    filter(n_traits >= min_traits_hotspot) %>%
    group_by(Group) %>%
    arrange(desc(n_traits), desc(n_assoc), desc(max_LOD), .by_group = TRUE) %>%
    slice_head(n = max_hotspots_per_group) %>%
    arrange(Chromosome, Hotspot_start, .by_group = TRUE) %>%
    ungroup()

  return(
    list(
      significant_snps = significant_snps,
      traits_by_group = traits_by_group,
      hotspots_by_group = hotspots_by_group
    )
  )
}

draw_circos_plot <- function(
    gwas_data,
    significant_snps,
    hotspots_by_group,
    chr_data,
    groups,
    color_map,
    lod_threshold,
    min_r2,
    max_r2,
    output_file
) {
  link_palette <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
    "#66a61e", "#e6ab02", "#a6761d", "#1f78b4",
    "#b2df8a", "#fb9a99"
  )

  group_link_colors <- setNames(
    rep(link_palette, length.out = length(groups)),
    groups
  )

  png(output_file, width = 3600, height = 4000, res = 400, bg = "white")
  par(mar = c(0.05, 0.05, 0.05, 0.05))

  circos.clear()
  circos.par(
    canvas.xlim = c(-1.02, 1.02),
    canvas.ylim = c(-1.02, 1.02),
    gap.after = c(rep(5, 7), 20),
    start.degree = 98,
    track.margin = c(0.001, 0.001),
    cell.padding = c(0, 0, 0, 0),
    circle.margin = 0.01
  )

  circos.genomicInitialize(
    chr_data,
    plotType = c("axis", "labels"),
    labels.cex = 0.7
  )

  if (nrow(significant_snps) > 0) {
    density_df <- significant_snps %>%
      transmute(
        chr = Chromosome,
        start = Position,
        end = Position + 1
      )

    circos.genomicDensity(
      density_df,
      col = "#5a6364",
      bg.col = "#f2f4f4",
      track.height = 0.05,
      border = "grey50"
    )
  }

  for (group_name in groups) {
    group_data <- gwas_data %>%
      filter(Group == group_name)

    max_lod_group <- max(group_data$LOD, na.rm = TRUE)
    if (max_lod_group < lod_threshold) {
      max_lod_group <- lod_threshold + 1
    }

    axis_ceiling <- ceiling(max_lod_group)
    if (axis_ceiling %% 2 != 0) {
      axis_ceiling <- axis_ceiling + 1
    }

    axis_ticks <- c(0, axis_ceiling / 2, axis_ceiling)
    axis_ticks <- axis_ticks[axis_ticks > 0]

    circos.track(
      ylim = c(0, axis_ceiling),
      panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        xlim <- CELL_META$xlim

        chromosome_data <- group_data %>%
          filter(Chromosome == chr)

        circos.lines(
          xlim,
          rep(lod_threshold, 2),
          lty = 3,
          col = "grey85",
          lwd = 0.7
        )

        if (chr == "Pd01") {
          circos.lines(
            rep(xlim[1], 2),
            c(0, axis_ceiling),
            col = "black",
            lwd = 1
          )

          for (tick_value in axis_ticks) {
            tick_pos_x <- xlim[1] - (xlim[2] - xlim[1]) * 0.012

            circos.lines(
              c(xlim[1], tick_pos_x),
              rep(tick_value, 2),
              col = "black",
              lwd = 1
            )

            circos.text(
              x = tick_pos_x,
              y = tick_value,
              labels = as.character(tick_value),
              sector.index = "Pd01",
              facing = "downward",
              adj = c(1.2, 0.4),
              cex = 0.4,
              font = 1
            )
          }
        }

        if (nrow(chromosome_data) > 0) {
          for (i in seq_len(nrow(chromosome_data))) {
            point_color <- chromosome_data$Color[i]

            if (chromosome_data$LOD[i] < lod_threshold) {
              point_color <- paste0(point_color, "35")
            }

            circos.points(
              x = chromosome_data$Position[i],
              y = chromosome_data$LOD[i],
              col = point_color,
              pch = 16,
              cex = chromosome_data$Cex_Size[i]
            )
          }
        }
      },
      track.height = 0.065,
      bg.col = "#fefefe",
      bg.border = "grey88"
    )

    circos.text(
      CELL_META$xlim[2] + 12000000,
      axis_ceiling / 2,
      group_name,
      sector.index = "Pd08",
      cex = 0.4,
      font = 2,
      adj = c(1, 0.4),
      facing = "downward"
    )
  }

  if (nrow(hotspots_by_group) > 0) {
    link_radius <- get_most_inside_radius() - 0.01

    for (group_name in groups) {
      group_hotspots <- hotspots_by_group %>%
        filter(Group == group_name) %>%
        arrange(Chromosome, Hotspot_start)

      n_hotspots <- nrow(group_hotspots)

      if (n_hotspots == 2) {
        circos.link(
          sector.index1 = group_hotspots$Chromosome[1],
          point1 = group_hotspots$Hotspot_mid[1],
          sector.index2 = group_hotspots$Chromosome[2],
          point2 = group_hotspots$Hotspot_mid[2],
          rou1 = link_radius,
          rou2 = link_radius,
          col = adjustcolor(group_link_colors[group_name], alpha.f = 0.35),
          border = adjustcolor(group_link_colors[group_name], alpha.f = 0.55),
          lwd = 1.4,
          h.ratio = 0.65
        )
      }

      if (n_hotspots >= 3) {
        for (i in seq_len(n_hotspots - 1)) {
          circos.link(
            sector.index1 = group_hotspots$Chromosome[i],
            point1 = group_hotspots$Hotspot_mid[i],
            sector.index2 = group_hotspots$Chromosome[i + 1],
            point2 = group_hotspots$Hotspot_mid[i + 1],
            rou1 = link_radius,
            rou2 = link_radius,
            col = adjustcolor(group_link_colors[group_name], alpha.f = 0.35),
            border = adjustcolor(group_link_colors[group_name], alpha.f = 0.55),
            lwd = 1.4,
            h.ratio = 0.65
          )
        }

        circos.link(
          sector.index1 = group_hotspots$Chromosome[n_hotspots],
          point1 = group_hotspots$Hotspot_mid[n_hotspots],
          sector.index2 = group_hotspots$Chromosome[1],
          point2 = group_hotspots$Hotspot_mid[1],
          rou1 = link_radius,
          rou2 = link_radius,
          col = adjustcolor(group_link_colors[group_name], alpha.f = 0.35),
          border = adjustcolor(group_link_colors[group_name], alpha.f = 0.55),
          lwd = 1.4,
          h.ratio = 0.65
        )
      }
    }
  }

  dev.off()

  return(group_link_colors)
}

draw_compact_legend <- function(
    gwas_data,
    groups,
    color_map,
    group_link_colors,
    min_r2,
    max_r2,
    output_file
) {
  png(output_file, width = 3600, height = 1450, res = 400, bg = "white")

  par(mar = c(0.02, 0.02, 0.02, 0.02))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  pos_x_cols <- seq(0.05, 0.84, length.out = 7)
  current_col <- 1
  y_start <- 0.93
  y_step <- 0.038

  max_traits_per_group <- max(sapply(groups, function(group_name) {
    length(unique(gwas_data$Trait[gwas_data$Group == group_name]))
  }))

  y_bottom_traits <- y_start - 0.04 - ((max_traits_per_group - 1) * y_step)
  y_link_legend <- y_bottom_traits - 0.05
  y_r2 <- y_link_legend - 0.08
  box_bottom <- y_r2 - 0.055

  rect(
    xleft = 0.03,
    ybottom = box_bottom,
    xright = 0.97,
    ytop = 0.98,
    col = "#fafafa",
    border = "grey75",
    lwd = 1.2
  )

  for (group_name in groups) {
    group_traits <- sort(unique(gwas_data$Trait[gwas_data$Group == group_name]))
    x_coord <- pos_x_cols[current_col]

    text(
      x_coord,
      y_start,
      labels = group_name,
      adj = 0,
      font = 2,
      cex = 0.62,
      col = "black"
    )

    segments(
      x_coord,
      y_start - 0.012,
      x_coord + 0.10,
      y_start - 0.012,
      col = "grey70",
      lwd = 0.8
    )

    y_item <- y_start - 0.04

    for (trait_name in group_traits) {
      points(
        x_coord + 0.005,
        y_item,
        col = color_map[trait_name],
        pch = 16,
        cex = 0.95
      )

      text(
        x_coord + 0.016,
        y_item,
        labels = trait_name,
        adj = 0,
        cex = 0.42,
        col = "grey20"
      )

      y_item <- y_item - y_step
    }

    current_col <- current_col + 1
  }

  hotspot_title <- "Hotspot links:"
  text(
    0.05,
    y_link_legend,
    labels = hotspot_title,
    font = 2,
    cex = 0.58,
    adj = 0
  )

  x_cursor <- 0.18

  for (group_name in groups) {
    segments(
      x0 = x_cursor,
      y0 = y_link_legend,
      x1 = x_cursor + 0.018,
      y1 = y_link_legend,
      col = adjustcolor(group_link_colors[group_name], alpha.f = 0.9),
      lwd = 2.4
    )

    text(
      x = x_cursor + 0.026,
      y = y_link_legend,
      labels = group_name,
      adj = 0,
      cex = 0.35,
      col = "grey20"
    )

    x_cursor <- x_cursor + 0.10
  }

  r2_multiplier <- if (max_r2 <= 1) 100 else 1

  r2_min <- max(1, round(min_r2 * r2_multiplier))
  r2_max <- round(max_r2 * r2_multiplier)
  r2_mid <- round((r2_min + r2_max) / 2)

  if (max_r2 == min_r2) {
    example_cex <- c(0.8, 0.8, 0.8)
  } else {
    source_r2_values <- c(r2_min, r2_mid, r2_max) / r2_multiplier
    example_cex <- 0.3 + 1.1 * ((source_r2_values - min_r2) / (max_r2 - min_r2))
  }

  r2_labels <- paste0(c(r2_min, r2_mid, r2_max), "%")

  text(
    0.05,
    y_r2,
    labels = "R2 scale (point size):",
    font = 2,
    cex = 0.58,
    adj = 0
  )

  x_cursor <- 0.25

  for (i in seq_along(r2_labels)) {
    points(
      x_cursor,
      y_r2,
      col = "#5a6364",
      pch = 16,
      cex = example_cex[i] * 1.35
    )

    text(
      x_cursor + 0.020,
      y_r2,
      labels = r2_labels[i],
      adj = 0,
      cex = 0.48,
      col = "grey30"
    )

    x_cursor <- x_cursor + 0.11
  }

  dev.off()
}

merge_vertical_pngs <- function(
    top_image,
    bottom_image,
    output_file,
    res = 400
) {
  top <- png::readPNG(top_image)
  bottom <- png::readPNG(bottom_image)

  top_h <- dim(top)[1]
  top_w <- dim(top)[2]
  bottom_h <- dim(bottom)[1]
  bottom_w <- dim(bottom)[2]

  if (top_w != bottom_w) {
    stop("Both images must have the same width.", call. = FALSE)
  }

  png(
    output_file,
    width = top_w,
    height = top_h + bottom_h,
    res = res,
    bg = "white"
  )

  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 1,
        heights = unit(c(top_h, bottom_h), "null")
      )
    )
  )

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.raster(
    top,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = FALSE
  )
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.raster(
    bottom,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = FALSE
  )
  popViewport()

  dev.off()
}


############################################################
# 3. Main execution
############################################################

gwas_data <- prepare_gwas_data(
  input_file = input_file,
  min_methods_per_trait = min_methods_per_trait
)

plot_data <- assign_trait_colors_and_sizes(gwas_data)
gwas_data <- plot_data$data

hotspot_results <- detect_hotspots(
  gwas_data = gwas_data,
  chr_data = chr_data,
  lod_threshold = lod_threshold,
  hotspot_window_size = hotspot_window_size,
  min_prop_traits_hotspot = min_prop_traits_hotspot,
  min_abs_traits_hotspot = min_abs_traits_hotspot,
  max_hotspots_per_group = max_hotspots_per_group
)

readr::write_csv(
  gwas_data,
  file.path(output_dir, "filtered_GWAS_for_circos.csv")
)

write.table(
  hotspot_results$traits_by_group,
  file.path(output_dir, "hotspot_thresholds_by_group.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  hotspot_results$hotspots_by_group,
  file.path(output_dir, "hotspots_by_group.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

circos_file <- file.path(output_dir, "temp_circos.png")
legend_file <- file.path(output_dir, "temp_legend.png")

group_link_colors <- draw_circos_plot(
  gwas_data = gwas_data,
  significant_snps = hotspot_results$significant_snps,
  hotspots_by_group = hotspot_results$hotspots_by_group,
  chr_data = chr_data,
  groups = plot_data$groups,
  color_map = plot_data$color_map,
  lod_threshold = lod_threshold,
  min_r2 = plot_data$min_r2,
  max_r2 = plot_data$max_r2,
  output_file = circos_file
)

draw_compact_legend(
  gwas_data = gwas_data,
  groups = plot_data$groups,
  color_map = plot_data$color_map,
  group_link_colors = group_link_colors,
  min_r2 = plot_data$min_r2,
  max_r2 = plot_data$max_r2,
  output_file = legend_file
)

merge_vertical_pngs(
  top_image = circos_file,
  bottom_image = legend_file,
  output_file = file.path(output_dir, "Circos_GWAS_Final_Report.png"),
  res = 400
)

merge_vertical_pngs(
  top_image = circos_file,
  bottom_image = legend_file,
  output_file = file.path(output_dir, "Circos_GWAS_Final_Report_paper.png"),
  res = 400
)

cat("\nCircos plot completed successfully.\n")
cat("Final report saved to:", file.path(output_dir, "Circos_GWAS_Final_Report.png"), "\n")

############################################################
# End of script
############################################################
