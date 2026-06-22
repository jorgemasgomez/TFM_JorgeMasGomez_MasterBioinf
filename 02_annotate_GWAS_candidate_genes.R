############################################################
# Script 02 - GWAS SNP Annotation and Candidate Gene Retrieval
#
# Purpose:
#   Annotate significant GWAS markers by:
#     1. Creating SNP windows around significant markers.
#     2. Finding nearby genes from a GFF3 annotation file.
#     3. Retrieving Prunus dulcis annotations from Ensembl Plants.
#     4. Retrieving Arabidopsis ortholog annotations.
#     5. Producing collapsed candidate-gene annotation tables.
#
# Expected inputs:
#   - merged_GWAS_results.csv or resultado_unido.csv
#   - Prunus dulcis GFF3 file
#
# Main outputs:
#   - GWAS_candidate_genes_raw.tsv
#   - GWAS_candidate_genes_collapsed.tsv
#   - Arabidopsis_annotations_collapsed.tsv
#   - GWAS_candidate_genes_final_annotated.tsv
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
load_package("tidyr")
load_package("GenomicRanges")
load_package("rtracklayer")
load_package("biomaRt")


############################################################
# 1. User configuration
############################################################

gwas_results_file <- "merged_GWAS_results.csv"
gff_file <- "anotaciones_texasv2_ncbi.gff"
output_dir <- "GWAS_annotation_outputs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Window around each SNP, in kilobases.
window_kb <- 25

# Column names in the GWAS output.
snp_id_col <- "RS."
trait_col <- "Trait.name"
chromosome_col <- "Chromosome"
position_col <- "Marker.position..bp."

# Chromosome name mapping for Prunus dulcis Texas v2 / NCBI.
chrom_map <- c(
  "Pd01" = "NC_047650.1",
  "Pd02" = "NC_047651.1",
  "Pd03" = "NC_047652.1",
  "Pd04" = "NC_047653.1",
  "Pd05" = "NC_047654.1",
  "Pd06" = "NC_047655.1",
  "Pd07" = "NC_047656.1",
  "Pd08" = "NC_047657.1"
)


############################################################
# 2. Helper functions
############################################################

collapse_values <- function(x) {
  x <- unique(na.omit(x))
  x <- x[x != ""]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = " | ")
}

standardise_gwas_snps <- function(
    gwas_results_file,
    chromosome_col,
    position_col,
    snp_id_col,
    trait_col,
    chrom_map
) {
  snps <- read.csv(
    gwas_results_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(chromosome_col, position_col, snp_id_col, trait_col)
  missing_cols <- setdiff(required_cols, colnames(snps))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required GWAS columns are missing: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  snps <- snps %>%
    mutate(
      Chromosome = paste0("Pd", sprintf("%02d", as.numeric(.data[[chromosome_col]]))),
      Chromosome_ncbi = chrom_map[Chromosome],
      SNP_position = as.numeric(.data[[position_col]]),
      SNP_ID = .data[[snp_id_col]],
      Trait = .data[[trait_col]]
    ) %>%
    filter(
      !is.na(Chromosome_ncbi),
      !is.na(SNP_position),
      !is.na(SNP_ID),
      !is.na(Trait)
    )

  return(snps)
}

create_snp_windows <- function(snps, window_kb) {
  snp_gr <- GenomicRanges::GRanges(
    seqnames = snps$Chromosome_ncbi,
    ranges = IRanges::IRanges(
      start = snps$SNP_position,
      width = 1
    ),
    SNP_ID = snps$SNP_ID,
    Trait = snps$Trait,
    Chromosome = snps$Chromosome
  )

  snp_windows <- resize(
    snp_gr,
    width = 2 * window_kb * 1000 + 1,
    fix = "center"
  )

  return(
    list(
      snp_gr = snp_gr,
      snp_windows = snp_windows
    )
  )
}

extract_genes_from_gff <- function(gff_file, snp_gr, snp_windows) {
  gff <- rtracklayer::import(gff_file)

  genes <- gff[gff$type == "gene"]
  mrna <- gff[gff$type == "mRNA"]

  overlaps <- GenomicRanges::findOverlaps(genes, snp_windows)

  if (length(overlaps) == 0) {
    stop("No genes were found within the selected SNP windows.", call. = FALSE)
  }

  genes_near <- genes[S4Vectors::queryHits(overlaps)]
  snp_hits <- snp_gr[S4Vectors::subjectHits(overlaps)]

  genes_info <- data.frame(
    gene_id = genes_near$gene,
    seqnames = as.character(seqnames(genes_near)),
    start = start(genes_near),
    end = end(genes_near),
    SNP_ID = S4Vectors::mcols(snp_hits)$SNP_ID,
    Trait = S4Vectors::mcols(snp_hits)$Trait,
    SNP_position = start(snp_hits),
    stringsAsFactors = FALSE
  )

  mrna_df <- as.data.frame(S4Vectors::mcols(mrna))

  if (all(c("gene", "product") %in% colnames(mrna_df))) {
    mrna_df <- mrna_df[, c("gene", "product")]
    colnames(mrna_df) <- c("gene_id", "annotation")

    genes_info <- merge(
      genes_info,
      mrna_df,
      by = "gene_id",
      all.x = TRUE
    )
  } else {
    warning(
      "The GFF mRNA attributes 'gene' and/or 'product' were not found. ",
      "GFF functional annotation will be skipped.",
      call. = FALSE
    )
    genes_info$annotation <- NA_character_
  }

  genes_info <- unique(genes_info)

  return(genes_info)
}

download_pdulcis_annotations <- function(entrez_ids) {
  mart_plants <- biomaRt::useMart(
    "plants_mart",
    host = "https://plants.ensembl.org"
  )

  mart_almond <- biomaRt::useDataset(
    "pdulcis_eg_gene",
    mart = mart_plants
  )

  cat("\nDownloading Prunus dulcis base annotations...\n")
  df_base <- biomaRt::getBM(
    filters = "entrezgene_id",
    attributes = c(
      "entrezgene_id",
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "description"
    ),
    values = entrez_ids,
    mart = mart_almond
  )

  cat("Downloading WikiGene annotations...\n")
  df_wiki <- biomaRt::getBM(
    filters = "entrezgene_id",
    attributes = c(
      "entrezgene_id",
      "wikigene_name",
      "wikigene_description"
    ),
    values = entrez_ids,
    mart = mart_almond
  )

  cat("Downloading GO, KEGG and UniProt annotations...\n")
  df_go_kegg <- biomaRt::getBM(
    filters = "entrezgene_id",
    attributes = c(
      "entrezgene_id",
      "go_id",
      "name_1006",
      "definition_1006",
      "kegg_enzyme",
      "uniprot_gn_id"
    ),
    values = entrez_ids,
    mart = mart_almond
  )

  cat("Downloading InterPro annotations...\n")
  df_interpro <- biomaRt::getBM(
    filters = "entrezgene_id",
    attributes = c(
      "entrezgene_id",
      "interpro",
      "interpro_short_description",
      "interpro_description"
    ),
    values = entrez_ids,
    mart = mart_almond
  )

  cat("Downloading Arabidopsis homologs...\n")
  df_homologs <- biomaRt::getBM(
    filters = "entrezgene_id",
    attributes = c(
      "ensembl_gene_id",
      "athaliana_eg_homolog_ensembl_gene",
      "athaliana_eg_homolog_associated_gene_name",
      "athaliana_eg_homolog_ensembl_peptide",
      "athaliana_eg_homolog_subtype",
      "athaliana_eg_homolog_orthology_type",
      "athaliana_eg_homolog_perc_id",
      "athaliana_eg_homolog_perc_id_r1",
      "athaliana_eg_homolog_wga_coverage",
      "athaliana_eg_homolog_orthology_confidence"
    ),
    values = entrez_ids,
    mart = mart_almond
  )

  df_ensembl <- merge(df_base, df_wiki, by = "entrezgene_id", all = TRUE)
  df_ensembl <- merge(df_ensembl, df_go_kegg, by = "entrezgene_id", all = TRUE)
  df_ensembl <- merge(df_ensembl, df_interpro, by = "entrezgene_id", all = TRUE)
  df_ensembl <- merge(df_ensembl, df_homologs, by = "ensembl_gene_id", all = TRUE)

  return(df_ensembl)
}

collapse_candidate_genes <- function(candidate_df, chrom_map) {
  candidate_df[candidate_df == ""] <- NA

  collapsed <- candidate_df %>%
    group_by(
      gene_id,
      seqnames,
      start,
      end,
      SNP_ID,
      Trait,
      SNP_position,
      athaliana_eg_homolog_ensembl_gene
    ) %>%
    summarise(
      across(everything(), collapse_values),
      .groups = "drop"
    )

  collapsed[collapsed == ""] <- NA

  reverse_chrom_map <- names(chrom_map)
  names(reverse_chrom_map) <- chrom_map

  collapsed <- collapsed %>%
    mutate(chromosome = reverse_chrom_map[seqnames]) %>%
    relocate(chromosome, .after = seqnames)

  return(collapsed)
}

download_arabidopsis_annotations <- function(candidate_df) {
  high_confidence_genes <- candidate_df %>%
    filter(
      !is.na(athaliana_eg_homolog_orthology_confidence),
      athaliana_eg_homolog_orthology_confidence == 1,
      !is.na(athaliana_eg_homolog_ensembl_gene),
      athaliana_eg_homolog_ensembl_gene != ""
    )

  arabidopsis_ids <- unique(high_confidence_genes$athaliana_eg_homolog_ensembl_gene)

  cat("\nHigh-confidence Arabidopsis homologs detected:", length(arabidopsis_ids), "\n")

  if (length(arabidopsis_ids) == 0) {
    warning("No high-confidence Arabidopsis homologs were detected.", call. = FALSE)
    return(data.frame())
  }

  mart_plants <- biomaRt::useMart(
    "plants_mart",
    host = "https://plants.ensembl.org"
  )

  mart_arabidopsis <- biomaRt::useDataset(
    "athaliana_eg_gene",
    mart = mart_plants
  )

  cat("Downloading Arabidopsis GO annotations...\n")
  annot_go <- biomaRt::getBM(
    filters = "ensembl_gene_id",
    attributes = c(
      "ensembl_gene_id",
      "go_id",
      "name_1006",
      "definition_1006"
    ),
    values = arabidopsis_ids,
    mart = mart_arabidopsis
  )

  cat("Downloading Arabidopsis InterPro and WikiGene annotations...\n")
  annot_inter_wiki <- biomaRt::getBM(
    filters = "ensembl_gene_id",
    attributes = c(
      "ensembl_gene_id",
      "interpro",
      "interpro_description",
      "wikigene_description"
    ),
    values = arabidopsis_ids,
    mart = mart_arabidopsis
  )

  annot_at <- merge(
    annot_go,
    annot_inter_wiki,
    by = "ensembl_gene_id",
    all = TRUE
  )

  colnames(annot_at) <- c(
    "athaliana_eg_homolog_ensembl_gene",
    "TAIR_go_id",
    "TAIR_go_name",
    "TAIR_go_definition",
    "TAIR_interpro",
    "TAIR_interpro_description",
    "TAIR_wikigene_description"
  )

  annot_at[annot_at == ""] <- NA

  annot_at_collapsed <- annot_at %>%
    group_by(athaliana_eg_homolog_ensembl_gene) %>%
    summarise(
      across(everything(), collapse_values),
      .groups = "drop"
    )

  annot_at_collapsed[annot_at_collapsed == ""] <- NA

  return(annot_at_collapsed)
}


############################################################
# 3. Main execution
############################################################

snps <- standardise_gwas_snps(
  gwas_results_file = gwas_results_file,
  chromosome_col = chromosome_col,
  position_col = position_col,
  snp_id_col = snp_id_col,
  trait_col = trait_col,
  chrom_map = chrom_map
)

snp_ranges <- create_snp_windows(
  snps = snps,
  window_kb = window_kb
)

genes_info <- extract_genes_from_gff(
  gff_file = gff_file,
  snp_gr = snp_ranges$snp_gr,
  snp_windows = snp_ranges$snp_windows
)

genes_info$entrezgene_id <- sub("^LOC", "", genes_info$gene_id)
entrez_ids <- unique(genes_info$entrezgene_id)

pdulcis_annotations <- download_pdulcis_annotations(entrez_ids)

candidate_genes_raw <- merge(
  genes_info,
  pdulcis_annotations,
  by = "entrezgene_id",
  all.x = TRUE
)

readr::write_tsv(
  candidate_genes_raw,
  file.path(output_dir, "GWAS_candidate_genes_raw.tsv"),
  na = ""
)

candidate_genes_collapsed <- collapse_candidate_genes(
  candidate_df = candidate_genes_raw,
  chrom_map = chrom_map
)

readr::write_tsv(
  candidate_genes_collapsed,
  file.path(output_dir, "GWAS_candidate_genes_collapsed.tsv"),
  na = ""
)

arabidopsis_annotations <- download_arabidopsis_annotations(candidate_genes_raw)

if (nrow(arabidopsis_annotations) > 0) {
  readr::write_tsv(
    arabidopsis_annotations,
    file.path(output_dir, "Arabidopsis_annotations_collapsed.tsv"),
    na = ""
  )

  final_annotated <- candidate_genes_collapsed %>%
    left_join(
      arabidopsis_annotations,
      by = "athaliana_eg_homolog_ensembl_gene"
    ) %>%
    group_by(gene_id, seqnames, start, end, SNP_ID, Trait, SNP_position) %>%
    summarise(
      across(everything(), collapse_values),
      .groups = "drop"
    )
} else {
  final_annotated <- candidate_genes_collapsed
}

final_annotated[final_annotated == ""] <- NA

readr::write_tsv(
  final_annotated,
  file.path(output_dir, "GWAS_candidate_genes_final_annotated.tsv"),
  na = ""
)

cat("\nAnnotation completed successfully.\n")
cat("Final output:", file.path(output_dir, "GWAS_candidate_genes_final_annotated.tsv"), "\n")
cat("Rows:", nrow(final_annotated), "\n")

############################################################
# End of script
############################################################
