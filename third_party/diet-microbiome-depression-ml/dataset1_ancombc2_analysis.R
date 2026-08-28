#!/usr/bin/env Rscript

# Dataset 1 post-hoc differential-abundance analysis with ANCOM-BC2
#
# The script fits:
#   1. an unadjusted model: depression group;
#   2. an age-adjusted model: depression group + age.
#
# Dataset 1 contains only male participants, so gender is not included.
#
# Required input files in the prepared-data directory:
#   metadata_ancombc_matched.txt
#   counts_ancombc_asv.txt
#   taxonomy_ancombc_asv.txt
#
# Usage:
# Rscript src/dataset1_ancombc2_analysis.R \
#   prepared_data/dataset1_ancombc2 \
#   results/dataset1/ancombc2

required_packages <- c("ANCOMBC")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install the required Bioconductor package before running this script."
  )
}

library(ANCOMBC)


# =============================================================================
# Paths
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

prepared_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path("prepared_data", "dataset1_ancombc2")
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("results", "dataset1", "ancombc2")
}

metadata_file <- file.path(
  prepared_dir,
  "metadata_ancombc_matched.txt"
)

counts_file <- file.path(
  prepared_dir,
  "counts_ancombc_asv.txt"
)

taxonomy_file <- file.path(
  prepared_dir,
  "taxonomy_ancombc_asv.txt"
)

required_files <- c(
  metadata_file,
  counts_file,
  taxonomy_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# Analysis settings
# =============================================================================

ANCOMBC_SEED <- 123
ALPHA <- 0.05
N_CORES <- 2

# Final SVM-prioritized taxa interpreted in the thesis.
# The ambiguous feature labelled only as "uncultured" is excluded.
PRIORITIZED_TAXA <- c(
  "Dielma",
  "Lachnospiraceae_ND3007_group",
  "Lachnospiraceae_FCS020_group",
  "Coprobacter",
  "Intestinibacter",
  "Oscillibacter",
  "Rhizobiales_unclassified",
  "Catabacter",
  "Anaerofustis",
  "Eisenbergiella"
)


# =============================================================================
# Helper functions
# =============================================================================

require_columns <- function(data, required, data_name) {
  missing <- setdiff(required, colnames(data))

  if (length(missing) > 0) {
    stop(
      data_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}


normalise_name <- function(value) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  value <- gsub("\\\\_", "_", value)
  value <- trimws(tolower(value))
  value <- gsub("[[:space:]-]+", "_", value)
  value <- gsub("_+", "_", value)
  value
}


parse_logical <- function(value) {
  if (length(value) == 0 || is.na(value)) {
    return(FALSE)
  }

  if (is.logical(value)) {
    return(isTRUE(value))
  }

  if (is.numeric(value)) {
    return(value != 0)
  }

  tolower(trimws(as.character(value))) %in% c(
    "true", "t", "1", "yes", "y"
  )
}


direction_from_lfc <- function(value) {
  if (is.na(value)) {
    return(NA_character_)
  }

  if (value > 0) {
    return("higher_in_depression")
  }

  if (value < 0) {
    return("higher_in_control")
  }

  "no_direction"
}


prepare_zero_table <- function(zero_ind) {
  if (is.null(zero_ind)) {
    return(data.frame())
  }

  zero_table <- as.data.frame(
    zero_ind,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"taxon" %in% colnames(zero_table)) {
    zero_table$taxon <- rownames(zero_table)
  }

  zero_table$taxon <- trimws(as.character(zero_table$taxon))
  rownames(zero_table) <- NULL

  zero_table[, c(
    "taxon",
    setdiff(colnames(zero_table), "taxon")
  ), drop = FALSE]
}


structural_zero_direction <- function(zero_table, taxon_id) {
  if (nrow(zero_table) == 0) {
    return(NA_character_)
  }

  row <- zero_table[
    zero_table$taxon == as.character(taxon_id),
    ,
    drop = FALSE
  ]

  if (nrow(row) == 0) {
    return(NA_character_)
  }

  depression_columns <- grep(
    "depression",
    colnames(zero_table),
    value = TRUE,
    ignore.case = TRUE
  )

  control_columns <- grep(
    "control",
    colnames(zero_table),
    value = TRUE,
    ignore.case = TRUE
  )

  depression_zero <- if (length(depression_columns) > 0) {
    parse_logical(row[[depression_columns[[1]]]][[1]])
  } else {
    FALSE
  }

  control_zero <- if (length(control_columns) > 0) {
    parse_logical(row[[control_columns[[1]]]][[1]])
  } else {
    FALSE
  }

  if (depression_zero && !control_zero) {
    return("higher_in_control")
  }

  if (control_zero && !depression_zero) {
    return("higher_in_depression")
  }

  NA_character_
}


standardise_result <- function(result) {
  result <- as.data.frame(
    result,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"taxon" %in% colnames(result)) {
    result$taxon <- rownames(result)
  }

  result$taxon <- trimws(as.character(result$taxon))
  rownames(result) <- NULL

  result[, c(
    "taxon",
    setdiff(colnames(result), "taxon")
  ), drop = FALSE]
}


add_taxonomy <- function(result, taxonomy_lookup) {
  result <- standardise_result(result)

  taxonomy_rows <- taxonomy_lookup[
    match(result$taxon, taxonomy_lookup$taxon),
    setdiff(colnames(taxonomy_lookup), "taxon"),
    drop = FALSE
  ]

  rownames(taxonomy_rows) <- NULL
  cbind(result, taxonomy_rows)
}


save_group_effect <- function(result_with_taxonomy, output_file) {
  effect_columns <- grep(
    "groupdepression",
    colnames(result_with_taxonomy),
    value = TRUE,
    ignore.case = FALSE
  )

  taxonomy_columns <- intersect(
    c(
      "taxon_name",
      "taxlevel",
      "daughterlevels",
      "total",
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    ),
    colnames(result_with_taxonomy)
  )

  columns_to_keep <- unique(c(
    "taxon",
    taxonomy_columns,
    effect_columns
  ))

  group_effect <- result_with_taxonomy[
    ,
    columns_to_keep,
    drop = FALSE
  ]

  if ("lfc_groupdepression" %in% colnames(group_effect)) {
    group_effect$direction <- vapply(
      group_effect$lfc_groupdepression,
      direction_from_lfc,
      character(1)
    )
  }

  write.csv(
    group_effect,
    output_file,
    row.names = FALSE
  )

  group_effect
}


run_ancombc2_model <- function(
  count_matrix,
  metadata,
  fixed_formula
) {
  set.seed(ANCOMBC_SEED)

  ancombc2(
    data = count_matrix,
    meta_data = metadata,
    fix_formula = fixed_formula,
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "group",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = ALPHA,
    n_cl = N_CORES,
    verbose = TRUE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    iter_control = list(
      tol = 1e-5,
      max_iter = 20,
      verbose = FALSE
    ),
    em_control = list(
      tol = 1e-5,
      max_iter = 100
    ),
    lme_control = NULL,
    mdfdr_control = list(
      fwer_ctrl_method = "BH",
      B = 100
    ),
    trend_control = NULL
  )
}


extract_effect_value <- function(row, column_name) {
  if (
    nrow(row) == 0 ||
    !column_name %in% colnames(row)
  ) {
    return(NA)
  }

  row[[column_name]][[1]]
}


extract_prioritized_result <- function(
  result,
  zero_table,
  taxon_id
) {
  result_row <- result[
    result$taxon == as.character(taxon_id),
    ,
    drop = FALSE
  ]

  zero_direction <- structural_zero_direction(
    zero_table,
    taxon_id
  )

  if (nrow(result_row) == 0) {
    return(list(
      lfc = NA_real_,
      se = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      q_value = NA_real_,
      differentially_abundant = NA,
      passed_sensitivity = NA,
      robust_significant = NA,
      direction = zero_direction,
      structural_zero = !is.na(zero_direction)
    ))
  }

  lfc <- extract_effect_value(
    result_row,
    "lfc_groupdepression"
  )

  list(
    lfc = lfc,
    se = extract_effect_value(
      result_row,
      "se_groupdepression"
    ),
    statistic = extract_effect_value(
      result_row,
      "W_groupdepression"
    ),
    p_value = extract_effect_value(
      result_row,
      "p_groupdepression"
    ),
    q_value = extract_effect_value(
      result_row,
      "q_groupdepression"
    ),
    differentially_abundant = extract_effect_value(
      result_row,
      "diff_groupdepression"
    ),
    passed_sensitivity = extract_effect_value(
      result_row,
      "passed_ss_groupdepression"
    ),
    robust_significant = extract_effect_value(
      result_row,
      "diff_robust_groupdepression"
    ),
    direction = direction_from_lfc(lfc),
    structural_zero = FALSE
  )
}


prefix_result <- function(values, prefix) {
  names(values) <- paste0(prefix, names(values))
  values
}


# =============================================================================
# Load and validate input data
# =============================================================================

metadata <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- read.delim(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxonomy <- read.delim(
  taxonomy_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

require_columns(
  metadata,
  c("sample_id", "group", "age"),
  "Metadata"
)

require_columns(
  counts,
  "ASV_ID",
  "Count table"
)

require_columns(
  taxonomy,
  "ASV_ID",
  "Taxonomy table"
)


# =============================================================================
# Prepare metadata
# =============================================================================

metadata$sample_id <- trimws(
  as.character(metadata$sample_id)
)

metadata$group <- trimws(
  tolower(as.character(metadata$group))
)

metadata$age <- suppressWarnings(
  as.numeric(metadata$age)
)

metadata$group <- factor(
  metadata$group,
  levels = c("control", "depression")
)

if (anyNA(metadata$sample_id) || any(metadata$sample_id == "")) {
  stop("Metadata contains missing or empty sample IDs.")
}

if (anyNA(metadata$group)) {
  stop(
    "Metadata group values must be exactly 'control' or 'depression'."
  )
}

if (anyNA(metadata$age)) {
  stop("Metadata contains missing or non-numeric age values.")
}

if (anyDuplicated(metadata$sample_id)) {
  stop("Metadata contains duplicated sample IDs.")
}

rownames(metadata) <- metadata$sample_id


# =============================================================================
# Prepare count matrix
# =============================================================================

counts$ASV_ID <- trimws(
  as.character(counts$ASV_ID)
)

if (anyNA(counts$ASV_ID) || any(counts$ASV_ID == "")) {
  stop("Count table contains missing or empty ASV IDs.")
}

if (anyDuplicated(counts$ASV_ID)) {
  stop("Count table contains duplicated ASV IDs.")
}

count_columns <- setdiff(
  colnames(counts),
  "ASV_ID"
)

count_data <- counts[
  ,
  count_columns,
  drop = FALSE
]

count_data[] <- lapply(
  count_data,
  function(value) {
    suppressWarnings(as.numeric(as.character(value)))
  }
)

if (anyNA(as.matrix(count_data))) {
  stop("Count table contains missing or non-numeric count values.")
}

if (any(as.matrix(count_data) < 0)) {
  stop("Count table contains negative values.")
}

count_matrix <- as.matrix(count_data)
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- counts$ASV_ID


# =============================================================================
# Align metadata and count data
# =============================================================================

shared_samples <- rownames(metadata)[
  rownames(metadata) %in% colnames(count_matrix)
]

if (length(shared_samples) == 0) {
  stop(
    "No matching sample IDs were found between metadata and count data."
  )
}

metadata_only <- setdiff(
  rownames(metadata),
  colnames(count_matrix)
)

counts_only <- setdiff(
  colnames(count_matrix),
  rownames(metadata)
)

write.csv(
  data.frame(sample_id = metadata_only),
  file.path(
    output_dir,
    "sample_ids_only_in_metadata.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(sample_id = counts_only),
  file.path(
    output_dir,
    "sample_ids_only_in_counts.csv"
  ),
  row.names = FALSE
)

metadata <- metadata[
  shared_samples,
  ,
  drop = FALSE
]

count_matrix <- count_matrix[
  ,
  shared_samples,
  drop = FALSE
]

if (!identical(
  rownames(metadata),
  colnames(count_matrix)
)) {
  stop(
    "Metadata rows and count-matrix columns are not aligned."
  )
}

if (length(unique(metadata$group)) < 2) {
  stop("Both depression groups must be present in the aligned data.")
}

write.csv(
  metadata,
  file.path(
    output_dir,
    "metadata_used_in_ancombc2.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    sample_id = colnames(count_matrix),
    library_size = colSums(count_matrix)
  ),
  file.path(
    output_dir,
    "library_sizes_used_in_ancombc2.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Prepare taxonomy lookup
# =============================================================================

taxonomy$ASV_ID <- trimws(
  as.character(taxonomy$ASV_ID)
)

if (anyDuplicated(taxonomy$ASV_ID)) {
  stop("Taxonomy table contains duplicated ASV IDs.")
}

taxonomy$taxon <- taxonomy$ASV_ID

taxonomy_columns <- intersect(
  c(
    "taxon_name",
    "taxlevel",
    "daughterlevels",
    "total",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  colnames(taxonomy)
)

taxonomy_lookup <- taxonomy[
  ,
  c("taxon", taxonomy_columns),
  drop = FALSE
]


# =============================================================================
# Fit ANCOM-BC2 models
# =============================================================================

cat(
  "Samples used:",
  nrow(metadata),
  "\n"
)

cat(
  "Taxa used:",
  nrow(count_matrix),
  "\n"
)

cat("Group counts:\n")
print(table(metadata$group))

cat("Age summary:\n")
print(summary(metadata$age))


# Unadjusted model
unadjusted_output <- run_ancombc2_model(
  count_matrix,
  metadata,
  fixed_formula = "group"
)

unadjusted_result <- standardise_result(
  unadjusted_output$res
)

unadjusted_with_taxonomy <- add_taxonomy(
  unadjusted_result,
  taxonomy_lookup
)

unadjusted_zero_table <- prepare_zero_table(
  unadjusted_output$zero_ind
)

write.csv(
  unadjusted_result,
  file.path(
    output_dir,
    "dataset1_ancombc2_unadjusted_full.csv"
  ),
  row.names = FALSE
)

write.csv(
  unadjusted_with_taxonomy,
  file.path(
    output_dir,
    "dataset1_ancombc2_unadjusted_with_taxonomy.csv"
  ),
  row.names = FALSE
)

unadjusted_group_effect <- save_group_effect(
  unadjusted_with_taxonomy,
  file.path(
    output_dir,
    "dataset1_ancombc2_unadjusted_group_effect.csv"
  )
)

if (nrow(unadjusted_zero_table) > 0) {
  write.csv(
    unadjusted_zero_table,
    file.path(
      output_dir,
      "dataset1_ancombc2_unadjusted_structural_zeros.csv"
    ),
    row.names = FALSE
  )
}


# Age-adjusted model
adjusted_output <- run_ancombc2_model(
  count_matrix,
  metadata,
  fixed_formula = "group + age"
)

adjusted_result <- standardise_result(
  adjusted_output$res
)

adjusted_with_taxonomy <- add_taxonomy(
  adjusted_result,
  taxonomy_lookup
)

adjusted_zero_table <- prepare_zero_table(
  adjusted_output$zero_ind
)

write.csv(
  adjusted_result,
  file.path(
    output_dir,
    "dataset1_ancombc2_age_adjusted_full.csv"
  ),
  row.names = FALSE
)

write.csv(
  adjusted_with_taxonomy,
  file.path(
    output_dir,
    "dataset1_ancombc2_age_adjusted_with_taxonomy.csv"
  ),
  row.names = FALSE
)

adjusted_group_effect <- save_group_effect(
  adjusted_with_taxonomy,
  file.path(
    output_dir,
    "dataset1_ancombc2_age_adjusted_group_effect.csv"
  )
)

if (nrow(adjusted_zero_table) > 0) {
  write.csv(
    adjusted_zero_table,
    file.path(
      output_dir,
      "dataset1_ancombc2_age_adjusted_structural_zeros.csv"
    ),
    row.names = FALSE
  )
}


# =============================================================================
# Extract results for the SVM-prioritized taxa
# =============================================================================

taxonomy_name <- if (
  "taxon_name" %in% colnames(taxonomy_lookup)
) {
  taxonomy_lookup$taxon_name
} else {
  rep("", nrow(taxonomy_lookup))
}

taxonomy_genus <- if (
  "Genus" %in% colnames(taxonomy_lookup)
) {
  taxonomy_lookup$Genus
} else {
  rep("", nrow(taxonomy_lookup))
}

taxonomy_lookup$taxon_name_normalised <- normalise_name(
  taxonomy_name
)

taxonomy_lookup$genus_normalised <- normalise_name(
  taxonomy_genus
)

prioritized_rows <- vector(
  "list",
  length(PRIORITIZED_TAXA)
)

for (index in seq_along(PRIORITIZED_TAXA)) {
  feature <- PRIORITIZED_TAXA[[index]]
  feature_normalised <- normalise_name(feature)

  matches <- taxonomy_lookup[
    taxonomy_lookup$taxon_name_normalised == feature_normalised |
      taxonomy_lookup$genus_normalised == feature_normalised,
    ,
    drop = FALSE
  ]

  if (nrow(matches) == 0) {
    warning(
      "No exact taxonomy match was found for: ",
      feature
    )

    prioritized_rows[[index]] <- data.frame(
      Feature = feature,
      taxon = NA_character_,
      taxonomy_match_count = 0,
      stringsAsFactors = FALSE
    )

    next
  }

  if (nrow(matches) > 1) {
    warning(
      nrow(matches),
      " exact taxonomy matches were found for ",
      feature,
      "; the first match is used."
    )
  }

  taxon_id <- matches$taxon[[1]]

  unadjusted_values <- extract_prioritized_result(
    unadjusted_group_effect,
    unadjusted_zero_table,
    taxon_id
  )

  adjusted_values <- extract_prioritized_result(
    adjusted_group_effect,
    adjusted_zero_table,
    taxon_id
  )

  row_values <- c(
    list(
      Feature = feature,
      taxon = taxon_id,
      taxonomy_match_count = nrow(matches)
    ),
    prefix_result(
      unadjusted_values,
      "unadjusted_"
    ),
    prefix_result(
      adjusted_values,
      "age_adjusted_"
    )
  )

  prioritized_rows[[index]] <- as.data.frame(
    row_values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

prioritized_summary <- do.call(
  rbind,
  prioritized_rows
)

write.csv(
  prioritized_summary,
  file.path(
    output_dir,
    "dataset1_ancombc2_svm_prioritized_taxa.csv"
  ),
  row.names = FALSE,
  na = ""
)


# =============================================================================
# Reproducibility information
# =============================================================================

session_file <- file.path(
  output_dir,
  "sessionInfo.txt"
)

sink(session_file)
print(sessionInfo())
sink()

cat("\nANCOM-BC2 analysis completed.\n")
cat("Results saved to:\n")
cat(output_dir, "\n")

