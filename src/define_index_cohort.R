## Define index cohort for a disease family
##
## For a given disease family, extract all patients who have at least one
## claim matching the index ICD-10 codes. Determine index date (first
## diagnosis), baseline/follow-up windows, and demographics.
##
## Usage:
##   Rscript src/define_index_cohort.R \
##     <disease_family> <families_tsv> <tekiyo_csv> \
##     <coo_path> <col_ids_path> <row_ids_path> \
##     <output_cohort_tsv>

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7)

disease_family <- args[1]
families_path  <- args[2]
tekiyo_path    <- args[3]
coo_path       <- args[4]
col_ids_path   <- args[5]
row_ids_path   <- args[6]
out_path       <- args[7]

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

cat(sprintf("=== define_index_cohort: %s ===\n", disease_family))

# --- Load disease family definition ---
families <- fread(families_path, sep = "\t")
fam_row <- families[disease_family == disease_family]
stopifnot(nrow(fam_row) == 1)

index_codes <- trimws(unlist(strsplit(fam_row$index_icd10, ",")))
cat(sprintf("Index ICD-10 codes: %s\n", paste(index_codes, collapse = ", ")))

# --- Load column (disease code) mapping ---
# Format: code,id,name (CSV with potential commas in name)
col_lines <- readLines(col_ids_path, encoding = "UTF-8")
col_lines <- sub("^﻿", "", col_lines)
col_split <- function(s) {
  pos <- gregexpr(",", s, fixed = TRUE)[[1]]
  if (identical(pos, -1L) || length(pos) < 2L) return(NULL)
  code <- substr(s, 1, pos[1] - 1)
  id   <- substr(s, pos[1] + 1, pos[2] - 1)
  c(code = trimws(code), id = as.integer(trimws(id)))
}
col_parsed <- do.call(rbind, lapply(col_lines, col_split))
col_map <- data.table(
  code = col_parsed[, "code"],
  col_id = as.integer(col_parsed[, "id"])
)

# --- Identify columns matching index disease codes (prefix match) ---
is_index_col <- function(code, index_codes) {
  any(vapply(index_codes, function(ic) startsWith(code, ic), logical(1)))
}
col_map[, is_index := vapply(code, is_index_col, logical(1),
                              index_codes = index_codes)]
index_col_ids <- col_map[is_index == TRUE, col_id]
cat(sprintf("Columns matching index codes: %d\n", length(index_col_ids)))

# --- Load row ID mapping ---
row_ids <- fread(row_ids_path, header = FALSE)
# Assume format: patient_id (one per line or two-column with row_number)
if (ncol(row_ids) == 1) {
  row_ids[, row_id := .I]
  setnames(row_ids, "V1", "kojin_id")
} else {
  setnames(row_ids, c("kojin_id", "row_id"))
}

# --- Scan COO matrix for patients with index disease codes ---
# COO format: row_index col_index (1-based or 0-based)
# For large files, stream and filter
cat("Scanning COO for patients with index disease ...\n")

# Adjust col_ids for 0-based if needed
index_col_set <- as.character(index_col_ids)

# Read COO in chunks to find rows (patients) that have index disease columns
# TODO: For production, use Julia streaming for large COO.zst files
# For now, use data.table::fread if the file is text-based
# If zst-compressed, decompress on the fly:
if (grepl("\\.zst$", coo_path)) {
  con <- pipe(sprintf("zstd -dc %s", shQuote(coo_path)), open = "r")
} else {
  con <- file(coo_path, open = "r")
}

index_rows <- integer(0)
chunk_size <- 10000000L
repeat {
  chunk <- readLines(con, n = chunk_size)
  if (length(chunk) == 0) break
  parts <- strsplit(trimws(chunk), "\\s+")
  rows <- as.integer(vapply(parts, `[`, character(1), 1))
  cols <- as.integer(vapply(parts, `[`, character(1), 2))
  hit <- rows[cols %in% index_col_ids]
  index_rows <- unique(c(index_rows, hit))
  cat(sprintf("  scanned %d lines, %d index patients so far\n",
              length(chunk), length(index_rows)))
}
close(con)

cat(sprintf("Total patients with index disease: %d\n", length(index_rows)))

# --- Load demographics (tekiyo) ---
tek <- fread(tekiyo_path)

# --- Build cohort table ---
cohort_rows <- row_ids[row_id %in% index_rows]
cohort <- merge(cohort_rows, tek, by = "kojin_id", all.x = TRUE)

# Add placeholder columns for index_date / baseline / follow-up
# TODO: Implement proper index date extraction from claims dates
cohort[, index_date := NA_character_]
cohort[, baseline_start := observable_start_ym]
cohort[, baseline_end := observable_end_ym]
cohort[, followup_start := NA_character_]
cohort[, followup_end := observable_end_ym]

# Demographics
cohort[, age_at_baseline := NA_real_]  # TODO: compute from birth_ym
cohort[, sex := sex_code]

# Death flag
if ("death_flag" %in% names(cohort)) {
  cohort[, event_death := as.integer(death_flag == 1)]
} else {
  cohort[, event_death := 0L]
}

# Time to event (years)
# TODO: Implement proper time calculation from dates
cohort[, time_years := NA_real_]

cat(sprintf("Cohort size: %d patients\n", nrow(cohort)))
cat(sprintf("  deaths: %d (%.1f%%)\n",
            sum(cohort$event_death, na.rm = TRUE),
            100 * mean(cohort$event_death, na.rm = TRUE)))

# --- Output ---
fwrite(cohort, out_path, sep = "\t")
cat(sprintf("Written: %s\n", out_path))
