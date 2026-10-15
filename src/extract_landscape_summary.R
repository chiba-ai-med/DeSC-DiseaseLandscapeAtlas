## Extract per-patient landscape summary
##
## Reads Landscaper output (Allstates, Basin, SubGraph, Energy) and
## Scores.mm to assign each patient a basin ID and energy value.
## Merges with cohort demographics for downstream Cox analysis.
##
## Usage:
##   Rscript src/extract_landscape_summary.R \
##     <cohort_tsv> <scores_mm> <allstates_tsv> <basin_tsv> \
##     <subgraph_tsv> <energy_tsv> <output_tsv>

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7)

cohort_path    <- args[1]
scores_mm_path <- args[2]
allstates_path <- args[3]
basin_path     <- args[4]
subgraph_path  <- args[5]
energy_path    <- args[6]
out_path       <- args[7]

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

cat("=== extract_landscape_summary ===\n")

# --- Load Allstates (pattern definitions: 2^k rows, k columns of -1/+1) ---
allstates <- as.matrix(fread(allstates_path, header = FALSE))
n_patterns <- nrow(allstates)
k <- ncol(allstates)
cat(sprintf("Allstates: %d patterns, k=%d\n", n_patterns, k))

# --- Load Basin (basin representative patterns) ---
basin <- fread(basin_path, header = FALSE)
cat(sprintf("Basins: %d\n", nrow(basin)))

# --- Load SubGraph (pattern -> basin assignment) ---
subgraph <- fread(subgraph_path, header = FALSE)
setnames(subgraph, c("pattern_id", "basin_id"))

# --- Load Energy (per-pattern energy) ---
energy <- fread(energy_path, header = FALSE)
setnames(energy, "energy")
energy[, pattern_id := .I]

# --- Read Scores.mm and assign each row to a pattern ---
# Matrix Market: sparse binary, non-zero = +1, zero = -1
cat("Reading Scores.mm ...\n")
scores <- readMM(scores_mm_path)
n_patients <- nrow(scores)
cat(sprintf("Patients: %d\n", n_patients))

# Convert each row to pattern ID via bit encoding
# Pattern: for each row, create +-1 vector from sparse, then match to allstates
assign_pattern <- function(row_sparse, allstates_mat) {
  k <- ncol(allstates_mat)
  x <- rep(-1L, k)
  nz <- which(row_sparse != 0)
  x[nz] <- 1L
  # Match against allstates
  for (pid in seq_len(nrow(allstates_mat))) {
    if (all(allstates_mat[pid, ] == x)) return(pid)
  }
  return(NA_integer_)
}

cat("Assigning patterns to patients ...\n")
pattern_ids <- integer(n_patients)
for (i in seq_len(n_patients)) {
  row_vec <- scores[i, , drop = FALSE]
  pattern_ids[i] <- assign_pattern(row_vec, allstates)
  if (i %% 10000 == 0) cat(sprintf("  %d / %d\n", i, n_patients))
}

# --- Build person-level table ---
person <- data.table(
  row_id = seq_len(n_patients),
  pattern_id = pattern_ids
)

# Merge basin assignment
person <- merge(person, subgraph, by = "pattern_id", all.x = TRUE)

# Merge energy
person <- merge(person, energy, by = "pattern_id", all.x = TRUE)

# --- Merge with cohort demographics ---
cohort <- fread(cohort_path)
# Assume cohort has row_id column for joining
if ("row_id" %in% names(cohort)) {
  # Remap: cohort row_ids are original, person row_ids are 1..n after filtering
  # TODO: proper remapping using row_ids.txt
  person <- cbind(person, cohort[seq_len(n_patients),
                                  .(kojin_id, age_at_baseline, sex,
                                    event_death, time_years)])
}

# --- Output ---
fwrite(person, out_path, sep = "\t")
cat(sprintf("Written: %s (%d rows)\n", out_path, nrow(person)))

# Summary stats
cat(sprintf("  Basins found: %s\n",
            paste(sort(unique(person$basin_id)), collapse = ", ")))
cat(sprintf("  Energy range: [%.3f, %.3f]\n",
            min(person$energy, na.rm = TRUE),
            max(person$energy, na.rm = TRUE)))
