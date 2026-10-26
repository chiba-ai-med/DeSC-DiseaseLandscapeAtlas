## Build state_id -> basin_id mapping table for one disease
##
## Reads StatusNetwork.tsv to identify truly isolated vertices (in-degree=0
## AND out-degree=0), then maps SubGraph.tsv rows to non-isolated state IDs.
## States with SubGraph=128 rows (full format, e.g. dementia/diabetes) are
## handled by direct index mapping.
##
## Output: state_basin_map.tsv with columns:
##   state_id  (1-indexed, 1..128)
##   basin_id  (integer; 0 = non-basin / isolated)
##   energy    (from E.tsv)
##   isolated  (1 if both in-degree and out-degree are 0, else 0)
##
## Usage:
##   Rscript src/build_state_basin_map.R <landscaper_dir> <output_tsv>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
ls_dir   <- args[1]
out_path <- args[2]

sn_path <- file.path(ls_dir, "StatusNetwork.tsv")
sg_path <- file.path(ls_dir, "SubGraph.tsv")
e_path  <- file.path(ls_dir, "E.tsv")

sn <- as.matrix(fread(sn_path, header = FALSE))
n_states <- nrow(sn)
stopifnot(ncol(sn) == n_states)

subgraph_raw <- fread(sg_path, header = FALSE)[[1]]
energies <- fread(e_path, header = FALSE)[[1]]
stopifnot(length(energies) == n_states)

out_deg <- rowSums(sn)
in_deg  <- colSums(sn)
isolated <- as.integer(out_deg == 0L & in_deg == 0L)
n_isolated <- sum(isolated)
non_isolated_ids <- which(isolated == 0L)

n_sg <- length(subgraph_raw)
if (n_sg == n_states) {
  basin_ids <- subgraph_raw
} else if (n_sg == length(non_isolated_ids)) {
  basin_ids <- rep(0L, n_states)
  basin_ids[non_isolated_ids] <- subgraph_raw
} else {
  stop(sprintf("SubGraph has %d rows, expected %d (full) or %d (non-isolated)",
               n_sg, n_states, length(non_isolated_ids)))
}

result <- data.table(
  state_id = seq_len(n_states),
  basin_id = basin_ids,
  energy   = energies,
  isolated = isolated
)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
fwrite(result, out_path, sep = "\t")
cat(sprintf("Written: %s (%d states, %d isolated, %d basins, SubGraph=%d rows)\n",
            out_path, n_states, n_isolated,
            length(unique(basin_ids[basin_ids > 0])), n_sg))
