## Tag cloud of characteristic comorbidity codes per basin
##
## Binary COO format:
##   Header: nrows (Int32), ncols (Int32)
##   Body: (col_id, row_id) pairs as Int32
##
## Each row gets a basin_id via Scores.mm binarization -> Allstates -> SubGraph.
## Enrichment = (freq in basin) / (freq outside basin).
##
## Usage:
##   Rscript src/tagcloud_per_basin.R <disease> <base_dir> <out_dir> <icd10_en_tsv>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
disease    <- args[1]
base_dir   <- args[2]
out_dir    <- args[3]
icd10_path <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Tag cloud per basin: %s ===\n", disease))

# --- Build row -> basin_id mapping ---
ls_dir <- file.path(base_dir, "landscaper", "7")
subgraph <- fread(file.path(ls_dir, "SubGraph.tsv"), header = FALSE)[[1]]

allstates <- as.matrix(fread(file.path(ls_dir, "Allstates.tsv"), header = FALSE))
n_patterns <- nrow(allstates)

as_keys <- integer(n_patterns)
for (p in seq_len(n_patterns)) {
  as_keys[p] <- sum(2L^(which(allstates[p, ] > 0) - 1L))
}

cat("Reading Scores.mm ...\n")
pca_dir <- file.path(base_dir, "pca", "7")
S <- readMM(file.path(pca_dir, "Scores.mm"))
n_rows <- nrow(S)

S_bin <- as(S > 0, "TsparseMatrix")
ij_dt <- data.table(row = S_bin@i + 1L, col = S_bin@j + 1L)
ij_dt[, bit := 2L^(col - 1L)]
key_dt <- ij_dt[, .(key = sum(bit)), by = row]
rm(S, S_bin, ij_dt)

row_keys <- integer(n_rows)
row_keys[key_dt$row] <- key_dt$key
rm(key_dt)

key_to_pid <- setNames(seq_len(n_patterns), as.character(as_keys))
row_basin <- integer(n_rows)
ukeys <- unique(row_keys)
for (k in ukeys) {
  pid <- key_to_pid[as.character(k)]
  if (!is.na(pid)) row_basin[row_keys == k] <- subgraph[pid]
}
rm(row_keys, ukeys, allstates)

basins <- sort(unique(row_basin[row_basin > 0L]))
cat(sprintf("  %d rows, basins: %s\n", n_rows, paste(basins, collapse = ",")))

# --- Load col_ids ---
col_ids_dt <- fread(file.path(base_dir, "cohort", "col_ids.txt"),
                    header = FALSE, col.names = c("code", "col_id"))

# --- Load ICD-10 English names ---
icd10_en <- fread(icd10_path, sep = "\t", header = TRUE,
                  select = c("icd10_code", "english_name"))
setnames(icd10_en, c("code", "en_name"))
col_ids_dt <- merge(col_ids_dt, icd10_en, by = "code", all.x = TRUE)
col_ids_dt[is.na(en_name) | en_name == "", en_name := code]
# Wrap long names
col_ids_dt[, en_label := vapply(en_name, function(s) {
  paste(strwrap(s, width = 28), collapse = "\n")
}, character(1))]
cat(sprintf("  %d/%d codes have English names\n",
            sum(col_ids_dt$en_name != col_ids_dt$code), nrow(col_ids_dt)))

# --- Read binary COO ---
cat("Reading binary COO ...\n")
coo_path <- file.path(base_dir, "cohort", "feature_matrix.coo.zst")
con <- pipe(sprintf("zstd -d -c '%s'", coo_path), "rb")
header <- readBin(con, "integer", n = 2, size = 4, endian = "little")
cat(sprintf("  header: nrows=%d, ncols=%d\n", header[1], header[2]))

raw <- readBin(con, "integer", n = 2e9, size = 4, endian = "little")
close(con)

n_vals <- length(raw)
cat(sprintf("  %d pairs read\n", n_vals %/% 2))

coo_dt <- data.table(
  row_id = raw[seq.int(1L, n_vals, 2L)],
  col_id = raw[seq.int(2L, n_vals, 2L)]
)
rm(raw); gc()

coo_dt[, basin_id := row_basin[row_id]]
coo_dt <- coo_dt[basin_id > 0L]
coo_dt <- merge(coo_dt, col_ids_dt[, .(col_id, code, en_label)],
                by = "col_id", all.x = FALSE)
cat(sprintf("  %d entries with basin and code\n", nrow(coo_dt)))

# --- Compute enrichment: freq_in_basin / freq_outside_basin ---
n_per_basin <- tabulate(row_basin, nbins = max(basins))

code_label <- unique(coo_dt[, .(code, en_label)])
basin_code <- coo_dt[, .(n_in = .N), by = .(basin_id, code)]
overall_code <- coo_dt[, .(n_total = .N), by = code]
basin_code <- merge(basin_code, overall_code, by = "code")
basin_code <- merge(basin_code, code_label, by = "code", all.x = TRUE)

for (b in basins) {
  nb <- n_per_basin[b]
  n_out <- sum(n_per_basin[basins]) - nb
  basin_code[basin_id == b, `:=`(
    freq_in  = n_in / nb,
    n_out    = n_total - n_in,
    freq_out = (n_total - n_in) / n_out
  )]
}
basin_code[, enrichment := freq_in / pmax(freq_out, 1e-10)]

# --- Generate tag clouds ---
for (b in basins) {
  nb <- n_per_basin[b]
  df <- basin_code[basin_id == b]

  top <- df[n_in >= 100][order(-enrichment)][1:min(40, .N)]

  if (nrow(top) < 5) {
    cat(sprintf("  Basin %d: too few enriched codes, skipping\n", b))
    next
  }

  n_top <- nrow(top)
  set.seed(42)
  n_cols <- ceiling(sqrt(n_top * 1.5))
  n_rows <- ceiling(n_top / n_cols)
  top[, x := ((seq_len(.N) - 1L) %% n_cols) + runif(.N, -0.3, 0.3)]
  top[, y := ((seq_len(.N) - 1L) %/% n_cols) + runif(.N, -0.2, 0.2)]

  enrich_range <- range(top$enrichment)
  top[, size_scaled := 2.5 + 5.5 * (enrichment - enrich_range[1]) /
        max(enrich_range[2] - enrich_range[1], 1e-6)]

  g <- ggplot(top, aes(x = x, y = y, label = en_label,
                        size = size_scaled, color = freq_in)) +
    geom_text(fontface = "bold", show.legend = TRUE,
              lineheight = 0.8) +
    scale_size_identity() +
    scale_color_viridis_c(option = "plasma", name = "Freq in\nbasin",
                          begin = 0.15, end = 0.95) +
    labs(subtitle = sprintf("%s / Basin %d (n=%s)",
                            disease, b, format(nb, big.mark = ","))) +
    theme_void(base_size = 11) +
    theme(plot.subtitle = element_text(face = "bold", size = 14),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))

  out_path <- file.path(out_dir, sprintf("tagcloud_basin%d.png", b))
  ggsave(out_path, g, width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  Written: %s\n", out_path))
}

# --- Save enrichment table ---
enrich_path <- file.path(out_dir, "basin_code_enrichment.tsv")
fwrite(basin_code[order(basin_id, -enrichment)], enrich_path, sep = "\t")
cat(sprintf("  Written: %s\n", enrich_path))

cat("Done.\n")
