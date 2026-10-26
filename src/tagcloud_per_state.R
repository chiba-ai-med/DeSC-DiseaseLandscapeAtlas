## Tag cloud of characteristic comorbidity codes per Allstates pattern
##
## Same as tagcloud_per_basin.R but generates one tag cloud per pattern
## (up to 128 = 2^7 patterns) instead of per basin.
##
## Usage:
##   Rscript src/tagcloud_per_state.R <disease> <base_dir> <out_dir> <icd10_en_tsv>

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
cat(sprintf("=== Tag cloud per state: %s ===\n", disease))

# --- Build row -> pattern_id mapping ---
ls_dir <- file.path(base_dir, "landscaper", "7")
allstates <- as.matrix(fread(file.path(ls_dir, "Allstates.tsv"), header = FALSE))
sbmap <- fread(file.path(ls_dir, "state_basin_map.tsv"))
subgraph <- sbmap$basin_id  # length 128, indexed by state_id
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
row_pattern <- integer(n_rows)
ukeys <- unique(row_keys)
for (k in ukeys) {
  pid <- key_to_pid[as.character(k)]
  if (!is.na(pid)) row_pattern[row_keys == k] <- pid
}
rm(row_keys, ukeys)

patterns <- sort(unique(row_pattern[row_pattern > 0L]))
cat(sprintf("  %d rows, %d patterns with data (of %d total)\n",
            n_rows, length(patterns), n_patterns))

# --- Load col_ids ---
col_ids_dt <- fread(file.path(base_dir, "cohort", "col_ids.txt"),
                    header = FALSE, col.names = c("code", "col_id"))

# --- Load ICD-10 English names ---
icd10_en <- fread(icd10_path, sep = "\t", header = TRUE,
                  select = c("icd10_code", "english_name"))
setnames(icd10_en, c("code", "en_name"))
col_ids_dt <- merge(col_ids_dt, icd10_en, by = "code", all.x = TRUE)
col_ids_dt[is.na(en_name) | en_name == "", en_name := code]
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

coo_dt[, pattern_id := row_pattern[row_id]]
coo_dt <- coo_dt[pattern_id > 0L]
coo_dt <- merge(coo_dt, col_ids_dt[, .(col_id, code, en_label)],
                by = "col_id", all.x = FALSE)
cat(sprintf("  %d entries with pattern and code\n", nrow(coo_dt)))

# --- Compute enrichment per pattern ---
n_per_pattern <- tabulate(row_pattern, nbins = n_patterns)

code_label <- unique(coo_dt[, .(code, en_label)])
pat_code <- coo_dt[, .(n_in = .N), by = .(pattern_id, code)]
overall_code <- coo_dt[, .(n_total = .N), by = code]
pat_code <- merge(pat_code, overall_code, by = "code")
pat_code <- merge(pat_code, code_label, by = "code", all.x = TRUE)

n_total_rows <- sum(n_per_pattern[patterns])
for (p in patterns) {
  np <- n_per_pattern[p]
  n_out <- n_total_rows - np
  pat_code[pattern_id == p, `:=`(
    freq_in  = n_in / np,
    n_out    = n_total - n_in,
    freq_out = (n_total - n_in) / pmax(n_out, 1L)
  )]
}
pat_code[, enrichment := freq_in / pmax(freq_out, 1e-10)]

# --- State label: pattern sign vector + basin ---
state_labels <- character(n_patterns)
for (p in seq_len(n_patterns)) {
  signs <- ifelse(allstates[p, ] > 0, "+", "-")
  state_labels[p] <- sprintf("S%d [%s] B%d",
                              p, paste(signs, collapse = ""), subgraph[p])
}

# --- Generate tag clouds ---
min_count <- 30
n_written <- 0L
for (p in patterns) {
  np <- n_per_pattern[p]
  if (np < 100) {
    next
  }

  df <- pat_code[pattern_id == p]
  top <- df[n_in >= min_count][order(-enrichment)][1:min(40, .N)]

  if (nrow(top) < 5) next

  n_top <- nrow(top)
  set.seed(42)
  n_cols_grid <- ceiling(sqrt(n_top * 1.5))
  top[, x := ((seq_len(.N) - 1L) %% n_cols_grid) + runif(.N, -0.3, 0.3)]
  top[, y := ((seq_len(.N) - 1L) %/% n_cols_grid) + runif(.N, -0.2, 0.2)]

  enrich_range <- range(top$enrichment)
  top[, size_scaled := 2.5 + 5.5 * (enrichment - enrich_range[1]) /
        max(enrich_range[2] - enrich_range[1], 1e-6)]

  g <- ggplot(top, aes(x = x, y = y, label = en_label,
                        size = size_scaled, color = freq_in)) +
    geom_text(fontface = "bold", show.legend = TRUE,
              lineheight = 0.8) +
    scale_size_identity() +
    scale_color_viridis_c(option = "plasma", name = "Freq in\nstate",
                          begin = 0.15, end = 0.95) +
    labs(subtitle = sprintf("%s / %s (n=%s)",
                            disease, state_labels[p],
                            format(np, big.mark = ","))) +
    theme_void(base_size = 11) +
    theme(plot.subtitle = element_text(face = "bold", size = 12),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))

  out_path <- file.path(out_dir, sprintf("tagcloud_state%d.png", p))
  ggsave(out_path, g, width = 14, height = 10, dpi = 200, bg = "white")
  n_written <- n_written + 1L
}
cat(sprintf("  %d tag clouds written\n", n_written))

# --- Save enrichment table ---
enrich_path <- file.path(out_dir, "state_code_enrichment.tsv")
pat_code[, state_label := state_labels[pattern_id]]
fwrite(pat_code[order(pattern_id, -enrichment)], enrich_path, sep = "\t")
cat(sprintf("  Written: %s\n", enrich_path))

cat("Done.\n")
