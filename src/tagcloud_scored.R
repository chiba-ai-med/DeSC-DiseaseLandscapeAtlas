## Tag cloud per basin using multiple scoring methods
##
## Reads basin_code_enrichment.tsv (precomputed) and scores codes by:
##   enrichment — freq_in / freq_out (original, rare-code-dominated)
##   tfidf      — TF-IDF: tf(basin,code) * idf(code)
##   pmi        — Pointwise mutual information with Laplace smoothing
##   logitdiff  — logit(freq_in) - logit(freq_overall)
##   raw_freq   — plain freq_in (baseline)
##
## Usage:
##   Rscript src/tagcloud_scored.R <disease> <method> <enrichment_tsv> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## Usage:
##   Rscript src/tagcloud_scored.R <disease> <method> <enrichment_tsv> <out_dir> [patient_scores_tsv]
args <- commandArgs(trailingOnly = TRUE)
disease    <- args[1]
method     <- args[2]
enrich_tsv <- args[3]
out_dir    <- args[4]
scores_tsv <- if (length(args) >= 5) args[5] else NULL

stopifnot(method %in% c("enrichment", "tfidf", "pmi", "logitdiff", "raw_freq"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("=== tagcloud_scored: %s / %s ===\n", disease, method))

## --- Read enrichment table ---
dt <- fread(enrich_tsv)
basins <- sort(unique(dt$basin_id))
codes  <- sort(unique(dt$code))
n_basins <- length(basins)
n_codes  <- length(codes)

cat(sprintf("  %d basins, %d codes\n", n_basins, n_codes))

## --- Build Basin x Code count matrix (n_in) ---
basin_n <- dt[, .(n_basin = sum(n_in)), by = basin_id]
setkey(basin_n, basin_id)

## --- Patient count per basin ---
if (!is.null(scores_tsv) && file.exists(scores_tsv)) {
  ps <- fread(scores_tsv, select = "basin_id")
  basin_pat <- ps[basin_id > 0, .(n_patients = .N), by = basin_id]
} else {
  basin_pat <- dt[freq_in > 0.01, .(n_patients = as.integer(median(n_in / freq_in))),
                  by = basin_id]
}
setkey(basin_pat, basin_id)

## --- Compute scores per method ---
score_dt <- copy(dt[, .(code, basin_id, n_in, n_total, en_label, freq_in)])

## Pre-join basin size to avoid vector-length mismatch
score_dt <- merge(score_dt, basin_n, by = "basin_id")
N_total <- sum(basin_n$n_basin)

if (method == "enrichment") {
  score_dt[, `:=`(
    freq_out_val = (n_total - n_in) / pmax(N_total - n_basin, 1)
  )]
  score_dt[, score := freq_in / pmax(freq_out_val, 1e-10)]

} else if (method == "raw_freq") {
  score_dt[, score := freq_in]

} else if (method == "tfidf") {
  code_df <- dt[n_in > 0, .(df = uniqueN(basin_id)), by = code]
  score_dt <- merge(score_dt, code_df, by = "code", all.x = TRUE)
  score_dt[is.na(df), df := 0]
  score_dt[, idf := log((n_basins + 1) / (df + 1)) + 1]
  score_dt[, tf := n_in / n_basin]
  score_dt[, score := tf * idf]

} else if (method == "pmi") {
  alpha <- 1
  score_dt[, score := log((n_in + alpha) * N_total /
                          ((n_basin + alpha) * (n_total + alpha)))]

} else if (method == "logitdiff") {
  logit <- function(p) log(p) - log1p(-p)
  alpha <- 1
  score_dt[, score := logit((n_in + alpha) / (n_basin + 2 * alpha)) -
                      logit((n_total + alpha) / (N_total + 2 * alpha))]
}

## --- Generate tag clouds ---
method_label <- c(
  enrichment = "Enrichment ratio",
  tfidf      = "TF-IDF",
  pmi        = "PMI",
  logitdiff  = "Log-odds difference",
  raw_freq   = "Raw frequency"
)

for (b in basins) {
  df <- score_dt[basin_id == b]
  nb <- basin_n[.(b), n_basin]

  # Filter: freq_in > 1% AND n_in >= 100
  df <- df[freq_in > 0.01 & n_in >= 100]
  top <- df[order(-score)][1:min(30, .N)]

  if (nrow(top) < 3) {
    cat(sprintf("  Basin %d: too few codes after filtering, skipping\n", b))
    next
  }

  n_top <- nrow(top)
  set.seed(42)
  n_cols <- ceiling(sqrt(n_top * 1.5))
  n_rows_g <- ceiling(n_top / n_cols)
  top[, x := ((seq_len(.N) - 1L) %% n_cols) + runif(.N, -0.3, 0.3)]
  top[, y := ((seq_len(.N) - 1L) %/% n_cols) + runif(.N, -0.2, 0.2)]

  score_range <- range(top$score)
  top[, size_scaled := 3 + 5 * (score - score_range[1]) /
        max(score_range[2] - score_range[1], 1e-6)]

  g <- ggplot(top, aes(x = x, y = y, label = en_label,
                        size = size_scaled, color = freq_in)) +
    geom_text(fontface = "bold", show.legend = TRUE, lineheight = 0.8) +
    scale_size_identity() +
    scale_color_viridis_c(option = "plasma", name = "Freq in\nbasin",
                          begin = 0.15, end = 0.95,
                          labels = scales::percent) +
    labs(title = sprintf("%s — Basin %d (n=%s) — %s",
                         disease, b, format(nb, big.mark = ","),
                         method_label[method]),
         subtitle = sprintf("Top %d codes (freq > 1%%)", n_top)) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))

  out_path <- file.path(out_dir, sprintf("tagcloud_%s_basin%d.png", method, b))
  ggsave(out_path, g, width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  Written: %s\n", out_path))
}

## --- Generate bar charts (quantitative) ---
n_show <- as.integer(Sys.getenv("BARCHART_TOP_N", "5"))

## Capitalize disease name for display
disease_label <- paste0(toupper(substring(disease, 1, 1)), substring(disease, 2))
disease_label <- gsub("_", " ", disease_label)

## Truncate long en_label (max 40 chars)
truncate_label <- function(s, maxlen = 40) {
  ifelse(nchar(s) > maxlen, paste0(substr(s, 1, maxlen - 1), "…"), s)
}

for (b in basins) {
  df <- score_dt[basin_id == b]
  nb <- basin_n[.(b), n_basin]

  df <- df[freq_in > 0.01 & n_in >= 100]
  top <- df[order(-score)][1:min(n_show, .N)]

  if (nrow(top) < 2) next

  np <- basin_pat[.(b), n_patients]
  top[, label := paste0(truncate_label(en_label), " (", code, ")")]
  top[, label := factor(label, levels = rev(label))]

  g <- ggplot(top, aes(x = label, y = score, fill = freq_in)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", freq_in * 100)),
              hjust = -0.1, size = 4, fontface = "bold") +
    coord_flip() +
    scale_fill_viridis_c(option = "plasma", name = "Freq",
                         begin = 0.15, end = 0.95,
                         labels = scales::percent) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = sprintf("%s — Basin %d (n=%s patients)", disease_label, b,
                         format(np, big.mark = ",")),
         x = NULL, y = method_label[method]) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size = 11),
          plot.background = element_rect(fill = "white", color = NA))

  fig_h <- (1.5 + 0.7 * nrow(top)) * 0.8
  out_path <- file.path(out_dir, sprintf("barchart_%s_basin%d.png", method, b))
  ggsave(out_path, g, width = 8, height = fig_h, dpi = 300, bg = "white")
  cat(sprintf("  Written: %s\n", out_path))
}

## --- Save scored table ---
score_path <- file.path(out_dir, sprintf("basin_code_%s.tsv", method))
fwrite(score_dt[order(basin_id, -score)], score_path, sep = "\t")
cat(sprintf("  Written: %s\n", score_path))

cat("Done.\n")
