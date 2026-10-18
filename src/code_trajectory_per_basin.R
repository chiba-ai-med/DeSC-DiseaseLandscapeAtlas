## ICD-10 code trajectory per basin before death/hospitalization
##
## For each basin, tracks how disease code composition changes over time
## before the event. Reads binary COO matrix to get raw comorbidity codes.
##
## Usage:
##   Rscript src/code_trajectory_per_basin.R <disease> <base_dir> <out_dir> \
##     <traj_tsv> <endpoint> <icd10_en_tsv>
##
##   endpoint: "death" or "hosp"
##   traj_tsv: trajectories_dead.tsv or trajectories_hospitalized.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
disease    <- args[1]
base_dir   <- args[2]
out_dir    <- args[3]
traj_path  <- args[4]
endpoint   <- args[5]
icd10_path <- args[6]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Code trajectory per basin (%s): %s ===\n", endpoint, disease))

period_col <- if (endpoint == "death") "periods_before_death" else "periods_before_hosp"

# --- 1. Load trajectory (already has kojin_id, window, basin_id, periods_before) ---
traj <- fread(traj_path)
traj[, kojin_id := as.character(kojin_id)]
cat(sprintf("  %d trajectory rows, %d patients\n", nrow(traj), uniqueN(traj$kojin_id)))

max_period <- min(19L, max(traj[[period_col]]))
traj <- traj[get(period_col) <= max_period]

# --- 2. Load person_map: new_row_id -> kojin_id ---
pm <- fread(file.path(base_dir, "cohort", "person_map.tsv"))
pm[, kojin_id := as.character(kojin_id)]

# --- 3. Load window lookup: old_row_id -> window ---
wl <- fread("data/rowid_window_lookup.tsv", header = FALSE,
            col.names = c("old_row_id", "window"))

# --- 4. Load row_ids: old_row_id -> new_row_id ---
ri <- fread(file.path(base_dir, "cohort", "row_ids.txt"), header = FALSE,
            col.names = c("old_row_id", "new_row_id"))

# --- 5. Build new_row_id -> window mapping ---
ri_wl <- merge(ri, wl, by = "old_row_id")

# --- 6. Build new_row_id -> kojin_id + window ---
row_info <- merge(ri_wl[, .(new_row_id, window)],
                  pm[, .(new_row_id, kojin_id)],
                  by = "new_row_id")
rm(pm, wl, ri, ri_wl); gc()

# --- 7. Filter to trajectory patients + windows ---
traj_key <- unique(traj[, .(kojin_id, window, basin_id, periods = get(period_col))])
row_traj <- merge(row_info, traj_key, by = c("kojin_id", "window"))
cat(sprintf("  %d COO rows matched to trajectory\n", nrow(row_traj)))
rm(row_info, traj_key); gc()

# --- 8. Load col_ids ---
col_ids_dt <- fread(file.path(base_dir, "cohort", "col_ids.txt"),
                    header = FALSE, col.names = c("code", "col_id"))

# --- 9. Load ICD-10 English names ---
icd10_en <- fread(icd10_path, sep = "\t", header = TRUE,
                  select = c("icd10_code", "english_name"))
setnames(icd10_en, c("code", "en_name"))
col_ids_dt <- merge(col_ids_dt, icd10_en, by = "code", all.x = TRUE)
col_ids_dt[is.na(en_name) | en_name == "", en_name := code]
col_ids_dt[, en_label := vapply(en_name, function(s) {
  paste(strwrap(s, width = 35), collapse = "\n")
}, character(1))]

# --- 10. Read binary COO ---
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

# --- 11. Filter COO to trajectory rows ---
target_rows <- row_traj[, .(new_row_id, basin_id, periods)]
setnames(target_rows, "new_row_id", "row_id")
setkey(target_rows, row_id)
setkey(coo_dt, row_id)

coo_traj <- coo_dt[target_rows, nomatch = 0L]
rm(coo_dt); gc()
cat(sprintf("  %d COO entries in trajectory windows\n", nrow(coo_traj)))

coo_traj <- merge(coo_traj, col_ids_dt[, .(col_id, code, en_label)],
                  by = "col_id", all.x = FALSE)

# --- 12. Count patients per basin x period (denominator) ---
n_patients_bp <- row_traj[, .(n_patients = uniqueN(kojin_id)),
                           by = .(basin_id, periods)]
rm(row_traj); gc()

# --- 13. Compute code frequency per basin x period ---
code_freq <- coo_traj[, .(n_codes = .N), by = .(basin_id, periods, code, en_label)]
code_freq <- merge(code_freq, n_patients_bp, by = c("basin_id", "periods"))
code_freq[, freq := n_codes / n_patients]
rm(coo_traj); gc()

# --- 14. Save full table ---
tsv_path <- file.path(out_dir,
  sprintf("basin_code_trajectory_%s.tsv", endpoint))
fwrite(code_freq[order(basin_id, periods, -freq)], tsv_path, sep = "\t")
cat(sprintf("  Written: %s\n", tsv_path))

# --- 15. Heatmap per basin (top 30 codes by overall frequency) ---
basins <- sort(unique(code_freq$basin_id[code_freq$basin_id > 0L]))

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

event_label <- if (endpoint == "death") "death" else "hospitalization"

for (b in basins) {
  df <- code_freq[basin_id == b]
  if (nrow(df) < 10) next

  top_codes <- df[, .(total = sum(n_codes)), by = .(code, en_label)]
  setorder(top_codes, -total)
  top30 <- top_codes[1:min(30, nrow(top_codes))]

  hm <- df[code %in% top30$code]
  hm[, en_label := factor(en_label, levels = rev(top30$en_label))]

  g <- ggplot(hm, aes(x = -periods, y = en_label, fill = freq)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "inferno", name = "Freq\n(per patient)",
                          begin = 0.05, end = 0.95) +
    scale_x_continuous(
      breaks = seq(-max_period, 0, by = 2),
      labels = function(x) ifelse(x == 0, event_label, as.character(-x))) +
    labs(title = sprintf("%s / Basin %d: ICD-10 code trajectory before %s",
                          disease, b, event_label),
         x = sprintf("Periods before %s", event_label),
         y = NULL) +
    theme_pub +
    theme(axis.text.y = element_text(size = 7))

  png_path <- file.path(out_dir,
    sprintf("basin_code_traj_%s_basin%d.png", endpoint, b))
  ggsave(png_path, g, width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  Written: %s\n", png_path))
}

cat("Done.\n")
