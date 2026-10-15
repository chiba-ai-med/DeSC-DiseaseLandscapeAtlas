## Systematic PC-energy correlation analysis across all disease families
##
## For each disease x PC: compute correlation with energy
## Identify the most energy-driving PC per disease, extract top loadings
##
## Usage:
##   Rscript src/pc_energy_analysis.R <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fids <- list.dirs("output", recursive = FALSE, full.names = FALSE)

all_scores <- list()
all_loadings <- list()

# Load ICD-10 English names
en_path <- "data/col_id_disease_name_en_small.txt"
icd_en <- list()
if (file.exists(en_path)) {
  en_lines <- readLines(en_path, encoding = "UTF-8")
  for (ln in en_lines) {
    parts <- strsplit(ln, ",")[[1]]
    if (length(parts) >= 2) {
      icd_en[[trimws(parts[1])]] <- trimws(paste(parts[-(1:2)], collapse = ","))
    }
  }
}

for (fid in fids) {
  scores_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
  eigvec_path <- sprintf("output/%s/all_comorbid/pca/7/Eigen_vectors.csv", fid)
  col_ids_path <- sprintf("output/%s/all_comorbid/cohort/col_ids.txt", fid)
  mort_path <- sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check.tsv", fid)

  if (!file.exists(scores_path) || !file.exists(eigvec_path)) next

  pt <- fread(scores_path)
  pt <- pt[is.finite(energy)]
  if (nrow(pt) < 100) next

  pt[, death_bin := as.integer(death == "dead")]

  # Mortality Cox result
  mort_hr <- NA_real_
  mort_ok <- NA
  if (file.exists(mort_path)) {
    m <- fread(mort_path)
    mort_hr <- m$HR[1]
    mort_ok <- m$direction_ok[1]
  }

  pc_cols <- paste0("PC", 1:7)
  pc_cols <- intersect(pc_cols, names(pt))

  for (pc in pc_cols) {
    vals <- pt[[pc]]
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) next

    cor_energy <- cor(vals, pt$energy, use = "pairwise.complete.obs")
    cor_death <- cor(vals, pt$death_bin, use = "pairwise.complete.obs")

    all_scores[[length(all_scores) + 1]] <- data.table(
      disease = fid, pc = pc,
      cor_energy = cor_energy,
      cor_death = cor_death,
      mort_HR = mort_hr,
      mort_pass = mort_ok,
      n_patients = nrow(pt)
    )
  }

  # --- Top loadings for each PC ---
  if (!file.exists(col_ids_path)) next
  eigvec <- tryCatch(as.matrix(read.csv(eigvec_path, header = FALSE)),
                     error = function(e) NULL)
  if (is.null(eigvec)) next

  col_lines <- readLines(col_ids_path, encoding = "UTF-8")
  code_names <- character(nrow(eigvec))
  for (ln in col_lines) {
    parts <- strsplit(ln, ",")[[1]]
    if (length(parts) >= 2) {
      code <- trimws(parts[1])
      id <- as.integer(trimws(parts[2]))
      if (!is.na(id) && id >= 1 && id <= length(code_names)) {
        code_names[id] <- code
      }
    }
  }

  n_top <- 10
  for (pc_idx in seq_along(pc_cols)) {
    if (pc_idx > ncol(eigvec)) break
    v <- eigvec[, pc_idx]
    top_idx <- order(-abs(v))[1:min(n_top, length(v))]
    for (r in seq_along(top_idx)) {
      ti <- top_idx[r]
      code <- code_names[ti]
      en_name <- if (!is.null(icd_en[[code]])) icd_en[[code]] else ""
      all_loadings[[length(all_loadings) + 1]] <- data.table(
        disease = fid, pc = pc_cols[pc_idx],
        rank = r,
        icd10 = code,
        disease_name_en = en_name,
        loading = v[ti]
      )
    }
  }

  cat(sprintf("  %s: done\n", fid))
}

scores_dt <- rbindlist(all_scores)
loadings_dt <- rbindlist(all_loadings)

# --- Best PC per disease (highest |cor_energy|) ---
scores_dt[, abs_cor_energy := abs(cor_energy)]
best_pc <- scores_dt[, .SD[which.max(abs_cor_energy)], by = disease]
best_pc <- best_pc[order(-abs_cor_energy)]

fwrite(scores_dt, file.path(out_dir, "pc_energy_correlation_all.tsv"), sep = "\t")
fwrite(best_pc, file.path(out_dir, "pc_energy_best.tsv"), sep = "\t")
fwrite(loadings_dt, file.path(out_dir, "pc_top_loadings_all.tsv"), sep = "\t")

cat("\n=== Disease ranking by |cor(best_PC, energy)| ===\n")
print(best_pc[, .(disease, pc, cor_energy, cor_death, mort_HR,
                   mort_pass, n_patients)], nrows = 60)

# --- Best PC loadings ---
best_loadings <- merge(best_pc[, .(disease, pc)], loadings_dt,
                       by = c("disease", "pc"))
fwrite(best_loadings, file.path(out_dir, "best_pc_loadings.tsv"), sep = "\t")

# --- Plot 1: |cor(PC, energy)| ranking ---
best_pc[, disease_f := factor(disease, levels = rev(disease))]
best_pc[, pass_label := ifelse(mort_pass, "PASS", "FAIL")]
g1 <- ggplot(best_pc, aes(x = abs_cor_energy, y = disease_f,
                            fill = pass_label)) +
  geom_col() +
  scale_fill_manual(values = c(PASS = "#1565C0", FAIL = "#C62828")) +
  labs(title = "Best PC-energy correlation by disease family",
       subtitle = "|cor(PC, energy)| — which PC drives energy most?",
       x = "|Correlation with energy|", y = NULL, fill = "Mortality") +
  theme_minimal(base_size = 10)
ggsave(file.path(out_dir, "pc_energy_ranking.png"),
       g1, width = 10, height = 14, dpi = 150)

# --- Plot 2: heatmap disease x PC ---
scores_dt[, disease_f := factor(disease, levels = rev(best_pc$disease))]
g2 <- ggplot(scores_dt, aes(x = pc, y = disease_f, fill = abs_cor_energy)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#1565C0") +
  labs(title = "PC x Disease: |Correlation with energy|",
       x = NULL, y = NULL, fill = "|cor|") +
  theme_minimal(base_size = 9)
ggsave(file.path(out_dir, "pc_energy_heatmap.png"),
       g2, width = 8, height = 14, dpi = 150)

# --- Plot 3: cor_energy vs cor_death scatter ---
g3 <- ggplot(scores_dt, aes(x = cor_energy, y = cor_death, color = pc)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_vline(xintercept = 0, linetype = 2, color = "gray40") +
  labs(title = "PC correlation with energy vs death (all diseases x PCs)",
       x = "cor(PC, energy)", y = "cor(PC, death)", color = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "pc_energy_vs_death.png"),
       g3, width = 9, height = 7, dpi = 150)

cat("\nDone.\n")
