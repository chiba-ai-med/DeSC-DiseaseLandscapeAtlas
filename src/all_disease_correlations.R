## Correlation matrix across all diseases: PC1-7, energy, age, sex, death,
## obs_years, hospitalized, n_admissions, kenshin_data_ari, (n_diagnoses for HT)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "summary"

fids <- list.dirs("output", recursive = FALSE, full.names = FALSE)

all_pts <- list()
for (fid in fids) {
  path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
  if (!file.exists(path)) next
  pt <- fread(path)
  pt[, disease := fid]

  ndiag_path <- sprintf("output/%s/all_comorbid/cohort/n_diagnoses.tsv", fid)
  if (file.exists(ndiag_path)) {
    nd <- fread(ndiag_path)
    pt <- merge(pt, nd, by = "kojin_id", all.x = TRUE)
  } else {
    pt[, n_diagnoses := NA_integer_]
  }

  all_pts[[length(all_pts) + 1]] <- pt
  cat(sprintf("  %s: %d patients\n", fid, nrow(pt)))
}

dt <- rbindlist(all_pts, fill = TRUE)
cat(sprintf("\nTotal: %d patients from %d diseases\n", nrow(dt), length(all_pts)))

dt[, death_bin := as.integer(death == "dead")]
dt[, sex_f := as.integer(sex == "female")]
dt[, kenshin := as.integer(kenshin_data_ari)]
dt[, hosp_bin := as.integer(hospitalized == "yes")]
dt[is.na(hosp_bin), hosp_bin := as.integer(hospitalized)]

# --- Correlation matrix (all diseases, without n_diagnoses) ---
vars_all <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
              "energy", "age", "sex_f", "death_bin", "obs_years",
              "hosp_bin", "n_admissions", "kenshin")
var_labels <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
                "Energy", "Age", "Sex(F)", "Death", "Obs years",
                "Hospitalized", "N admissions", "Kenshin")

mat <- as.matrix(dt[, ..vars_all])
cor_mat <- cor(mat, use = "pairwise.complete.obs")
rownames(cor_mat) <- var_labels
colnames(cor_mat) <- var_labels

cat("\n=== Correlation matrix (all 51 diseases combined) ===\n")
print(round(cor_mat, 3))

cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("Var1", "Var2", "cor"))

g1 <- ggplot(cor_dt, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 2.5) +
  scale_fill_gradient2(low = "#C62828", mid = "white", high = "#1565C0",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation matrix (all 51 diseases combined)",
       x = NULL, y = NULL, fill = "cor") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "all_disease_correlation_matrix.png"),
       g1, width = 10, height = 9, dpi = 150)

# --- With n_diagnoses (hypertension only) ---
dt_ht <- dt[disease == "hypertension" & !is.na(n_diagnoses)]
vars_ht <- c(vars_all, "n_diagnoses")
var_labels_ht <- c(var_labels, "N diagnoses")

mat_ht <- as.matrix(dt_ht[, ..vars_ht])
cor_ht <- cor(mat_ht, use = "pairwise.complete.obs")
rownames(cor_ht) <- var_labels_ht
colnames(cor_ht) <- var_labels_ht

cat("\n=== Correlation matrix (hypertension, with n_diagnoses) ===\n")
print(round(cor_ht, 3))

cor_ht_dt <- as.data.table(as.table(cor_ht))
setnames(cor_ht_dt, c("Var1", "Var2", "cor"))

g2 <- ggplot(cor_ht_dt, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 2.5) +
  scale_fill_gradient2(low = "#C62828", mid = "white", high = "#1565C0",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation matrix (hypertension, with n_diagnoses)",
       x = NULL, y = NULL, fill = "cor") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "hypertension_correlation_matrix.png"),
       g2, width = 11, height = 10, dpi = 150)

# --- Per-disease correlations ---
cat("\n=== Per-disease correlations (key pairs) ===\n")
key_pairs <- list(
  c("energy", "age"),
  c("energy", "sex_f"),
  c("energy", "death_bin"),
  c("energy", "hosp_bin"),
  c("energy", "n_admissions"),
  c("energy", "kenshin"),
  c("energy", "obs_years"),
  c("age", "death_bin"),
  c("PC1", "age"),
  c("PC1", "energy")
)

per_disease_rows <- list()
for (fid in sort(unique(dt$disease))) {
  sub <- dt[disease == fid]
  for (pair in key_pairs) {
    v1 <- sub[[pair[1]]]
    v2 <- sub[[pair[2]]]
    r <- cor(v1, v2, use = "pairwise.complete.obs")
    per_disease_rows[[length(per_disease_rows) + 1]] <- data.table(
      disease = fid, var1 = pair[1], var2 = pair[2], cor = r
    )
  }
}
per_disease_dt <- rbindlist(per_disease_rows)
per_disease_wide <- dcast(per_disease_dt, disease ~ var1 + var2, value.var = "cor")
fwrite(per_disease_wide, file.path(out_dir, "per_disease_key_correlations.tsv"), sep = "\t")

cat("Key correlations by disease:\n")
print(per_disease_dt[var1 == "energy" & var2 == "kenshin"][order(cor)],
      nrows = 60)

# --- Pair plot (all diseases combined, sampled) ---
cat("\nGenerating pair plot (sampled 100k patients)...\n")
set.seed(42)
samp <- dt[sample(.N, min(100000, .N))]

plot_vars <- c("PC1", "PC2", "PC3", "energy", "age", "sex_f",
               "death_bin", "kenshin", "hosp_bin", "n_admissions", "obs_years")
plot_labels <- c("PC1", "PC2", "PC3", "Energy", "Age", "Sex(F)",
                 "Death", "Kenshin", "Hosp", "N adm", "Obs yrs")

# Use GGally if available, otherwise simple cor heatmap is enough
pair_data <- samp[, ..plot_vars]
setnames(pair_data, plot_labels)

# Scatterplot matrix with correlations
pair_cor <- cor(pair_data, use = "pairwise.complete.obs")
pair_cor_dt <- as.data.table(as.table(pair_cor))
setnames(pair_cor_dt, c("Var1", "Var2", "cor"))

# Simple hex-binned pair plots for key combinations
plot_pairs <- list(
  c("Energy", "Age"),
  c("Energy", "Death"),
  c("Energy", "Kenshin"),
  c("Energy", "N adm"),
  c("Energy", "PC1"),
  c("PC1", "Age"),
  c("PC1", "Kenshin"),
  c("Age", "Death"),
  c("Kenshin", "Death")
)

plots <- list()
for (p in plot_pairs) {
  g <- ggplot(pair_data, aes(x = .data[[p[1]]], y = .data[[p[2]]])) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(trans = "log10") +
    labs(x = p[1], y = p[2]) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "none")
  plots[[length(plots) + 1]] <- g
}

for (i in seq_along(plots)) {
  fname <- sprintf("pair_%s_vs_%s.png", plot_pairs[[i]][1], plot_pairs[[i]][2])
  fname <- gsub("[() ]", "", fname)
  ggsave(file.path(out_dir, fname), plots[[i]], width = 5, height = 4, dpi = 150)
}

# --- Save tables ---
write.csv(round(cor_mat, 4), file.path(out_dir, "all_disease_correlation_matrix.csv"))
write.csv(round(cor_ht, 4), file.path(out_dir, "hypertension_correlation_matrix.csv"))

cat("\nDone.\n")
