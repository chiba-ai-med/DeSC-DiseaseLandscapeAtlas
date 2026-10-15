## Hypertension deep dive: what comorbidities distinguish low vs high energy?
##
## Split patients into energy tertiles, compare ICD-10 code prevalence
## Show that low energy (severe) patients have downstream organ damage

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

out_dir <- "summary/hypertension_deep_dive"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fid <- "hypertension"

# --- Load patient scores ---
pt <- fread(sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid))
pt <- pt[is.finite(energy)]
pt[, death_bin := as.integer(death == "dead")]

# Energy tertiles
pt[, E_tertile := cut(energy, breaks = quantile(energy, c(0, 1/3, 2/3, 1)),
                       labels = c("Low", "Mid", "High"), include.lowest = TRUE)]

cat("=== Energy tertile summary ===\n")
print(pt[, .(n = .N, death_rate = mean(death_bin),
             mean_age = mean(age), pct_female = mean(sex == "female")),
         by = E_tertile][order(E_tertile)])

# --- Load COO matrix to get per-patient comorbidity vectors ---
coo_path <- sprintf("output/%s/all_comorbid/cohort/coo.tsv.zst", fid)
col_ids_path <- sprintf("output/%s/all_comorbid/cohort/col_ids.txt", fid)
row_ids_path <- sprintf("output/%s/all_comorbid/cohort/row_ids.txt", fid)
person_map_path <- sprintf("output/%s/all_comorbid/cohort/person_map.tsv", fid)

# Read col_ids (ICD-10 code -> column index)
col_lines <- readLines(col_ids_path, encoding = "UTF-8")
col_map <- data.table(
  icd10 = character(length(col_lines)),
  col_id = integer(length(col_lines))
)
for (i in seq_along(col_lines)) {
  parts <- strsplit(col_lines[i], ",")[[1]]
  col_map[i, icd10 := trimws(parts[1])]
  col_map[i, col_id := as.integer(trimws(parts[2]))]
}

# Read person_map to link new_row_id -> kojin_id
pm <- fread(person_map_path)

# For each patient (using latest period), get their comorbidity profile
# But COO is too large to load in R. Use patient_scores which already has kojin_id.
# Alternative: use PCA scores directly - the eigenvectors tell us what each PC means.

# More efficient approach: for each ICD-10 code, compute prevalence in each tertile
# using the COO matrix column sums per tertile group.

# Actually, let's use the binarized PCA scores approach:
# patient PC scores * eigenvectors = approximate reconstruction of comorbidities
# But this loses sparsity.

# Simplest correct approach: read eigenvectors, reconstruct importance per ICD code
# via PC scores weighted by eigenvector loadings.

eigvec_path <- sprintf("output/%s/all_comorbid/pca/7/Eigen_vectors.csv", fid)
eigvec <- as.matrix(read.csv(eigvec_path, header = FALSE))

# For each patient, approximate comorbidity profile = PC_scores %*% t(eigvec)
# Then average by tertile to get "mean comorbidity profile"

pc_cols <- paste0("PC", 1:7)
pc_mat <- as.matrix(pt[, ..pc_cols])

# Mean PC scores per tertile
mean_pcs <- pt[, lapply(.SD, mean), .SDcols = pc_cols, by = E_tertile]
cat("\n=== Mean PC scores by energy tertile ===\n")
print(mean_pcs)

# Reconstruct mean comorbidity profile per tertile
# profile = mean_PC %*% t(eigvec)  [n_tertiles x n_codes]
mean_pc_mat <- as.matrix(mean_pcs[, ..pc_cols])
profiles <- mean_pc_mat %*% t(eigvec)
colnames(profiles) <- col_map$icd10

# Difference: Low - High
diff_profile <- profiles[1, ] - profiles[3, ]  # Low minus High
diff_dt <- data.table(
  icd10 = names(diff_profile),
  diff_score = as.numeric(diff_profile),
  low_score = as.numeric(profiles[1, ]),
  high_score = as.numeric(profiles[3, ])
)
diff_dt[, abs_diff := abs(diff_score)]
diff_dt <- diff_dt[order(-abs_diff)]

# Add disease names
icd_labels <- list(
  I10 = "Essential hypertension",
  I509 = "Heart failure",
  I209 = "Ischemic heart disease",
  I639 = "Cerebral infarction",
  I693 = "Sequelae of cerebral infarction",
  N189 = "CKD",
  N19 = "Renal failure",
  E14 = "Diabetes",
  E119 = "Type 2 diabetes",
  E785 = "Hyperlipidemia",
  J449 = "COPD",
  G309 = "Alzheimer's",
  F009 = "Dementia",
  K210 = "GERD",
  K590 = "Constipation/IBS",
  H522 = "Astigmatism",
  G470 = "Insomnia",
  M5456 = "Low back pain",
  M8199 = "Osteoporosis",
  K295 = "Chronic gastritis",
  M171 = "Knee osteoarthritis",
  L309 = "Dermatitis",
  J189 = "Pneumonia",
  K053 = "Chronic periodontitis",
  K029 = "Dental caries",
  K081 = "Loss of teeth",
  T888 = "Surgical complications",
  J304 = "Allergic rhinitis",
  I489 = "Atrial fibrillation",
  C61 = "Prostate cancer",
  N40 = "BPH",
  D509 = "Iron deficiency anemia",
  I251 = "Atherosclerotic heart disease",
  I259 = "Chronic IHD",
  I500 = "Congestive heart failure",
  I639 = "Cerebral infarction",
  I64 = "Stroke NOS",
  I679 = "Cerebrovascular disease",
  N180 = "End-stage renal disease",
  N183 = "CKD stage 3",
  N184 = "CKD stage 4",
  N185 = "CKD stage 5",
  H269 = "Cataract",
  Z961 = "Intraocular lens"
)
diff_dt[, label := sapply(icd10, function(x) {
  if (x %in% names(icd_labels)) paste0(x, " (", icd_labels[[x]], ")")
  else x
})]

fwrite(diff_dt, file.path(out_dir, "comorbidity_diff_low_vs_high.tsv"), sep = "\t")

cat("\n=== Top 30 comorbidities enriched in LOW energy (severe) ===\n")
print(diff_dt[diff_score > 0][1:30, .(icd10, label, diff_score)])

cat("\n=== Top 30 comorbidities enriched in HIGH energy (mild) ===\n")
print(diff_dt[diff_score < 0][order(diff_score)][1:30, .(icd10, label, diff_score)])

# --- Plot: top enriched in low vs high ---
top_n <- 25
top_low <- diff_dt[diff_score > 0][1:top_n]
top_high <- diff_dt[diff_score < 0][order(diff_score)][1:top_n]
plot_dt <- rbind(
  top_low[, .(label, diff_score, direction = "Enriched in LOW energy (severe)")],
  top_high[, .(label, diff_score, direction = "Enriched in HIGH energy (mild)")]
)
plot_dt[, label := factor(label, levels = rev(label))]

g1 <- ggplot(plot_dt, aes(x = diff_score, y = label, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Enriched in LOW energy (severe)" = "#C62828",
                                "Enriched in HIGH energy (mild)" = "#1565C0")) +
  geom_vline(xintercept = 0, color = "gray20") +
  labs(title = "Hypertension: comorbidities distinguishing low vs high energy patients",
       subtitle = "PCA-reconstructed comorbidity difference (Low energy - High energy)",
       x = "Difference score (Low - High)", y = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "comorbidity_low_vs_high.png"),
       g1, width = 11, height = 12, dpi = 150)

# --- Check specific organ damage codes ---
organ_damage <- c("I509", "I209", "I639", "I693", "N189", "N19",
                  "I489", "I500", "I251", "I259", "I679",
                  "N180", "N183", "N184", "N185", "I64",
                  "G309", "F009", "J189")
cat("\n=== Organ damage codes: Low vs High energy ===\n")
od_dt <- diff_dt[icd10 %in% organ_damage][order(-diff_score)]
print(od_dt[, .(icd10, label, diff_score, low_score, high_score)])

cat("\nDone.\n")
