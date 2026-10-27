suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(grid)
})
base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)

survfit_to_dt <- function(sf) {
  strata_names <- if (!is.null(sf$strata)) {
    rep(names(sf$strata), sf$strata)
  } else { rep("All", length(sf$time)) }
  data.table(time = sf$time, surv = sf$surv,
             upper = sf$upper, lower = sf$lower,
             strata = strata_names)
}

T1 <- list(
  Homogeneous = c("breast_cancer","colorectal_cancer","interstitial_lung_disease",
    "ischemic_heart_disease","liver_disease","lung_cancer","metastatic_cancer",
    "obesity","stomach_cancer","stroke"),
  Reversed = c("autoimmune_connective_tissue","cirrhosis_fibrosis","copd","gout",
    "hip_fracture","pancreatic_cancer","tuberculosis"),
  Symmetric = c("ckd","hypertension","malnutrition"),
  Convergent = c("atrial_fibrillation","dementia","heart_failure","sleep_disorders"),
  Progressive = c("anemia","aspiration_pneumonia","asthma","diabetes",
    "hematologic_malignancy","inflammatory_bowel_disease","malignancy",
    "osteoarthritis","osteoporosis","peptic_ulcer"),
  "Structured reversed" = c("acute_renal_failure","anxiety",
    "benign_prostatic_hyperplasia","delirium","depression","dyslipidemia",
    "epilepsy","liver_cancer","parkinson","pneumonia","prostate_cancer",
    "rheumatoid_arthritis","schizophrenia","sepsis","thyroid_disorder",
    "urinary_tract_infection","viral_hepatitis")
)
classify <- function(disease) {
  for (pat in names(T1)) if (disease %in% T1[[pat]]) return(pat)
  NA_character_
}
classify_v <- Vectorize(classify)

name_map <- c(
  copd="COPD", ckd="CKD",
  acute_renal_failure="AKI", aspiration_pneumonia="Asp. pneumonia",
  atrial_fibrillation="AF", autoimmune_connective_tissue="Autoimmune CT",
  benign_prostatic_hyperplasia="BPH", breast_cancer="Breast ca.",
  cirrhosis_fibrosis="Cirrhosis", colorectal_cancer="Colorectal ca.",
  hematologic_malignancy="Hem. malig.", hip_fracture="Hip fracture",
  inflammatory_bowel_disease="IBD", interstitial_lung_disease="ILD",
  ischemic_heart_disease="IHD", liver_cancer="Liver ca.",
  lung_cancer="Lung ca.", metastatic_cancer="Metastatic ca.",
  osteoarthritis="OA", pancreatic_cancer="Pancreatic ca.",
  peptic_ulcer="Peptic ulcer", prostate_cancer="Prostate ca.",
  rheumatoid_arthritis="RA", sleep_disorders="Sleep dis.",
  stomach_cancer="Stomach ca.", thyroid_disorder="Thyroid dis.",
  urinary_tract_infection="UTI", viral_hepatitis="Viral hep."
)
nice <- function(d) ifelse(d %in% names(name_map), name_map[d], gsub("_"," ",d))

# ============================================================
# Item 2a: Sepsis mortality KM — redraw from existing cohort
# ============================================================
cat("=== FigS8: Sepsis mortality KM ===\n")
sep <- fread("output/sepsis/all_comorbid/landmark/landmark_cohort.tsv",
             select = c("basin_id","shibou_flg_riyouka","event_death","fu_time"))
sep_mort <- sep[shibou_flg_riyouka == 1L]
sep_mort[, basin_label := paste0("Basin ", basin_id)]

for (b in sort(unique(sep_mort[basin_id > 0]$basin_id))) {
  sub <- sep_mort[basin_id == b]
  cat(sprintf("  B%d: n=%d, deaths=%d, death_rate=%.1f%%\n",
      b, nrow(sub), sum(sub$event_death), mean(sub$event_death)*100))
}

fit_sm <- survfit(Surv(fu_time, event_death) ~ basin_label, data = sep_mort)
km_sm <- survfit_to_dt(fit_sm)
km_sm[, strata := gsub("basin_label=", "", strata)]

sd_sm <- tryCatch({
  sd <- survdiff(Surv(fu_time, event_death) ~ basin_label, data = sep_mort)
  1 - pchisq(sd$chisq, length(sd$n) - 1)
}, error = function(e) NA)

g_sm <- ggplot(km_sm, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 0.8) +
  labs(title = sprintf("Sepsis - Mortality KM (n = %s)",
                       format(nrow(sep_mort), big.mark = ",")),
       subtitle = if (is.na(sd_sm)) "" else if (sd_sm < 0.001) "p < 0.001" else sprintf("p = %.3f", sd_sm),
       x = "Years from landmark", y = "Survival probability", color = "Basin") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + coord_cartesian(xlim = c(0, 8))
ggsave("plot/Figures/Supplementary/FigS8_sepsis_km_mortality.png", g_sm, width = 8, height = 6, dpi = 200)
cat("  Written: FigS8_sepsis_km_mortality.png\n")
rm(sep, sep_mort); gc(verbose = FALSE)

# ============================================================
# Item 2b: Thyroid composite — rebuild from existing panels
# ============================================================
cat("\n=== Fig06: Thyroid composite ===\n")
thy_panels <- c(
  "plot/Figures/Fig06/Fig06a_thyroid_basin3_vs_4.png",
  "plot/Figures/Fig06/Fig06b_thyroid_km_inter_basin.png",
  "plot/Figures/Fig06/Fig06c_thyroid_disorder_discon.png",
  "plot/Figures/Fig06/Fig06d_thyroid_disorder_sankey.png"
)
# Check existence
for (p in thy_panels) cat(sprintf("  %s: %s\n", basename(p), file.exists(p)))

# Since we can't use png package, save panel b as updated (already done in 5553ec8)
cat("  KM panel already updated. Composite requires png package — skipping assembly.\n")

# ============================================================
# Item 2c: Diabetes composite — rebuild from existing panels
# ============================================================
cat("\n=== FigS7: Diabetes composite ===\n")
dia_panels <- c(
  "plot/Figures/Supplementary/Fig06a_diabetes_basin1_vs_3.png",
  "plot/Figures/Supplementary/Fig06b_diabetes_km_inter_basin.png",
  "plot/Figures/Supplementary/Fig06c_diabetes_discon.png",
  "plot/Figures/Supplementary/Fig06d_diabetes_sankey.png"
)
for (p in dia_panels) cat(sprintf("  %s: %s\n", basename(p), file.exists(p)))
cat("  KM panel already updated. Composite requires png package — skipping assembly.\n")

# ============================================================
# Item 3: Recompute asymmetry from deceased-only transition data
# ============================================================
cat("\n=== FigS4: Recompute asymmetry ===\n")

trans_files <- Sys.glob("plot/*/basin_transition_matrix.tsv")
cat(sprintf("  Found %d deceased-only transition files\n", length(trans_files)))

all_asym <- list()
for (f in trans_files) {
  d <- basename(dirname(f))
  tm <- fread(f)
  # Only cross-basin transitions
  cross <- tm[from != to]
  if (nrow(cross) == 0) next

  n_basins <- length(unique(c(cross$from, cross$to)))
  total_steps <- sum(tm$N)
  same_basin <- sum(tm[from == to]$N)
  diff_basin <- sum(cross$N)
  move_pct <- diff_basin / total_steps * 100

  # Compute pairwise asymmetry
  pairs <- unique(t(apply(cbind(cross$from, cross$to), 1, sort)))
  asyms <- numeric()
  for (r in seq_len(nrow(pairs))) {
    b1 <- pairs[r, 1]; b2 <- pairs[r, 2]
    r12 <- cross[from == b1 & to == b2]$prob
    r21 <- cross[from == b2 & to == b1]$prob
    if (length(r12) == 1 && length(r21) == 1 && (r12 + r21) > 0) {
      asyms <- c(asyms, abs(r12 - r21) / (r12 + r21))
    }
  }
  mean_asym <- if (length(asyms) > 0) mean(asyms) else NA_real_

  all_asym[[d]] <- data.table(
    disease = d, n_basins = n_basins,
    total_steps = total_steps, same_basin = same_basin,
    diff_basin = diff_basin, per_step_move_pct = move_pct,
    mean_asymmetry = mean_asym, n_pairs = length(asyms))
}
asym_dt <- rbindlist(all_asym)
asym_dt[, pattern := classify_v(disease)]

# Report key diseases
cat("  Key asymmetry values (deceased-only):\n")
for (d in c("dementia","diabetes","anemia","thyroid_disorder","sepsis")) {
  row <- asym_dt[disease == d]
  if (nrow(row) > 0)
    cat(sprintf("    %s: asym=%.4f, move=%.2f%%\n", d, row$mean_asymmetry, row$per_step_move_pct))
}

fwrite(asym_dt, "summary/transition_asymmetry_summary.tsv", sep = "\t")
cat("  Updated transition_asymmetry_summary.tsv\n")

# Also update asymmetry_corrected.tsv
fwrite(asym_dt, "plot/Figures/asymmetry_corrected.tsv", sep = "\t")
cat("  Updated asymmetry_corrected.tsv\n")

# ============================================================
# Redraw FigS4 pattern space scatter
# ============================================================
cat("\n=== FigS4: Pattern space scatter ===\n")
multi <- asym_dt[!pattern %in% c("Homogeneous", "Reversed") & !is.na(pattern)]
multi[, label := nice(disease)]
cat(sprintf("  Multi-basin diseases: %d\n", nrow(multi)))

pattern_colors <- c(
  "Progressive" = "#E41A1C", "Convergent" = "#377EB8",
  "Symmetric" = "#4DAF4A", "Structured reversed" = "#FF7F00"
)
pattern_shapes <- c(
  "Progressive" = 16, "Convergent" = 17,
  "Symmetric" = 15, "Structured reversed" = 18
)

highlight <- c("dementia","diabetes","thyroid_disorder","sepsis","anemia",
  "heart_failure","ckd","peptic_ulcer","aspiration_pneumonia",
  "hematologic_malignancy","inflammatory_bowel_disease",
  "acute_renal_failure","osteoporosis","parkinson")
multi[, show_label := disease %in% highlight]

png("plot/Figures/Supplementary/FigS_pattern_space.png",
    width = 4200, height = 3000, res = 300)
p4 <- ggplot(multi, aes(x = per_step_move_pct, y = mean_asymmetry,
                         color = pattern, shape = pattern)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text(data = multi[show_label == TRUE],
            aes(label = label), size = 3.2, vjust = -0.8, show.legend = FALSE) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_shape_manual(values = pattern_shapes, name = "Pattern") +
  labs(title = "ELA Behavioral Patterns (deceased-patient transitions)",
       subtitle = sprintf("%d multi-basin diseases", nrow(multi)),
       x = "Per-step transition rate (%)",
       y = "Transition asymmetry index") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = c(0.85, 0.25),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank()
  )
print(p4)
dev.off()
cat("  Written: FigS_pattern_space.png\n")

cat("\nDone.\n")
