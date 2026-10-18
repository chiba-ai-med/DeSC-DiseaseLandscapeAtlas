## Comprehensive Basin exploration across 51 diseases
## Finds surprising/unexpected patterns for specialist consultation
##
## Usage: Rscript src/explore_basin_findings.R

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Comprehensive Basin Exploration ===\n\n")

diseases <- list.dirs("output", recursive = FALSE, full.names = FALSE)
diseases <- diseases[diseases != "hf_cancer"]

out_dir <- "summary/basin_exploration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Basin-level outcomes: death rate, hosp rate, size, energy
# ============================================================
cat("--- 1. Basin-level outcomes ---\n")
basin_stats <- list()
for (dis in diseases) {
  path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", dis)
  if (!file.exists(path)) next
  dt <- fread(path)
  if (nrow(dt) == 0 || !"basin_id" %in% names(dt)) next

  bs <- dt[, .(
    n_patients = .N,
    frac = .N / nrow(dt),
    death_rate = mean(death == "dead", na.rm = TRUE),
    hosp_rate = mean(hospitalized == "yes", na.rm = TRUE),
    mean_energy = mean(energy, na.rm = TRUE),
    mean_age = mean(age, na.rm = TRUE),
    n_basins = uniqueN(basin_id)
  ), by = basin_id]
  bs[, disease := dis]
  basin_stats[[dis]] <- bs
}
basin_all <- rbindlist(basin_stats)
fwrite(basin_all, file.path(out_dir, "basin_outcomes.tsv"), sep = "\t")

# --- Find basins with extreme death rates ---
cat("\n=== HIGH MORTALITY BASINS (death_rate > 0.15) ===\n")
high_mort <- basin_all[death_rate > 0.15 & n_patients > 100]
setorder(high_mort, -death_rate)
for (i in seq_len(min(30, nrow(high_mort)))) {
  r <- high_mort[i]
  cat(sprintf("  %s Basin%d: death=%.1f%%, n=%d (%.0f%% of cohort), age=%.1f, energy=%.3f\n",
      r$disease, r$basin_id, r$death_rate*100, r$n_patients, r$frac*100, r$mean_age, r$mean_energy))
}

# --- Find basins where death rate differs greatly between basins ---
cat("\n=== DISEASES WITH LARGEST INTER-BASIN MORTALITY GAP ===\n")
gap <- basin_all[, .(
  max_death = max(death_rate),
  min_death = min(death_rate),
  gap = max(death_rate) - min(death_rate),
  n_basins = uniqueN(basin_id)
), by = disease]
setorder(gap, -gap)
for (i in seq_len(min(20, nrow(gap)))) {
  r <- gap[i]
  cat(sprintf("  %s: %.1f%% - %.1f%% = gap %.1f%% (%d basins)\n",
      r$disease, r$max_death*100, r$min_death*100, r$gap*100, r$n_basins))
}

# --- Tiny but lethal basins ---
cat("\n=== SMALL BUT LETHAL BASINS (< 5% of patients, death > 10%) ===\n")
tiny_lethal <- basin_all[frac < 0.05 & death_rate > 0.10 & n_patients > 50]
setorder(tiny_lethal, -death_rate)
for (i in seq_len(min(20, nrow(tiny_lethal)))) {
  r <- tiny_lethal[i]
  cat(sprintf("  %s Basin%d: death=%.1f%%, n=%d (%.1f%% of cohort), age=%.1f\n",
      r$disease, r$basin_id, r$death_rate*100, r$n_patients, r$frac*100, r$mean_age))
}

# ============================================================
# 2. Cross-disease code enrichment patterns
# ============================================================
cat("\n--- 2. Cross-disease code enrichment ---\n")
enrich_list <- list()
for (dis in diseases) {
  path <- sprintf("plot/%s/basin_code_enrichment.tsv", dis)
  if (!file.exists(path)) next
  dt <- tryCatch(fread(path, quote = "\""), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) next
  dt[, disease := dis]
  # Keep top enriched codes per basin (enrichment > 2 and freq_in > 0.01)
  dt_top <- dt[enrichment > 2 & freq_in > 0.01]
  enrich_list[[dis]] <- dt_top[, .(disease, code, basin_id, enrichment, freq_in, freq_out, n_in)]
}
enrich_all <- rbindlist(enrich_list, fill = TRUE)
cat(sprintf("  %d enriched code-basin-disease combinations\n", nrow(enrich_all)))

# --- Codes enriched in unexpected diseases ---
# Find codes that appear as highly enriched in diseases where you wouldn't expect them
# e.g., psychiatric codes enriched in cancer, cancer codes enriched in psychiatric diseases
cat("\n=== CODES ENRICHED IN POTENTIALLY UNEXPECTED CONTEXTS ===\n")

# Define disease-code category mismatches
psych_codes <- "^F[0-4]"  # mental/behavioral
cancer_codes <- "^C"
cardio_codes <- "^I"
resp_codes <- "^J"
neuro_codes <- "^G"
renal_codes <- "^N1[7-9]"
metabolic_codes <- "^E[0-1]"

psych_diseases <- c("schizophrenia", "depression", "anxiety", "delirium", "sleep_disorders")
cancer_diseases <- c("lung_cancer", "stomach_cancer", "colorectal_cancer", "breast_cancer",
                      "liver_cancer", "pancreatic_cancer", "prostate_cancer",
                      "hematologic_malignancy", "malignancy", "metastatic_cancer")
cardio_diseases <- c("heart_failure", "ischemic_heart_disease", "atrial_fibrillation")
renal_diseases <- c("ckd", "acute_renal_failure")

# Psych codes in non-psych diseases
psych_in_nonpsych <- enrich_all[grepl(psych_codes, code) & !disease %in% psych_diseases & enrichment > 3]
setorder(psych_in_nonpsych, -enrichment)
if (nrow(psych_in_nonpsych) > 0) {
  cat("\n  Psychiatric codes (F0-F4) enriched in non-psychiatric diseases:\n")
  for (i in seq_len(min(15, nrow(psych_in_nonpsych)))) {
    r <- psych_in_nonpsych[i]
    cat(sprintf("    %s Basin%d: %s (enrichment=%.1fx, freq=%.1f%%)\n",
        r$disease, r$basin_id, r$code, r$enrichment, r$freq_in*100))
  }
}

# Cancer codes in non-cancer diseases
cancer_in_noncancer <- enrich_all[grepl(cancer_codes, code) & !disease %in% cancer_diseases & enrichment > 3]
setorder(cancer_in_noncancer, -enrichment)
if (nrow(cancer_in_noncancer) > 0) {
  cat("\n  Cancer codes (C) enriched in non-cancer diseases:\n")
  for (i in seq_len(min(15, nrow(cancer_in_noncancer)))) {
    r <- cancer_in_noncancer[i]
    cat(sprintf("    %s Basin%d: %s (enrichment=%.1fx, freq=%.1f%%)\n",
        r$disease, r$basin_id, r$code, r$enrichment, r$freq_in*100))
  }
}

# Renal codes in non-renal diseases
renal_in_nonrenal <- enrich_all[grepl(renal_codes, code) & !disease %in% renal_diseases & enrichment > 3]
setorder(renal_in_nonrenal, -enrichment)
if (nrow(renal_in_nonrenal) > 0) {
  cat("\n  Renal codes (N17-N19) enriched in non-renal diseases:\n")
  for (i in seq_len(min(15, nrow(renal_in_nonrenal)))) {
    r <- renal_in_nonrenal[i]
    cat(sprintf("    %s Basin%d: %s (enrichment=%.1fx, freq=%.1f%%)\n",
        r$disease, r$basin_id, r$code, r$enrichment, r$freq_in*100))
  }
}

# ============================================================
# 3. Basin similarity across diseases
# ============================================================
cat("\n--- 3. Cross-disease basin similarity ---\n")

# For each disease-basin, get top 20 enriched codes (by enrichment)
basin_signatures <- enrich_all[, .SD[order(-enrichment)][1:min(20, .N)], by = .(disease, basin_id)]
basin_signatures[, basin_key := paste(disease, basin_id, sep = "_B")]

# Compute Jaccard similarity between all basin pairs across diseases
keys <- unique(basin_signatures$basin_key)
# Only compare high-mortality basins to keep it manageable
high_mort_keys <- basin_all[death_rate > 0.08 & n_patients > 100, paste(disease, basin_id, sep = "_B")]
keys_hm <- intersect(keys, high_mort_keys)

cat(sprintf("  Comparing %d high-mortality basins\n", length(keys_hm)))

if (length(keys_hm) > 2 && length(keys_hm) <= 200) {
  sigs <- split(basin_signatures[basin_key %in% keys_hm]$code,
                basin_signatures[basin_key %in% keys_hm]$basin_key)

  sim_list <- list()
  for (i in seq_along(keys_hm)) {
    for (j in seq_len(i-1)) {
      a <- sigs[[keys_hm[i]]]
      b <- sigs[[keys_hm[j]]]
      if (is.null(a) || is.null(b)) next
      jaccard <- length(intersect(a, b)) / length(union(a, b))
      if (jaccard > 0.2) {
        sim_list[[length(sim_list)+1]] <- data.table(
          basin_a = keys_hm[i], basin_b = keys_hm[j],
          jaccard = jaccard, shared = length(intersect(a, b)))
      }
    }
  }

  if (length(sim_list) > 0) {
    sim_dt <- rbindlist(sim_list)
    setorder(sim_dt, -jaccard)
    # Filter out same-disease pairs
    sim_dt[, dis_a := sub("_B.*", "", basin_a)]
    sim_dt[, dis_b := sub("_B.*", "", basin_b)]
    cross_dis <- sim_dt[dis_a != dis_b]

    cat("\n=== SIMILAR BASINS ACROSS DIFFERENT DISEASES (Jaccard > 0.2) ===\n")
    for (i in seq_len(min(30, nrow(cross_dis)))) {
      r <- cross_dis[i]
      shared_codes <- intersect(sigs[[r$basin_a]], sigs[[r$basin_b]])
      cat(sprintf("  %s ~ %s: Jaccard=%.2f (%d shared codes: %s)\n",
          r$basin_a, r$basin_b, r$jaccard, r$shared,
          paste(head(shared_codes, 5), collapse = ", ")))
    }
    fwrite(sim_dt, file.path(out_dir, "basin_similarity.tsv"), sep = "\t")
  }
}

# ============================================================
# 4. Age paradox: basins where young patients die more
# ============================================================
cat("\n--- 4. Age paradox ---\n")
cat("\n=== BASINS WITH HIGH MORTALITY DESPITE LOW MEAN AGE ===\n")
age_paradox <- basin_all[death_rate > 0.05 & mean_age < 60 & n_patients > 100]
setorder(age_paradox, mean_age)
for (i in seq_len(min(15, nrow(age_paradox)))) {
  r <- age_paradox[i]
  cat(sprintf("  %s Basin%d: age=%.1f, death=%.1f%%, n=%d\n",
      r$disease, r$basin_id, r$mean_age, r$death_rate*100, r$n_patients))
}

# ============================================================
# 5. Hospitalization without death: basins with high hosp but low death
# ============================================================
cat("\n--- 5. High hospitalization, low death basins ---\n")
cat("\n=== BASINS WITH HIGH HOSP (>30%) BUT LOW DEATH (<5%) ===\n")
hosp_not_death <- basin_all[hosp_rate > 0.30 & death_rate < 0.05 & n_patients > 100]
setorder(hosp_not_death, -hosp_rate)
for (i in seq_len(min(15, nrow(hosp_not_death)))) {
  r <- hosp_not_death[i]
  cat(sprintf("  %s Basin%d: hosp=%.1f%%, death=%.1f%%, n=%d, age=%.1f\n",
      r$disease, r$basin_id, r$hosp_rate*100, r$death_rate*100, r$n_patients, r$mean_age))
}

# ============================================================
# Summary
# ============================================================
cat("\n\n=== SUMMARY ===\n")
cat(sprintf("  Diseases analyzed: %d\n", uniqueN(basin_all$disease)))
cat(sprintf("  Total basin-disease combinations: %d\n", nrow(basin_all)))
cat(sprintf("  High mortality basins (>15%%): %d\n", nrow(basin_all[death_rate > 0.15 & n_patients > 100])))
cat(sprintf("  Cross-disease enrichment entries: %d\n", nrow(enrich_all)))

cat("\nDone.\n")
