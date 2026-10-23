## Atlas overview: 3-variable visualization
## death_rate x rows/patient, colored by per-basin majority HR<1
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

diseases <- fread("config/disease_families.tsv")[family_id != "hf_cancer", family_id]

# Build per-disease summary from per-basin Cox results
build_summary <- function(diseases, outcome_suffix = "") {
  cox_file <- if (outcome_suffix == "")
    "energy_sanity_check_per_basin.tsv"
  else
    sprintf("energy_sanity_check_per_basin_%s.tsv", outcome_suffix)

  res <- rbindlist(lapply(diseases, function(d) {
    f <- sprintf("output/%s/all_comorbid/cox/7/%s", d, cox_file)
    if (!file.exists(f)) return(NULL)
    pb <- fread(f)
    n_basins <- nrow(pb)
    n_lt1 <- sum(pb$direction_ok)
    median_hr <- median(pb$HR)
    data.table(disease = d, n_basins = n_basins, n_lt1 = n_lt1,
               majority_ok = n_lt1 >= n_basins / 2,
               median_hr = median_hr)
  }))
  res
}

mort <- build_summary(diseases, "")
hosp <- build_summary(diseases, "hosp")

# Get rows/patient and death rate
sm <- fread("summary/atlas_summary.tsv")
sm[, rows_per_patient := NA_real_]
sm[, death_rate := NA_real_]
for (i in seq_len(nrow(sm))) {
  fid <- sm$disease_family[i]
  stats <- sprintf("output/%s/all_comorbid/cohort/matrix_stats.tsv", fid)
  pt_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
  if (file.exists(stats)) {
    st <- fread(stats, header = FALSE, sep = "\t")
    n_rows <- as.numeric(st[V1 == "n_patients", V2])
    n_pat <- sm$n_patients[i]
    if (n_pat > 0) sm$rows_per_patient[i] <- n_rows / n_pat
  }
  if (file.exists(pt_path)) {
    pt <- fread(pt_path, select = c("death"))
    sm$death_rate[i] <- mean(pt$death == "dead", na.rm = TRUE)
  }
}
meta <- sm[, .(disease = disease_family, rows_per_patient, death_rate)]

plot_atlas <- function(pb, meta, title_outcome, outfile) {
  df <- merge(pb, meta, by = "disease")
  df[, status := ifelse(majority_ok, "HR < 1", "HR ≥ 1")]
  df[, label := gsub("_", " ", disease)]
  df[, hr_effect := abs(log(median_hr))]

  n_ok <- sum(df$majority_ok)
  n_total <- nrow(df)

  g <- ggplot(df, aes(x = rows_per_patient, y = death_rate,
                       color = status, size = hr_effect, label = label)) +
    geom_point(alpha = 0.7) +
    geom_text_repel(size = 3.5, max.overlaps = 50, show.legend = FALSE,
                    segment.color = "grey60", segment.size = 0.3) +
    scale_color_manual(values = c("HR < 1" = "#2196F3", "HR ≥ 1" = "#F44336")) +
    scale_size_continuous(range = c(2, 10), name = "|log(HR)|") +
    labs(title = sprintf("51 Disease Families: %s (HR < 1 in %d/%d, per-basin majority)",
                          title_outcome, n_ok, n_total),
         x = "Rows per patient (observation density)", y = "Death rate",
         color = "Energy Cox") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 18, face = "bold"),
          legend.position = "right")

  ggsave(outfile, g, width = 14, height = 10, dpi = 300)
  cat(sprintf("%s: HR<1 in %d/%d (per-basin majority)\n", title_outcome, n_ok, n_total))
}

plot_atlas(mort, meta, "Mortality",
           "plot/Figures/Fig03_atlas_overview.png")
plot_atlas(hosp, meta, "Hospitalization",
           "plot/Figures/FigS_atlas_hospitalization.png")

cat("Done.\n")
