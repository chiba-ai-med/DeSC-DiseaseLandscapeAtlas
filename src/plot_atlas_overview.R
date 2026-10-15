## Atlas overview: 3-variable visualization
## death_rate x rows/patient x HR, colored by PASS/FAIL
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Load summary
sm <- fread("summary/atlas_summary.tsv")

# Get rows/patient and death rate from patient scores
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

sm[, status := ifelse(direction_ok, "PASS", "FAIL")]
sm[, label := disease_family]

# --- Plot 1: HR vs death_rate ---
g1 <- ggplot(sm, aes(x = death_rate, y = HR, color = status, label = label)) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray40") +
  geom_point(size = 3) +
  geom_text(size = 2.2, vjust = -1, hjust = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c(PASS = "#2196F3", FAIL = "#F44336")) +
  labs(title = "Energy HR vs Death Rate by Disease Family",
       x = "Death rate in cohort", y = "HR per +1 SD Energy",
       color = NULL) +
  theme_minimal(base_size = 12)
ggsave("plot/atlas_hr_vs_deathrate.png", g1, width = 12, height = 8, dpi = 150)

# --- Plot 2: HR vs rows/patient ---
g2 <- ggplot(sm, aes(x = rows_per_patient, y = HR, color = status, label = label)) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray40") +
  geom_point(size = 3) +
  geom_text(size = 2.2, vjust = -1, hjust = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c(PASS = "#2196F3", FAIL = "#F44336")) +
  labs(title = "Energy HR vs Observation Density by Disease Family",
       x = "Rows per patient (observation density)", y = "HR per +1 SD Energy",
       color = NULL) +
  theme_minimal(base_size = 12)
ggsave("plot/atlas_hr_vs_rowsperpt.png", g2, width = 12, height = 8, dpi = 150)

# --- Plot 3: death_rate vs rows/patient, sized by |log(HR)| ---
sm[, hr_effect := abs(log(HR))]
g3 <- ggplot(sm, aes(x = rows_per_patient, y = death_rate,
                       color = status, size = hr_effect, label = label)) +
  geom_point(alpha = 0.7) +
  geom_text(size = 2.2, vjust = -1.5, show.legend = FALSE) +
  scale_color_manual(values = c(PASS = "#2196F3", FAIL = "#F44336")) +
  scale_size_continuous(range = c(1, 8), name = "|log(HR)|") +
  labs(title = "Disease Family Atlas: Observation Density x Death Rate x Energy HR",
       x = "Rows per patient", y = "Death rate", color = NULL) +
  theme_minimal(base_size = 12)
ggsave("plot/atlas_3var.png", g3, width = 13, height = 9, dpi = 150)

# --- Plot 4: HR bar chart sorted ---
sm_sorted <- sm[order(HR)]
sm_sorted[, label := factor(label, levels = label)]
g4 <- ggplot(sm_sorted, aes(x = HR, y = label, fill = status)) +
  geom_col() +
  geom_vline(xintercept = 1, linetype = 2, color = "gray20") +
  scale_fill_manual(values = c(PASS = "#2196F3", FAIL = "#F44336")) +
  labs(title = "Energy Cox Sanity Check: All 51 Disease Families",
       subtitle = "HR < 1 = low energy = high risk (PASS)",
       x = "HR per +1 SD Energy", y = NULL, fill = NULL) +
  theme_minimal(base_size = 10)
ggsave("plot/atlas_hr_barplot.png", g4, width = 10, height = 14, dpi = 150)

cat("Done. Plots in plot/\n")
