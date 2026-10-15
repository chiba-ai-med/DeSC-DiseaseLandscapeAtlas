## Compare all-cause vs disease-specific hospitalization Cox HR
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

fids <- list.dirs("output", recursive = FALSE, full.names = FALSE)

rows <- list()
for (fid in fids) {
  all_f <- sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check_hosp.tsv", fid)
  spec_f <- sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check_hosp_specific.tsv", fid)
  if (!file.exists(spec_f)) next
  if (file.exists(all_f)) {
    a <- fread(all_f)
    rows[[length(rows) + 1]] <- data.table(
      disease = fid, type = "All-cause", HR = a$HR,
      HR_lower = a$HR_lower, HR_upper = a$HR_upper,
      status = ifelse(a$direction_ok, "PASS", "FAIL"))
  }
  s <- fread(spec_f)
  rows[[length(rows) + 1]] <- data.table(
    disease = fid, type = "Disease-specific", HR = s$HR,
    HR_lower = s$HR_lower, HR_upper = s$HR_upper,
    status = ifelse(s$direction_ok, "PASS", "FAIL"))
}

dt <- rbindlist(rows)
dt[, type := factor(type, levels = c("All-cause", "Disease-specific"))]

# Sort by disease-specific HR
spec_order <- dt[type == "Disease-specific"][order(HR), disease]
dt[, disease := factor(disease, levels = spec_order)]

# --- Bar chart ---
g1 <- ggplot(dt, aes(x = HR, y = disease, fill = interaction(type, status))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray20") +
  scale_fill_manual(
    values = c("All-cause.PASS" = "#90CAF9", "All-cause.FAIL" = "#EF9A9A",
               "Disease-specific.PASS" = "#1565C0", "Disease-specific.FAIL" = "#C62828"),
    labels = c("All-cause PASS", "All-cause FAIL",
               "Specific PASS", "Specific FAIL")) +
  labs(title = "Hospitalization Cox: All-cause vs Disease-specific (primary diagnosis)",
       subtitle = sprintf("PASS rate: All-cause %d/51 -> Disease-specific %d/51",
                          nrow(dt[type == "All-cause" & status == "PASS"]),
                          nrow(dt[type == "Disease-specific" & status == "PASS"])),
       x = "HR per +1 SD Energy", y = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave("plot/hosp_comparison_barplot.png", g1, width = 11, height = 16, dpi = 150)

# --- Summary table ---
summary_dt <- dcast(dt, disease ~ type, value.var = c("HR", "status"))
setnames(summary_dt, c("disease",
  "HR_allcause", "HR_specific",
  "status_allcause", "status_specific"))
summary_dt <- summary_dt[order(HR_specific)]
fwrite(summary_dt, "summary/hosp_comparison.tsv", sep = "\t")

cat("Done.\n")
cat(sprintf("  plot/hosp_comparison_barplot.png\n"))
cat(sprintf("  summary/hosp_comparison.tsv\n"))
