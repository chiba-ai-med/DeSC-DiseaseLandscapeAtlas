suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

tsv <- file.path(base, "plot/diabetes/tagcloud_scored/logitdiff/basin_code_logitdiff.tsv")
dt <- fread(tsv)

ps <- fread(file.path(base, "output/diabetes/all_comorbid/patient/patient_scores.tsv"),
            select = "basin_id")
bpat <- ps[basin_id > 0, .(n_patients = .N), by = basin_id]

n1 <- bpat[basin_id == 1, n_patients]
n2 <- bpat[basin_id == 2, n_patients]

top_n <- 10
b1 <- dt[basin_id == 1 & freq_in > 0.01][order(-score)][1:top_n]
b2 <- dt[basin_id == 2 & freq_in > 0.01][order(-score)][1:top_n]

truncate_label <- function(s, maxlen = 40) {
  ifelse(nchar(s) > maxlen, paste0(substr(s, 1, maxlen - 1), "…"), s)
}

b1[, label := paste0(truncate_label(en_label), " (", code, ")")]
b2[, label := paste0(truncate_label(en_label), " (", code, ")")]

b1[, basin := sprintf("Basin 1  HR=0.84 (n=%s)", format(n1, big.mark = ","))]
b2[, basin := sprintf("Basin 2  HR=1.17 (n=%s)", format(n2, big.mark = ","))]

b1[, label := factor(label, levels = rev(label))]
b2[, label := factor(label, levels = rev(label))]

combined <- rbindlist(list(b1, b2))
combined[, basin := factor(basin, levels = c(
  sprintf("Basin 1  HR=0.84 (n=%s)", format(n1, big.mark = ",")),
  sprintf("Basin 2  HR=1.17 (n=%s)", format(n2, big.mark = ","))))]

sz <- 16

g <- ggplot(combined, aes(x = label, y = score, fill = freq_in)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", freq_in * 100)),
            hjust = -0.1, size = sz / 3.5, fontface = "bold") +
  coord_flip() +
  facet_wrap(~basin, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", name = "Freq",
                       begin = 0.15, end = 0.95,
                       labels = scales::percent) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Diabetes - Top 10 comorbidity codes by basin",
       x = NULL, y = "Log-odds ratio") +
  theme_bw(base_size = sz) +
  theme(plot.title = element_text(face = "bold", size = sz),
        axis.text.y = element_text(size = sz * 0.9),
        axis.text.x = element_text(size = sz),
        axis.title.x = element_text(size = sz),
        strip.text = element_text(size = sz * 0.75, face = "bold"),
        legend.text = element_text(size = sz - 2),
        legend.title = element_text(size = sz - 2),
        plot.background = element_rect(fill = "white", color = NA))

fig_h <- (1.5 + 0.7 * top_n) * 0.8
out <- file.path(base, "plot/Figures/Fig06/Fig06a_diabetes_basin1_vs_2.png")
ggsave(out, g, width = 16 * 0.8, height = fig_h, dpi = 300, bg = "white")
cat(sprintf("Written: %s\n", out))
