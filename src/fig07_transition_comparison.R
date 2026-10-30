suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

diseases <- c(
  "dementia"      = "Dementia\n(Symmetric)",
  "heart_failure" = "Heart Failure\n(Convergent)",
  "diabetes"      = "Diabetes\n(Progressive)"
)

all_dt <- list()
for (d in names(diseases)) {
  tsv <- file.path(base, sprintf("plot/%s/basin_transition_matrix.tsv", d))
  dt <- fread(tsv)
  # Keep only basin-to-basin (exclude No basin)
  dt <- dt[grepl("^Basin", from) & grepl("^Basin", to)]
  dt[, disease := diseases[d]]
  all_dt[[d]] <- dt
}
combined <- rbindlist(all_dt)
combined[, disease := factor(disease, levels = unname(diseases))]

sz <- 16

g <- ggplot(combined, aes(x = to, y = from, fill = prob)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.1f%%", prob * 100)),
            size = sz * 1.2 / .pt, fontface = "bold",
            color = ifelse(combined$prob > 0.5, "white", "black")) +
  scale_y_discrete(limits = rev) +
  facet_wrap(~disease, scales = "free", nrow = 1) +
  scale_fill_gradient2(low = "#F7F7F7", mid = "#6BAED6", high = "#08306B",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Transition\nprobability",
                       labels = scales::percent) +
  labs(title = "Basin transition probabilities across ELA patterns",
       x = "To", y = "From") +
  theme_minimal(base_size = sz) +
  theme(
    plot.title = element_text(face = "bold", size = sz * 1.5),
    axis.text = element_text(size = sz * 1.1),
    axis.title = element_text(size = sz * 1.2),
    strip.text = element_text(size = sz * 1.5, face = "bold"),
    legend.text = element_text(size = sz * 0.9),
    legend.title = element_text(size = sz * 0.9),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

out <- file.path(base, "plot/Figures/Fig07/Fig07_comparison.png")
ggsave(out, g, width = 22 * 0.729, height = 5.5, dpi = 300, bg = "white")
cat(sprintf("Written: %s\n", out))
