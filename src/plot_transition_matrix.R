library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_transition_matrix.R <transition.tsv> <disease_label> <pattern_label> <output.png>")
}

trans_file <- args[1]
disease_label <- args[2]
pattern_label <- args[3]
out_png <- args[4]

dt <- fread(trans_file)
dt <- dt[from_basin > 0 & to_basin > 0]

basins <- sort(unique(c(dt$from_basin, dt$to_basin)))
nb <- length(basins)

mat <- matrix(0, nb, nb, dimnames = list(basins, basins))
for (i in seq_len(nrow(dt))) {
  r <- as.character(dt$from_basin[i])
  cc <- as.character(dt$to_basin[i])
  mat[r, cc] <- dt$rate[i]
}

plot_df <- expand.grid(From = paste("Basin", basins), To = paste("Basin", basins))
plot_df$rate <- NA_real_
for (i in seq_len(nrow(plot_df))) {
  ri <- which(paste("Basin", basins) == plot_df$From[i])
  ci <- which(paste("Basin", basins) == plot_df$To[i])
  plot_df$rate[i] <- mat[ri, ci]
}

plot_df$From <- factor(plot_df$From, levels = rev(paste("Basin", basins)))
plot_df$To <- factor(plot_df$To, levels = paste("Basin", basins))

plot_df$pct <- sprintf("%.1f%%", plot_df$rate * 100)

plot_df$sq_size <- sqrt(pmax(plot_df$rate, 0.005))
max_sq <- max(plot_df$sq_size)
plot_df$sq_size <- plot_df$sq_size / max_sq * 0.95

plot_df$text_col <- ifelse(plot_df$rate > 0.15, "white", "grey30")

base_blue <- "#4878A8"

p <- ggplot(plot_df, aes(x = To, y = From)) +
  geom_tile(fill = "#E8E8F0", colour = "white", linewidth = 2) +
  geom_tile(aes(width = sq_size, height = sq_size, fill = rate),
            colour = NA) +
  scale_fill_gradient(low = "#B0C4DE", high = base_blue, guide = "none",
                      limits = c(0, 1)) +
  geom_text(aes(label = pct, colour = text_col), size = if (nb <= 2) 8 else 5,
            fontface = "bold") +
  scale_colour_identity() +
  labs(title = disease_label, subtitle = pattern_label, x = "To", y = "From") +
  theme_minimal(base_size = if (nb <= 2) 18 else 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = if (nb <= 2) 24 else 18),
    plot.subtitle = element_text(hjust = 0.5, colour = "red", face = "italic",
                                 size = if (nb <= 2) 18 else 14),
    axis.title = element_text(size = if (nb <= 2) 16 else 12),
    axis.text = element_text(size = if (nb <= 2) 14 else 11),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  coord_fixed()

w <- if (nb <= 2) 8 else 10
h <- if (nb <= 2) 8 else 8
ggsave(out_png, p, width = w, height = h, dpi = 300, bg = "white")
cat("Saved:", out_png, "\n")
