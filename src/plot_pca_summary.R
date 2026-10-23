library(data.table)
library(ggplot2)

outdir <- "plot/Figures"

diseases <- fread("config/disease_families.tsv")[["family_id"]]
diseases <- diseases[diseases != "hf_cancer"]

total_var_df <- fread("/tmp/total_variance_51.tsv")

## --- Compute PC variances and correlations ---
var_list <- vector("list", length(diseases))
cor_list <- vector("list", length(diseases))

for (i in seq_along(diseases)) {
  d <- diseases[i]
  f <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", d)
  if (!file.exists(f)) next
  tv_row <- total_var_df[disease == d]
  if (nrow(tv_row) == 0) next
  tv <- tv_row$total_var

  dt <- fread(f, select = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","age","sex"))

  vars <- vapply(dt[, paste0("PC",1:7), with=FALSE], var, numeric(1), na.rm = TRUE)
  var_list[[i]] <- data.frame(
    disease = d, PC = 1:7,
    true_cum_pct = cumsum(as.numeric(vars)) / tv * 100,
    stringsAsFactors = FALSE
  )

  age_r <- vapply(paste0("PC",1:7), function(pc) cor(dt[[pc]], dt$age, use="complete.obs"), numeric(1))
  sex_num <- ifelse(dt$sex == "male", 1, 0)
  sex_r <- vapply(paste0("PC",1:7), function(pc) cor(dt[[pc]], sex_num, use="complete.obs"), numeric(1))
  cor_list[[i]] <- data.frame(
    disease = d, PC = 1:7,
    cor_age = as.numeric(age_r),
    cor_sex = as.numeric(sex_r),
    stringsAsFactors = FALSE
  )
}

var_df <- do.call(rbind, var_list[!vapply(var_list, is.null, logical(1))])
cor_df <- do.call(rbind, cor_list[!vapply(cor_list, is.null, logical(1))])

n_diseases <- length(unique(var_df$disease))
med7 <- median(var_df$true_cum_pct[var_df$PC == 7])

## --- Cumulative variance plot (true % of total) ---
p_var <- ggplot(var_df, aes(x = PC, y = true_cum_pct, group = disease)) +
  geom_line(alpha = 0.15, colour = "grey40") +
  geom_point(alpha = 0.15, size = 0.8, colour = "grey40") +
  stat_summary(aes(group = 1), fun = median, geom = "line",
               colour = "red", linewidth = 1.2) +
  stat_summary(aes(group = 1), fun = median, geom = "point",
               colour = "red", size = 3) +
  scale_x_continuous(breaks = 1:7, labels = paste0("PC", 1:7)) +
  scale_y_continuous(limits = c(0, max(var_df$true_cum_pct) * 1.1),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Cumulative variance explained (%)",
       title = sprintf("Cumulative variance explained (%d diseases)", n_diseases),
       subtitle = sprintf("Median at PC7: %.1f%%", med7)) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14))

ggsave(file.path(outdir, "Fig02c_cumulative_variance.png"), p_var,
       width = 7, height = 6, dpi = 300)
cat(sprintf("Cumulative variance plot: %d diseases, median PC7=%.1f%%\n", n_diseases, med7))

## --- PC-correlation heatmap (vertical, Age left + Sex right) ---
sort_order <- cor_df[cor_df$PC == 1, c("disease","cor_age")]
sort_order <- sort_order[order(sort_order$cor_age), ]
cor_df$disease <- factor(cor_df$disease, levels = sort_order$disease)

age_long <- cor_df[, c("disease","PC","cor_age")]
age_long$panel <- "Age"
names(age_long)[3] <- "r"

sex_long <- cor_df[, c("disease","PC","cor_sex")]
sex_long$panel <- "Sex"
names(sex_long)[3] <- "r"

hm_df <- rbind(age_long, sex_long)
hm_df$panel <- factor(hm_df$panel, levels = c("Age", "Sex"))

p_hm <- ggplot(hm_df, aes(x = factor(PC), y = disease, fill = r)) +
  geom_tile() +
  facet_wrap(~panel) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-0.55, 0.55),
                       oob = scales::squish, name = "r") +
  scale_x_discrete(labels = paste0("PC", 1:7)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 20, face = "bold"),
        panel.grid = element_blank())

ggsave(file.path(outdir, "Fig02d_pc_correlation_vertical.png"), p_hm,
       width = 12, height = 18, dpi = 300)

cat("Done.\n")
