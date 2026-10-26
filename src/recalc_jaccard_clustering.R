library(data.table)

setwd("/home/koki/dev/DeSC-DiseaseLandscapeAtlas")
diseases <- fread("config/disease_families.tsv")$family_id
diseases <- diseases[diseases != "hf_cancer"]

# Read top-30 codes per basin per disease
top30 <- list()
for (d in diseases) {
  f <- file.path("plot", d, "basin_code_enrichment.tsv")
  if (!file.exists(f)) { cat("MISSING:", d, "\n"); next }
  dt <- fread(f)
  for (b in unique(dt$basin_id)) {
    sub <- dt[basin_id == b][order(-freq_in)][1:min(30, .N)]
    key <- paste0(d, ":", b)
    top30[[key]] <- sub$code
  }
}
cat(length(top30), "basin code sets loaded\n")

# Basin-level Jaccard
keys <- names(top30)
rows <- list()
for (i in seq_along(keys)) {
  for (j in seq_len(i - 1)) {
    s1 <- top30[[keys[i]]]
    s2 <- top30[[keys[j]]]
    inter <- length(intersect(s1, s2))
    uni <- length(union(s1, s2))
    jac <- inter / uni
    parts_i <- strsplit(keys[i], ":")[[1]]
    parts_j <- strsplit(keys[j], ":")[[1]]
    rows[[length(rows) + 1]] <- data.table(
      disease1 = parts_i[1], basin1 = as.integer(parts_i[2]),
      disease2 = parts_j[1], basin2 = as.integer(parts_j[2]),
      intersection = inter, union = uni, jaccard = jac
    )
  }
}
jac_basin <- rbindlist(rows)
fwrite(jac_basin, "summary/jaccard_all_basins.tsv", sep = "\t")
cat("Basin-level Jaccard:", nrow(jac_basin), "pairs\n")

# Disease-level max Jaccard (exclude schizophrenia as in manuscript)
diseases_50 <- diseases[diseases != "schizophrenia"]
n <- length(diseases_50)
dmat <- matrix(0, n, n, dimnames = list(diseases_50, diseases_50))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) { dmat[i, j] <- 1; next }
    keys_i <- grep(paste0("^", diseases_50[i], ":"), names(top30), value = TRUE)
    keys_j <- grep(paste0("^", diseases_50[j], ":"), names(top30), value = TRUE)
    max_jac <- 0
    for (ki in keys_i) {
      for (kj in keys_j) {
        s1 <- top30[[ki]]; s2 <- top30[[kj]]
        jac <- length(intersect(s1, s2)) / length(union(s1, s2))
        if (jac > max_jac) max_jac <- jac
      }
    }
    dmat[i, j] <- max_jac
  }
}

# Clustering
dist_mat <- as.dist(1 - dmat)
hc <- hclust(dist_mat, method = "average")

# Optimal leaf ordering
if (requireNamespace("seriation", quietly = TRUE)) {
  library(seriation)
  olo <- seriate(dist_mat, method = "OLO", control = list(hclust = hc))
  hc <- olo[[1]]$hclust
  if (is.null(hc)) hc <- hclust(dist_mat, method = "average")
}

# Cut tree at 2 groups
grp2 <- cutree(hc, k = 2)
for (g in 1:2) {
  members <- names(grp2[grp2 == g])
  cat(sprintf("\nGroup %d (%d diseases): %s\n", g, length(members), paste(members, collapse = ", ")))
}

# Save disease-level max Jaccard matrix
jac_disease <- data.table(disease = diseases_50)
for (d in diseases_50) jac_disease[[d]] <- dmat[, d]
fwrite(jac_disease, "summary/jaccard_disease_maxjac.tsv", sep = "\t")

# Save dendrogram order for Figure 4
cat("\nDendrogram leaf order:\n")
cat(paste(hc$labels[hc$order], collapse = "\n"), "\n")

# Generate Figure 4 PNG
png("plot/Figures/Fig04/Fig04_jaccard_clustering.png", width = 1200, height = 800, res = 150)
par(mar = c(4, 2, 2, 12))
plot(as.dendrogram(hc), horiz = TRUE, main = "Disease clustering by basin code similarity",
     xlab = "1 - max Jaccard", cex = 0.6)
dev.off()
cat("\nFigure 4 saved.\n")

# Generate Figure S2 (2-block annotation)
acute_frailty <- names(grp2[grp2 == grp2["acute_renal_failure"]])
chronic <- names(grp2[grp2 != grp2["acute_renal_failure"]])
cat(sprintf("\nAcute/Frailty (%d): %s\n", length(acute_frailty), paste(acute_frailty, collapse = ", ")))
cat(sprintf("Chronic/Lifestyle (%d): %s\n", length(chronic), paste(chronic, collapse = ", ")))

png("plot/Figures/Supplementary/FigS_jaccard_2blocks.png", width = 1400, height = 900, res = 150)
par(mar = c(4, 2, 2, 14))
dend <- as.dendrogram(hc)
leaf_cols <- ifelse(hc$labels[hc$order] %in% acute_frailty, "red", "steelblue")
plot(dend, horiz = TRUE, main = "Jaccard clustering with Acute/Frailty annotation",
     xlab = "1 - max Jaccard", cex = 0.5)
dev.off()
cat("Figure S2 saved.\n")
