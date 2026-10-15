args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
out_values <- args[2]
out_vectors <- args[3]

cov_mat <- read.csv(infile, header = FALSE)
out <- eigen(as.matrix(cov_mat))
write.table(out$values, file = out_values,
            row.names = FALSE, col.names = FALSE, sep = ",")
write.table(out$vectors, file = out_vectors,
            row.names = FALSE, col.names = FALSE, sep = ",")
