setwd("path-to-dir")
library(data.table)
library(ggplot2)

# Singl-cell data can be obtained from the following sources:
# Emont, M.P., Jacobs, C., Essene, A.L. et al. A single-cell atlas of human and mouse white adipose tissue. Nature 603, 926–933 (2022). https://doi.org/10.1038/s41586-022-04518-2
# Wirka, R.C., Wagh, D., Paik, D.T. et al. Atheroprotective roles of smooth muscle cell phenotypic modulation and the TCF21 disease gene as revealed by single-cell analysis. Nat Med 25, 1280–1289 (2019). https://doi.org/10.1038/s41591-019-0512-5

genename <- "COL6A3"

# Read metadata and UMAP data
metainfo <- as.data.frame(fread("path-to-metadata-file"))
metainfo <- metainfo[-1, ]
rownames(metainfo) <- metainfo$NAME

umap <- as.data.frame(fread("path-to-umap-file"))
umap <- umap[-1, ]
rownames(umap) <- umap$NAME

metainfo <- metainfo[rownames(umap), ]

# Read data from the first dataset
barcode1 <- as.data.frame(fread("path-to-barcode-file-1", header = FALSE))
gene1 <- as.data.frame(fread("path-to-gene-file-1", header = FALSE))
exp1 <- as.data.frame(fread("path-to-expression-file-1"))  # This file was zipped to save space

# Read data from the second dataset
barcode2 <- as.data.frame(fread("path-to-barcode-file-2", header = FALSE))
gene2 <- as.data.frame(fread("path-to-gene-file-2", header = FALSE))
exp2 <- as.data.frame(fread("path-to-expression-file-2"))  # This file was zipped to save space

# Select expression data for the gene of interest
idx1 <- which(exp1$`%%MatrixMarket` == which(gene1$V2 == genename))
selectexp1 <- exp1[idx1, ]

idx2 <- which(exp2$`%%MatrixMarket` == which(gene2$V2 == genename))
selectexp2 <- exp2[idx2, ]

# Set row names and combine datasets
rownames(selectexp1) <- barcode1$V1[selectexp1$matrix]
rownames(selectexp2) <- barcode2$V1[selectexp2$matrix]
selectexp <- rbind.data.frame(selectexp1, selectexp2)

# Prepare plotting data
pltdat <- metainfo
pltdat$gene_count <- 0
pltdat[rownames(selectexp), ]$gene_count <- selectexp$coordinate
pltdat <- pltdat[rownames(umap), ]
pltdat$umap_x <- as.numeric(umap$X)
pltdat$umap_y <- as.numeric(umap$Y)

# Simplify cluster names
pltdat$cluster <- factor(pltdat$cluster)

# Simplify the ggplot commands
# First plot: Clusters on UMAP
ggplot(pltdat, aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP dimension 1", y = "UMAP dimension 2", color = "Cluster") +
  theme_minimal() -> plt1
ggsave("path-to-dir/WAT_cluster.pdf", plt1, height = 8, width = 9.5)

# Second plot: Gene expression on UMAP
ggplot(pltdat, aes(x = umap_x, y = umap_y, color = gene_count)) +
  geom_point(size = 0.5) +
  labs(
    x = "UMAP dimension 1",
    y = "UMAP dimension 2",
    color = paste(genename, "Expression")
  ) +
  theme_minimal() +
  scale_color_viridis_c() -> plt2
ggsave("path-to-dir/WAT_gene_expression.pdf", plt2, height = 8, width = 10)

# Statistical tests for enrichment
permutefun <- function(replicate, cell, label, expression) {
  nexp_cell <- sum(label == cell & expression > 0)
  nexp_cell_permute <- numeric(replicate)
  for (ite in 1:replicate) {
    permute_expression <- sample(expression, length(expression), replace = FALSE)
    nexp_cell_permute[ite] <- sum(label == cell & permute_expression > 0)
  }
  p_value <- (sum(nexp_cell_permute >= nexp_cell) + 1) / (replicate + 1)  # Continuity correction
  cat(
    cell, "| expression in", nexp_cell, "/", sum(label == cell),
    "| prop. expression", round(nexp_cell / sum(label == cell), 3),
    "| mean exp level", round(mean(expression[label == cell]), 3), "vs. overall",
    round(mean(expression), 3), "| p-value =", p_value, "\n"
  )
}

for (celltype in sort(unique(pltdat$cluster))) {
  permutefun(1000, celltype, pltdat$cluster, pltdat$gene_count)
}
