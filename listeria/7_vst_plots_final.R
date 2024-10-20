library(ggplot2)
library(ggpubr)
library(DESeq2)
library(umap)
library(dendextend)
library(pvclust)
library(corrplot)
library(ggrepel)
library(viridis)
library(khroma)
library(gridExtra)

# set seed
set.seed(42)

################ prepare input data

# load raw counts for protein-coding genes
gene_counts <- read.delim(file = "gene_counts_pc.tsv",
                        header = T, sep = "\t",
                        stringsAsFactors = F,
                        check.names = F)
head(gene_counts)

# separate metadata and counts
counts_tab <- gene_counts[, grepl(pattern = "_trial", x = colnames(gene_counts))]
rownames(counts_tab) <- gene_counts$locus_tag
head(counts_tab)

meta_tab <- gene_counts[, !(grepl(pattern = "_trial", x = colnames(gene_counts)))]
head(meta_tab)

# make data frame of conditions
cond <- data.frame(sample = colnames(counts_tab))

# make trial column
cond$trial <- unlist(lapply(strsplit(x = cond$sample, split = "_"), "[[", 2))
cond
table(cond$trial)

# make group column
cond$group <- unlist(lapply(strsplit(x = cond$sample, split = "_"), "[[", 3))
cond
table(cond$group)

# make ultrasound column
cond$ultrasound <- ifelse(test = grepl(pattern = "ultrasound", x = cond$group),
                          yes = "ultrasound", no = "not_ultrasound")
cond
table(cond$ultrasound)

# make chlorine column
cond$chlorine <- ifelse(test = grepl(pattern = "chlorine", x = cond$group),
                          yes = "chlorine", no = "not_chlorine")
cond
table(cond$chlorine)

# make water column
cond$water <- ifelse(test = grepl(pattern = "water", x = cond$group),
                        yes = "water", no = "not_water")
cond
table(cond$water)

# make sequencing depth column
table(colnames(counts_tab)%in%cond$sample)
table(cond$sample%in%colnames(counts_tab))
table(colnames(counts_tab)==cond$sample)
table(cond$sample==colnames(counts_tab))

cond$seqdepth <- colSums(counts_tab)
cond


# note: deseq need column names here to match rownames of count table.
# they must also be in the same order
rownames(cond) <- cond$sample
table(colnames(counts_tab)==rownames(cond))
table(rownames(cond)==colnames(counts_tab))

################ get variance stabilized estimate

# create DESeq data set object
dds <- DESeqDataSetFromMatrix(countData = counts_tab,
                              colData = cond,
                              design= ~ 1) # initialize with no design

# run DESeq 
dds <- DESeq(dds)

# apply variance stabilizing transformation (vst)
vsd <- vst(dds)

# store vst-transformed counts as variable foo
foo <- vsd@assays@data[[1]]

# sanity check: column names of vst-transformed counts match sample order in conditions data frame
table(colnames(foo)==cond$sample)
table(cond$sample==colnames(foo))

################ hierarchical clustering of vst-transformed counts

# hierarchical clustering of vst-transformed counts
h <- pvclust(data = foo, method.hclust = "average", method.dist = "euclidean", nboot = 1000)
# plot basic dendrogram
plot(h)

# make it nicer
# par(mar = c(bottom, left, top, right))
pdf(file = "sample_hclust.pdf", width = 8.5, height = 11)
par(mar = c(15, 3, 1, 1))
h.dend <-
  as.dendrogram(h) %>%
  set("branches_lwd", 3) %>% 
  pvclust_show_signif(h) %>% 
  set("by_labels_branches_col",
      value = cond$sample[which(cond$group=="water")],
      TF_value = "#66CCEE") %>%
  set("by_labels_branches_col",
      value = cond$sample[which(cond$group=="water-ultrasound")],
      TF_value = "#AA3377") %>%
  set("by_labels_branches_col",
      value = cond$sample[which(cond$group=="chlorine-ultrasound")],
      TF_value = "#CCBB44") %>%
  set("by_labels_branches_col",
      value = cond$sample[which(cond$group=="chlorine")],
      TF_value = "#228833") %>%
  plot(horiz = F)

h %>% text()
dev.off()

# test if linkage method affects dendrogram
dend1 <- as.dendrogram(hclust(d = dist(t(foo)), method = "average"))
dend2 <- as.dendrogram(hclust(d = dist(t(foo)), method = "complete"))
dend3 <- as.dendrogram(hclust(d = dist(t(foo)), method = "centroid"))
dend4 <- as.dendrogram(hclust(d = dist(t(foo)), method = "median"))
dend5 <- as.dendrogram(hclust(d = dist(t(foo)), method = "single"))
dend6 <- as.dendrogram(hclust(d = dist(t(foo)), method = "ward.D"))
dend7 <- as.dendrogram(hclust(d = dist(t(foo)), method = "ward.D2"))
dend8 <- as.dendrogram(hclust(d = dist(t(foo)), method = "mcquitty"))

dend1to8 <- dendlist("Average" = dend1, "Complete" = dend2,
                     "Centroid" = dend3, "Median" = dend4,
                     "Single" = dend5, "Ward" = dend6,
                     "Ward.D2" = dend7, "McQuitty" = dend8)

pdf(file = "sample_hclust_corrplot.pdf", width = 8.5, height = 8.5)
par(mfrow=c(2,2))
corrplot(cor.dendlist(dend1to8), method = "pie", type = "lower")
corrplot(cor.dendlist(dend1to8), method = "ellipse", type = "lower")
corrplot(cor.dendlist(dend1to8), method = "number", type = "lower")
corrplot(cor.dendlist(dend1to8), method = "color", type = "lower")
dev.off()

################ UMAP of of vst-transformed counts

# run UMAP on vst-transformed counts
out <- umap(t(foo))

# add UMAP results to condition data frame
table(rownames(out$layout)==cond$sample)
table(cond$sample==rownames(out$layout))
cond$umap_x <- out$layout[,1]
cond$umap_y <- out$layout[,2]

# plot UMAP results, colored by group
bright <- c("water" = "#66CCEE", "water-ultrasound" = "#AA3377",
            "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")

umap.group.hull <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = group)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = group, fill = group), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = bright) + 
  scale_color_manual(values = bright)

umap.group <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = group)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = bright)

umap.group.text <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = group, label = sample)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = bright) +
  geom_text_repel()
umap.group.text

pdf(file = "umap_group.pdf", width = 11, height = 5)
grid.arrange(umap.group.hull, umap.group, nrow = 1)
dev.off()

# plot UMAP results, colored by trial
light <- color("light")

umap.trial.hull <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = trial)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = trial, fill = trial), geom = "polygon", alpha = 0.1) + 
  theme_bw() + 
  scale_fill_manual(values = light(8)) + 
  scale_color_manual(values = light(8))

umap.trial <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = trial)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_fill_manual(values = light(8)) + 
  scale_color_manual(values = light(8))

pdf(file = "umap_trial.pdf", width = 11, height = 5)
grid.arrange(umap.trial.hull, umap.trial, nrow = 1)
dev.off()

# plot UMAP results, colored by ultrasound
umap.ultrasound.hull <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = ultrasound)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = ultrasound, fill = ultrasound), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("gray", "deeppink3")) + 
  scale_color_manual(values = c("gray", "deeppink3"))

umap.ultrasound <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = ultrasound)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_fill_manual(values = c("gray", "deeppink3")) + 
  scale_color_manual(values = c("gray", "deeppink3"))

pdf(file = "umap_ultrasound.pdf", width = 11, height = 5)
grid.arrange(umap.ultrasound.hull, umap.ultrasound, nrow = 1)
dev.off()

# plot UMAP results, colored by chlorine
umap.chlorine.hull <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = chlorine)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = chlorine, fill = chlorine), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("turquoise3", "gray")) + 
  scale_color_manual(values = c("turquoise3", "gray"))

umap.chlorine <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = chlorine)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_fill_manual(values = c("turquoise3", "gray")) + 
  scale_color_manual(values = c("turquoise3", "gray"))

pdf(file = "umap_chlorine.pdf", width = 11, height = 5)
grid.arrange(umap.chlorine.hull, umap.chlorine, nrow = 1)
dev.off()

# plot UMAP results, colored by water
umap.water.hull <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = water)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = water, fill = water), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("gray", "dodgerblue")) + 
  scale_color_manual(values = c("gray", "dodgerblue"))

umap.water <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = water)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_fill_manual(values = c("gray", "dodgerblue")) + 
  scale_color_manual(values = c("gray", "dodgerblue"))

pdf(file = "umap_water.pdf", width = 11, height = 5)
grid.arrange(umap.water.hull, umap.water, nrow = 1)
dev.off()

# sequencing depth
pdf(file = "umap_seqdepth.pdf", width = 8.5, height = 5)
ggplot(data = cond, aes(x = umap_x, y = umap_y, color = log10(seqdepth), label = sample)) + 
  geom_point(size = 3) +
  theme_bw() + 
  scale_color_viridis() + 
  geom_text_repel()
dev.off()

################ PCA of of vst-transformed counts
pca <- prcomp(t(foo), scale = F, center = T)
pca3col <- as.data.frame(cbind(pca$x[,1],pca$x[,2],pca$x[,3]))

table(rownames(pca3col)==cond$sample)
table(cond$sample==rownames(pca3col))
cond$PC1 <- pca3col$V1
cond$PC2 <- pca3col$V2
cond$PC3 <- pca3col$V3

# plot PCA, color by group
pca.group.hull <- ggplot(data = cond, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = group, fill = group), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = bright) + 
  scale_color_manual(values = bright)

pca.group <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = group)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = bright)

pca.group.text <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = group, label = sample)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = bright) +
  geom_text_repel()
pca.group.text

pdf(file = "pca_group.pdf", width = 11, height = 5)
grid.arrange(pca.group.hull, pca.group, nrow = 1)
dev.off()

# plot PCA, color by trial
pca.trial.hull <- ggplot(data = cond, aes(x = PC1, y = PC2, color = trial)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = trial, fill = trial), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = light(8)) + 
  scale_color_manual(values = light(8))

pca.trial <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = trial)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = light(8))

pdf(file = "pca_trial.pdf", width = 11, height = 5)
grid.arrange(pca.trial.hull, pca.trial, nrow = 1)
dev.off()

# plot PCA, color by ultrasound
pca.ultrasound.hull <- ggplot(data = cond, aes(x = PC1, y = PC2, color = ultrasound)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = ultrasound, fill = ultrasound), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("gray", "deeppink3")) + 
  scale_color_manual(values = c("gray", "deeppink3"))

pca.ultrasound <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = ultrasound)) + 
  geom_point() +
  theme_bw() +
  scale_fill_manual(values = c("gray", "deeppink3")) + 
  scale_color_manual(values = c("gray", "deeppink3"))

pdf(file = "pca_ultrasound.pdf", width = 11, height = 5)
grid.arrange(pca.ultrasound.hull, pca.ultrasound, nrow = 1)
dev.off()

# plot PCA, color by chlorine
pca.chlorine.hull <- ggplot(data = cond, aes(x = PC1, y = PC2, color = chlorine)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = chlorine, fill = chlorine), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("turquoise3", "gray")) + 
  scale_color_manual(values = c("turquoise3", "gray"))

pca.chlorine <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = chlorine)) + 
  geom_point() +
  theme_bw() +
  scale_fill_manual(values = c("turquoise3", "gray")) + 
  scale_color_manual(values = c("turquoise3", "gray"))

pdf(file = "pca_chlorine.pdf", width = 11, height = 5)
grid.arrange(pca.chlorine.hull, pca.chlorine, nrow = 1)
dev.off()

# plot PCA, color by water
pca.water.hull <- ggplot(data = cond, aes(x = PC1, y = PC2, color = water)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = water, fill = water), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("gray", "dodgerblue")) + 
  scale_color_manual(values = c("gray", "dodgerblue"))

pca.water <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = water)) + 
  geom_point() +
  theme_bw() +
  scale_fill_manual(values = c("gray", "dodgerblue")) + 
  scale_color_manual(values = c("gray", "dodgerblue"))

pdf(file = "pca_water.pdf", width = 11, height = 5)
grid.arrange(pca.water.hull, pca.water, nrow = 1)
dev.off()

# sequencing depth
pdf(file = "pca_seqdepth.pdf", width = 11, height = 5)
ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = log10(seqdepth), label = sample)) + 
  geom_point() +
  theme_bw() + 
  scale_color_viridis() + 
  geom_text_repel()
dev.off()
