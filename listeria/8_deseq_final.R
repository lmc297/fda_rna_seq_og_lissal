library(ggplot2)
library(DESeq2)
library(umap)
library(ggpubr)
library(viridis)
library(ggrepel)
library(ggVennDiagram)
library(pheatmap)

################ prepare input data
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

################ differential expression testing

# create DESeq object from count matrix, this time with group in our design
dds <- DESeqDataSetFromMatrix(countData = counts_tab,
                              colData = cond,
                              design= ~ group) 

# run DESeq
dds <- DESeq(dds)

################ water-ultrasound vs water

# test water-ultrasound vs water
res <- results(dds, c("group", "water-ultrasound", "water"))

# reorder based on p-value
res <- res[order(res$pvalue),] 

# add name column with gene IDs
res$name <- rownames(res)

# save results as res_wu (water-ultrasound)
res_wu <- res

# add siggene column
# adjusted p-value < 0.05 & log2fc >= 1 == Upregulated
# adjusted p-value < 0.05 & log2fc <= -1 == Downregulated
# Everything else == Not significant
res$siggene <- ifelse(test = (res$padj<0.05 & res$log2FoldChange>=1.0), yes = "Upregulated",
                    no = ifelse(test = (res$padj<0.05 & res$log2FoldChange<=-1.0), yes = "Downregulated",
                                no = "Not significant"))
table(res$siggene)

# get quantiles of -log10(p-values)
quantile(x = -log10(res$pvalue),
         probs = c(0.75, 0.90, 0.95, 0.99), na.rm = T)

# plot histogram of -log10(p-values)
pdf(file = "hist_ultrasound_vs_water.pdf", width = 8, height = 8)
ggplot(data = as.data.frame(res), mapping = aes(x = -log10(pvalue))) + 
  geom_histogram() +
  theme_bw()
dev.off()

# get most significant p-values
supersig <- rownames(res)[which(-log10(res$pvalue)>=72.13)]
supersig

# plot volcano plot
pdf(file = "volcano_ultrasound_vs_water.pdf", width = 8, height = 8)
ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(pvalue), color = siggene)) + 
  geom_point() + 
  theme_bw() +
  geom_vline(xintercept = 1, color = "red3", linetype = "dashed") + 
  geom_vline(xintercept = -1, color = "blue3", linetype = "dashed") +
  scale_color_manual(values = c("blue3", "gray54", "red3")) +
  geom_hline(yintercept = 72.13, linetype = "dashed") +
  geom_text_repel(aes(label = ifelse(name%in%supersig, name, "")))
dev.off()

# save results to TSV
write.table(x = res, file = "deseq_ultrasound_vs_water.tsv",
            append = F, quote = F,
            sep = "\t", row.names = T,
            col.names = T)

################ chlorine vs water
res <- results(dds, c("group","chlorine", "water")) 
res <- res[order(res$pvalue),] 
res$name <- rownames(res)
res_wc <- res
res
res$siggene <- ifelse(test = (res$padj<0.05 & res$log2FoldChange>=1.0), yes = "Upregulated",
                      no = ifelse(test = (res$padj<0.05 & res$log2FoldChange<=-1.0), yes = "Downregulated",
                                  no = "Not significant"))
table(res$siggene)
quantile(x = -log10(res$pvalue),
         probs = c(0.75, 0.90, 0.95, 0.99), na.rm = T)

pdf(file = "hist_chlorine_vs_water.pdf", width = 8, height = 8)
ggplot(data = as.data.frame(res), mapping = aes(x = -log10(pvalue))) + 
  geom_histogram() +
  theme_bw()
dev.off()

supersig <- rownames(res)[which(-log10(res$pvalue)>=8.96)]
supersig

pdf(file = "volcano_chlorine_vs_water.pdf", width = 8, height = 8)
ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(pvalue), color = siggene)) + 
  geom_point() + 
  theme_bw() +
  geom_vline(xintercept = 1, color = "red3", linetype = "dashed") + 
  geom_vline(xintercept = -1, color = "blue3", linetype = "dashed") +
  scale_color_manual(values = c("blue3", "gray54", "red3")) +
  geom_hline(yintercept = 8.96, linetype = "dashed") +
  geom_text_repel(aes(label = ifelse(name%in%supersig, name, "")))
dev.off()

write.table(x = res, file = "deseq_chlorine_vs_water.tsv",
            append = F, quote = F,
            sep = "\t", row.names = T,
            col.names = T)

################ chlorine-ultrasound vs water

res <- results(dds, c("group", "chlorine-ultrasound", "water")) 
res <- res[order(res$pvalue),] 
res$name <- rownames(res)
res_wcu <- res
res

res$siggene <- ifelse(test = (res$padj<0.05 & res$log2FoldChange>=1.0), yes = "Upregulated",
                      no = ifelse(test = (res$padj<0.05 & res$log2FoldChange<=-1.0), yes = "Downregulated",
                                  no = "Not significant"))
table(res$siggene)
quantile(x = -log10(res$pvalue),
         probs = c(0.75, 0.90, 0.95, 0.99), na.rm = T)

pdf(file = "hist_chlorine-ultrasound_vs_water.pdf", width = 8, height = 8)
ggplot(data = as.data.frame(res), mapping = aes(x = -log10(pvalue))) + 
  geom_histogram() + 
  theme_bw()
dev.off()


supersig <- rownames(res)[which(-log10(res$pvalue)>=68.09)]
supersig

pdf(file = "volcano_chlorine-ultrasound_vs_water.pdf", width = 8, height = 8)
ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(pvalue), color = siggene)) + 
  geom_point() + 
  theme_bw() +
  geom_vline(xintercept = 1, color = "red3", linetype = "dashed") + 
  geom_vline(xintercept = -1, color = "blue3", linetype = "dashed") +
  scale_color_manual(values = c("blue3", "gray54", "red3")) +
  geom_hline(yintercept = 68.09, linetype = "dashed") +
  geom_text_repel(aes(label = ifelse(name%in%supersig, name, "")))
dev.off()

write.table(x = res, file = "deseq_chlorine-ultrasound_vs_water.tsv",
            append = F, quote = F,
            sep = "\t", row.names = T,
            col.names = T)

################ venn diagram of significant genes

# create a function to get "top genes" for each condition
# top genes have adjusted p-value < 0.05 AND a log2 fold change >= 1 or <= -1
topgenes <- function(res.df){
  yeet <- rownames(res.df)[which(res.df$padj<0.05 & abs(res.df$log2FoldChange)>=1.0)]
  return(yeet)
}

# make a list of top genes for each condition
list4venn <- list(
  Chlorine = topgenes(res_wc),
  Ultrasound = topgenes(res_wu),
  ChlorineUltrasound = topgenes(res_wcu)
)

# plot venn diagram of top genes
pdf(file = "topgenes_venn.pdf", width = 8, height = 8)
ggVennDiagram(list4venn, label_alpha = 0, color = rep("black", 3)) + 
  scale_fill_gradient(low="grey90",high = "deeppink3") + 
  scale_color_manual(values = rep("black", 3))
dev.off()

# sanity checks (see if we get the same results as the venn)
table(topgenes(res_wc)%in%topgenes(res_wu))
table(topgenes(res_wc)%in%topgenes(res_wcu))
table(topgenes(res_wcu)%in%topgenes(res_wu))


################ heatmap of most significant genes

# get the 50 most significant genes
supertop <- rbind(res_wc, res_wu, res_wcu)
supertop <- supertop[order(supertop$pvalue, decreasing = F),]
head(supertop)

all_topgenes <- supertop[1:50,]
dim(all_topgenes)
head(all_topgenes)

# scale and center raw counts
counts.scaled <- scale(x = t(counts_tab), center = T, scale = T)
head(counts.scaled[,1:5])

# plot heatmap of most significant genes
# note that batch effects not removed; interpret with caution
row.metadata <- unlist(lapply(strsplit(x = rownames(counts.scaled), split = "_"), "[[", 3))
row.metadata <- as.data.frame(gsub(pattern = "-", replacement = "_", x = row.metadata))
rownames(row.metadata) <- rownames(counts.scaled)
colnames(row.metadata) <- "Group"
row.metadata

ann_colors <- list(Group = c(water="#66CCEE", water_ultrasound="#AA3377",
             chlorine="#228833", chlorine_ultrasound="#CCBB44"))


# save heatmap as PDF
pdf(file = "heatmap_top50genes.pdf", width = 11, height = 8)
pheatmap(mat = counts.scaled[,which(colnames(counts.scaled)%in%all_topgenes$name)],
         fontsize_col = 5, fontsize_row = 8,
         annotation_row = row.metadata,
         annotation_colors = ann_colors)
dev.off()
