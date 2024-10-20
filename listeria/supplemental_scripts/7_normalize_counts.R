library(DESeq2)
library(edgeR)

# inspired by BactSeq: https://github.com/adamd3/BactSeq/blob/main/bin/normalise_counts.R

# load counts for protein-coding genes
gene_counts <- read.delim(file = "gene_counts_pc.tsv",
                         header = T, sep = "\t",
                         check.names = F,
                         stringsAsFactors = F)
head(gene_counts)

# separate metadata and counts
counts_tab <- gene_counts[, grepl(pattern = "_trial", x = colnames(gene_counts))]
rownames(counts_tab) <- gene_counts$locus_tag
head(counts_tab)

meta_tab <- gene_counts[, !(grepl(pattern = "_trial", x = colnames(gene_counts)))]
head(meta_tab)

################## normalize counts with DESeq

# create DESeq data set from count matrix 
colData <- data.frame(sample_name = colnames(counts_tab))
colData

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_tab),
  colData = colData, design = ~1 # initialize with no design
)

# estimate size factors
dds <- estimateSizeFactors(dds)

# provide counts scaled by size factors
deseq_norm <- counts(dds, normalized = T)

# add 1 to scaled counts and take log2
deseq_df <- as.data.frame(log2(deseq_norm + 1))

# save deseq_df, which contains log2 scaled counts
write.table(x = deseq_df,
            file = "deseq_counts.tsv",
            col.names = T,
            row.names = T,
            sep = "\t", quote = F)

################## normalize counts with edgeR/CPM

# create DGEList object from count table
y <- DGEList(counts = counts_tab)

# normalize library sizes via edgeR TMM method
y <- calcNormFactors(y, method = "TMM")

# compute counts per million (cpm)
# cpm normalized for lib size but not gene length
cpm_df <- as.data.frame(edgeR::cpm(y, log = F))
head(cpm_df)

# save cpm normalized counts as cpm_counts.tsv
write.table(x = cpm_df,
            file = "cpm_counts.tsv",
            col.names = T,
            row.names = T,
            sep = "\t",
            quote = F)

################## normalize counts with edgeR/RPKM

# create DGEList object with gene lengths included
y <- DGEList(counts = counts_tab,
             genes = data.frame(gene.length = as.numeric(meta_tab$gene_length)))

# normalize library sizes using TMM
y <- calcNormFactors(y, method = "TMM")

# compute reads per kilobase per million (rpkm).
# rpkm normalized for lib size + gene length
rpkm_df <- as.data.frame(edgeR::rpkm(y, log = F))

# save rpkm normalized counts as rpkm_counts.tsv
write.table(x = rpkm_df, file = "rpkm_counts.tsv",
  col.names = T, row.names = T,
  sep = "\t", quote = F)
