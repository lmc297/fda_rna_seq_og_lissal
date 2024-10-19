library(ggplot2)
library(stringr)
library(reshape2)
library(ape)
library(khroma)
library(Rsubread)

setwd(dir = "~/Google Drive/joelle_fda/rnaseq_101724")

############# Read in reference genome annotation GFF

ref_annot <- ape::read.gff(file = "reference_genome/listeria_GCA_016306305.1.gff",
                           na.strings = c(".", "?"), GFF3 = T)
head(ref_annot)
table(ref_annot$type)
ref_annot <- subset(ref_annot, type == "gene")

gene_attr <- stringr::str_split(ref_annot$attributes, ";")

locus_tags <- unlist(lapply(gene_attr, function(x) {
  x[grepl("locus_tag", x)]
}))

gene_biotypes <- unlist(lapply(gene_attr, function(x) {
  x[grepl("gene_biotype", x)]
}))

common_gene_names <- unlist(lapply(gene_attr, function(x) {
  x <- x[grepl("gene=", x)]
  x[identical(x, character(0))] <- ""
  x
}))

gene_lengths <- (ref_annot$end - ref_annot$start) + 1

ref_gene_df <- data.frame(
  locus_tag = locus_tags,
  biotype = gene_biotypes,
  gene_name = common_gene_names,
  gene_length = gene_lengths
)
ref_gene_df$locus_tag <- gsub("locus_tag=", "", ref_gene_df$locus_tag)
ref_gene_df$biotype <- gsub("gene_biotype=", "", ref_gene_df$biotype)
ref_gene_df$gene_name <- gsub("gene=", "", ref_gene_df$gene_name)

head(ref_gene_df)
table(ref_gene_df$biotype)

write.table(x = ref_gene_df, file = "reference_genome/ref_gene_df.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F)

############# Read in BAM files

bamfiles <- Sys.glob(paths = "trimmed_reads_stranded/*.bam")

gene_counts <- Rsubread::featureCounts(
  files = bamfiles,
  annot.ext = "reference_genome/listeria_GCA_016306305.1.gff",
  isGTFAnnotationFile = T,
  GTF.featureType = "gene",
  GTF.attrType = "locus_tag",
  nthreads = 1,
  countMultiMappingReads = T,
  fraction = T, # assign fractional counts to multimappers
  isPairedEnd = T,
  strandSpecific = 2
)

stats_mat <- as.data.frame(gene_counts$stat)
colnames(stats_mat) <- gsub(pattern = ".bam", replacement = "", x = colnames(stats_mat))
head(stats_mat)

counts_mat <- as.data.frame(gene_counts$counts)
head(counts_mat)
colnames(counts_mat) <- gsub(pattern = ".bam", replacement = "", x = colnames(counts_mat))
head(counts_mat)
table(rownames(counts_mat)%in%ref_gene_df$locus_tag)
table(ref_gene_df$locus_tag%in%rownames(counts_mat))
table(rownames(counts_mat)==ref_gene_df$locus_tag)
table(ref_gene_df$locus_tag==rownames(counts_mat))
counts_mat <- cbind(ref_gene_df, counts_mat)
head(counts_mat)

write.table(x = counts_mat, file = "gene_counts.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F)

gene_counts_pc <- counts_mat[counts_mat$biotype == "protein_coding", ]
head(gene_counts_pc)
write.table(x = gene_counts_pc, file = "gene_counts_pc.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F)

all_biotypes <- unique(counts_mat$biotype)
all_biotypes

biotype_counts <- data.frame(do.call(
  cbind,
  lapply(all_biotypes, function(biotype) {
    colSums(gene_counts$counts[counts_mat$biotype == biotype, , drop = FALSE])
  })
))
colnames(biotype_counts) <- all_biotypes
rownames(biotype_counts) <- gsub(pattern = ".bam", replacement = "", x = rownames(biotype_counts))
head(biotype_counts)

counts_summary <- data.frame(rownames(biotype_counts),
                             unlist(lapply(strsplit(x = rownames(biotype_counts), split = "_"), "[[", 3)),
                             unlist(lapply(strsplit(x = rownames(biotype_counts), split = "_"), "[[", 2)))
colnames(counts_summary) <- c("sample", "group", "rep_no")
counts_summary
table(counts_summary$sample==rownames(biotype_counts))
table(rownames(biotype_counts)==counts_summary$sample)
counts_summary <- cbind(counts_summary, biotype_counts)
counts_summary

countfiles <- Sys.glob(paths = "trimmed_reads_stranded/*.counts")
fvec <- c()
cvec <- c()
for (i in 1:length(countfiles)){
  fname <- countfiles[i]
  fname <- gsub(pattern = "trimmed_reads_stranded/", replacement = "", x = fname)
  fname <- gsub(pattern = ".counts", replacement = "", x = fname)
  fvec <- c(fvec, fname)
  cf <- read.delim(file = countfiles[i],
                   header = F, sep = "\t",
                   stringsAsFactors = F,
                   check.names = F)
  cvec <- c(cvec, cf$V1)
}
merged_total_counts <- data.frame(fvec, cvec)
colnames(merged_total_counts) <- c("sample", "mapped")
merged_total_counts

table(counts_summary$sample%in%merged_total_counts$sample)
table(merged_total_counts$sample%in%counts_summary$sample)
table(counts_summary$sample==merged_total_counts$sample)
table(merged_total_counts$sample==counts_summary$sample)
#counts_summary <- merge(counts_summary, merged_total_counts, by = "sample")
merged_total_counts <- merged_total_counts[match(counts_summary$sample, merged_total_counts$sample),]
table(counts_summary$sample==merged_total_counts$sample)
table(merged_total_counts$sample==counts_summary$sample)
counts_summary$mapped <- merged_total_counts$mapped
counts_summary
counts_summary$missing <- counts_summary$mapped - rowSums(counts_summary[,4:(ncol(counts_summary)-1)])
counts_summary
counts_summary$other <- counts_summary$mapped - counts_summary$protein_coding
counts_summary

write.table(x = counts_summary, file = "counts_summary.tsv",
  col.names = T, row.names = F,
  sep = "\t", quote = F)


counts_melt <- reshape2::melt(
  counts_summary,
  id.vars = c("sample"),
  measure.vars = c(
    "protein_coding",
    "other"
    #all_biotypes
  )
)

counts_melt

p1 <- ggplot(
  counts_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("# read pairs") + 
  theme_bw() + 
  scale_fill_manual(values = c("turquoise3", "deeppink3")) +
  scale_color_manual(values = c("turquoise3", "deeppink3"))
p1

propCols <- counts_summary[c(
  "protein_coding",
  "other"
  # all_biotypes
)] / counts_summary$mapped

propCols$sample <- counts_summary$sample

# rowSums(propCols) ## each row should sum to 1
prop_melt <- reshape2::melt(
  propCols,
  id.vars = c("sample"),
  measure.vars = c(
    "protein_coding",
    "other"
    # all_biotypes
  )
)

p2 <- ggplot(
  prop_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("Proportion of reads") + 
  theme_bw() + 
  scale_fill_manual(values = c("turquoise3", "deeppink3")) +
  scale_color_manual(values = c("turquoise3", "deeppink3"))
p2

pdf(file = "bar_counts_total.pdf", width = 8.5, height = 6)
p1
dev.off()

pdf(file = "bar_counts_proportion.pdf", width = 8.5, height = 6)
p2
dev.off()


counts_melt2 <- reshape2::melt(
  counts_summary,
  id.vars = c("sample"),
  measure.vars = c(
    all_biotypes,
    "missing"
    #all_biotypes[all_biotypes%in%colnames(counts_summary.expanded)]
  )
)

counts_melt2
table(counts_melt2$variable)
vibrant <- color("muted")

p3 <- ggplot(
  counts_melt2,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("# read pairs") + 
  theme_bw() + 
  scale_fill_manual(values = vibrant(n = 8)) +
  scale_color_manual(values = vibrant(n = 8))
p3

propCols.expanded <- counts_summary[c(
  all_biotypes, 
  "missing"
)] / counts_summary$mapped

propCols.expanded$sample <- counts_summary$sample

# rowSums(propCols) ## each row should sum to 1
prop_melt.expanded <- reshape2::melt(
  propCols.expanded,
  id.vars = c("sample"),
  measure.vars = c(
    all_biotypes,
    "missing"
  )
)

p4 <- ggplot(
  prop_melt.expanded,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("Proportion of reads") + 
  theme_bw() + 
  scale_fill_manual(values = vibrant(n = 8)) +
  scale_color_manual(values = vibrant(n = 8))
p4

pdf(file = "bar_counts_total_expanded.pdf", width = 8.5, height = 6)
p3
dev.off()

pdf(file = "bar_counts_proportion_expanded.pdf", width = 8.5, height = 6)
p4
dev.off()

############# stats

stats_mat = setNames(data.frame(t(stats_mat[,-1])), stats_mat[,1])
stats_mat
stats_mat <- stats_mat[,which(colSums(stats_mat)>0)]
stats_mat
stats_mat$sample <- rownames(stats_mat)
stats_mat

stats_melt <- reshape2::melt(
  stats_mat,
  id.vars = c("sample"),
  measure.vars = c(
    colnames(stats_mat[1:(ncol(stats_mat)-1)])
    #all_biotypes
  )
)

stats_melt
table(stats_melt$variable)
vibrant <- color("vibrant")

p5 <- ggplot(
  stats_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("# read pairs") + 
  theme_bw() + 
  scale_fill_manual(values = vibrant(n = 4)) +
  scale_color_manual(values = vibrant(n = 4))
p5

statCols <- stats_mat[c(
  colnames(stats_mat)[1:(ncol(stats_mat)-1)]
)] / rowSums(stats_mat[,1:(ncol(stats_mat)-1)])

statCols$sample <- stats_mat$sample
statCols

# rowSums(propCols) ## each row should sum to 1
propStat_melt <- reshape2::melt(
  statCols,
  id.vars = c("sample"),
  measure.vars = c(
    colnames(statCols[1:(ncol(statCols)-1)])
  )
)
propStat_melt

p6 <- ggplot(
  propStat_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("Proportion of reads") + 
  theme_bw() + 
  scale_fill_manual(values = vibrant(n = 4)) +
  scale_color_manual(values = vibrant(n = 4))
p6

pdf(file = "bar_stats_total.pdf", width = 8.5, height = 6)
p5
dev.off()

pdf(file = "bar_stats_proportion.pdf", width = 8.5, height = 6)
p6
dev.off()
