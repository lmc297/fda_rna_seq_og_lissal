library(ggplot2)
library(clusterProfiler)

# uses clusterProfiler: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html

###################### set seed

set.seed(42)

###################### define functions: kegg pathway enrichment

# define run.kegg.pathway function, to perform kegg pathway enrichment
run.kegg.pathway <- function(genes.og, updown, fc, p.cutoff, p.correction, treatment, term2gene){
  
    # get significant genes, which are < user-supplied p-value cutoff
    genes <- genes.og[which(genes.og$padj < p.cutoff),]
    
    # if users want upregulated genes:
    # take significant genes with fold change >= user-supplied fold-change cutoff
    if (updown == "up"){
        sig.genes <- rownames(genes[which(genes$log2FoldChange >= fc),])
    }
    
    # if users want downregulated genes:
    # take significant genes with fold change <= user-supplied fold-change cutoff
    else if (updown == "down"){
        sig.genes <- rownames(genes[which(genes$log2FoldChange <= fc),])
    }
    
    # throw an error if user doesn't specify up or down
    else{
        stop("Please select 'up' or 'down' to query over-represented pathways among
              up- andn down-regulated genes, respectively.")
    }
    
    # KEGG pathway over-representation analysis
    # use enricher
    kk <- enricher(gene = sig.genes,
                     TERM2GENE = term2gene,
                     pvalueCutoff = p.cutoff,
                     pAdjustMethod = p.correction,
                     minGSSize = 1, maxGSSize = 1e9)

    # KEGG pathway gene set enrichment analysis
    # rank genes by log2 fold change from highest to lowest
    genes.og <- genes.og[order(genes.og$log2FoldChange, decreasing = T),]
    allgenes.ranked <- genes.og$log2FoldChange
    names(allgenes.ranked) <- rownames(genes.og)
    print(allgenes.ranked)[1:5]
    
    # use GSEA
    kk2 <- GSEA(geneList = allgenes.ranked,
                   TERM2GENE = term2gene,
                   minGSSize = 1,
                   maxGSSize = 1e9,
                   eps = 0,
                   pvalueCutoff = p.cutoff,
                   pAdjustMethod = p.correction,
                   verbose = F)

    # get results with adjusted p<0.05
    kk.final <- kk@result[which(kk@result$p.adjust < 0.05),]
    kk2.final <- kk2@result[which(kk2@result$p.adjust < 0.05),]

    # add treatment columns to final data frames
    kk.final$treatment <- rep(treatment, nrow(kk.final))
    kk2.final$treatment <- rep(treatment, nrow(kk2.final))

    # return results for both analyses
    return(list("pathway.overrep" = kk.final,
                "pathway.gsea" = kk2.final))
}

###################### define functions: kegg module enrichment

# define run.kegg.module function to perform kegg module enrichment
run.kegg.module <- function(genes.og, updown, fc, p.cutoff, p.correction, treatment, term2gene){
  
    # get genes that have p-value < user-specified cutoff
    genes <- genes.og[which(genes.og$padj < p.cutoff),]
    
    # get significant up or down-regulated genes (or exit if up or down not specified)
    if (updown == "up"){
        sig.genes <- rownames(genes[which(genes$log2FoldChange >= fc),])
    }
    else if (updown == "down"){
        sig.genes <- rownames(genes[which(genes$log2FoldChange <= fc),])
    }
    else{
        stop("Please select 'up' or 'down' to query over-represented pathways among
              up- andn down-regulated genes, respectively.")
    }

    # KEGG module over-representation analysis
    # uses enricher
    mkk <- enricher(gene = sig.genes,
                    TERM2GENE = term2gene,
                    pvalueCutoff = p.cutoff,
                    pAdjustMethod = p.correction,
                    minGSSize = 1,
                    maxGSSize = 1e9)

    # KEGG module gene set enrichment analysis
    # rank genes by log2 fold change from highest to lowest
    genes.og <- genes.og[order(genes.og$log2FoldChange, decreasing = T),]
    allgenes.ranked <- genes.og$log2FoldChange
    names(allgenes.ranked) <- rownames(genes.og)
    print(allgenes.ranked)[1:5]
    
    # use GSEA function
    mkk2 <- GSEA(geneList = allgenes.ranked,
                 TERM2GENE = term2gene,
                 minGSSize = 1,
                 maxGSSize = 1,
                 eps = 0,
                 pvalueCutoff = p.cutoff,
                 pAdjustMethod = p.correction,
                 verbose = F)

    # get significant results with adjusted p<0.05
    mkk.final <- mkk@result[which(mkk@result$p.adjust < 0.05),]
    mkk2.final <- mkk2@result[which(mkk2@result$p.adjust < 0.05),]

    # add treatment column to final data frame
    mkk.final$treatment <- rep(treatment, nrow(mkk.final))
    mkk2.final$treatment <- rep(treatment, nrow(mkk2.final))

    # return final data frames
    return(list("module.overrep" = mkk.final,
                "module.gsea" = mkk2.final))
}

###################### load/prepare kegg annotations

# load kegg pathways TSV
# this file contains two columns (no header)
# column 1 = kegg pathway term
# column 2 = gene name
T2G.pathway <- read.delim(file = "../eggnog/kegg_pathways.tsv",
                          header = F, sep = "\t",
                          stringsAsFactors = F,
                          check.names = F)
head(T2G.pathway)

# add column names to kegg pathways
colnames(T2G.pathway) <- c("TermID", "geneID")
head(T2G.pathway)

# load kegg modules TSV
# this file contains two columns (no header)
# column 1 = kegg module term
# column 2 = gene name
T2G.module <- read.delim(file = "../eggnog/kegg_modules.tsv",
                          header = F, sep = "\t",
                          stringsAsFactors = F,
                          check.names = F)
head(T2G.module)

# add column names to kegg modules
colnames(T2G.module) <- c("TermID", "geneID")
head(T2G.module)

###################### kegg over-representation analysis: ultrasound vs water

# load annotated DESeq data
deseq.ultra <- read.delim(file = "../10_topgo_cds/annot_deseq_ultrasound_vs_water.tsv",
                          header = T, sep = "\t",
                          stringsAsFactors = F,
                          check.names = F)
head(deseq.ultra)

# remove rows (genes) where log2 fold change is NA
deseq.ultra <- deseq.ultra[!(is.na(deseq.ultra$log2FoldChange)),]

# perform kegg pathway enrichment analysis (upregulated genes)
deseq.ultra.up.pathway <- run.kegg.pathway(genes.og = deseq.ultra,
                                       term2gene = T2G.pathway,
                                       updown = "up", fc = 1.0,
                                       p.cutoff = 0.05,
                                       p.correction = "fdr",
                                       treatment = "ultrasound")

dim(deseq.ultra.up.pathway$pathway.overrep)
dim(deseq.ultra.up.pathway$pathway.gsea)

# perform kegg pathway enrichment analysis (downregulated genes)
deseq.ultra.down.pathway <- run.kegg.pathway(genes.og = deseq.ultra,
                                           term2gene = T2G.pathway,
                                           updown = "down", fc = -1.0,
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "ultrasound")

dim(deseq.ultra.down.pathway$pathway.overrep)
dim(deseq.ultra.down.pathway$pathway.gsea)

# perform kegg module enrichment analysis (upregulated genes)
deseq.ultra.up.module <- run.kegg.module(genes.og = deseq.ultra,
                                           term2gene = T2G.module,
                                           updown = "up", fc = 1.0,
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "ultrasound")

dim(deseq.ultra.up.module$module.overrep)
dim(deseq.ultra.up.module$module.gsea)

# perform kegg module enrichment analysis (downregulated genes)
deseq.ultra.down.module <- run.kegg.module(genes.og = deseq.ultra,
                                             term2gene = T2G.module,
                                             updown = "down", fc = -1.0,
                                             p.cutoff = 0.05,
                                             p.correction = "fdr",
                                             treatment = "ultrasound")

dim(deseq.ultra.down.module$module.overrep)
dim(deseq.ultra.down.module$module.gsea)

###################### kegg over-representation analysis: chlorine vs water

deseq.chlor <- read.delim(file = "../10_topgo_cds/annot_deseq_chlorine_vs_water.tsv",
                          header = T, sep = "\t",
                          stringsAsFactors = F,
                          check.names = F)
head(deseq.chlor)

deseq.chlor <- deseq.chlor[!(is.na(deseq.chlor$log2FoldChange)),]
head(deseq.chlor)

deseq.chlor.up.pathway <- run.kegg.pathway(genes.og = deseq.chlor,
                                           term2gene = T2G.pathway,
                                           updown = "up", fc = 1.0,
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "chlorine")

dim(deseq.chlor.up.pathway$pathway.overrep)
dim(deseq.chlor.up.pathway$pathway.gsea)

deseq.chlor.down.pathway <- run.kegg.pathway(genes.og = deseq.chlor,
                                             term2gene = T2G.pathway,
                                             updown = "down", fc = -1.0,
                                             p.cutoff = 0.05,
                                             p.correction = "fdr",
                                             treatment = "chlorine")

dim(deseq.chlor.down.pathway$pathway.overrep)
dim(deseq.chlor.down.pathway$pathway.gsea)

deseq.chlor.up.module <- run.kegg.module(genes.og = deseq.chlor,
                                         term2gene = T2G.module,
                                         updown = "up", fc = 1.0,
                                         p.cutoff = 0.05,
                                         p.correction = "fdr",
                                         treatment = "chlorine")

dim(deseq.chlor.up.module$module.overrep)
dim(deseq.chlor.up.module$module.gsea)

deseq.chlor.down.module <- run.kegg.module(genes.og = deseq.chlor,
                                           term2gene = T2G.module,
                                           updown = "down", fc = -1.0,
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "chlorine")

dim(deseq.chlor.down.module$module.overrep)
dim(deseq.chlor.down.module$module.gsea)

###################### kegg over-representation analysis: ultrasound-chlorine vs water

deseq.chlorultra <- read.delim(file = "../10_topgo_cds/annot_deseq_chlorine-ultrasound_vs_water.tsv",
                          header = T, sep = "\t",
                          stringsAsFactors = F,
                          check.names = F)
head(deseq.chlorultra)

deseq.chlorultra <- deseq.chlorultra[!(is.na(deseq.chlorultra$log2FoldChange)),]
head(deseq.chlorultra)

deseq.chlorultra.up.pathway <- run.kegg.pathway(genes.og = deseq.chlorultra,
                                           term2gene = T2G.pathway,
                                           updown = "up", fc = 1.0,
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "chlorine-ultrasound")

dim(deseq.chlorultra.up.pathway$pathway.overrep)
dim(deseq.chlorultra.up.pathway$pathway.gsea)

deseq.chlorultra.down.pathway <- run.kegg.pathway(genes.og = deseq.chlorultra,
                                             term2gene = T2G.pathway,
                                             updown = "down", fc = -1.0,
                                             p.cutoff = 0.05,
                                             p.correction = "fdr",
                                             treatment = "chlorine-ultrasound")

dim(deseq.chlorultra.down.pathway$pathway.overrep)
dim(deseq.chlorultra.down.pathway$pathway.gsea)

deseq.chlorultra.up.module <- run.kegg.module(genes.og = deseq.chlorultra,
                                         term2gene = T2G.module,
                                         updown = "up", fc = 1.0,
                                         p.cutoff = 0.05,
                                         p.correction = "fdr",
                                         treatment = "chlorine-ultrasound")

dim(deseq.chlorultra.up.module$module.overrep)
dim(deseq.chlorultra.up.module$module.gsea)

deseq.chlorultra.down.module <- run.kegg.module(genes.og = deseq.chlorultra,
                                           term2gene = T2G.module,
                                           updown = "down", fc = -1.0,
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "chlorine-ultrasound")

dim(deseq.chlorultra.down.module$module.overrep)
dim(deseq.chlorultra.down.module$module.gsea)

###################### create and save final TSV files

# make data frame of enriched kegg pathways for upregulated genes
final.pathway.up <- rbind(deseq.ultra.up.pathway$pathway.overrep,
                          deseq.chlor.up.pathway$pathway.overrep,
                          deseq.chlorultra.up.pathway$pathway.overrep)

# save TSV file
write.table(x = final.pathway.up, file = "listeria_kegg_pathway_up.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# make data frame of enriched kegg pathways for downregulated genes
final.pathway.down <- rbind(deseq.ultra.down.pathway$pathway.overrep,
                          deseq.chlor.down.pathway$pathway.overrep,
                          deseq.chlorultra.down.pathway$pathway.overrep)

# save TSV file
write.table(x = final.pathway.down, file = "listeria_kegg_pathway_down.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final GSEA data frame for kegg pathways
final.pathway.gsea <- rbind(deseq.ultra.up.pathway$pathway.gsea,
                            deseq.chlor.up.pathway$pathway.gsea,
                            deseq.chlorultra.up.pathway$pathway.gsea)

# save TSV file
write.table(x = final.pathway.gsea, file = "listeria_kegg_pathway_gsea.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final data frame with kegg module enrichment results for upregulated genes
final.module.up <- rbind(deseq.ultra.up.module$module.overrep,
                         deseq.chlor.up.module$module.overrep,
                         deseq.chlorultra.up.module$module.overrep)

# save TSV file
write.table(x = final.module.up, file = "listeria_kegg_module_up.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final data frame with kegg module enrichment results for downregulated genes
final.module.down <- rbind(deseq.ultra.down.module$module.overrep,
                           deseq.chlor.down.module$module.overrep,
                           deseq.chlorultra.down.module$module.overrep)

# save TSV file
write.table(x = final.module.down, file = "listeria_kegg_module_down.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final GSEA data frame for kegg modules
final.module.gsea <- rbind(deseq.ultra.up.module$module.gsea,
                           deseq.chlor.up.module$module.gsea,
                           deseq.chlorultra.up.module$module.gsea)

# save TSV file
write.table(x = final.module.gsea, file = "listeria_kegg_module_gsea.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

###################### create and save plots

# for kegg pathways, take -log10 of the raw p-values
final.pathway.up$log10p <- -log10(final.pathway.up$pvalue)
final.pathway.down$log10p <- log10(final.pathway.down$pvalue)

# merge overrepresented pathways
final.pathway.overrep <- rbind(final.pathway.up, final.pathway.down)
final.pathway.overrep$direction <- c(rep("Up-regulated", nrow(final.pathway.up)),
                                         rep("Down-regulated", nrow(final.pathway.down)))
head(final.pathway.overrep)

# create KEGG_term column
final.pathway.overrep$KEGG_term <- final.pathway.overrep$Description

# overview of results
table(final.pathway.overrep$direction)
summary(final.pathway.overrep$log10p)

# make bar chart of enriched kegg pathways
pdf(file = "plot_pathway_overrep.pdf",width = 8.5, height = 11)
plot.pathway.overrep <- ggplot(final.pathway.overrep, aes(x=reorder(KEGG_term, log10p), y=log10p, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Pathway") +
    ylab("-log10 P-Value") +
    scale_fill_manual(values = c("#33BBEE", "#EE3377", "yellow")) +
    theme_bw(base_size=12) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
plot.pathway.overrep
dev.off()

# for kegg modules, take -log10 of the raw p-values
final.module.up$log10p <- -log10(final.module.up$pvalue)
final.module.down$log10p <- log10(final.module.down$pvalue)

# merge overrepresented modules
final.module.overrep <- rbind(final.module.up, final.module.down)
final.module.overrep$direction <- c(rep("Up-regulated", nrow(final.module.up)),
                                     rep("Down-regulated", nrow(final.module.down)))
head(final.module.overrep)

# create KEGG_term column
final.module.overrep$KEGG_term <- final.module.overrep$Description

# overview of results
table(final.module.overrep$direction)
summary(final.module.overrep$log10p)

# make bar chart of enriched kegg modules
pdf(file = "plot_module_overrep.pdf",width = 15, height = 11)
plot.module.overrep <- ggplot(final.module.overrep, aes(x=reorder(KEGG_term, log10p), y=log10p, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG module") +
    ylab("-log10 P-Value") +
    scale_fill_manual(values = c("#33BBEE", "#EE3377", "yellow")) +
    theme_bw(base_size=12) +
guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
plot.module.overrep
dev.off()

# add KEGG_term column to final.pathway.gsea
final.pathway.gsea$KEGG_term <- final.pathway.gsea$Description

# plot GSEA results
pdf(file = "plot_pathway_gsea.pdf",width = 15, height = 11)
plot.pathway.gsea <- ggplot(final.pathway.gsea, aes(x=reorder(KEGG_term, NES), y=NES, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Pathway") +
    ylab("-log10 P-Value") +
    scale_fill_manual(values = c("#33BBEE", "#EE3377", "yellow")) +
    theme_bw(base_size=12) +
guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
plot.pathway.gsea
dev.off()
