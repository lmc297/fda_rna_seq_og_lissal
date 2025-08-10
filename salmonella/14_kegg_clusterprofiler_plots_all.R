library(ggplot2)
library(clusterProfiler)

# uses clusterProfiler: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html

###################### set seed

set.seed(42)

###################### load final TSV files

# load data frame of enriched kegg pathways for upregulated genes
final.pathway.up <- read.delim(file = "kegg_pathway_up/salmonella_kegg_pathway_up.tsv",
                               header = T, sep = "\t",
                               stringsAsFactors = F, check.names = F)
head(final.pathway.up)
tail(final.pathway.up)

# load data frame of enriched kegg pathways for downregulated genes
final.pathway.down <- read.delim(file = "kegg_pathway_down/salmonella_kegg_pathway_down.tsv",
                                 header = T, sep = "\t",
                                 stringsAsFactors = F, check.names = F)
head(final.pathway.down)
tail(final.pathway.down)

# load final GSEA data frame for kegg pathways
final.pathway.gsea <- read.delim(file = "kegg_pathway_gsea/salmonella_kegg_pathway_gsea.tsv",
                                 header = T, sep = "\t",
                                 stringsAsFactors = F, check.names = F)
head(final.pathway.gsea)
tail(final.pathway.gsea)

# load final data frame with kegg module enrichment results for upregulated genes
final.module.up <- read.delim(file = "kegg_module_up/salmonella_kegg_module_up.tsv",
                              header = T, sep = "\t",
                              stringsAsFactors = F, check.names = F)
head(final.module.up)
tail(final.module.up)

# load final data frame with kegg module enrichment results for downregulated genes
final.module.down <- read.delim(file = "kegg_module_down/salmonella_kegg_module_down.tsv",
                                header = T, sep = "\t",
                                stringsAsFactors = F, check.names = F)
head(final.module.down)
tail(final.module.down)

# load final GSEA data frame for kegg modules
final.module.gsea <- read.delim(file = "kegg_module_gsea/salmonella_kegg_module_gsea.tsv",
                                header = T, sep = "\t",
                                stringsAsFactors = F, check.names = F)
head(final.module.gsea)
tail(final.module.gsea)

###################### KEGG pathway enrichment

# for kegg pathways, take -log10 of the raw p-values
final.pathway.up$log10p <- -log10(final.pathway.up$pvalue)
final.pathway.down$log10p <- log10(final.pathway.down$pvalue)

# merge overrepresented pathways
final.pathway.overrep <- rbind(final.pathway.up, final.pathway.down)
final.pathway.overrep$direction <- c(rep("Up-regulated", nrow(final.pathway.up)),
                                         rep("Down-regulated", nrow(final.pathway.down)))
head(final.pathway.overrep)

# create KEGG_term column
final.pathway.overrep$KEGG_term <- gsub(pattern = " - Salmonella enterica subsp. enterica serovar Typhimurium LT2",
                                        replacement = "", x = final.pathway.overrep$Description)
head(final.pathway.overrep)

# overview of results
table(final.pathway.overrep$direction)
summary(final.pathway.overrep$log10p)

# make bar chart of enriched kegg pathways
pdf(file = "plot_pathway_overrep_all.pdf",width = 11, height = 8.5)
plot.pathway.overrep <- ggplot(final.pathway.overrep, aes(x=reorder(KEGG_term, log10p), y=log10p, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Pathway") +
    ylab("(-)log10 P-Value") +
    scale_y_continuous(limits = c(-5, 8), breaks = seq(-5, 8, by = 1)) +
  scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
  theme_bw(base_size=12) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
plot.pathway.overrep
dev.off()

###################### KEGG module enrichment


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
pdf(file = "plot_module_overrep_all.pdf",width = 11, height = 8.5)
plot.module.overrep <- ggplot(final.module.overrep, aes(x=reorder(KEGG_term, log10p), y=log10p, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Module") +
    ylab("(-)log10 P-Value") +
    scale_y_continuous(limits = c(-4, 5), breaks = seq(-4, 5, by = 1)) +
  scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
  theme_bw(base_size=12) +
guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
plot.module.overrep
dev.off()

###################### KEGG pathway gsea

# add KEGG_term column to final.pathway.gsea
final.pathway.gsea$KEGG_term <- gsub(pattern = " - Salmonella enterica subsp. enterica serovar Typhimurium LT2",
                                     replacement = "", x = final.pathway.gsea$Description)

# overview of results
summary(final.pathway.gsea$pvalue)
summary(final.pathway.gsea$p.adjust)
table(final.pathway.gsea$treatment)
summary(final.pathway.gsea$NES)

# plot GSEA results
pdf(file = "plot_pathway_gsea_all.pdf",width = 11, height = 8.5)
plot.pathway.gsea <- ggplot(final.pathway.gsea, aes(x=reorder(KEGG_term, NES), y=NES, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Pathway") +
    ylab("NES") +
    scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 1)) +
  scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
  theme_bw(base_size=12) +
guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
plot.pathway.gsea
dev.off()

###################### KEGG module gsea


# add KEGG_term column to final.module.gsea
final.module.gsea$KEGG_term <- final.module.gsea$Description
head(final.module.gsea)

# overview of results
summary(final.module.gsea$pvalue)
summary(final.module.gsea$p.adjust)
table(final.module.gsea$treatment)
summary(final.module.gsea$NES)

# plot GSEA results
pdf(file = "plot_module_gsea_all.pdf",width = 15, height = 6)
plot.module.gsea <- ggplot(final.module.gsea, aes(x=reorder(KEGG_term, NES), y=NES, fill = treatment)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("KEGG Module") +
  ylab("NES") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
  theme_bw(base_size=12) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
plot.module.gsea
dev.off()

# make final merged plots: pathway overrepresentation and GSEA
pdf(file = "SuppFigS21_salmonella_kegg_pathway.pdf", width = 30, height = 11)
grid.arrange(plot.pathway.overrep, plot.pathway.gsea, nrow = 1)
dev.off()

# make final merged plots: module overrepresentation and GSEA
pdf(file = "SuppFigS22_salmonella_kegg_module.pdf", width = 30, height = 11)
grid.arrange(plot.module.overrep, plot.module.gsea, nrow = 1)
dev.off()



