library(ggplot2)
library(KEGGREST)
library(readxl)

###################### set seed

set.seed(42)

###################### load final TSV files

# load data frame of enriched kegg pathways for upregulated genes
final.pathway.up <- read.delim(file = "kegg_pathway_up/listeria_kegg_pathway_up.tsv",
                               header = T, sep = "\t",
                               stringsAsFactors = F, check.names = F)
head(final.pathway.up)
tail(final.pathway.up)

# load data frame of enriched kegg pathways for downregulated genes
final.pathway.down <- read.delim(file = "kegg_pathway_down/listeria_kegg_pathway_down.tsv",
                                 header = T, sep = "\t",
                                 stringsAsFactors = F, check.names = F)
head(final.pathway.down)
tail(final.pathway.down)

# load final GSEA data frame for kegg pathways
final.pathway.gsea <- read.delim(file = "kegg_pathway_gsea/listeria_kegg_pathway_gsea.tsv",
                                 header = T, sep = "\t",
                                 stringsAsFactors = F, check.names = F)
head(final.pathway.gsea)
tail(final.pathway.gsea)

# load final data frame with kegg module enrichment results for upregulated genes
final.module.up <- read.delim(file = "kegg_module_up/listeria_kegg_module_up.tsv",
                              header = T, sep = "\t",
                              stringsAsFactors = F, check.names = F)
head(final.module.up)
tail(final.module.up)

# load final data frame with kegg module enrichment results for downregulated genes
final.module.down <- read.delim(file = "kegg_module_down/listeria_kegg_module_down.tsv",
                                header = T, sep = "\t",
                                stringsAsFactors = F, check.names = F)
head(final.module.down)
tail(final.module.down)

# load final GSEA data frame for kegg modules
final.module.gsea <- read.delim(file = "kegg_module_gsea/listeria_kegg_module_gsea.tsv",
                                header = T, sep = "\t",
                                stringsAsFactors = F, check.names = F)
head(final.module.gsea)
tail(final.module.gsea)

###################### plot kegg pathways

# for kegg pathways, take -log10 of the raw p-values
final.pathway.up$log10p <- -log10(final.pathway.up$pvalue)
final.pathway.down$log10p <- log10(final.pathway.down$pvalue)

# merge overrepresented pathways
final.pathway.overrep <- rbind(final.pathway.up, final.pathway.down)
final.pathway.overrep$direction <- c(rep("Up-regulated", nrow(final.pathway.up)),
                                         rep("Down-regulated", nrow(final.pathway.down)))
head(final.pathway.overrep)

# get KEGG pathway names
KEGGREST_Name <- c()
for (i in 1:length(final.pathway.overrep$ID)){
  pathway.id <- final.pathway.overrep$ID[i]
  pathway.get <- keggGet(dbentries = pathway.id)
  KEGGREST_Name <- c(KEGGREST_Name, pathway.get[[1]]$NAME)
}

KEGGREST_Name
final.pathway.overrep$KEGGREST_Name <- KEGGREST_Name

# create KEGG_term column
final.pathway.overrep$KEGG_term <- paste(final.pathway.overrep$Description, final.pathway.overrep$KEGGREST_Name, sep = ", ")
head(final.pathway.overrep)

# overview of results
table(final.pathway.overrep$direction)
summary(final.pathway.overrep$log10p)
table(final.pathway.overrep$treatment)

# make bar chart of enriched kegg pathways
pdf(file = "plot_pathway_overrep_all.pdf",width = 11, height = 8.5)
plot.pathway.overrep <- ggplot(final.pathway.overrep, aes(x=reorder(KEGG_term, log10p), y=log10p, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Pathway") +
    ylab("(-)log10 P-Value") +
    scale_y_continuous(limits = c(-10, 18), breaks = seq(-10, 18, by = 2)) +
    scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
    theme_bw(base_size=12) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
plot.pathway.overrep
dev.off()

###################### plot kegg modules

# load kegg modules
# obtained from old version of kegg (12 December 2017)
# many modules assigned by eggnog mapper are not in the current version of kegg
module.dict <- read_excel(path = "module_html/kegg_module_wayback_12122017.xlsx", col_names = F)
head(module.dict)
module.dict <- data.frame(module.dict$...1, module.dict$...18)
colnames(module.dict) <- c("Module", "Name")
head(module.dict)

# for kegg modules, take -log10 of the raw p-values
final.module.up$log10p <- -log10(final.module.up$pvalue)
final.module.down$log10p <- log10(final.module.down$pvalue)

# merge overrepresented modules
final.module.overrep <- rbind(final.module.up, final.module.down)
final.module.overrep$direction <- c(rep("Up-regulated", nrow(final.module.up)),
                                     rep("Down-regulated", nrow(final.module.down)))
head(final.module.overrep)

# get KEGG module names
KEGGREST_Name <- c()
for (i in 1:length(final.module.overrep$ID)){
  module.id <- final.module.overrep$ID[i]
  if (module.id%in%module.dict$Module){
  module.get <- module.dict[which(module.dict$Module==module.id), "Name"]}
  else{
    module.get <- "NA"}
  KEGGREST_Name <- c(KEGGREST_Name, module.get)
}

final.module.overrep$KEGGREST_Name <- KEGGREST_Name
final.module.overrep$KEGGREST_Name

# create KEGG_term column
final.module.overrep$KEGG_term <- paste(final.module.overrep$Description, final.module.overrep$KEGGREST_Name, sep = ", ")
head(final.module.overrep)

# overview of results
table(final.module.overrep$direction)
summary(final.module.overrep$log10p)
table(final.module.overrep$treatment)

# make bar chart of enriched kegg modules
pdf(file = "plot_module_overrep_all.pdf",width = 11, height = 8.5)
plot.module.overrep <- ggplot(final.module.overrep, aes(x=reorder(KEGG_term, log10p), y=log10p, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Module") +
    ylab("(-)log10 P-Value") +
    scale_y_continuous(limits = c(-5, 9), breaks = seq(-5, 9, by = 1)) +
  scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
  theme_bw(base_size=12) +
guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
plot.module.overrep
dev.off()

###################### plot kegg pathway GSEA

# get KEGG pathway names
KEGGREST_Name <- c()
for (i in 1:length(final.pathway.gsea$ID)){
  pathway.id <- final.pathway.gsea$ID[i]
  print(pathway.id)
  if (pathway.id=="ko01130"|pathway.id=="map01130"){
    # from https://web.archive.org/web/20171215045609/https://www.genome.jp/kegg/pathway.html
    # December 15, 2017
    # 01130Biosynthesis of antibiotics
    pathway.get <- "Biosynthesis of antibiotics"}
  else{
  pathway.get <- keggGet(dbentries = pathway.id)
  pathway.get <- pathway.get[[1]]$NAME}
  KEGGREST_Name <- c(KEGGREST_Name, pathway.get)
}

final.pathway.gsea$KEGGREST_Name <- KEGGREST_Name
final.pathway.gsea$KEGGREST_Name

# create KEGG_term column
final.pathway.gsea$KEGG_term <- paste(final.pathway.gsea$Description, final.pathway.gsea$KEGGREST_Name, sep = ", ")
head(final.pathway.gsea)

# overview of results
summary(final.pathway.gsea$pvalue)
summary(final.pathway.gsea$p.adjust)
table(final.pathway.gsea$treatment)
summary(final.pathway.gsea$NES)

# plot pathway GSEA results
pdf(file = "plot_pathway_gsea_all.pdf",width = 11, height = 8.5)
plot.pathway.gsea <- ggplot(final.pathway.gsea, aes(x=reorder(KEGG_term, NES), y=NES, fill = treatment)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("KEGG Pathway") +
    ylab("NES") +
    scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
  theme_bw(base_size=12) +
guides(colour=guide_legend(override.aes=list(size=2.5))) +
    theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
plot.pathway.gsea
dev.off()

###################### plot kegg module GSEA

# get KEGG module names
KEGGREST_Name <- c()
for (i in 1:length(final.module.gsea$ID)){
  module.id <- final.module.gsea$ID[i]
  if (module.id%in%module.dict$Module){
    module.get <- module.dict[which(module.dict$Module==module.id), "Name"]}
  else{
    module.get <- "NA"}
  KEGGREST_Name <- c(KEGGREST_Name, module.get)
}

final.module.gsea$KEGGREST_Name <- KEGGREST_Name
final.module.gsea$KEGGREST_Name

# create KEGG_term column
final.module.gsea$KEGG_term <- paste(final.module.gsea$Description, final.module.gsea$KEGGREST_Name, sep = ", ")
head(final.module.gsea)

# overview of results
summary(final.module.gsea$pvalue)
summary(final.module.gsea$p.adjust)
table(final.module.gsea$treatment)
summary(final.module.gsea$NES)

# plot module GSEA results
pdf(file = "plot_module_gsea_all.pdf",width = 11, height = 8.5)
plot.module.gsea <- ggplot(final.module.gsea, aes(x=reorder(KEGG_term, NES), y=NES, fill = treatment)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("KEGG Module") +
  ylab("NES") +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  scale_fill_manual(values = c("ultrasound" = "#AA3377", "chlorine-ultrasound" = "#CCBB44", "chlorine" = "#228833")) +
  theme_bw(base_size=12) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
plot.module.gsea
dev.off()

# make final merged plots: pathway overrepresentation and GSEA
pdf(file = "SuppFigS19_listeria_kegg_pathway.pdf", width = 30, height = 11)
grid.arrange(plot.pathway.overrep, plot.pathway.gsea, nrow = 1)
dev.off()
