rm(list=ls())
library(dplyr)
library(tibble)
library(xml2)
library(clusterProfiler)
library(data.table)
library(ggplot2)
library("RColorBrewer")
library("gplots")
library(stringr)
library(tidyr)
library(BioCircos)
library(EnhancedVolcano)
library("DESeq2")
library(hrbrthemes)
library(viridis)
library("pvca")
library("Biobase")
library("rtracklayer")
library("dplyr")
library("comprehenr")
library("genefilter")
library(reshape2)
library(ComplexHeatmap)
dev.off()

outputPrefix <- "excerpt_RNA_merged_hDRG"
object_excerpt <- load("./AgghDRGResults/exceRpt_smallRNAQuants_ReadCounts.RData")
barplot_ncRNA = read.table("./AgghDRGResults/exceRpt_biotypeCounts.txt")
barplot_ncRNA$biotype = rownames(barplot_ncRNA)

vector_biotypes = c("miRNA","piRNA","tRNA","rRNA","snoRNA","snRNA", "lincRNA","protein_coding","misc_RNA")

# make the vector for the new biotypes
empty_vector = c()
for(i in barplot_ncRNA$biotype){
  if(i %in% vector_biotypes){
    empty_vector <- append(empty_vector, i)
  }
  else{
    empty_vector <- append(empty_vector, "other")
  }
}

colnames_bar = c( "ReferenceID",
                  "hDRG_1",
                  "hDRG_2",
                  "hDRG_3",
                  "hDRG_4",
                  "hDRG_5",
                  "hDRG_6",
                  "hDRG_7",
                  "biotype")


# greate the biotype graph
barplot_ncRNA$biotype = empty_vector
colnames(barplot_ncRNA) = colnames_bar[2:9]
# make an aggregation of the plots to get the time course
barplot_agg = gather(barplot_ncRNA, Sample, Counts, 1:7) 

#agg data based on sample and biotype
barplot_bar = barplot_agg %>%
  group_by(Sample, biotype) %>%
  summarise(percentage = sum(Counts)) %>%
  mutate(freq = percentage/ sum(percentage))

# creating the ggplot object
group.colors <- c(miRNA = "#262424", 
                  lincRNA = "#eb1b0c",
                  rRNA ="#4a3d3d",
                  snoRNA = "#dcde57",
                  snRNA = "#3dd12a",
                  misc_RNA = "blue",
                  tRNA = "#aa34cf",
                  protein_coding = "#878787",
                  other = "lightgreen",
                  piRNA = "skyblue")

ggplot(barplot_bar, aes(x = Sample, y = freq, fill = biotype)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=group.colors) +
  xlab("Samples")+ 
  ylab("Frequency")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        legend.position="bottom",
        legend.box.background = element_rect(colour = "black"),
  ) +
  facet_grid(. ~ "Frequency per ncRNA type")

#make a pie chart for the different ncRNA types aggregated on the Day
pieplot_chart <- barplot_bar
pieplot_chart$time <- str_split_fixed(pieplot_chart$Sample,"_",2)[,1]


# group by the biotype and per timepoint for the creation of the pie-chart
pieplot_chart_trial <- pieplot_chart %>% 
  dplyr::group_by(time, biotype) %>% 
  dplyr::summarise(pie = sum(percentage)) %>%
  dplyr::mutate(freq = round(100* pie/ sum(pie),2))

# define the positions of the label within the pie chart
df2 <- pieplot_chart_trial %>% 
  mutate(csum = rev(cumsum(rev(freq))), 
         pos = freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), freq/2, pos))


#drawing the pie charts with custom 
ggplot(data=as.data.frame(pieplot_chart_trial), aes(x="", y=freq, color = biotype, fill = biotype)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  coord_polar("y") + 
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(freq, "%")),
                   size = 3, nudge_x = 1, show.legend = FALSE,
                   color = "white"
  ) +
  guides(fill = guide_legend(title = "biotype")) +
  facet_wrap(.~ time, ncol = 3)+
  theme(legend.position="bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        strip.background = element_rect(colour="black", fill = "grey"))



trna_expression <- read.csv("tRNA_files_expression.csv")
trna_expression <- trna_expression[,2:9]
# reduce heterogenity by merging and aggregating counts belonging to the
# same isoacceptor
rownames(trna_expression) = trna_expression$X
trna_expression$trna_acceptor = str_split_fixed(rownames(trna_expression), "-(?=[^-]+$)",2)[,1] 
trna_expression = trna_expression %>% group_by(trna_acceptor) %>% summarise_if(is.numeric,sum)
trna_expression = as.data.frame(trna_expression)
colnames(trna_expression) = colnames_bar[1:8]
rownames(trna_expression) = trna_expression$ReferenceID
trna_expression = trna_expression[2:8]
trna_expression$biotype = "tRNA"


####################################################################
biotype_gen <- as.data.frame(as.matrix(str_split_fixed(rownames(exprs.gencode),":",2)))
exprs.gencode <- as.data.frame(exprs.gencode)
exprs.miRNA <- as.data.frame(exprs.miRNA)
exprs.piRNA <- as.data.frame(exprs.piRNA)
exprs.piRNA$biotype = "piRNA"
exprs.miRNA$biotype = "miRNA"
exprs.gencode$biotype = biotype_gen$V2
exprs.gencode = exprs.gencode[exprs.gencode$biotype != "miRNA",]
colnames(exprs.piRNA) = colnames_bar[2:9]
colnames(exprs.miRNA) = colnames_bar[2:9]
colnames(exprs.gencode)= colnames_bar[2:9]
cts <- rbind(exprs.gencode,exprs.miRNA, exprs.piRNA, trna_expression)

# make an annotation data table
annotation = cts
annotation$mean = rowMeans(annotation[1:7], na.rm = FALSE)
annotation =  annotation[annotation$mean > 1,]
annotation$ReferenceId = rownames(annotation)
annotation = annotation[,c("ReferenceId","biotype")] 

# remove the biotype again for differential gene expression analysis
cts = subset(cts, select = -c(biotype))
cts = cts[rownames(cts) %in% annotation$ReferenceId,]
cts <- mutate_all(cts, function(x) as.integer(as.character(x)))


#build the deseq2 object iwth metadata and design
cell_line = rep(c("hDRG"), each = 7)
sampleTable <- data.frame(sampleName = colnames(cts), cellline = cell_line)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sampleTable, 
                              design = ~1)
dds <- estimateSizeFactors(dds)
colData(dds)<- sampleTable

# get the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- merge(normalized_counts, annotation, by.x = 0, by.y = 0)
normalized_counts$mean <- rowMeans(normalized_counts[,2:8], na.rm=T)
normalized_counts$symbol <- str_split_fixed(normalized_counts$Row.names, ":", 2)[,1]
empty_vector = c()
for(i in normalized_counts$biotype){
  if(i %in% vector_biotypes){
    empty_vector <- append(empty_vector, i)
  }
  else{
    empty_vector <- append(empty_vector, "other")
  }
}  
normalized_counts$biotype = empty_vector


##
# Check hDRG expression patterns
# Most significant ncRNAs per biotype
##
# make a groupby in for biotype and Row.names
data = normalized_counts %>%
  dplyr::arrange(desc(mean)) %>%
  select(biotype, symbol, mean) %>%
  dplyr::group_by(biotype, symbol) %>%
  summarise(top_mean = sum(mean)) %>%
  top_n(20, top_mean) %>%
  arrange(desc(top_mean))

data$symbol <- str_split_fixed(data$symbol,"gb", 2)[,1]
# draw the top expressed ncRNA per aggregation

library(ggh4x)
ggplot(data, aes(x = reorder(symbol, -top_mean), y = top_mean))+
  geom_col(fill="forest green")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(biotype), nrow = 4, ncol = 4, scales = "free")



# get the variance stabilized counts
vsd <- varianceStabilizingTransformation(dds, blind=T)


# make a sample to sample PCA and clustering 
pca <- plotPCA(vsd, returnData = T, intgroup = "sampleName")
percentVar <- round(100 * attr(pca, "percentVar"))
head(pca)

ggplot_pca <- function(data, percentVar){
  p1 <- ggplot(data, aes(PC1, PC2, color = sampleName))+
    geom_point(size=3)+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))+
    stat_ellipse()+
    facet_grid(. ~ "hDRG ncRNA distribution")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p1)
}
p1 <- ggplot_pca(pca, percentVar)
p1

## 
# Hierachical Clustering of DRG samples
##
vsd_cor = cor(assay(vsd))
library(pheatmap)

### Plot heatmap
pheatmap(vsd_cor, annotation = sampleTable, color=colorRampPalette(c("purple", "black", "orange"))(75))
#######
#Compare with iPSC-derives senory neuronshttp://127.0.0.1:17195/graphics/cf46ff5e-ed97-4d37-a468-a7433db98152.png
# load vsd file from the ncRNA analysis of iPSC-derived sensory neurons
iPSC_vsd = read.csv("excerpt_RNA_-vsd.csv", header = T, row.names = 1)
iPSC_normalized = read.csv("excerpt_RNA_merged_trnanormalized_expression_counts.csv", header = T)
iPSC_normalized$mean = rowMeans(iPSC_normalized[,c("D26_1",
                                                   "D26_2",
                                                   "D26_3",
                                                   "D26_4",
                                                   "D26_5",
                                                   "D26_6",
                                                   "D26_7",
                                                   "D26_8",
                                                   "D26_9",
                                                   "D36_1",
                                                   "D36_2",
                                                   "D36_3",
                                                   "D36_4",
                                                   "D36_5",
                                                   "D36_6",
                                                   "D36_7",
                                                   "D36_8",
                                                   "D36_9")])
iPSC_normalized = merge(iPSC_normalized, annotation, by.x = "Row.names", by.y = 0)

empty_vector = c()
for(i in iPSC_normalized$biotype){
  if(i %in% vector_biotypes){
    empty_vector <- append(empty_vector, i)
  }
  else{
    empty_vector <- append(empty_vector, "other")
  }
}  
iPSC_normalized$biotype = empty_vector

#####
# get the top 500 expressed ncRNAs in iPSC and DRG
####
# draw a venn diagram
library(VennDiagram)
top_500_hDRG = normalized_counts[order(normalized_counts$mean, decreasing = T),][1:500,]
top_500_iPSC = iPSC_normalized[order(iPSC_normalized$mean, decreasing = T),][1:500,]
intersection = intersect(top_500_iPSC$Row.names, top_500_hDRG$Row.names)

p1 <- venn.diagram(
  x = list(top_500_iPSC$Row.names, top_500_hDRG$Row.names),
  category.names = c("iDN Nociceptor" , "hDRGs"),
  lwd = 2,
  lty = 'blank',
  fill = c("#56B4E9", "#009E73"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055),
  filename = NULL
)

grid.draw(p1)

dev.off()
# rearrange the iPSC derived sensory neuron data
data_iPSC = iPSC_normalized %>%
  dplyr::arrange(desc(mean)) %>%
  select(biotype, Row.names, mean) %>%
  dplyr::group_by(biotype, Row.names) %>%
  summarise(top_mean = sum(mean)) %>%
  top_n(20, top_mean) %>%
  arrange(desc(top_mean))

#retrieve the intersection between the top iPSC-derived sensory neuron ncRNAs and hDRGS ncRNAs
data_iPSC$Row.names <- str_split_fixed(data_iPSC$Row.names, ":", 2)[,1]
data_iPSC$Row.names <- str_split_fixed(data_iPSC$Row.names,"gb", 2)[,1]




ggplot(data_iPSC, aes(x = reorder(Row.names, -top_mean), y = top_mean))+
  geom_col(fill="purple")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(biotype), nrow = 4, ncol = 4, scales = "free")

# calculate the groups ike for hDRGs


DRG_vsd = vsd[rownames(vsd) %in% rownames(iPSC_vsd)]

# retrieve the correlation matrix
matrix_DRG = as.data.frame(assay(DRG_vsd))
matrix_all = merge(matrix_DRG, iPSC_vsd, by= 0, all.x = T)
rownames(matrix_all) = matrix_all$Row.names
matrix_all = matrix_all[2:62]
order_number <- str_split_fixed(colnames(matrix_all),"_", 2)[,1][8:61]
ordering <- factor(order_number, levels = unique(order_number))
colnames(cormat) <- ordering

# plot the correlation matrix
cormat <- round(cor(matrix_all),2)
cormat <- cormat[1:7,8:61]
colnames(cormat) <- ordering
cormat <- cormat[,order(colnames(cormat))]

annotation <- data.frame(Var1 = colnames(cormat), Var2 = ordering)
annotation

pheatmap::pheatmap(cormat,
                   cluster_cols=F,
                   clusterin_distance_rows = "correlation",
                   annotation_legend = F,
                   color = hcl.colors(50, "Lisbon"),
                   border_color="grey",
                   show_colnames = F,
                   fontsize_number = 7,
                   gaps_col = c(9,18,27,36,45)
)
