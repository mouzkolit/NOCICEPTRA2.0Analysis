library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(dplyr)
library(tibble)
library(xml2)
library(clusterProfiler)
library(DOSE)
library(data.table)
library(ggplot2)
library("RColorBrewer")
library("gplots")
library(org.Hs.eg.db)
library("biomaRt")
library(stringr)

outputPrefix = "human_drg"

# load the gtf file 
txdb <- makeTxDbFromGFF("gencode.v39.annotation.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- tx2gene[, 2:1]
tx2gene$TXNAME <- str_split_fixed(tx2gene$TXNAME,"\\.",3)[,1]
tx2gene$GENEID <- str_split_fixed(tx2gene$GENEID,"\\.",3)[,1]

# provide gene level estimation of counts using salmon
samples <- read.csv("DRG_samples.csv", header = TRUE, check.names = F)
files <- file.path(paste(samples$samples, ".gz_quant", sep = ""), "quant.sf")
names(files) <- paste0(samples$samples)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)


#give sample names 
#please add here the name of the samples

sampleNames <- samples$samples
condition = samples$condition
sampleTable <- data.frame(sampleName = sampleNames, condition = condition)

# # run the differential expression analysis 
dds <- DESeqDataSetFromTximport(txi.salmon, colData = sampleTable, design = ~1)
#dds <- DESeq(dds, test = "LRT", reduced=~1)

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

vsd <- varianceStabilizingTransformation(dds, blind=F)
vsd_counts = as.data.frame(assay(vsd))

mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",         
                            mart    = useMart("ENSEMBL_MART_ENSEMBL"))      
a_sig <- vsd_counts
b_sig <- getBM(attributes=c("ensembl_gene_id",
                            "external_gene_name",
                            "gene_biotype"
                            
),
filters = "ensembl_gene_id",
values=rownames(a_sig),
mart=mart)

# retrieve only the lncRNAs from the analysis
resdata_final <- merge(a_sig, b_sig, by.x = 0, by.y = "ensembl_gene_id")
annotation_slice = resdata_final[resdata_final$gene_biotype == "lncRNA",]
rownames(annotation_slice) = annotation_slice$Row.names
annotation_slice = annotation_slice[,2:7]

# read the iPSC lncRNA table to construct a correlation matrix
# 
lnc_iPSC = read.csv("Salmon_lnc_counts_stabilized_counts.csv")
rownames(lnc_iPSC) = lnc_iPSC$Row.names
lnc_iPSC = lnc_iPSC[,5:58]
columns = c(rep("Day00", each = 9), 
            rep("Day05", each = 9),
            rep("Day09", each = 9),
            rep("Day16", each = 9),
            rep("Day26", each = 9),
            rep("Day36", each = 9))

colnames(lnc_iPSC) <- columns


merged_human_ipsc = merge(lnc_iPSC, annotation_slice, by = 0)
rownames(merged_human_ipsc) = merged_human_ipsc$Row.names
merged_human_ipsc = merged_human_ipsc[,2:61]

library(reshape2)
# creating correlation matrix
corr_mat <- cor(merged_human_ipsc)
corr_mat_trial <- corr_mat[1:54,55:60]
# reduce the size of correlation matrix
melted_corr_mat <- melt(corr_mat_trial)
colnames(melted_corr_mat) <- c("iPSC_Samples", "DRG_Samples", "correlation")
# head(melted_corr_mat)

library(ggplot2)
ggplot(data = melted_corr_mat, aes(x=iPSC_Samples, y=DRG_Samples,fill=correlation)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_gradient(low = "#fde725",
                      high = "#440154",
                      guide = "colorbar")+
  coord_fixed(ratio = 2)
