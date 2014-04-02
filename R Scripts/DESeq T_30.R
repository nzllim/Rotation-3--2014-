# Main Pipeline Using DESeq
# (WT vs sxy), at time = 30.
# by Nathaniel (April 2014)

library(DESeq)

# Reset Workspace
rm(list = ls())

# Loads RNA-Seq Raw Fragment Count
raw.reads <- read.table('Data/fragment_counts.txt', header = TRUE, sep = '\t')

# Loads Design Matrix for limma
frag.matrix <- read.table('Data/FragmentMatrix.txt', header = TRUE, sep = '\t')

# Loads Truncated GeneList (made from rc.bed)
gene.list <- read.table('Data/GeneID.txt', header = TRUE)
rownames(gene.list) <- gene.list$ID

# Loads crp-S Target List
crp.list <- read.table('Data/crp.txt', header = TRUE)

# Associating GeneID as KEY identifier for all datasets
rownames(gene.list) <- gene.list$ID
rownames(raw.reads) <- gene.list$ID

# Creating Data Subset (with Restructuring)
# Wild Type and Sxy (K and S); Time = 30
small.meta <- subset(frag.matrix, frag.matrix$Type %in% c('K', 'S') & frag.matrix$Timepoint == '30')
small.meta$Type <- factor(small.meta$Type)
small.raw <- raw.reads[,colnames(raw.reads) %in% small.meta$Sample]

# Unloading Unneeded Files
rm(raw.reads)
rm(frag.matrix)

# Running Dataset through DESeq pipeline
deseq.dat <- newCountDataSet(small.raw, small.meta$Type)
deseq.dat <- estimateSizeFactors(deseq.dat)
deseq.dat <- estimateDispersions(deseq.dat)
deseq.output <- nbinomTest(deseq.dat, 'K', 'S')
rownames(deseq.output) <- deseq.output$id

# Throwing out Hits
# **Current Null Hypothesis: TypeS (i.e. WT against sxy at time = 30)
top.hits <- rownames(deseq.output[deseq.output$padj < 0.05, ])
output.genes <- gene.list[gene.list$ID %in% top.hits,]

# Creating crp-S Subset List
crps.list <- subset(crp.list, prom %in% c('crps', 'sxy'))
rownames(crps.list) <- NULL

# Displaying the top hits
output.genes

# Comparing if all crp-S entries fit the RNA-Seq Data
nrow(crps.list)
nrow(output.genes[output.genes$ID %in% crps.list$orf,])