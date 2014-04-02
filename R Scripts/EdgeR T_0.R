# Main Pipeline Using EdgeR
# (WT vs sxy), at time = 0.
# by Nathaniel (April 2014)

library(limma)
library(edgeR)

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
# Wild Type and Sxy (K and S); Time = 0
small.meta <- subset(frag.matrix, frag.matrix$Type %in% c('K', 'S') & frag.matrix$Timepoint == '0')
small.meta$Type <- factor(small.meta$Type)
small.raw <- raw.reads[,colnames(raw.reads) %in% small.meta$Sample]

# Unloading Unneeded Files
rm(raw.reads)
rm(frag.matrix)

# Running Dataset through EdgeR pipeline
# 1. Creating DGEList from Trimmed Data
dge.list <- DGEList(counts = small.raw, group = small.meta$Timepoint)

# 2. Creating Design Matrix
sm.design <- model.matrix(~Type, small.meta)

# 3. ???
dge.common.list <- estimateGLMCommonDisp(dge.list, sm.design, verbose = FALSE)
dge.trend.list <- estimateGLMTrendedDisp(dge.common.list, sm.design)
dge.tag.list <- estimateGLMTagwiseDisp(dge.trend.list, sm.design)

# 4. Fitting and Throwing Out Hits
# **Current Null Hypothesis: TypeS (i.e. WT against sxy at time = 30)
dge.fit <- glmFit(dge.tag.list, sm.design)
dge.lrt <- glmLRT(dge.fit, coef = 'TypeS')
dge.top <- topTags(dge.lrt, n = Inf)
top.hits <- rownames(dge.top$table[dge.top$table$FDR < 0.05,])
dge.cpm <- cpm(dge.tag.list)[top.hits,]
output.genes <- gene.list[gene.list$ID %in% top.hits,]
rownames(output.genes) <- NULL

# Creating crp-S Subset List
crps.list <- subset(crp.list, prom %in% c('crps', 'sxy'))
rownames(crps.list) <- NULL

# Displaying the top hits
output.genes

# Comparing if all crp-S entries fit the RNA-Seq Data
nrow(crps.list)
nrow(output.genes[output.genes$ID %in% crps.list$orf,])
output.genes[output.genes$ID %in% crps.list$orf,]