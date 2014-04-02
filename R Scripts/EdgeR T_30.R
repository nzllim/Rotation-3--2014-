# Main Pipeline Using EdgeR
# (WT vs sxy), at time = 30.
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
# Wild Type and Sxy (K and S); Time = 30
small.meta <- subset(frag.matrix, frag.matrix$Type %in% c('K', 'S') & frag.matrix$Timepoint == '30')
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

# =================================================
# Progress 2 Starts

# Displaying Expression Differential Summary
summary(dge.summ <- decideTestsDGE(dge.lrt, p = 0.05, adjust = 'BH'))

# Searching for UP-regulated genes
dge.up <- dge.lrt$table[(dge.summ == 1),]
output.up.genes <- gene.list[gene.list$ID %in% rownames(dge.up),]
show(output.up.genes)
show(data.frame(ID = rownames(dge.up), logFC = dge.up$logFC))
nrow(dge.up[rownames(dge.up) %in% crps.list$orf,])

# Searching for DOWN-regulated genes
dge.down <- dge.lrt$table[(dge.summ == -1),]
output.down.genes <- gene.list[gene.list$ID %in% rownames(dge.down),]
show(output.down.genes)
show(data.frame(ID = rownames(dge.down), logFC = dge.down$logFC))
nrow(dge.down[rownames(dge.down) %in% crps.list$orf,])

# Searching for the non-differential crp-S gene
output.odd.crps <- crps.list[-which(crps.list$orf %in% top.hits),]
output.odd.crps <- gene.list[gene.list$ID %in% output.odd.crps$orf,]
show(output.odd.crps)