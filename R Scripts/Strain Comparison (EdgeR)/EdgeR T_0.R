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

# Loads Truncated GeneList
gene.list <- read.table('Data/GeneID.txt', header = TRUE)

# Loads crp-S Target List
crp.list <- read.table('Data/crp.txt', header = TRUE)

# Associating GeneID as KEY identifier for all datasets
rownames(gene.list) <- gene.list$ID
rownames(raw.reads) <- gene.list$ID

# Creating Data Subset (with Restructuring)
# Wild Type and Sxy (K and S); Time = 0
small.meta <- subset(frag.matrix, frag.matrix$Type %in% c('K', 'S') & frag.matrix$Timepoint %in% '0')
small.meta$Type <- factor(small.meta$Type)
small.raw <- raw.reads[,colnames(raw.reads) %in% small.meta$Sample]

# Unloading Unneeded Files
rm(raw.reads)
rm(frag.matrix)

# Running Dataset through EdgeR pipeline
# 1. Creating DGEList from Trimmed Data
dge.list <- DGEList(counts = small.raw, group = small.meta$Type)

# 2. Creating Design Matrix
sm.design <- model.matrix(~Type, small.meta)

# 3. Continuing EdgeR Pipeline
dge.common.list <- estimateGLMCommonDisp(dge.list, sm.design)
dge.trend.list <- estimateGLMTrendedDisp(dge.common.list, sm.design)
dge.tag.list <- estimateGLMTagwiseDisp(dge.trend.list, sm.design)

# 4. Fitting and Throwing Out Hits
# **Current Null Hypothesis: TypeS (i.e. WT against sxy at time = 0)
dge.fit <- glmFit(dge.tag.list, sm.design)
dge.lrt <- glmLRT(dge.fit, coef = 'TypeS')
dge.top <- topTags(dge.lrt, n = Inf)
dge.cpm <- cpm(dge.tag.list)

final.top.table <- dge.top$table[dge.top$table$FDR < 0.05, ]
final.top.cpm <- dge.cpm[rownames(final.top.table),]
final.top.genes <- gene.list[gene.list$ID %in% rownames(final.top.table), ]

# Creating crp-S Subset List
crps.list <- subset(crp.list, prom %in% c('crps', 'sxy'))
rownames(crps.list) <- NULL

# Creating crp-N Subset List
crpn.list <- subset(crp.list, prom %in% 'crpn')
rownames(crpn.list) <- NULL

# ====================================
# Analysis Output Begins Here
# ====================================

# Number of differentially expressed genes:
nrow(final.top.genes)

# Number of down-regulated genes:
dge.down <- final.top.table[final.top.table$logFC < 0, ]
nrow(dge.down)

# Number of up-regulated genes:
dge.up <- final.top.table[final.top.table$logFC > 0, ]
nrow(dge.up)

# crp-S genes:
# Number of differentially expressed genes:
crps.diff <- final.top.table[rownames(final.top.table) %in% crps.list$orf, ]
nrow(crps.diff)

# Number of non-differentially expressed genes:
tmp.crps.non <- crps.list[-which(crps.list$orf %in% rownames(crps.diff)), ]
crps.non <- dge.top$table[as.character(tmp.crps.non$orf), ]
rm(tmp.crps.non)
nrow(crps.non)

# Number of down-regulated genes:
crps.down <- crps.diff[crps.diff$logFC < 0, ]
nrow(crps.down)

# Number of up-regulated genes:
crps.up <- crps.diff[crps.diff$logFC > 0, ]
nrow(crps.up)

# crp-N genes:
# Number of differentially expressed genes:
crpn.diff <- final.top.table[rownames(final.top.table) %in% crpn.list$orf, ]
nrow(crpn.diff)

# Number of non-differentially expressed genes:
tmp.crpn.non <- crpn.list[-which(crpn.list$of %in% rownames(crpn.diff)), ]
crpn.non <- dge.top$table[as.character(tmp.crpn.non$orf), ]
rm(tmp.crpn.non)
nrow(crpn.non)

# Number of down-regulated genes:
crpn.down <- crpn.diff[crpn.diff$logFC < 0, ]
nrow(crpn.down)

# Number of up-regulated genes:
crpn.up <- crpn.diff[crpn.diff$logFC > 0, ]
nrow(crpn.up)

# Displaying all differentially expressed genes
show(final.top.genes)

# Displaying all down-regulated genes
dge.down.genes <- gene.list[gene.list$ID %in% rownames(dge.down), ]
show(dge.down.genes)

# Displaying all up-regulated genes
dge.up.genes <- gene.list[gene.list$ID %in% rownames(dge.up), ]
show(dge.up.genes)

# Creating Clean Tables for exporting:
# Complete Normalized Counts
export.dge.cpm <- data.frame(gene.list[rownames(dge.cpm),], dge.cpm)
write.table(export.dge.cpm, file = 'EdgeR T0 Normalized Counts.tsv', sep = '\t', row.names = FALSE, col.names = TRUE)

# Differentially Expressed Genes (with Norm Counts)
export.final.top <- data.frame(final.top.genes, final.top.table[rownames(final.top.genes),], final.top.cpm[rownames(final.top.genes),])
write.table(export.final.top, file = 'EdgeR T0 Top Hits.tsv', sep = '\t', row.names = FALSE, col.names = TRUE)