# estimate enhancer activity cell type specificity

# load libraries
library(RColorBrewer)
library(gplots)
library(edgeR)

# load data
enh <- read.table("STARR_starrs_peaks.enhancer_activity.eFDR05.ann.high_motif.scATAC_ACRs.bed")
acrs <- read.table("maize_scATAC_atlas_ACR_celltype_CPM.V5.txt")
motifs <- read.table("TFBS_peaks.motifs.ENRICHED.enhancer_activity.bed")

# subset for representative leaf-derived clusters
keep <- c("bulliform.2.26",
          "bundle_sheath.2.16",
          "ground_meristem.7.69",
          "guard_cell.7.74",
          "guard_mother_cell.7.71",
          "L1_SAM.4.46",
          "leaf_provascular.7.67",
          "mesophyll.2.14",
          "parenchyma.10.90",
          "protoderm.7.72",
          "stomatal_precursor.7.75",
          "subsidiary.7.68")
all.acrs <- acrs

# rescale acrs
acrs <- cpm(acrs, log=F)
acrs <- acrs[,keep]

# subset enhancers by genomic feature
enh <- subset(enh, enh$V8=="intergenic")

# get overlapping regions from the scATAC matrix
enh$ids <- paste(enh$V11,enh$V12,enh$V13,sep="_")
enh <- enh[order(enh$V11, decreasing=T),]
enh <- enh[!duplicated(enh$ids),]
shared <- intersect(enh$ids, rownames(acrs))
rownames(enh) <- enh$ids
enh <- enh[shared,]
enh$starrIDs <- paste(enh$V1, enh$V2, enh$V3,sep="_")

# filter motifs
motifs$starrIDs <- paste(motifs$V5, motifs$V6, motifs$V7, sep="_")
motifs <- motifs[motifs$starrIDs %in% unique(enh$starrIDs),]

# normalize acrs
#acrs.z <- as.matrix(t(scale(t(acrs))))
#acrs <- t(apply(acrs.z, 1, pnorm))
acrs <- t(apply(acrs, 1, function(x){x/max(x)}))

# iterate over each cell type
cts <- colnames(acrs)
outs <- lapply(cts, function(x){
    access <- acrs[rownames(enh),x]
    names(access) <- enh$starrIDs
    motif.scores <- access[motifs$starrIDs] * as.numeric(as.character(motifs$V15))
    mtf <- data.frame(motif=motifs$V4, score=motif.scores)
    aves <- aggregate(score~motif, data=mtf, FUN=mean)
    score <- aves$score
    names(score) <- aves$motif
    return(score)
})
outs <- do.call(cbind, outs)
colnames(outs) <- cts
vars <- apply(outs, 1, var)
outs <- outs[vars > 0,]
z <- as.matrix(t(scale(t(outs))))

# check correlation among motifs
cors <- cor(t(z))

#
maxes <- apply(z, 1, max)
#z <- z[maxes > 1.5,]
#z[z > 3] <- 3
#z[z < -3] <- -3

# get acr z-scores before reducing data set
# acrs.z <- t(as.matrix(scale(t(acrs))))
# acrs <- acrs[shared,]
# acrs.z <- acrs.z[shared,]
# 
# # rescale ACR accessibility for each cell type by maximum accessibility
# acrs <- t(apply(acrs.z, 1, pnorm))
# 
# # scale accessibility scores by enhancer activity
# acr.act <- t(t(as.matrix(acrs)) %*% diag(x=enh$V11))
# 
# # ave
# ct.score <- colMeans(acr.act)
# ct.score <- (ct.score - mean(ct.score)) / sd(ct.score)
# ct.score <- ct.score[order(ct.score, decreasing=T)]
# 
# # plot
# pdf("Celltype_STARR_score.pdf", width=12, height=8)
# barplot(ct.score, las=2, col=ifelse(ct.score > 1, "firebrick3", ifelse(ct.score < -1, "dodgerblue3", "grey75")),
#         cex.names=0.5)
# dev.off()
# 
# # convert to z-scores
# z <- as.matrix(t(scale(t(acrs))))
# 
# # cluster columns
co <- hclust(dist(t(outs)))$order
# 
# # reorder rows
z <- z[,co]
row.o <- apply(z, 1, which.max)
z <- z[order(row.o, decreasing=F),]
# 
# # cap
# z.og <- z
z[z < -3] <- -3
z[z > 3] <- 3

# get family
tfs <- data.frame(do.call(rbind, strsplit(rownames(z), "\\.")))
cols2 <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tfs$X1)))
tfs$cols2 <- cols2[factor(tfs$X1)]

# visualize
pdf("celltype_starr_bias.pdf", width=10, height=10)
heatmap.2(z, scale="none", trace='none', 
          RowSideColors=tfs$cols,
          col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
          useRaster=T, Colv=F, Rowv=F, dendrogram="none", margins=c(9,9))
dev.off()