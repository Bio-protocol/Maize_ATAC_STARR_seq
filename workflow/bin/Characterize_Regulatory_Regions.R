# Analyze regulatory regions 

# load libraries
library(vioplot)
library(dplyr)
library(MASS)
library(RColorBrewer)
library(scales)

# load data
starr <- read.table("STARR_merged_peaks.enhancer_activity.eFDR05.ann.bed")
control <- read.table("STARR_CONTROL.enhancer_activity.ann.bed")

# select random control regions to match the filtered STARR peaks
control <- control[sample(nrow(starr)),]

# rename columns for clarity (frac_RR_motif = fraction of regulatory region covered by motifs)
starr[,7:15] <- NULL
control[,7:15] <- NULL
colnames(starr)[4:7] <- c("activity", "motif_counts", "frac_RR_motif", "gene_distance")
colnames(control)[4:7] <- c("activity", "motif_counts", "frac_RR_motif", "gene_distance")

# classify
starr$class <- ifelse((starr$gene_distance < 0 & starr$gene_distance > -200), "TSS", 
                     ifelse(starr$gene_distance < -200 & starr$gene_distance > -2000, "promoter", 
                            ifelse(starr$gene_distance > 0 & starr$gene_distance < 1000, "TTS", 
                                   ifelse(starr$gene_distance == 0, "genic", "intergenic"))))
control$class <- ifelse((control$gene_distance < 0 & control$gene_distance > -200), "TSS", 
                     ifelse(control$gene_distance < -200 & control$gene_distance > -2000, "promoter", 
                            ifelse(control$gene_distance > 0 & control$gene_distance < 1000, "TTS", 
                                   ifelse(control$gene_distance == 0, "genic", "intergenic"))))

# plot distribution
pdf("STARR_peak_control_genomic_distribution.pdf", width=10, height=5)
layout(matrix(c(1:2), nrow=1))
pie(table(starr$class))
pie(table(control$class))
dev.off()

# estimate regulatory region size (in log10 scale)
starr$size <- log10(starr$V3-starr$V2)
control$size <- log10(control$V3-control$V2)

# compare sizes between peaks and controls (sanity check)
pdf("STARR_peak_control_sizes.pdf", width=5, height=6)
vioplot(starr$size, control$size, 
        ylab="Interval size (log10)",
        col=c("dodgerblue", "grey75"),
        names=c(paste0("STARR peaks \n (n=",nrow(starr),")"),
                paste0("Control regions \n (n=",nrow(control),")")))
dev.off()

# compare motif counts
pval <- wilcox.test(starr$motif_counts, control$motif_counts)$p.value
pval <- ifelse(pval==0, 2.2e-16, pval)
mean.peak <- mean(starr$motif_counts)
mean.cont <- mean(control$motif_counts)

# find 95% quantile for control motif count 
upper.threshold <- quantile(control$motif_counts, 0.95)

# plot
pdf("STARR_peak_control_motif_counts.pdf", width=5, height=6)
vioplot(log1p(starr$motif_counts), log1p(control$motif_counts), 
        ylab="log2(Motif counts + 1)",
        col=c("dodgerblue", "grey75"),
        names=c(paste0("STARR peaks \n (n=",nrow(starr),")"),
                paste0("Control regions \n (n=",nrow(control),")")),
        ylim=c(0,8),
        areaEqual=T,
        h=0.25)
mtext(paste0("Wilcoxon Rank Sum P-value = ", signif(pval, digits=3)))
text(1, 7.5, labels=paste0("Mean = ", signif(mean.peak, digits=3)))
text(2, 7.5, labels=paste0("Mean = ", signif(mean.cont, digits=3)))
points(1, log1p(upper.threshold), col="red", pch="-")
points(2, log1p(upper.threshold), col="red", pch="-")
dev.off()

# split STARR regions by motif counts based on 95% quantile control dist
starr$group <- ifelse(starr$motif_counts >= upper.threshold, "high", "low")
pval <- kruskal.test(starr$activity, starr$group)$p.value
pdf("STARR_peak_activity_vs_group.pdf", width=5, height=6)
vioplot(starr$activity~starr$group,
        ylab="Enhancer activity",
        col=c("dodgerblue4", "dodgerblue"),
        names=c(paste0("Motif-enriched \n STARR peaks \n (n=",nrow(starr[starr$group=="high",]),")"),
                paste0("Motif-depleted \n STARR peaks \n (n=",nrow(starr[starr$group=="low",]),")")),
        areaEqual=F,
        xlab="",
        h=0.25)
mtext(paste0("Kruskal-Wallis rank sum P-value = ", signif(pval, digits=3)))
dev.off()

# compare STARR region size
pval <- kruskal.test(starr$size, starr$group)$p.value
pval <- ifelse(pval==0, 2.2e-16, pval)
pdf("STARR_peak_size_vs_group.pdf", width=5, height=6)
vioplot(starr$size~starr$group,
        ylab="Interval size (log10)",
        col=c("dodgerblue4", "dodgerblue"),
        names=c(paste0("Motif-enriched \n STARR peaks \n (n=",nrow(starr[starr$group=="high",]),")"),
                paste0("Motif-depleted \n STARR peaks \n (n=",nrow(starr[starr$group=="low",]),")")),
        areaEqual=F,
        xlab="",
        h=0.25)
mtext(paste0("Kruskal-Wallis rank sum P-value = ", signif(pval, digits=3)))
dev.off()

# compare motif coverage
pval <- kruskal.test(starr$frac_RR_motif, starr$group)$p.value
pval <- ifelse(pval==0, 2.2e-16, pval)
pdf("STARR_peak_motif_coverage_vs_group.pdf", width=5, height=6)
vioplot(starr$frac_RR_motif~starr$group,
        ylab="Fraction motif coverage",
        col=c("dodgerblue4", "dodgerblue"),
        names=c(paste0("Motif-enriched \n STARR peaks \n (n=",nrow(starr[starr$group=="high",]),")"),
                paste0("Motif-depleted \n STARR peaks \n (n=",nrow(starr[starr$group=="low",]),")")),
        areaEqual=F,
        xlab="")
mtext(paste0("Kruskal-Wallis rank sum P-value = ", signif(pval, digits=3)))
dev.off()

# split by group
starr.me <- subset(starr, starr$group=="high")
starr.md <- subset(starr, starr$group=="low")
write.table(starr.me, file="STARR_starrs_peaks.enhancer_activity.eFDR05.ann.high_motif.bed",
            quote=F, row.names=F, col.names=F, sep="\t")
write.table(starr.md, file="STARR_starrs_peaks.enhancer_activity.eFDR05.ann.low_motif.bed",
            quote=F, row.names=F, col.names=F, sep="\t")

