# load libraries
library(scales)

# load data
starr <- read.table("STARR_merged_peaks.enhancer_activity.bed")
con <- read.table("STARR_CONTROL.enhancer_activity.bed")

# set missing to 0
starr$V4[starr$V4=='.'] <- 0
con$V4[con$V4=='.'] <- 0

# convert to numeric
starr$V4 <- as.numeric(as.character(starr$V4))
con$V4 <- as.numeric(as.character(con$V4))

# get empirical thresholds
fdr <- 0.05
threshold <- quantile(con$V4, (1-fdr))

# filter STARR regulatory regions
filtered <- subset(starr, starr$V4 >= threshold)

# estimate fraction of retained regions
frac <- signif(nrow(filtered)/nrow(starr), digits=4)

# set up multipanel plot area
pdf("Density_eFDR_STARR_Peak_Filtering.pdf", width=5, height=5)

# plot control/observed enhancer activities for STARR peaks with duplicates
den.starr <- density(starr$V4)
den.con <- density(con$V4)
plot(NA, 
     xlab="Enhancer Activity",
     ylab="Density",
     xlim=c(range(range(den.starr$x), range(den.con $x))),
     ylim=c(range(range(den.starr$y), range(den.con$y))))
polygon(x=c(min(den.starr$x), den.starr$x, max(den.starr$x)),
        y=c(0, den.starr$y, 0), col=alpha("darkorchid4", 0.5), border=NA)
polygon(x=c(min(den.con$x), den.con$x, max(den.con$x)),
        y=c(0, den.con$y, 0), col=alpha("grey80", 0.5), border=NA)
abline(v=threshold, col="red", lty=2)
mtext(paste0("STARR peaks pass = ",frac," (", nrow(filtered), "/", nrow(starr),")"))

legend("right", legend=c("STARR Peaks", "Control Peaks", paste0("eFDR = ", fdr)), col=c("darkorchid4", "grey75", "red"), border=c(NA, NA, "red"), pch=c(16, 16, NA), lty=c(NA, NA, 2))

# close device
dev.off()

# save filtered STARR regulatory regions
write.table(filtered, file="STARR_merged_peaks.enhancer_activity.eFDR05.bed", quote=F, row.names=F, col.names=F, sep="\t")

