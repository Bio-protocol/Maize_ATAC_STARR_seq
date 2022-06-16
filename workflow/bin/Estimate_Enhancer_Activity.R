# estimate enhancer activity #

# load data
a <- read.table("B73_maize_mRNA_DNA.activity.raw.bed")

# reformat
rownames(a) <- paste(a$V1,a$V2,a$V3,sep="_")
a[,1:3] <- NULL
a <- as.matrix(a)
colnames(a) <- c("mRNA", "DNA")
a <- a[rowSums(a)!=0,]
a <- a + 1

# normalize
a <- a %*% diag(x=1e6/colSums(a))
colnames(a) <- c("mRNA", "DNA")
a <- as.data.frame(a)


# estimate enhancer activity
a$enhancer_activity <- log2(a$mRNA/a$DNA)

# reformat output
rownames(a) <- gsub("scaf_","scaf", rownames(a))
df <- data.frame(do.call(rbind, strsplit(rownames(a), "_")))
df$X1 <- gsub("scaf", "scaf_", as.character(df$X1))
mrna <- df
dna <- df
df$X4 <- a$enhancer_activity
mrna$X4 <- a$mRNA
dna$X4 <- a$DNA

# cap negative activity at 0
df$X4 <- ifelse(df$X4 < 0, 0, df$X4)

# save enhancer activity BEDGRAPH file to disk
write.table(df, file="B73_maize.enhancer_activity.bdg",quote=F, row.names=F, col.names=F, sep="\t")
write.table(mrna, file="B73_maize.mRNA.bdg",quote=F, row.names=F, col.names=F, sep="\t")
write.table(dna, file="B73_maize.DNA.bdg",quote=F, row.names=F, col.names=F, sep="\t")

