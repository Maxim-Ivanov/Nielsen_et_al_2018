# Load the required libraries:
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(SummarizedExperiment)
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
library(edgeR)
library(reshape2)
library(ggplot2)

# Load the RangedSummarizedExperiment object with CAGEfightR output (see the 02-Calling_TSS_with_CAGEfightR.R pipeline):
rse <- readRDS("rse_FACT.RData")
mc <- mcols(rse)
table(mc$txType)
table(mc$wt_fact)

# Figure 5A (stacked barplot of TSS annotation frequencies):
fact_only <- rse[mc$wt_fact=="fact_only"]
basal <- rse[mc$wt_fact %in% c("both", "wt_only")]
fmc <- mcols(fact_only)
bmc <- mcols(basal)
tab <- cbind(table(bmc$txType), table(fmc$txType))
colnames(tab) <- c("Basal_TSS", "FACT_TSS")
write.table(tab, "Table_S3_counts.txt", sep="\t", quote=F)
ptab <- prop.table(tab, 2)
ptab <- round(ptab, 4)
write.table(ptab, "Table_S3_fractions.txt", sep="\t", quote=F)
melt_ptab <- melt(ptab)
melt_ptab$Var1 <- factor(melt_ptab$Var1, levels=c("promoter", "fiveUTR", "proximal", "exon", "intron", "threeUTR", "antisense", "intergenic"))
names(melt_ptab) <- c("Annotation", "Category", "Value")
melt_ptab$Category <- factor(melt_ptab$Category, levels = c("FACT_TSS", "Basal_TSS"))
pdf("Fig_5A.pdf")
ggplot(melt_ptab, aes(x=Category, y=Value, fill=Annotation)) + geom_bar(stat='identity', width=0.75) + scale_fill_brewer(palette='Paired', direction=-1) + labs(y='Fraction of TCs', x='TSS category') + scale_y_reverse() + coord_flip() + theme(aspect.ratio=0.5)
dev.off()

# Figures 5B and 5C (violin plots):
lcpm <- cpm(assay(rse), log=T)
df <- cbind(as.data.frame(lcpm), txType=mcols(rse)$txType, wt_fact=mcols(rse)$wt_fact)
df$txType <- factor(df$txType, levels=c("promoter", "fiveUTR", "proximal", "exon", "intron", "threeUTR", "antisense", "intergenic"))
dfm <- melt(df, id=c("txType", "wt_fact"))
dfm_basal <- subset(dfm, wt_fact %in% c("both", "wt_only"))
pdf("Fig_5B.pdf")
ggplot(dfm_basal[dfm_basal$value>-4,], aes(x=txType, y=value, fill=txType)) + geom_violin() + scale_fill_brewer(palette='Paired', direction=-1) + geom_boxplot(width=0.1, outlier.colour=NA, fill="white") + labs(y='Tag counts (log2CPM)', x='TSS category') + theme(aspect.ratio=0.5) + ylim(-3, 13)
dev.off()
dfm_fact <- subset(dfm, wt_fact=="fact_only")
pdf("Fig_5C.pdf")
ggplot(dfm_fact[dfm_fact$value>-4,], aes(x=txType, y=value, fill=txType)) + geom_violin() + scale_fill_brewer(palette='Paired', direction=-1) + geom_boxplot(width=0.1, outlier.colour=NA, fill="white") + labs(y='Tag counts (log2CPM)', x='TSS category') + theme(aspect.ratio=0.5) + ylim(-3, 13)
dev.off()

# Analyze intersections with CAGE data from Tokizawa:
tk <- read.table("Tokizawa_S3.txt", sep="\t", header=T)
tk_1bp <- GRanges(seqnames=tk$chr, ranges=IRanges(tk$summit-1, end=tk$summit+1), strand=tk$strand)
mcols(tk_1bp) <- tk[,c(4:9)]
# Basal TSS:
mean(basal %over% tk_1bp) # 0.5402674
basal_over <- basal[basal %over% tk_1bp]
basal_out <- basal[basal %outside% tk_1bp]
basal_whole_over_1bp <- round(prop.table(cbind(table(mcols(basal_over)$txType), table(mcols(basal_out)$txType)), 1), 4)
colnames(basal_whole_over_1bp) <- c("basal_over", "basal_outside")
# FACT-specific TSS:	
mean(fact_only %over% tk_1bp) # 0.2250587
fact_over <- fact_only[fact_only %over% tk_1bp]
fact_out <- fact_only[fact_only %outside% tk_1bp]
fact_whole_over_1bp <- round(prop.table(cbind(table(mcols(fact_over)$txType), table(mcols(fact_out)$txType)), 1), 4)
colnames(fact_whole_over_1bp) <- c("fact_over", "fact_outside")
# Combine results:
whole_over_1bp <- cbind(basal_whole_over_1bp, fact_whole_over_1bp)[,c(1,3)]
write.table(whole_over_1bp, "Table_S2.txt", sep="\t", quote=F, row.names=F)
whole_melt <- melt(whole_over_1bp)
levels(whole_melt$Var1) <- c("promoter", "fiveUTR", "proximal", "exon", "intron", "threeUTR", "antisense", "intergenic")
levels(whole_melt$Var2) <- c("Basal", "FACT")
pdf("Fig_S5E.pdf")
ggplot(whole_melt, aes(x=Var1, y=value, fill=Var2)) + geom_bar(stat="identity", position="dodge"); dev.off()

# Figure S5B (the reproducibility between technical replicates);
wt <- lcpm[,c(1,2)]
spt16 <- lcpm[,c(3,4)]
ssrp1 <- lcpm[,c(5,6)]
# Skip pairs with zero values:
wt_nz <- wt[wt[,1]>-4 & wt[,2]>-4,]
spt16_nz <- spt16[spt16[,1]>-4 & spt16[,2]>-4,]
ssrp1_nz <- ssrp1[ssrp1[,1]>-4 & ssrp1[,2]>-4,]
cor1 <- round(cor(wt_nz[,1], wt_nz[,2]), 4) # 0.9556
cor2 <- round(cor(spt16_nz[,1], spt16_nz[,2]), 4) # 0.9687
cor3 <- round(cor(ssrp1_nz[,1], ssrp1_nz[,2]), 4) # 0.945
pdf("Fig_S5B_wt.pdf")
smoothScatter(wt_nz, nbin=256, xlab="wt_rep1 (log2CPM)", ylab="wt_rep2 (log2CPM)", main=paste0("r = ", cor1)); dev.off()
pdf("Fig_S5B_spt16.pdf")
smoothScatter(spt16_nz, nbin=256, xlab="wt_rep1 (log2CPM)", ylab="wt_rep2 (log2CPM)", main=paste0("r = ", cor2)); dev.off()
pdf("Fig_S5B_ssrp1.pdf")
smoothScatter(ssrp1_nz, nbin=256, xlab="wt_rep1 (log2CPM)", ylab="wt_rep2 (log2CPM)", main=paste0("r = ", cor3)); dev.off()

# Stats for Figure 5D:
coding_fact <- rse[mcols(rse)$txType %in% c("exon", "intron", "antisense") & mcols(rse)$wt_fact=="fact_only"] # n = 11555
mcols(coding_fact)$txType <- droplevels(mcols(coding_fact)$txType)
mcols(coding_fact)$fact <- droplevels(mcols(coding_fact)$fact)
tbl <- table(mcols(coding_fact)$txType, mcols(coding_fact)$fact)
write.table(tbl, "Fig_5D.txt", sep="\t", quote=F)

# Analyze the distribution of exonic FACT-specific TSS along the exons:
source("extractSummits.R")
fact_exon <- rse[mc$wt_fact=="fact_only" & mc$txType=="exon"]
summits <- extractSummits(fact_exon)		# n = 10045
all_exons <- reduce(exons(txdb))			# n = 153940
summits <- summits[summits %over% all_exons] # n = 9908
overlap_matrix <- findOverlaps(summits, all_exons, type="within")
exons <- all_exons[subjectHits(overlap_matrix),]
df <- cbind(summit = as.data.frame(summits)$start, as.data.frame(exons, row.names=NULL)[,c("start", "end", "width", "strand")])
rel_pos_fw <- (df$summit - df$start) / df$width
rel_pos_rev <- (df$end - df$summit) / df$width
df$rel_pos <- ifelse(df$strand == "+", rel_pos_fw, rel_pos_rev)
pdf("Fig_S5G.pdf")
hist(df$rel_pos, breaks=50, main="Relative Position of FACT-only TSS within exons"); dev.off()
