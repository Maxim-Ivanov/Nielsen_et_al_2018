# Load the required libraries:
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(SummarizedExperiment)
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

# Collect stats for Figure 5D:
coding_fact <- rse[mcols(rse)$txType %in% c("exon", "intron", "antisense") & mcols(rse)$wt_fact=="fact_only"] # n = 11555
mcols(coding_fact)$txType <- droplevels(mcols(coding_fact)$txType)
mcols(coding_fact)$fact <- droplevels(mcols(coding_fact)$fact)
tbl <- table(mcols(coding_fact)$txType, mcols(coding_fact)$fact)
write.table(tbl, "Fig_5D.txt", sep="\t", quote=F)

# Figure 5G (boxplot of relative promoter expression level of genes with fact-specific TSS):
# i) extract FACT-specific exon TSS:
fact_exon <- rse[mcols(rse)$txType=="exon" & mcols(rse)$wt_fact=="fact_only"] # n=10045
# ii) find genes with FACT-specific exon TSS:
fact_exon_genes <- unique(mcols(fact_exon)$geneID) # n=5604
# iii) extract promoter TSS:
prom_tss <- rse[mcols(rse)$txType=="promoter"] # n=30487
# iv) extract WT and SSRP1 promoter expression values (CPM-normalized):
wt_prom_expr <- rowSums(cpm(assay(prom_tss))[, c(1, 2)])
ssrp_prom_expr <- rowSums(cmp(assay(prom_tss))[, c(5, 6)])
# v) summarize the promoter expression values within each gene (avoid zero values by adding pseudocount = 1):
wt_prom_sum <- tapply(wt_prom_expr, INDEX=mcols(prom_tss)$geneID, FUN=sum) + 1
ssrp_prom_sum <- tapply(ssrp_prom_expr, INDEX=mcols(prom_tss)$geneID, FUN=sum) + 1
# vi) extract names of genes with quantified promoter expression:
all_genes <- unique(mcols(prom_tss)$geneID) # n=19561
# vii) combine the results:
df <- data.frame("wt_prom_expr"=wt_prom_sum, "ssrp_prom_expr"=ssrp_prom_sum, "ssrp_wt_diff"=ssrp_prom_sum/wt_prom_sum)
df$fact_exon_tss <- as.factor(ifelse(rownames(df) %in% fact_exon_genes, "Yes", "No"))
# viii) calculate p-value for the difference of medians:
wilcox.test(ssrp_wt_diff ~ fact_exon_tss, data=df) # p-value < 2.2e-16
median(df$ssrp_wt_diff[df$fact_exon_tss=="Yes"]) # 1.077855 (slightly increased in SSRP compared to WT)
median(df$ssrp_wt_diff[df$fact_exon_tss=="No"]) # 1.028625
# ix) draw the boxplot:
p <- ggplot(df, aes(x=fact_exon_tss, y=ssrp_wt_diff)) + geom_boxplot(outlier.colour=NA) + ylim(NA, 2) + geom_hline(yintercept=1, colour="red")
ggsave("Boxplot.pdf", plot=p)

# Figure S5B (reproducibility between technical replicates of TSS-Seq):
counts <- as.data.frame(cpm(assays(rse)$counts)) # CPM-normalize tag counts
ps = 0.1 # add pseudocount to avoid taking log2 of zero
df_wt <- subset(counts, s1_wt != 0 | s2_wt != 0, c(s1_wt, s2_wt)) + ps # skip TSS clusters which were not detected in both replicates
df_spt16 <- subset(counts, s3_spt16 != 0 | s4_spt16 != 0, c(s3_spt16, s4_spt16)) + ps
df_ssrp1 <- subset(counts, s5_ssrp1 != 0 | s6_ssrp1 != 0, c(s5_ssrp1, s6_ssrp1)) + ps
cor_wt <- cor(df_wt[, 1], df_wt[, 2])
cor_spt16 <- cor(df_spt16[, 1], df_spt16[, 2])
cor_ssrp1 <- cor(df_ssrp1[, 1], df_ssrp1[, 2])
p1 <- ggplot(df_wt, aes(x = s1_wt, y = s2_wt)) + geom_point(alpha = 0.1) + 
  scale_x_continuous(trans = 'log2', limits = c(1, NA)) + scale_y_continuous(trans = 'log2', limits = c(1, NA)) + 
  ggtitle(paste0("Wild type (Pearson r = ", round(cor_wt, 4), ")")) + labs(x = "Rep1 (log2 CPM)", y = "Rep2 (log2 CPM)")
ggsave("TC_counts_rep1_vs_rep2_WT_log2.pdf", plot = p1)
p2 <- ggplot(df_spt16, aes(x = s3_spt16, y = s4_spt16)) + geom_point(alpha = 0.1) + 
  scale_x_continuous(trans = 'log2', limits = c(1, NA)) + scale_y_continuous(trans = 'log2', limits = c(1, NA)) + 
  ggtitle(paste0("spt16-1 (Pearson r = ", round(cor_spt16, 4), ")")) + labs(x = "Rep1 (log2 CPM)", y = "Rep2 (log2 CPM)")
ggsave("TC_counts_rep1_vs_rep2_SPT16_log2.pdf", plot = p2)
p3 <- ggplot(df_ssrp1, aes(x = s5_ssrp1, y = s6_ssrp1)) + geom_point(alpha = 0.1) + 
  scale_x_continuous(trans = 'log2', limits = c(1, NA)) + scale_y_continuous(trans = 'log2', limits = c(1, NA)) + 
  ggtitle(paste0("ssrp1-2 (Pearson r = ", round(cor_ssrp1, 4), ")")) + labs(x = "Rep1 (log2 CPM)", y = "Rep2 (log2 CPM)")
ggsave("TC_counts_rep1_vs_rep2_SSRP1_log2.pdf", plot = p3)

# Figure S5C (intersections with CAGE data from Tokizawa 2017 - PMID 28214361):
tk <- read.table("Tokizawa_S3.txt", sep="\t", header=T)
tk_1bp <- GRanges(seqnames=tk$chr, ranges=IRanges(tk$summit-1, end=tk$summit+1), strand=tk$strand)
mcols(tk_1bp) <- tk[,c(4:9)]
# Basal TSS:
mean(basal %over% tk_1bp) # 0.5402674
basal_over <- basal[basal %over% tk_1bp]
basal_out <- basal[basal %outside% tk_1bp]
basal_whole_over_1bp <- round(prop.table(cbind(table(mcols(basal_over)$txType), table(mcols(basal_out)$txType)), 1), 4)
colnames(basal_whole_over_1bp) <- c("basal_over", "basal_outside")
# fact-specific TSS:	
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

# Figure S5E (distribution of exonic FACT-specific TSS along the exons):

extractSummits <- function(gr) {
  stopifnot("thick" %in% names(mcols(gr)))
  swapped <- gr
  ranges(swapped) <- mcols(gr)$thick
  return(swapped)
}

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
pdf("Fig_S5E.pdf")
hist(df$rel_pos, breaks=50, main="Relative Position of FACT-only TSS within exons"); dev.off()
