# Load the required libraries
library(CAGEfightR)						# version 0.99.0
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BiocParallel)
register(MulticoreParam(4), default=T) # SnowParam() on Windows
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)
source("assignTxType_custom.R")

# Load the TSS-Seq BigWig files (see the 01-Alignment_of_5Cap-Seq_data.sh pipeline):
setwd("/path/to/bigwig/files")
bw_plus_filenames <- list.files(".", pattern="cov2_fw.bw$")
bw_minus_filenames <- list.files(".", pattern="cov2_rev.bw$")
bw_plus <- BigWigFileList(bw_plus_filenames)
bw_minus <- BigWigFileList(bw_minus_filenames)
sample_names <- sub('_cov2_fw.bw', '', bw_plus_filenames)
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names

# Make the design matrix:
design_matrix <- data.frame("Name"=sample_names, "BigWigPlus"=bw_plus_filenames, "BigWigMinus"=bw_minus_filenames, 
row.names=sample_names, genotype=factor(rep(c("wt", "spt16", "ssrp1"), each=2), levels=c("wt", "spt16", "ssrp1")))

# Quantify all tag clusters (TCs):
ctss <- quantifyCTSSs(plusStrand=bw_plus, minusStrand=bw_minus, design=design_matrix, genome=seqinfo(Athaliana))

# Call candidate TSS:
tss <- quickTSSs(ctss)

# Call candidate enhancers:
enhancers <- quickEnhancers(ctss)

# Annotate TSS and enhancers by genomic features (observe that a custom assignTxType() function is used):
rowRanges(tss)$txType <- suppressWarnings(assignTxType_custom(rowRanges(tss), txdb=txdb, asFactor=TRUE))
rowRanges(enhancers)$txType <- suppressWarnings(assignTxType_custom(rowRanges(enhancers), txdb=txdb, asFactor=TRUE))

# Remove enhancers overlapping known transcripts:
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))

# Combine candidate TSS and enhancers into a single RangedSummarizedExperiment object:
rowRanges(tss)$clusterType <- "TSS"
rowRanges(enhancers)$clusterType <- "enhancer"
rse <- combineClusters(tss, enhancers, removeIfOverlapping="object1")

# Remove low expressed TCs:
rse <- subsetBySupport(rse, inputAssay = "counts", outputColumn = "support", unexpressed = 0, minSamples = 1) # n = 96232

# Annotate TCs by gene IDs:
rse <- suppressWarnings(assignGeneID(rse, geneModels=txdb))

# Annotate TCs by gene names:
tair_ann <- import.gff3("Arabidopsis_thaliana.TAIR10.26.gff3")
mapping <- data.frame("geneID" = sub("gene:", "", tair_ann$ID), "name"=tair_ann$external_name)
tmp <- left_join(x=tibble(geneID=rowRanges(rse)$geneID), y=mapping, by="geneID", na_matches = "never")
mcols(rse) <- DataFrame(mcols(rse), tmp[,-1])
rm(tmp)

# For DE calling consider only strong peaks (TPM >= 1 in at least 2 samples):
rse2 <- rse[rowSums(cpm(assay(rse)) >= 1) >= 2]	# n = 25964

# Use DESeq2:
dds <- DESeqDataSet(rse2, design = ~ genotype)
dds <- DESeq(dds)

# Extract DE results:
spt16 <- results(dds, contrast = c("genotype", "spt16", "wt"))
summary(spt16)		# 664 up, 870 down
ssrp1 <- results(dds, contrast = c("genotype", "ssrp1", "wt"))
summary(ssrp1)		# 1183 up, 1777 down

# Combine DE results:
df1 <- as.data.frame(spt16)[,c(2,6)]
colnames(df1) <- c("log2FC_spt16", "padj_spt16")
df2 <- as.data.frame(ssrp1)[,c(2,6)]
colnames(df2) <- c("log2FC_ssrp1", "padj_ssrp1")
m <- cbind(df1, df2)
m[is.na(m)] <- 1			# replace padj=NA with padj=1

# Expand DE results to the original TC number:
m$key <- rownames(m)
orig <- data.frame("key"=names(rowRanges(rse)))
n <- left_join(orig, m, by=c("key"), all.x=TRUE)

# Add DE results to the RSE object:
mcols(rse) <- cbind(mcols(rse), n[-1])

##### Some extra calculations (not required in the general case)

# Add decisions on the sample specificity of called TSS:
mcols(rse)$nz_wt <- rowSums(assay(rse)[,c(1,2)]>0)
mcols(rse)$nz_spt16 <- rowSums(assay(rse)[,c(3,4)]>0)
mcols(rse)$nz_ssrp1 <- rowSums(assay(rse)[,c(5,6)]>0)
mc <- mcols(rse)
wt_only <- (mc$nz_spt16 + mc$nz_ssrp1)==0
fact_only <- mc$nz_wt==0
mcols(rse)$wt_fact <- "both"
mcols(rse)$wt_fact[wt_only] <- "wt_only"
mcols(rse)$wt_fact[fact_only] <- "fact_only"

no <- (mc$nz_spt16 + mc$nz_ssrp1) <= 1		
spt16_only <- mc$nz_spt16==2 & mc$nz_ssrp1==0
ssrp1_only <- mc$nz_spt16==0 & mc$nz_ssrp1==2
mcols(rse)$fact <- "both"
mcols(rse)$fact[no] <- "no"
mcols(rse)$fact[spt16_only] <- "spt16_only"
mcols(rse)$fact[ssrp1_only] <- "ssrp1_only"

table(mcols(rse)$wt_fact, mcols(rse)$fact)

# Add fold changes and also the explicit decisions on differential expression:
cm <- cpm(assay(rse))
mcols(rse)$spt16_FC <- rowMeans(cm[,c(3,4)]) / rowMeans(cm[,c(1,2)])
mcols(rse)$ssrp1_FC <- rowMeans(cm[,c(5,6)]) / rowMeans(cm[,c(1,2)])
mcols(rse)$spt16_dir <- ifelse(mcols(rse)$spt16_FC > 1, "up", "down")
mcols(rse)$ssrp1_dir <- ifelse(mcols(rse)$ssrp1_FC > 1, "up", "down")

mcols(rse)$spt16_de <- "nonDE"
spt16_up <- !is.na(mc$padj_spt16) & mc$padj_spt16<=0.1 & mc$log2FC_spt16>0
mcols(rse)$spt16_de[spt16_up] <- "up"
spt16_down <- !is.na(mc$padj_spt16) & mc$padj_spt16<=0.1 & mc$log2FC_spt16<0
mcols(rse)$spt16_de[spt16_down] <- "down"

mcols(rse)$ssrp1_de <- "nonDE"
ssrp1_up <- !is.na(mc$padj_ssrp1) & mc$padj_ssrp1<=0.1 & mc$log2FC_ssrp1>0
mcols(rse)$ssrp1_de[ssrp1_up] <- "up"
ssrp1_down <- !is.na(mc$padj_ssrp1) & mc$padj_ssrp1<=0.1 & mc$log2FC_ssrp1<0
mcols(rse)$ssrp1_de[ssrp1_down] <- "down"

# Convert vectors to factors:
mcols(rse)$name <- as.character(mcols(rse)$name)
mcols(rse)$wt_fact <- as.factor(mcols(rse)$wt_fact)
mcols(rse)$fact <- as.factor(mcols(rse)$fact)
mcols(rse)$spt16_dir <- as.factor(mcols(rse)$spt16_dir)
mcols(rse)$ssrp1_dir <- as.factor(mcols(rse)$ssrp1_dir)
mcols(rse)$spt16_de <- as.factor(mcols(rse)$spt16_de)
mcols(rse)$ssrp1_de <- as.factor(mcols(rse)$ssrp1_de)

# Save the final RangedSummarizedExperiment object:
saveRDS(rse, "rse_FACT.RData")
