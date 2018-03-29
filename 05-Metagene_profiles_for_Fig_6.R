# Load the required libraries:
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(SummarizedExperiment)
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
library(rtracklayer)
library(ggplot2)

# Load the required custom code:
source("batchReadTrackData.R")
source("metageneMatrix.R")
source("randomPositions.R")
source("drawMetagenePlot.R")
source("metagenePipeline.R")

# Load the RangedSummarizedExperiment object with CAGEfightR output (see the 02-Calling_TSS_with_CAGEfightR.R pipeline):
rse <- readRDS("rse_FACT.RData")

# Extract coordinates of exonic TSS summits:
fact_exon <- extractSummits(rse[mc$wt_fact=="fact_only" & mc$txType=="exon"])
basal_exon <- extractSummits(rse[mc$txType == "exon" & mc$wt_fact == "both"])
# Extract promoter TSS summits (to be used as positive control):
prom <- extractSummits(rse[mc$txType == "promoter"])
# Generate random exonic positions (negative control):
all_exons <- reduce(exons(txdb))
random_exon <- randomPositions(all_exons, 10000)
# Generate random genomic positions (negative control):
chroms <- as.data.frame(seqinfo(Athaliana))
chr_gr <- GRanges(seqnames=rownames(chroms), ranges=IRanges(c(1), end=chroms$seqlengths), strand=c("*"))
random_genome <- randomPositions(chr_gr, 10000, strand.at.random=TRUE)

# Combine all these coordinates into a single GRangesList object:
grl <- list("FACT-only_exonic_TSS" = fact_exon, "Basal_exonic_TSS" = basal_exon, "Random_exonic_positions" = random_exon, "Random_genomic_positions" = random_genome, "Promoter_TSS" = prom)

# Load ChIP-Seq data from Luo et al., 2012 (PMID 22962860):
luo_dir <- "/path/to/bedgraph/files"
luo_filenames <- list.files(luo_dir, pattern="*bg.gz$")
luo_data <- suppressWarnings(batchReadTrackData(luo_filenames, luo_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))

# Load ChIP-Seq data from Inagaki et al., 2017 (PMID 28100676):
inagaki_dir <- "/path/to/bedgraph/files"
inagaki_filenames <- list.files(inagaki_dir, pattern="*bdg.gz$")
inagaki_data <- suppressWarnings(batchReadTrackData(inagaki_filenames, inagaki_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))

# Load ChIP-seq data from Liu et al., 2016 (PMID 27225844):
liu_weigel_dir <- "/path/to/bedgraph/files"
liu_weigel_filenames <- list.files(liu_weigel_dir, pattern="*bdg.gz$")
liu_weigel_data <- suppressWarnings(batchReadTrackData(liu_weigel_filenames, liu_weigel_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))

# Load GRO-Seq data from Liu et al., 2018 (PMID 29379150):
liu_jacobsen_dir <- "/path/to/bedgraph/files"
liu_jacobsen_filenames <- list.files(liu_jacobsen_dir, pattern="fw_rev.bedgraph.gz$")
liu_jacobsen_data <- suppressWarnings(batchReadTrackData(liu_jacobsen_filenames, liu_jacobsen_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))

# Finally draw the metagene plots (400 bp windows were centered around the TSS summits and random positions):
r1 <- suppressWarnings(metagenePipeline(signal=luo_data, intervals=grl, skip="_treat_pileup.bg.gz", out_dir="./Luo_2012_ChIP-Seq", expand=400))
r2 <- suppressWarnings(metagenePipeline(signal=inagaki_data, intervals=grl, skip="_treat_pileup.bdg.gz", out_dir="./Inagaki_2017_ChIP-Seq", expand=400))
r3 <- suppressWarnings(metagenePipeline(signal=liu_weigel_data, intervals=grl, skip="_treat_pileup.bdg.gz", out_dir="./Liu_Weigel_2016_ChIP-Seq", expand=400))
r4 <- suppressWarnings(metagenePipeline(signal=liu_jacobsen_data, intervals=grl, skip="_fw_rev.bedgraph.gz", out_dir="./Liu_Jacobsen_2018_GRO-Seq", expand=400))
