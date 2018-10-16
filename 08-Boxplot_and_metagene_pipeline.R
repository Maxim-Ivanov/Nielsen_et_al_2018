# Load the required libraries:
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(SummarizedExperiment)
library(rtracklayer)
library(ggplot2)
library(ggpubr)

# Define custom functions:
generateWindows <- function(gr, width) {
  w <- suppressWarnings(resize(gr, width = width, fix = "center"))
  w <- w[!findOutOfBounds(w)]
  return(w)
}

convertListToDataFrame <- function(x, y, width = NULL) {
  df <- as.data.frame(do.call(cbind, x))
  df$group <- y
  if (!is.null(width) & is.numeric(width)) {
    df$pos <- seq(-width/2, length.out = width)
  }
  return(df)
}

extractSummits <- function(gr) {
  stopifnot("thick" %in% names(mcols(gr)))
  swapped <- gr
  ranges(swapped) <- mcols(gr)$thick
  return(swapped)
}

findOutOfBounds <- function(gr) {
  lookup <- as.list(seqlengths(seqinfo(gr)))
  lengths <- lookup[as.character(seqnames(gr))]
  out <- (start(gr) < 0) | (end(gr) > lengths)
  return(out)
}

metageneMatrix <- function(signal, intervals, na.as.zeros = FALSE, skip.zeros = FALSE, skip.outliers = 0.995) {
  zero_cnt <- outlier_cnt <- 0
  # Input intervals (GRanges) are expected to have the same width:
  if (length(unique(width(intervals))) != 1) {
    tbl <- as.data.frame(table(width(intervals)))
    common_width <- as.integer(tbl[which.max(tbl$Freq), 1])
    message("metageneMatrix(): Intervals have different length! Skipping ", sum(width(intervals) != common_width), " intervals;")
    intervals <- intervals[width(intervals) == common_width]
  }
  mlen <- width(intervals)[[1]]
  # Input signal (GRanges) may have any strandness:
  interval_strands <- list("+", "-", "*")
  signal_strands <- list(c("+", "*"), c("-", "*"), c("+", "-", "*"))
  matlist <- vector("list", 3)
  for (i in 1:3) {
    curr_signal <- signal[strand(signal) %in% signal_strands[[i]]]
    curr_intervals <- intervals[strand(intervals) %in% interval_strands[[i]]]
    rlelist <- coverage(signal, weight="score")
    rlelist_split <- rlelist[intervals]
    numlist <- as(rlelist_split, "NumericList")
    numlist <- revElements(numlist, strand(intervals) == "-")
    matlist[[i]] <- do.call(rbind, numlist)
  }
  mat <- do.call(rbind, matlist)
  gene_cov <- rowSums(mat, na.rm = TRUE)
  if (isTRUE(skip.zeros)) {
    zeros <- gene_cov == 0
    zero_cnt <- sum(zeros)
    if (zero_cnt > 0) {
      mat <- mat[!zeros, ]
      gene_cov <- rowSums(mat, na.rm = TRUE)
    }
  }
  if (is.numeric(skip.outliers)) {
    q <- quantile(gene_cov, skip.outliers)
    outliers <- gene_cov > q
    mat <- mat[!outliers, ]
    outlier_cnt <- sum(outliers)
  }
  message("\t\tProcessed ", length(intervals), " intervals, skipped ", zero_cnt, " zeros and ", outlier_cnt, " outliers;"); flush.console()
  return(mat)
}

# Load the RangedSummarizedExperiment object containing the CAGEfightR output (see the 02-Calling_TSS_with_CAGEfightR.R pipeline):
rse <- readRDS("rse_FACT.RData")

# Extract coordinates of exonic TSS summits:
fact_exon <- extractSummits(rowRanges(rse[mcols(rse)$wt_fact=="fact_only" & mcols(rse)$txType=="exon"]))
basal_exon <- extractSummits(rowRanges(rse[mcols(rse)$txType == "exon" & mcols(rse)$wt_fact == "both"]))
# Extract promoter TSS summits (to be used as positive control):
prom <- extractSummits(rowRanges(rse[mcols(rse)$txType == "promoter"]))

# Generate random positions within exons of genes which contain either basal or FACT-specific exonic TSS:
## i) Get names of these genes:
fact_genes <- unique(mcols(fact_exon)$geneID)
basal_genes <- unique(mcols(basal_exon)$geneID)
## ii) Extract exonic intervals belonging to these genes:
ebg <- exonsBy(txdb, by = "gene") # all Arabidopsis exons grouped by gene
exons_fact <- reduce(unlist(ebg[names(ebg) %in% fact_genes]))
exons_basal <- reduce(unlist(ebg[names(ebg) %in% basal_genes]))
## Generate regularly spaced positions throughout the whole genome:
chrom_sizes <- seqlengths(seqinfo(txdb))
tiles <- unlist(tileGenome(chrom_sizes, tilewidth = 10))
pos <- resize(tiles, width = 1, fix = "start")
## Assign strand info to these positions:
pos_fw <- pos_rev <- pos
strand(pos_fw) <- "+"
strand(pos_rev) <- "-"
pos_all <- sort(c(pos_fw, pos_rev))
## Retain only positions within exonic intervals in genes of interest:
fact_pos <- pos_all[pos_all %over% exons_fact]
basal_pos <- pos_all[pos_all %over% exons_basal]
## Choose random 10K positions of each kind:
fact_ctrl <- fact_pos[sample(1:length(fact_pos), size = 10000)]
basal_ctrl <- basal_pos[sample(1:length(basal_pos), size = 10000)]

# Combine all intervals:
all_intervals <- list("Promoter TSS" = granges(prom), "Matched control for basal TSS" = basal_ctrl, "Matched control for FACT-only TSS" = fact_ctrl,
                      "Basal exonic TSS" = granges(basal_exon), "FACT-only exonic TSS" = granges(fact_exon))

# Define the colors:
my_colors <- c("#00BF7D", "#E76BF3", "#00B0F6", "#F8766D", "#A3A500")
names(my_colors) <- names(all_intervals)

# Define comparisons for boxplots:
my_comparisons <- list(c("Basal exonic TSS", "FACT-only exonic TSS"), c("Matched control for basal TSS", "Matched control for FACT-only TSS"), 
                       c("Basal exonic TSS", "Matched control for basal TSS"), c("FACT-only exonic TSS", "Matched control for FACT-only TSS"))

# Define the names of subdirectories with Bedgraph files corresponding to each study:
studies <- list("Luo2013", "Inagaki2017", "Liu2016", "Liu2018", "Chen2017", "Gomez2018", "Nasrallah2018", "Torres2018", "Yelagandula2014", 
                "Zhou2017", "VanDijk2010", "Wang2015", "Zhu2015", "Cortijo2017", "Zhang2015", "Dai2017")

# The main loop:
boxplot_window_width <- 20
metagene_window_width <- 400
ggw <- ggh <- 20; ggu <- "cm" # dimensions for ggplot

for (i in seq_along(studies)) {
  study <- studies[[i]]
  filenames <- list.files(path = study, pattern="bedgraph.gz$")
  for (j in seq_along(filenames)) {
    filename <- filenames[[j]]
    message(study, "\t", filename); flush.console()
    name <- sub(".bedgraph.gz", "", filename)
    data <- suppressWarnings(trim(import(file.path(study, filename), format = "bedGraph", seqinfo = seqinfo(txdb))))
    # In case if Bedgraph file contains stranded data (e.g. GRO-Seq):
    if (any(score(data) < 0)) {
      strand(data) <- ifelse(score(data) >= 0, "+", "-")
      score(data) <- abs(score(data))
    }
    message("\tData loaded;"); flush.console()
    # Calculations for the boxplot:
    bxp_windows <- GRangesList(lapply(all_intervals, generateWindows, width = boxplot_window_width))
    message("\tmetageneMatrix is working (this can be slow)..."); flush.console()
    bxp_matlist <- lapply(bxp_windows, metageneMatrix, signal = data)
    bxp_veclist <- lapply(bxp_matlist, function(mat) { list("avg_bp" = apply(mat, 1, mean)) })
    bxp_dflist <- mapply(bxp_veclist, names(bxp_veclist), FUN = convertListToDataFrame, SIMPLIFY = FALSE)
    bxp_df <- do.call(rbind, bxp_dflist)
    # Find the maximal value on boxplot (skipping outliers):
    bxp_stats <- aggregate(bxp_df$avg_bp, by = list(bxp_df$group), FUN = function(x) { boxplot.stats(x)$stats[5]})
    ymax <- max(bxp_stats[, 2])
    # Draw the boxplot:
    ylabel <- paste0("Signal averaged within ", boxplot_window_width, " bp windows")
    xlabel <- "Groups of intervals"
    bxp <- ggboxplot(bxp_df, x = "group", y = "avg_bp", fill = "group", palette = my_colors, outlier.colour = NA, xlab = xlabel, ylab = ylabel, title = name, notch = TRUE) +
      stat_compare_means(comparisons = my_comparisons, label.y = c(ymax * 1.1, ymax * 1.1, ymax * 1.2, ymax * 1.3), tip.length = 0) + 
      coord_cartesian(ylim = c(0, ymax * 1.3)) + theme(axis.text.x = element_blank(), legend.position = "right", legend.title = element_blank())
    ggsave(filename = file.path(study, paste0(name, "_boxplot.png")), plot = bxp, width = ggw, height = ggh, units = ggu)
    ggsave(filename = file.path(study, paste0(name, "_boxplot.pdf")), plot = bxp, width = ggw, height = ggh, units = ggu)
    message("\tBoxplot saved;"); flush.console()
    # Calculations for the metagene plot:
    mtg_windows <- GRangesList(lapply(all_intervals, generateWindows, width = metagene_window_width))
    message("\tmetageneMatrix is working (this can be slow)..."); flush.console()
    mtg_matlist <- lapply(mtg_windows, metageneMatrix, signal = data)
    mtg_veclist <- lapply(mtg_matlist, function(mat) { 
      avg <- apply(mat, 2, mean)
      sem <- apply(mat, 2, sd) / sqrt(nrow(mat))
      lower <- avg - 2 * sem
      upper <- avg + 2 * sem
      return(list("avg" = avg, "lower" = lower, "upper" = upper))
    })
    mtg_dflist <- mapply(mtg_veclist, names(mtg_veclist), FUN = convertListToDataFrame, SIMPLIFY = FALSE, width = metagene_window_width)
    mtg_df <- do.call(rbind, mtg_dflist)
    # Draw the metagene plot:
    mtg <- ggplot(mtg_df, aes(x = pos, y = avg, color = group)) + geom_line(size = 1) + 
      geom_ribbon(aes(x = pos, ymin = lower, ymax = upper, fill = group), alpha = 0.2) + ylim(0, NA) + 
      geom_vline(xintercept = 0) + ggtitle(name) + xlab("Relative position from peak summit") + 
      ylab("Average ChIP-Seq signal (with 95% CI)") + theme_bw() + scale_colour_manual(values = my_colors) + 
      scale_fill_manual(values = my_colors)
    ggsave(filename = file.path(study, paste0(name, "_metagene.png")), plot = mtg, width = ggw, height = ggh, units = ggu)
    ggsave(filename = file.path(study, paste0(name, "_metagene.pdf")), plot = mtg, width = ggw, height = ggh, units = ggu)
    message("\tMetagene saved;\n"); flush.console()
  }
}

