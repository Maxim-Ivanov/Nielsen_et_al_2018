extractSummits <- function(gr) {
  if ("thick" %in% names(mcols(gr))) {
    summits <- GRanges(seqnames = seqnames(gr), ranges = mcols(gr)$thick, strand = strand(gr), seqinfo = seqinfo(gr))
  } else {
    summits <- GRanges(seqnames = seqnames(gr), ranges = ranges(gr), strand = strand(gr), seqinfo = seqinfo(gr))
  }
  summits
}