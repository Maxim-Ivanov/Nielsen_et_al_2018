# Given a set of genomic intervals, this function generates positions which are randomly distributed along covered portion of the genome;

randomPositions <- function(gr, n, strand.at.random=FALSE) {
  merged_gr <- reduce(gr)
  widths <- width(merged_gr)
  total_coverage <- sum(widths)
  starts <- cumsum(widths) - widths
  random_nums <- round(runif(n, 1, total_coverage))
  random_positions <- lapply(random_nums, function(x) {
    distances <- x - starts
    relative_position <- min(distances[distances >= 0])
    index <- distances == relative_position
    interval_start = start(merged_gr)[index]
    abs_position <- interval_start + relative_position
    seqname <- as.character(seqnames(merged_gr)[index])
    if (isTRUE(strand.at.random)) {
      strand <- sample(c("+", "-"), size=1)
    } else {
      strand <- as.character(strand(merged_gr)[index])
    }
    return(c(seqname, abs_position, strand))
  })
  mat <- do.call(rbind, random_positions)
  seqnames <- as.character(mat[,1])
  coords <- as.numeric(mat[,2])
  strands <- as.character(mat[,3])
  out_gr <- GRanges(seqnames=seqnames, ranges=IRanges(coords, end=coords), strand=strands, seqinfo=seqinfo(gr))
  out_gr
}