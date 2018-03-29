# The function batchReadTrackData() iterates over a list of track
# filenames (such as Bedgraph or Bigwig) and calls the import()
# function from the rtracklayer library. A GRangesList object is
# returned;
# If the input track files lack the built-in strand info (encoded
# as positive and negative values of the signal), then they are
# expected to have the same strandness. In this case the strand
# info is passed by the "strand" argument to the function call;
# It is highly recommended to pass a valid seqinfo object to the
# function call;

batchReadTrackData <- function(filenames, dir=".", format=NULL, strand="*", seqinfo=NULL) {
  require(rtracklayer)
  output <- vector("list", length(filenames))
  for (i in seq_along(filenames)) {
    fn <- filenames[i]
    cat(fn, "\n"); flush.console()
    if (is.null(format)) {
      gr <- trim(import(file.path(dir, fn), seqinfo=seqinfo))
    } else {
      gr <- trim(import(file.path(dir, fn), format=format, seqinfo=seqinfo))
    }
    if (any(score(gr) < 0)) {
      strand(gr) <- ifelse(score(gr) >= 0, "+", "-")
      score(gr) <- abs(score(gr))
    } else {
      strand(gr) <- strand
    }
    output[[i]] <- gr
    names(output)[[i]] <- basename(fn)
  }
  output <- GRangesList(output)
  return(output)
}