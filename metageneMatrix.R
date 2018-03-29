# The intervals are expected to be stranded;
# The signal can be either stranded (e.g. NET-Seq) or unstranded (e.g. ChIP-Seq);
# In the antisense mode (antisenseMode=TRUE) the intervals are flipped to the opposite strand;
#
# If at least one interval has a different length, then intervals can be either: a) scaled (extended or shrinked) to a certain common length (scale = TRUE), or b) remain unscaled (scale = FALSE).
# In the latter case (scale = FALSE), longer intervals are trimmed and shorter intervals are filled with NA values. The intervals can be anchored either at the start (default), or at the end (anchor = c("start", "end"));
# In both cases it is required to determine the length of the output matrix. The matrix.length argument can be either an arbitrary number, or "max" (width of the longest interval), or "min" (width of the shortest interval), otherwise the median interval length is used;
#
# na.as.zeros parameter allows to substitute NA values with zeros. By default it is FALSE, thus the output matrix may contain NA values (given that scale = FALSE). This has to be taken into account when applying statistical functions to the output matrix;
#
# skip.zeros parameter allows to remove intervals with zero signal (not a single tag along the whole interval);
#
# skip.outliers parameter allows to remove intervals with the highest signal (by default, the top 0.5% intervals are considered as potential outliers; the threshold can be changed arbitrarily). Use skip.outliers=FALSE to suppress the removal of outliers;
#
# equal.weights allows to normalize the integral signal over each interval to 1 (i.e. the signal magnitude is not taken into account when computing the matrix). This may be useful when the signal is expected to consist of sharp discrete peaks (e.g. 5Cap-Seq TSS), and their relative positions within intervals are of interest (whereas their heights are considered not important);

calcMatrixLength <- function(inter, m.length) {
  max_width <- max(width(inter))
  min_width <- min(width(inter))
  if (max_width == min_width) {
    out <- max_width
  } else {
    if (is.numeric(m.length)) {
      out <- m.length
    } else if (m.length == "max") {
      out <- max_width
    } else if (m.length == "min") {
      out <- min_width
    } else if (m.length == "median") {
      out <- round(median(width(inter)))
    } else {
      cat("Check matrix.length agrument!\n")
      out <- max_width
    }
  }
  return(out)
}

expandOrShrink <- function(x, mlen) {
  if (length(x) == mlen) {
    out <- x
  } else {
    expanded <- rep(x, each = mlen)
    breaks <- rep(seq(1, mlen), each = length(x))
    averaged <- tapply(expanded, breaks, mean)
    out <- as.numeric(averaged)
  }
  return(out)
}

trimOrFill <- function(x, mlen, anchor, na.as.zeros) {
  if (length(x) == mlen) {
    out <- x
  } else if (length(x) > mlen) {
    if (anchor == "start") {
      out <- x[1:mlen]
    } else {
      out <- x[(length(x)-mlen+1):length(x)]
    }
  } else {
    if (anchor == "start") {
      out <- c(x, rep(NA, (mlen-length(x))))
    } else {
      out <- c(rep(NA, (mlen-length(x))), x)
    }
    if (isTRUE(na.as.zeros)) {
      out[is.na(out)] <- 0
    }
  }
  return(out)
}

calcMatrix <- function(signal, intervals, strand=NA, max_w, min_w, scale, mlen, anchor, na.as.zeros) {
  if (!is.na(strand)) {
    intervals <- intervals[strand(intervals)==strand]
    signal <- signal[strand(signal)==strand]
  }
  rlelist <- coverage(signal, weight="score")
  rlelist_split <- rlelist[intervals]
  numlist <- as(rlelist_split, "NumericList")
  numlist <- revElements(numlist, strand(intervals) == "-")
  if (max_w != min_w) {
    if (isTRUE(scale)) {
      cat("Expanding and shrinking...\n")
      numlist <- lapply(numlist, expandOrShrink, mlen=mlen)
    } else {
      cat("Trimming and filling...\n")
      numlist <- lapply(numlist, trimOrFill, mlen=mlen, anchor=anchor, na.as.zeros=na.as.zeros)
    }
  }
  mat <- do.call(rbind, numlist)
  return(mat)
}

metageneMatrix <- function(signal, intervals, scale = TRUE, matrix.length = "max", anchor = "start", na.as.zeros = FALSE, skip.zeros = TRUE, skip.outliers = 0.995, equal.weights = FALSE, antisenseMode = FALSE) {
  require(GenomicRanges)
  mlen <- calcMatrixLength(inter = intervals, m.length = matrix.length)
  max_w <- max(width(intervals)); min_w <- min(width(intervals))
  if (any(strand(intervals)=="*")) {
    message("Intervals contain unstranded records (will be considered as forward);")
    strand(intervals) <- ifelse(strand(intervals)=="*", "+", strand(intervals))
  }
  if (isTRUE(antisenseMode)) {
    strand(intervals) <- ifelse(strand(intervals)=="+", "-", "+")
  }
  if (any(strand(signal)=="*")) {
    message("Signal contains unstranded records;")
    mat <- calcMatrix(signal, intervals, max_w=max_w, min_w=min_w, scale=scale, mlen=mlen, anchor=anchor, na.as.zeros=na.as.zeros)
  } else {
    mat_fw <- calcMatrix(signal, intervals, strand="+", max_w=max_w, min_w=min_w, scale=scale, mlen=mlen, anchor=anchor, na.as.zeros=na.as.zeros)
    mat_rev <- calcMatrix(signal, intervals, strand="-", max_w=max_w, min_w=min_w, scale=scale, mlen=mlen, anchor=anchor, na.as.zeros=na.as.zeros)
    mat <- rbind(mat_fw, mat_rev)
  }
  if (isTRUE(antisenseMode)) {
    mat <- mat[, ncol(mat):1, drop=FALSE]
  }
  gene_cov <- rowSums(mat, na.rm=TRUE)
  if (isTRUE(skip.zeros)) {
    zeros <- gene_cov == 0
    if (sum(zeros) > 0) {
      mat <- mat[!zeros, ]
      gene_cov <- rowSums(mat, na.rm=TRUE)
      cat("Skipped", sum(zeros), "intervals with zero signal;\n")
    }
  }
  if (is.numeric(skip.outliers)) {
    q <- quantile(gene_cov, skip.outliers)
    outliers <- gene_cov > q
    mat <- mat[!outliers, ]
    cat("Skipped", sum(outliers), "potential outliers;\n")
  }
  if (isTRUE(equal.weights)) {
    mat <- t(apply(mat, 1, function(x) { x / sum(x, na.rm=TRUE) }))
    cat("All intervals were assigned equal weights;")
  }
  return(mat)
}
