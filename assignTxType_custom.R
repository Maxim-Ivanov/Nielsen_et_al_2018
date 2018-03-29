# The original assignTxType() function from CAGEfightR 0.2.0 was modified to remove the excessive categories (PROMPTS and CDS) in the annotation hierarchy. In particular, the CDS category was removed because it consumes the vast majority of exonic and intronic TSS;

assignTxType_custom <- function (gr, txdb, tssUpstream = 100, tssDownstream = 100, proximalUpstream = 1000, asFactor = FALSE) {
  message("Running a custom version of assignTxType (CDS and PROMPT annotations were removed)...")
  Promoters <- trim(promoters(txdb, upstream = tssUpstream, downstream = tssDownstream))
  Proximal <- trim(promoters(txdb, upstream = proximalUpstream, downstream = 0))
  FivePrimeUTRs <- fiveUTRsByTranscript(txdb)
  ThreePrimeUTRs <- threeUTRsByTranscript(txdb)
  Introns <- intronsByTranscript(txdb)
  Exons <- exons(txdb)
  Antis <- transcripts(txdb)
  strand(Antis) <- ifelse(strand(Antis) == "+", "-", "+")
  featureType <- rep("intergenic", length(gr))
  featureType <- ifelse(overlapsAny(query = gr, subject = Antis), "antisense", featureType)
  featureType <- ifelse(overlapsAny(query = gr, subject = Introns), "intron", featureType)
  featureType <- ifelse(overlapsAny(query = gr, subject = Exons), "exon", featureType)
  featureType <- ifelse(overlapsAny(query = gr, subject = ThreePrimeUTRs), "threeUTR", featureType)
  featureType <- ifelse(overlapsAny(query = gr, subject = FivePrimeUTRs), "fiveUTR", featureType)
  featureType <- ifelse(overlapsAny(query = gr, subject = Proximal), "proximal", featureType)
  featureType <- ifelse(overlapsAny(query = gr, subject = Promoters), "promoter", featureType)
  if (asFactor) {
    featureType <- factor(featureType, levels = c("intergenic", "proximal", "promoter", "fiveUTR", "intron", "exon", "threeUTR", "antisense"))
  }
  featureType
}