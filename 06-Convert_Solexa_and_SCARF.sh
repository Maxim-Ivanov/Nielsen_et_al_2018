### Convert the ChIP-Seq data from van Dijk et al., 2010 (PMID 21050490) from Solexa and SCARF formats to the regular FASTQ ###

# Download the raw reads from NCBI GEO:
# H3K4me1 GSM296187
# H3K4me2 GSM296188
# H3K4me3 GSM296189

# Extract reads from archives:
for file in *tar.gz; do tar zxvf $file; done

# Convert samples H3K4me1 and H3K4me3 (Solexa format).
# (requires solexa2fastq.pl script which was downloaded from 
# http://seqanswers.com/forums/archive/index.php/t-282.html):

for file in *seq.txt; do fc=$( echo $file | cut -f 1 -d "_") && echo $file $fc && ~/solexa2fastq.pl $file | sed "s/@ILunknown_unknown/@HWI-EAS258_${fc}/" > ${file/seq.txt/fastq}; done

cat *fastq | gzip > vanDijk2010_${hist_mod}.fastq.gz

# Convert sample H3K4me2 (SCARF format).
# (requires fq_all2std.pl script which is a part of the MAQ 
# distribution and can be downloaded from maq.sourceforge.net):

for file in *sequence.txt; do echo $file && perl ~/fq_all2std.pl scarf2std $file > ${file/txt/fastq}; done
