### Change the format of GRO-Seq data from Liu et al., 2018 (PMID 29379150) ###

# Download the following samples (as "plus" and "minus" BigWig files) from NCBI GEO:
# wt_rep1 (GSM2667813)
# wt_rep2 (GSM2667814)

# Convert BigWig to bedGraph (required kentUtils):
for file in *bw; do echo $file && bigWigToBedGraph $file ${file/.bw/.bg}; done

# Modify the chromosome names:
for file in *bg; do echo $file && sed 's/chr6/chrMt/;s/chr7/chrPt/;s/chr//' $file > ${file/.bg/.bedgraph}; done

# Combine "plus" and "minus" files:
for file1 in *fw.bedgraph; do file2=${file1/fw/rev} && outfile=${file1/fw.bedgraph/fw_rev.bg.gz} && echo $file1 "+" $file2 "=" $outfile && awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | cat $file1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $outfile; done