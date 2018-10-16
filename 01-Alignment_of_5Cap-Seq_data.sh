# Quality and adapter trimming (observe the custom adapter 
# sequence):
for file in *fastq.gz; do echo $file && trim_galore --adapter "ATCTCGTATGCCG" $file; done

# Trim UMIs (first 8 nt of each read), append them to Fastq 
# headers (UMI-Tools required):
for file in *trimmed.fq.gz; do echo $file && umi_tools extract --stdin=${file} --bc-pattern=NNNNNNNN --stdout=${file/.fq.gz/_UMI.fq.gz}; done

# Align to TAIR10 using STAR:
for file in ./*UMI.fq.gz; do echo $file && STAR --genomeDir /index/tair10/star --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType EndToEnd --readFilesCommand zcat; done

# Sort SAM files and convert to BAM:
for file in *sam; do echo $file && samtools view -hu $file | samtools sort - -o ${file/.sam/_sorted.bam} && rm $file; done

# Filter out rRNA, tRNA and sn/snoRNA reads:
for file in *sorted.bam; do echo $file && bedtools intersect -v -abam $file -b Araport11_rRNA_tRNA_snRNA_snoRNA.bed > ${file/.bam/_filt.bam}; done

# Filter out multimapper reads:
for file in *filt.bam; do echo $file && samtools view -h -q 10 $file -o ${file/.bam/_mapq10.bam}; done

# Deduplicate on UMIs:
for file in *mapq10.bam; do echo $file && umi_tools dedup -I $file -S ${file/.bam/_dedup.bam}; done

# Make stranded Bedgraph files (only the first base of each read 
# is considered):
for str in "+" "-"; do
  echo $str
  [ "$str" = "+" ] && n="fw" || n="rev"
  for file in *dedup_sorted.bam; do
    echo $file && bedtools genomecov -ibam $file -bg -5 -strand $str | sort -k 1,1 -k 2,2n > ${file/.bam/}_${n}.bg
  done
done

# Remove singletons (positions covered by a single tag):
for file in *bg.gz; do echo $file && zcat $file | awk '{if ($4>1) print $0}' > ${file/.bg.gz/_cov2.bg}; done

# Ensure that all positions have 1 bp width (this is required for 
# further processing by CAGEfightR):
python3 Expand_bedGraph_to_single_base_resolution.py .

# Convert expanded Bedgraph files to Bigwig (requires kentUtils).
# These Bigwigs will be used as input for CAGEfightR:
for file in *cov2.bg; do echo $file && bedGraphToBigWig $file TAIR10.chrom.sizes ${file/bg/bw}; done

# Also make Bedgraph files for visualization in genomic browsers. # First, merge forward and reverse Bedgraph files corresponding 
# to the same sample (the strand info is encoded by the sign 
# of values in column 4):
for file1 in *fw.bg; do
  file2=${file1/fw/rev}
  outfile=${file1/fw.bg/fw_rev.bg.gz}
  echo $file1 "+" $file2 "=" $outfile
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | \
cat $file1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph \
color=0,100,200 altColor=200,100,0' | gzip > $outfile
done

# Then normalize the merged Bedgraph files to 1M tags:
for file in *fw_rev.bg.gz; do
  norm=$( zcat $file | sed '/^[#t]/d' | awk 'BEGIN{SUM=0}\
{SUM+=sqrt($4^2)*($3-$2)}END{print SUM / 1000000}' )
  echo $file $norm
  zcat $file | awk -v norm=$norm 'BEGIN{OFS="\t"}\
{if ($0~/^[#t]/) print $0; else print $1, $2, $3, $4 / norm}' | \
gzip > ${file/.bg.gz/_norm1M.bg.gz}
done