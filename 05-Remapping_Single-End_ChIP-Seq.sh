### Pipeline for Single-End ChIP-Seq data ###

# Trim adapters and align to TAIR10:
for file in *fastq.gz; do echo $file && STAR --genomeDir /index/tair10_star --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} --outSAMmultNmax 1 --alignEndsType Local --readFilesCommand zcat -- clip3pAdapterSeq "AGATCGGAAGAGC" --alignIntronMax 1; done

# Clean up in the working directory:
rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Convert to coordinate-sorted BAM (remove MAPQ values below 10):
for file in *sam; do echo $file && samtools view -hq 10 $file | samtools sort - -o ${file/.sam/.bam} && rm $file; done

# Extend reads to d/2 using MACS2:
for file in *bam; do echo $file && macs2 -t $file -n ${file/.bam/} -g 1.35e+08 -m 3,50 --half-ext --bdg > ${file/.bam/_log.txt} 2>&1; done

# Clean up:
rm *model.r *pvalue.bdg *qvalue.bdg *control_lambda.bdg *bed *xls *encodePeak *txt
for file in *pileup.bdg; do mv $file ${file/_treat_pileup.bdg/.bedgraph}; done

# Gzip Bedgraph files:
for file in *bedgraph; do echo $file && gzip $file; done
