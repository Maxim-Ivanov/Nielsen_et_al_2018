### Pipeline for Paired-End ChIP-Seq data ###

# Trim adapters and align to TAIR10:
for f1 in *_1.fastq.gz; do f2=${f1/_1.f/_2.f} && echo $f1 $f2 && STAR --genomeDir /index/tair10_star --readFilesIn $f1 $f2 --runThreadN 4 --outFileNamePrefix ${f1/1.fastq.gz/} --outSAMmultNmax 1 --alignEndsType Local --readFilesCommand zcat -- clip3pAdapterSeq "AGATCGGAAGAGC" --alignIntronMax 1; done

# Clean up in the working directory:
rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Convert to coordinate-sorted BAM (remove MAPQ values below 10 and singletons):
for file in *sam; do echo $file && samtools view -hq 10 -f 2 $file | samtools sort - -o ${file/.sam/.bam} && rm $file; done

# Filter forward reads in proper pairs, then add the positive TLEN value to POS (= insert start) to get the end position of the insert
for file in *bam; do echo $file && samtools view -f 2 -F 16 $file | awk '{print $3"\t"$4"\t"$4+$9}' > ${file/.bam/.bed}; done

# Only for MNase-Seq data (Dai et al., 2017):
# (additional filtering for insert size between 140 and 155 bp)
for file in Dai2017*bed; do awk '{W=$3-$2; if (W >= 140 && W <= 155) print}' $file > temp.bed && rm $file && mv temp.bed $file; done

# Convert BED files to Bedgraph:
for file in *insert.bed; do echo $file && bedtools genomecov -i $file -g $tair10ChrNameLength.txt -bg | sort -k1,1 -k2,2n | gzip > ${file/_insert.bed/.bedgraph.gz}; done


