#!/bin/bash

cd ~/Transcriptome_NGS

#first step --------------------------------------------------------------
#Do fastqc for all .fastq files
fastqc *.fastq

#Second step --------------------------------------------------------------
#to run the trimmomatic loop on all the .fastq files
for R1 in *.fastq ; do R2="${R1%_R1.fastq}_R2.fastq"; trimmomatic PE $R1 $R2 ${R1%.fastq}_trimmed.fq.gz ${R1%.fastq}_dr>

#Third step --------------------------------------------------------------
#Do fastqc for all .fastq files again to check if everyting is alright
fastqc *_trimmed.fq.gz

#Forth step --------------------------------------------------------------
#Do loop for pesudoalignment using kallisto against human index
for R1 in *R1_trimmed.fq.gz ;do R2="${R1%_R1_trimmed.fq.gz}_R2_trimmed.fq.gz" ; kallisto quant -i human2.idx -o "${R1%_>

#Fifth step  --------------------------------------------------------------
#Joining the matrices
output_file="merged_counts.tsv"
first_sample=$(ls SRR1000_folder/abundance.tsv | head -1)
cut -f 1 "$first_sample" > "$output_file"
for file in *_folder/abundance.tsv; do
sample=$(basename $(dirname "$file"))
cut -f 4 "$file" > tmp_counts.tsv
paste "$output_file" tmp_counts.tsv > tmp_merged.tsv
mv tmp_merged.tsv ./"$output_file";
done
rm tmp_counts.tsv
echo "✅ Merged counts saved in $output_file"
