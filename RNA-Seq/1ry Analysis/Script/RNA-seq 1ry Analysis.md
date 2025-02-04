# RNA-Seq Primary Analysis Pipeline
***
![RNA-seq Primary Analysis Pipeline pipeline template](https://github.com/user-attachments/assets/5172d0b1-9763-4598-bfe0-7865d79c9bf2)

## ***Zero Step*** (Samples Download) --------------------------------------
#### ***[Note:]*** 
You can download samples using SRAtoolkit from the command line, but in my case I downloaded the samples using SRA run selector tool on NCBI (SRA) database
```{bash}
mkdir ~/Transcriptome_NGS 
cd ~/Transcriptome_NGS
```
## ***First Step*** (Quality Check I) --------------------------------------
### Check the quality of all the samples (using FastQC)
#### ***[Note:]*** 
_All my samples were downloaded in the directory created at the very beginning_
```{bash}
fastqc *.fastq
```

## ***Second Step*** (QC) --------------------------------------------------
### Run trimmomatic tool to loop on all the .fastq files
```{bash}
for R1 in *.fastq ; do R2="${R1%_R1.fastq}_R2.fastq"; trimmomatic PE $R1 $R2 ${R1%.fastq}_trimmed.fq.gz ${R1%.fastq}_drop.fq.gz ${R2%.fastq}_trimmed.fq.gz ${R2%.fastq}_drop.fq.gz ILLUMINACLIP:./adaptor.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100; done>
```

## ***Third Step*** (Quality Check II) -------------------------------------
### Do fastqc for all .fastq files again to check if everything is alright
```{bash}
fastqc *_trimmed.fq.gz
```
#### ***[Note:]*** 
If the quality is acceptable then proceed to next step (i.e. no adapter sequences or low quality base pairs),
If not then you need to adjust the trimmomatic parameters according to your sample and do QC again 

## ***Forth Step*** (Alignment & Quantification)----------------------------
### Pesudoalignment (splice-aware alignment) using Kallisto against the human index
```{bash}
for R1 in *R1_trimmed.fq.gz ;do R2="${R1%_R1_trimmed.fq.gz}_R2_trimmed.fq.gz" ; kallisto quant -i human2.idx -o "${R1%_R1_trimmed.fq.gz}_folder" ./$R1 ./$R2 ;done 
```
*_(**The output will be a folder with the sample name containing an abundance.tsv file which is the estimated counts for that exact sample**)_
#### ***[Note:]*** 
The Fasta file supplied can be either in plaintext or gzipped format. Note: Do not supply the genome Fasta file; the Fasta file must be a transcriptome Fasta.
#### ***[Note:]*** 
*kallisto is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads. It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment. Pseudoalignment of reads preserves the key information needed for quantification, and kallisto is therefore not only fast, but also as accurate as existing quantification tools. In fact, because the pseudoalignment procedure is robust to errors in the reads, in many benchmarks kallisto significantly outperforms existing tools*


## ***Fifth Step*** (Merging the Counts Files) -----------------------------
### Joining the matrices
#### 1️⃣ Define Output File
```{bash}
output_file="merged_counts.tsv"
```
#### 2️⃣ Get the First Sample File
```{bash}
first_sample=$(ls SRR1000_folder/abundance.tsv | head -1)
```
#### 3️⃣ Extract target_id Column
```{bash}
cut -f 1 "$first_sample" > "$output_file"
```
#### 4️⃣ Loop Through All abundance.tsv Files
```{bash}
for file in *_folder/abundance.tsv; do
sample=$(basename $(dirname "$file"))                   # Extract Sample Name
cut -f 4 "$file" > tmp_counts.tsv                       # Extract Count Column from Each File
paste "$output_file" tmp_counts.tsv > tmp_merged.tsv    # Merge This Sample's Counts with the Main File
mv tmp_merged.tsv ./"$output_file";
done
rm tmp_counts.tsv                                       # Clean Up Temporary File
echo "✅ Merged counts saved in $output_file"
