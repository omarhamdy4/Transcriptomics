# RNA-Seq Primary Analysis Pipeline
***
![RNA-seq Primary Analysis Pipeline pipeline template](https://github.com/user-attachments/assets/3e921b76-6028-488b-92f5-dc2640a2e16a)

## ***Zero Step*** (Samples Download) --------------------------------------
#### ***[Note:]*** 
You can download samples using SRAtoolkit from the command line, but in my case I downloaded the samples using SRA run selector tool on NCBI (SRA) database
```{bash}
mkdir ~/Transcriptome_NGS 
cd ~/Transcriptome_NGS
```
## ***First Step*** (Quality Check I) --------------------------------------
### Check the quality of all the samples (using FastQC)
```{bash}
Do fastqc for all .fastq files
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
*_(**The output will be a folder with the sample name containing an abundance.tsv file in it for that exact sample**)_

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
