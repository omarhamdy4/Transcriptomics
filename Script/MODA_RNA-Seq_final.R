###########################################
#  - Multi Omics Course                   # 
#  - Gene Expression Profiling of Prostate#
#  - TCGA-PRAD project - RNAseq data      # 
#  - Assignment 4                         #    
#  - 2024- 1-24                           #
#  - Copyright: @Mohamed Hamed            #
###########################################
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

####load libraries ########  
library(readr)
library(org.Hs.eg.db)
library("vsn")
library(DESeq2)


#change the path parameters according to the location of your data in your laptop.
data.path="TCGA_data/TCGA_PRAD_STAR"
pheno.path="TCGA_data/gdc_sample_sheet.2022-12-09.tsv"

##### load the mRNA-Seq data #####
files <- list.files(path=data.path,recursive=T, pattern = "tsv")

# read the first file for the first time
file=files[1]
filepath=file.path(data.path,files[1])
colnames= unlist( strsplit ( readLines(filepath ,n=2) [2] ,"\t"))

temp <- read.table(filepath, header=F,skip = 6)
names(temp)=colnames
View(temp)
genetype_mapper=temp[,c(2,3)]

#create a storing object exp to save the whole readcounts of each file TPMs read in an iteration
exp=temp[temp$gene_type == "protein_coding" ,c(1,2)]

for(i in 1: length(files))
{
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  filepath=file.path(data.path,files[i])
  
  # read the next file  
  temp <- read.table(filepath, header=F,skip = 6)
  temp=temp [ temp[,3] == "protein_coding" , ]
  
  ## change the colname of the tpm with sample name 
  exp=cbind(exp,temp[,4])
  colnames(exp)[dim(exp)[2]]=file.id
}
View(exp)

#########################  Assignment #########################
## load the data using the two other approaches in the following link
# a professional example to follow is here using two ways 
#( reading from files , TCGAbiolink package) >>>  https://www.biostars.org/p/9500223/

### General Reminder
# don't forget to master the tidyvers packages https://www.tidyverse.org/packages/


# check duplciation of of gene symbols?  
x=duplicated(exp$gene_name)  
sum(x)

### yes .. why ? transcripts?  solutions : aggregation
exp.data=exp[ , 3:dim(exp)[2]]
exp.data=apply(exp.data,2, as.numeric)

#### remove  duplication by aggregation
exp.data.agg= aggregate(exp.data, list(exp$gene_name),FUN=mean)
genes=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[-1]
exp.data.agg=apply(exp.data.agg,2, as.numeric)
exp.data.agg=round(exp.data.agg)
rownames(exp.data.agg)=genes

file.ids=colnames(exp.data.agg)

###### load the mrna sample sheets  # sample sheets
pheno <- read_delim(pheno.path,"\t", escape_double = FALSE, trim_ws = TRUE)
View(pheno)
table(pheno$`Sample Type`)
names(pheno)=sub (" ", "_",names(pheno)) # rename the names of the pheno to add "_" instead of " "

#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 

file.ids.pheno=pheno$File_ID
index.files=match(file.ids,file.ids.pheno)
colnames(exp.data.agg)=pheno$Sample_ID [index.files]

save(exp.data.agg,genetype_mapper,  pheno, file="MODA_RNA_Seq.RDATA")




#### Exploratory analysis + filtration process : plz do it  ####################

# box plot , histogram , PCA

## remove the least 10 % variable genes ### plz do it here



#### Do differential analysis using Deseq2 package as it works on readcounts  ##########
## or use Limma package instead >> you can follow this tutorial here https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
##. read this link : https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

table(pheno$Sample_Type)

cond1="Solid Tissue Normal" 
cond2="Primary Tumor"

dds = DESeqDataSetFromMatrix( countData = exp.data.agg, colData = pheno , design = ~ Sample_Type)
dds.run = DESeq(dds)
### direct results or specifying the contrast (to make a res object based on two specific conditions/treatment)
#res=results(dds.run)
res=results(dds.run, contrast = c("Sample_Type",cond1 ,cond2) )

# remove nulls
res=res[complete.cases(res), ]
summary(res)


res.df=as.data.frame(res)
write.table(res.df, file = "res.txt")

plotMA(res, ylim=c(-1,1)) 
summary (res)

res.degs=res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>log2(2) ,  ]
res.degs=res.degs[order(res.degs$padj), ]

degs.genes=rownames(res.degs)

# or export them to a gene list to start the gene set enrichment analysis.
write.table(degs.genes, file="degs_prostate.txt", quote = F, col.names = F, row.names = F)

#### get the normalized and logged transformed values of all exp data
#using the the variance stabilizing transformation. vsn package

ntd=normTransform(dds)
exp.norm= assay(ntd)


### get the normalized expression levels of the degs ###################
degs.exp=exp.norm[degs.genes, ]

## 1- creating a heatmap for the top 100 DEG genes ## use always ComplexHeatMap
#  get the expression profiles of the top 100 degs only and create heatmap


### 2- creating  2D PCA for all degs and look how it segregates the normal and tumor samples 


### 3-creating  3D PCA  all degs and look how it segregates the normal and tumor samples 


### 4- creating  Volcano Plot using the (advanced methods)


### 5- DO enrichment analysis using DAVID, GeneTrail, Enrichr 


### 6- visualize the enrichment analysis results using GOPLOT Bioconductor package 


### 7- Run L1000CDS2 tool to characterize the small molecules that could reverse the disease signature 


### 8- creating  a schematic figure for study design (Bonus) for the prostate dataset(TCGA -PRAD)

