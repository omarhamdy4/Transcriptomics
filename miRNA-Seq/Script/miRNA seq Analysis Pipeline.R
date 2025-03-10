library(readr)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(factoextra)
library(ggfortify)
library(ggrepel)

setwd("E:/1.Fresh Grad/02_EgComBio2023/MODA miRNA")
data.path="E:/1.Fresh Grad/02_EgComBio2023/MODA miRNA/TCGA_BRCA/miRNA_samples_Transcriptome_profile"
pheno.path="E:/1.Fresh Grad/02_EgComBio2023/MODA miRNA/TCGA_BRCA"
file_pattern = "mirbase21.mirnas.quantification.txt"
##### load the miRNA-Seq data #####
files <- list.files(path=data.path,recursive=T, pattern = file_pattern)
# read the first file for the first time
file=files[1]
filepath=file.path(data.path,files[1])

temp <- read.table(filepath, header=T)
exp=as.data.frame(temp[,c(1)])
for(i in 1: length(files)){
    ## refer to the next file (note that we start from index 2, bec we already read the first file)
    file=files[i]
    file.id=strsplit(file,"/")[[1]][1]
    filepath=file.path(data.path,files[i])
    # read the next file  
    temp <- read.table(filepath, header=T)
    ## change the colname of the tpm with sample name 
    exp=cbind(exp,temp[,2])
    colnames(exp)[dim(exp)[2]]=file.id}
# Load Metadata
pheno <- read.delim("TCGA_BRCA/gdc_sample_sheet.2025-02-09.tsv")
table(pheno$`Sample.Type`)

sum(duplicated(rownames(exp)))  

mirnas.names=exp[,1]
exp=exp[,-1]
exp=apply (exp, 2,as.integer)
rownames(exp)=mirnas.names

# Changing column names of matrix from file ID to sample ID
exp2 <- exp
File.ID <- colnames(exp2)
exp2 <- rbind(exp2,File.ID)
tpheno <- pheno %>% data.frame() %>% select("File.ID","Sample.ID")
exp_final <- inner_join(data.frame(t(exp2)),tpheno, by="File.ID") %>% t()
sample.iD <- exp_final["Sample.ID",]
colnames(exp_final) <- sample.iD
exp_final <- exp_final[1:1881,]
miRNA.ID <- rownames(exp_final)
exp_final=apply (exp_final, 2,as.integer)
rownames(exp_final)=miRNA.ID

#Removing Zero variance miRNAs
dim(exp_final)
varrow <- apply(exp_final, 1, var, na.rm=TRUE)
cons_var <- (varrow == 0 | is.na(varrow))
exp_final <- exp_final[!cons_var,]
dim(exp_final)


## remove the low counts data(recomended by deseq2)
exp_final=exp_final[which(rowSums(exp_final) > 10),]
#check all samples are in metadata
all(colnames(exp_final) %in% pheno$Sample.ID)
#check all samples are in same order as metadata
all(colnames(exp_final) == pheno$Sample.ID)
#Re-order the matrix to be the same as in metadata
exp_final=exp_final[,pheno$Sample.ID]

table(pheno$Sample_Type)

dds = DESeqDataSetFromMatrix( countData = exp_final, colData = pheno , design = ~ Sample.Type)
dds.run = DESeq(dds)
### direct results or specifying the contrast (to make a res object based on two specific conditions/treatment)
#res=results(dds.run)
res=results(dds.run, contrast = c("Sample.Type","Primary Tumor" ,"Solid Tissue Normal") )
# remove nulls (recomended by deseq2)
res=res[complete.cases(res), ] 
summary(res)
plotMA(res)

res.df=as.data.frame(res)
res.dems=res.df[res.df$padj< 0.05,]

res.dems=res.dems[order(res.dems$padj), ]
dems=gsub("\\.","-",rownames(res.dems))
write.table(dems, file="dems_breast.txt", quote = F, col.names = T, row.names = F)

res.dems2=res.df[c(res.df$padj< 0.05&abs(res.df$log2FoldChange)>2),]
dems2=gsub("\\.","-",rownames(res.dems2))
write.table(dems2, file="dems_breast_diff.exp.txt", quote = F, col.names = F, row.names = F)
write.csv(res.dems2,"Dems_2.csv")


#### get the normalized and logged transformed values of all exp data
ntd=normTransform(dds)
exp.norm= assay(ntd)

### get the normalized expression levels of the dems ###################
dems.exp=exp.norm[dems,]
save(exp.norm,  pheno, file="MODA_miRNA_Seq.RDATA")
#load("MODA_miRNA_Seq.RDATA")

exp.norm.scaled=scale(exp.norm,center = TRUE, scale = TRUE)


boxplot(exp.norm.scaled[, 1:100], col = "purple",
        border = "blue",
        horizontal = FALSE,
        notch = FALSE,
        ylim=c(-1,max(exp.norm.scaled[, 1:100])))

df<- stack(as.data.frame(exp.norm.scaled)[, 1:30])
ggplot(df, aes(x= values, fill=ind))+geom_density(alpha=0.5)+theme_classic()+scale_x_continuous(limits=c(-2,max(exp.norm.scaled[, 1:100])))

options(repr.plot.width=10,repr.plot.height=8)
par(mar = c(8,5,2,2),mfrow=c(1,3))
boxplot(exp[1:100,1:100], main="Before log2" ,horizontal=T, names=colnames(exp)[1:100],las=2,
        col = "lightgreen")
boxplot(exp.norm[1:100,1:100], main="After log2" ,horizontal=T, names=colnames(exp.norm)[1:100],
        las=2,col = "lightgreen")
boxplot(exp.norm.scaled[1:100,1:100], main="After log2 + Scaled " ,horizontal=T,names=colnames(exp.norm.scaled)[1:100],
        las=2,col = "lightgreen",ylim=c(-1,max(exp.norm.scaled[, 1:100])))

# Create column annotations for the heatmap
pheno <- mutate(pheno,
                Condition=case_when(
                    Sample.Type == "Primary Tumor" ~ "Breast_Cancer",
                    Sample.Type =="Solid Tissue Normal" ~ "Healthy_Control",
                    TRUE ~ "Metastatic" #this replaces the non-mentioned other values to be x
                ))
pheno.sorted=pheno[order(pheno$Condition), ]
sam_condition <- pheno.sorted$Condition
column_annot <- HeatmapAnnotation(
    Condition = sam_condition,
    col = list(Condition = c("Breast_Cancer" = "purple", "Healthy_Control" = "violet","Metastatic"="pink")))

# Generate the heatmap for DEMs
Heatmap(
    matrix = dems.exp,
    top_annotation = column_annot,
    row_title = "DEMs",
    column_title = "Samples",
    cluster_rows = TRUE,
    cluster_columns = F,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 2),
    column_names_rot = 45,
    column_names_centered = TRUE)



options(repr.plot.width=10,repr.plot.height=8)
par(mfrow=c(1,2))
pca.2 <- prcomp(t(log2(exp.norm + 1)), scale. = T)
autoplot(pca.2, data = pheno, colour = 'Condition',frame = Z,label = F, label.size = 3,shape="Condition")+ theme_classic()
#only DEMs
pca.3 <- prcomp(t(log2(dems.exp + 1)), scale. = T)
autoplot(pca.3, data = pheno, colour = 'Condition',frame = F,label = F, label.size = 3,shape="Condition")+ theme_classic()




res_df <- data.frame(res)
res_df$gene_symbol <- rownames(res_df)
theme_set(theme_classic(base_size = 20) +
              theme(
                  axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                  axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                  plot.title = element_text(hjust = 0.5)))
res_df$diffexpressed <- "NA"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "DOWN"
#head(res_df[order(res_df$padj) & res_df$diffexpressed == 'DOWN', ])
#for down miRNAs
res_df$delabel <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$log2FoldChange), "gene_symbol"],5), res_df$gene_symbol, NA)
#for UP miRNAs
res_df$delabel2 <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$log2FoldChange,decreasing=T), "gene_symbol"], 5), res_df$gene_symbol, NA)

ggplot(res_df, aes(x = log2FoldChange, y = -log(padj), col = res_df$diffexpressed),label=delabel) +
    geom_point(alpha = 0.4) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    labs(
        title = "Volcano Plot",
        x = "log Fold change",
        y = "-log10 adjusted Pvalue"
    ) +
    geom_vline(xintercept = c(1, -1), col="black",linetype = "dashed") +
    geom_hline(yintercept = -log(0.05), col="black",linetype = "dashed")+
    geom_point(size = 2)+
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), # to set the colours of our variable
                       labels = c("Downregulated", "Not significant", "Upregulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    coord_cartesian(ylim = c(0, 250), xlim = c(-8, 8)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Regulation', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adj.p-value")) + 
    scale_x_continuous(breaks = seq(-10, 10, 2))+ # to customise the breaks in the x axis
    ggtitle('Breast Cancer miRNAs volcano plot cancer vs healthy patients')+  # Plot title 
    geom_text_repel(label=res_df$delabel,max.overlaps = Inf)+ # To show all labels 
    geom_text_repel(label=res_df$delabel2,max.overlaps = Inf)





#ROC Curve
x <- select(pheno, "Condition")
colnames(x) <- "disease_state"
x <- mutate(x,
            value=case_when(
                disease_state =="Breast_Cancer" ~ "1",
                disease_state =="Healthy_Control" ~ "0",
                TRUE ~ NA #this replaces the non-mentioned other values to be x
            ))
x$value <- as.numeric(x$value)
Nx2 <- t(exp.norm.scaled)
CallMe <- function(q) {
    x$ex <- Nx2[,q]
    # Remove rows with missing values in x$value or x$ex
    x <- na.omit(x)
    glm.fit=glm(value~ex,data=x, family=binomial)
    # Debugging: Check lengths
    print(paste("Length of x$value:", length(x$value)))
    print(paste("Length of fitted values:", length(glm.fit$fitted.values)))
    library(pROC)
    par(pty="s")
    roc(x$value, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE,
        xlab="100-Specificity",ylab="Sensitivity", col="#377eb8",lwd=4, print.auc=TRUE, print.auc.x=45)    
    legend("bottomright" ,legend=c(q), col=c("#377eb8"), lwd=4)
}
CallMe("hsa.mir.940")



metastatic_samples <- unlist(c("TCGA-E2-A15E-06A","TCGA-E2-A15K-06A","TCGA-AC-A6IX-06A","TCGA-BH-A1ES-06A","TCGA-BH-A1FE-06A","TCGA-E2-A15A-06A","TCGA-BH-A18V-06A"))
#remove the metastatic samples
Nx3 <-exp.norm.scaled[,!colnames(exp.norm.scaled)%in% metastatic_samples]
Nx3 <- t(Nx3)
CallMeLOOP7<- function(file_path) {
    dataz <- read.csv(file_path, header = TRUE)
    mirnas <- dataz[, 1]
    auc_values <- vector(length = length(mirnas))  # Initialize vector to store AUC values
    for (i in seq_along(mirnas)) {
        q <- mirnas[i]
        x$ex <- Nx3[,q]
        x <- na.omit(x)
        glm.fit <- glm(value ~ ex, data = x, family = binomial)
        library(pROC)
        auc <- auc(x$value, glm.fit$fitted.values)
        auc_values[i] <- auc
        print(paste("mir",i,"named",q,"is finished with AUC =",auc_values[i]))
    }
    dataz$AUC <- auc_values  # Add AUC values as a new column
    write.csv(dataz, file = file_path, row.names = FALSE)  # Overwrite the original file with the updated data
}

CallMeLOOP7("E:/1.Fresh Grad/02_EgComBio2023/MODA miRNA/Dems_2.csv")

