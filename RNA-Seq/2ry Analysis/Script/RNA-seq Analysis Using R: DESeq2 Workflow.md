# RNA-seq Analysis Using R: DESeq2 Workflow
***
>#### This notebook demonstrates a step-by-step workflow for analyzing RNA-seq data using DESeq2. The analysis includes preprocessing the expression matrix, handling metadata, identifying differentially expressed genes (DEGs), and Plotting important figures such as volcano plots, Heatmap, and PCA.
## The Workflow goes as follow:
![Script Workflow template biorender](https://github.com/user-attachments/assets/54607715-5edc-40dc-8a71-ce8b3c55d34e)
***
## 1- Loading libraries
```{r}
library(readr)
library(org.Hs.eg.db)
library(vsn)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(rgl)
```
***
## 2- Loading the data and metadata
```{r}
read_counts <- read.delim("E:/1.Fresh Grad/02_EgComBio2023/MODA RNA_seq/GSE275290_raw_counts.tsv")
metadata <- read.csv("E:/1.Fresh Grad/02_EgComBio2023/MODA RNA_seq/Phenotable.csv")
```
***
## 3- Exploratory Data Analysis: Box plot and Histogram
```{r Exploration}
boxplot(log2(read_counts[,-1]+1), main="Exploratory Box plot", ylab="log2_count", las=2)
hist(as.matrix(log2(read_counts[,-1]+1)), main="Exploratory Histogram", xlab="Sample", ylab="log2_count", breaks=50)
```
***
## 4- Gene Annotation
#### Convert gene IDs to gene symbols using org.Hs.eg.db.
```{r Annotation}
gene_ids <- as.character(read_counts$GeneID)  # Ensure gene IDs are characters
gene_symbol <- mapIds(
  org.Hs.eg.db, keys = gene_ids, keytype = "ENTREZID",
  column = "SYMBOL", multiVals = "first")
gene_df <- data.frame(
  GeneID = names(gene_symbol),
  Gene_symbol = as.vector(gene_symbol),
  stringsAsFactors = FALSE)
data <- merge(read_counts, gene_df, by="GeneID", all.x=TRUE)
data <- data %>%
  dplyr::select(Gene_symbol, everything(), -GeneID)
```
***
## 5- Preprocessing
###    (a) metadata
```{r Preprocesing I}
meta <- metadata %>%
  select(Sample.Name, treatment) %>%
  rename(sampleid = Sample.Name, Condition = treatment)

meta <- meta %>%
  mutate(Condition = case_when(
    Condition == "FBZ for 48 hours" ~ "FBZ",
    Condition == "DMSO; Vehicle" ~ "Control",
    TRUE ~ Condition))
```
###    (b) Preprocess expression matrix
#### Clean the data (NAs, Duplicates, Zero variance genes)
#### A- Removing NAs
```{r Preprocesing II}
sum(is.na(data))
data <- na.omit(data)
```
#### B- Handling Duplicates
```{r}
sum(duplicated(data$Gene_symbol))
exp_data <- data %>% 
  dplyr::select(-Gene_symbol)
is.numeric(exp_data)  # Check if the data is numeric
exp_data_agg <- aggregate(exp_data, by=list(data$Gene_symbol), FUN=mean)
sum(duplicated(exp_data_agg$Group.1))  # Check for duplicated gene symbols after aggregation
row.names(exp_data_agg) <- exp_data_agg$Group.1
exp_data_agg <- exp_data_agg[,-1]
```
#### C- Removing Zero variance genes
```{r}
varrow <- apply(exp_data_agg, 1, var, na.rm=TRUE)
cons_var <- (varrow == 0 | is.na(varrow))
exp_data_agg <- exp_data_agg[!cons_var,]
```
***
## 6- Follow the exploratory analysis (PCA)
```{r PCA}
pca <- prcomp(t(log2(exp_data_agg + 1)), scale. = TRUE)
autoplot(pca, data = meta, colour = 'Condition',frame = T,label = T, label.size = 3,shape="Condition")
```
***
## 7- Differential expression analysis (DEseq2)
#### Prepare the data
```{r DEseq2}
exp <- exp_data_agg
all(colnames(exp) %in% meta$sampleid)
exp <- exp[, meta$sampleid]
exp <- round(exp)
meta$Condition <- factor(meta$Condition, levels = c("Control", "FBZ"))
class(meta$Condition)  # Check the class of the Condition column
```
#### Create DESeqDataSet object and Filter out low-count genes
```{r Follow DEseq2}
dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = meta,
  design = ~ Condition)
dds <- dds[rowSums(counts(dds)) >= 10,]
```
#### Run DESeq2 analysis
```{r}
dds_run <- DESeq(dds)
# Extract results for the contrast FBZ vs Control
res <- results(dds_run, contrast=c("Condition", "FBZ", "Control"), alpha=0.05)
res <- res[complete.cases(res),]  # Remove rows with NA values
summary(res)  # Summarize the results
```
***
## 8- Data normalization (VST)
```{r Normalization}
vsn_norm <- counts(dds, normalized=FALSE) %>% vsn2 %>% exprs
library(hexbin)
meanSdPlot(vsn_norm)
write.csv(vsn_norm, "vsn_norm.csv")
```
***
### 9- DEGs extraction
```{r DEGs}
degs <- res[res$padj < 0.05 & res$log2FoldChange > 1,]
degs.genes <- rownames(degs)
degs.exp <- vsn_norm[degs.genes,]
write.csv(degs.exp, "degs.exp.csv")
```
***
## 10- Downstream plots
###    (a) Heatmap for top 100 DEGs
```{r Heatmap}
deg_100 <- degs[order(degs$padj, -abs(degs$log2FoldChange)),]   #sort to get top 100
deg_100 <- degs[1:100,]                                         #select only them
deg_100_exp <- vsn_norm[rownames(deg_100),]                     #get their normalized expression values
sam_condition <- meta$Condition
column_annot <- HeatmapAnnotation(
  Condition = sam_condition,
  col = list(Condition = c("Control" = "limegreen", "FBZ" = "orange")))
Heatmap(
  matrix = deg_100_exp,
  top_annotation = column_annot,
  row_title = "Top 100 DEGs",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 4),
  column_names_rot = 45,
  column_names_centered = TRUE
)
```
###    (b) PCA on DEGs
```{r PCA again}
degs_normalized <- vsn_norm[rownames(vsn_norm) %in% rownames(degs),]
pca_result <- prcomp(t(degs_normalized), scale. = TRUE)
pca_data <-pca_result$x
pca_result$x
# Plot 2D PCA
autoplot(pca_result, data = meta, colour = 'Condition',frame = T,label = T, label.size = 3,shape="Condition")
# Plot 3D PCA
mycolors <- ifelse(meta$Condition == "FBZ", "red", "lightgreen")
plot3d(pca_result$x[, 1:3], col = mycolors, size = 12, type = "s", main = "3D PCA Plot")
```
###    (c) Volcano plot
```{r Volcano plot}
res_df <- data.frame(res)
res_df$gene_symbol <- rownames(res_df)
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)))
res_df$diffexpressed <- "NA"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$pvalue < 0.05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$pvalue < 0.05] <- "DOWN"
head(res_df[order(res_df$padj) & res_df$diffexpressed == 'DOWN', ])
res_df$delabel <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$padj), "gene_symbol"], 30), res_df$gene_symbol, NA)

ggplot(res_df, aes(x = log2FoldChange, y = -log(padj), col = res_df$diffexpressed),label=delabel) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(
    title = "Volcano Plot",
    x = "log Fold change",
    y = "-log10 adjusted Pvalue"
  ) +
  geom_vline(xintercept = c(1, -1), col="black",linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col="black",linetype = "dashed")+
    geom_point(size = 2)+
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
coord_cartesian(ylim = c(0, 150), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Regulation', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2))+ # to customise the breaks in the x axis
ggtitle('Ovarian cells volcano plot cancer vs healthy patients')+  # Plot title 
geom_text_repel(label=res_df$delabel,max.overlaps = Inf) # To show all labels 

```
