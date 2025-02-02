# RNA Microarray Analysis Pipeline 
>***
![Script workflow template](https://github.com/user-attachments/assets/0d29a885-c60a-4475-a9b2-a46618e53b81)
## 1- Library loading

> ------------------------------------------------------------------------

```{r setup, include=FALSE}
library(dplyr)
library(tidyverse)
library(affy)
library(affyPLM)
library(hgu133plus2.db)
library(GEOquery)
library(mouse4302.db)
library(limma)
library(ComplexHeatmap)
library(ggplot2)
library(ggfortify)
library(rgl)        #For 3D PCA plot
library(ggrepel)
library(RColorBrewer)

se
```

## 2- Reading Samples
> ------------------------------------------------------------------------

```{r}
setwd("E:/1.Fresh Grad/02_EgComBio2023/MODA microarray")
celFilesDirectory="GSE175844"
celFilsID= "GSE175844"
affyData = ReadAffy(celfile.path = celFilesDirectory, cdfname = "mouse4302cdf")
affyData
phenotable = read.csv("Metadata_GSE175844.csv")
```

```{r Explore the affy object}
head(sampleNames(affyData))
head(featureNames(affyData))
annotation(affyData)          # vip
View(exprs(affyData))   #Check the values if they were large or not i.e. need normalization or not
range(exprs(affyData))
```

## 3- Pre-processing
> ------------------------------------------------------------------------

### 3- (a) Expression matrix

```{r}
eset = threestep(affyData, 
                 background.method = "IdealMM",
                 normalize.method = "quantile", 
                 summary.method = "average.log")
range(exprs(eset))
data = exprs(eset)
```

### 3- (b) Phenodata

```{r}
phenotable = phenotable %>% 
  dplyr::mutate(sampleType = str_extract(Sample_source_name_ch1, "control|KLF4"))
```

## 4- Visualize after Pre-processing
> ------------------------------------------------------------------------

```{r}
par(mfrow = c(2, 2))
hist(affyData, main="Histogram of Raw data")        # Raw data histogram
hist(eset, main ="Histogram after 1ry analysis")    # Preprocessed data histogram
boxplot(affyData, main = "boxplot of Raw data", col = seq(1:ncol(eset)), las = 2)
boxplot(eset, main = "boxplot after 1ry analysis", col = seq(1:ncol(eset)), las = 2)
```

## 5- Gene Mapping
> ------------------------------------------------------------------------

```{r}
probe_id = rownames(data)
datax = cbind(probe_id, data)

ls("package:mouse4302.db")
mapper = mouse4302SYMBOL
map.df = as.data.frame(mapper)
data2 = merge(data,map.df, by ="probe_id" ,all.x=T)
data2 = data2[-1]

```

## 6- Quality Control (QC)
> ------------------------------------------------------------------------

```{r}
# remove all NAs
dim(data2)      
nulls = is.na(data2$symbol)
sum(nulls)          
data2 = data2[!nulls, ]
dim(data2)          

# Remove all duplication
duplicated = duplicated(data2$symbol)
sum(duplicated)          # 22274
# Solution is agreggation
exp.data = data2 %>% 
  dplyr::select(-last_col())
exp.data = apply(exp.data, 2, as.numeric)
mode(exp.data)
exp.data.agg = aggregate(exp.data, by = list(data2$symbol), FUN = mean)
dim(exp.data.agg)
rownames(exp.data.agg) = exp.data.agg$Group.1
exp.data.agg = exp.data.agg[-1]

# Remove genes with zero variance
varrow <- apply(exp.data.agg, 1, var, na.rm=TRUE)
cons_var <- (varrow == 0 | is.na(varrow))
exp.data.agg <- exp.data.agg[!cons_var,]


#Renaming the samples (colnames)
names(exp.data.agg)
names(exp.data.agg) = sapply(strsplit(names(exp.data.agg), "_"), `[`, 1)
```



```{r Temporary save}
exp = exp.data.agg
save(exp, phenotable, file ="GSE175844.RDATA")
load("GSE175844.RDATA")
```

## 7- Differential Expression Analysis (Limma)
> ------------------------------------------------------------------------

```{r}
all(colnames(exp) == phenotable$Sample_geo_accession)
unique(phenotable$sampleType)
groups = factor(phenotable$sampleType, levels = c("KLF4", "control"))
design = model.matrix(~0 + groups)
colnames(design) = levels(groups)
design

fit = lmFit(exp, design)
contrast.matrix = makeContrasts(KLF4VsCTRL = KLF4 - control,levels = design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)
# Now extract the degs 
degs_KLF4 = topTable(fit2, coef = "KLF4VsCTRL", number = Inf, adjust = "BH")
```

## 9- Downstream plots
>***

### 9- (a) Heatmap for top 100 DEGs

```{r Heatmap Style 1}
deg_100 <- degs_KLF4[order(degs_KLF4$P.Value, -abs(degs_KLF4$logFC)),]   #sort to get top 100
deg_100 <- degs_KLF4[1:100,]                                         #select only them
deg_100_exp <- as.matrix(exp[rownames(deg_100),])                     #get their normalized expression values

# Create column annotations for the heatmap
sam_condition <- phenotable$sampleType
column_annot <- HeatmapAnnotation(
  Condition = sam_condition,
  col = list(Condition = c("KLF4" = "purple", "control" = "violet")))

# Generate the heatmap
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
```{r Heatmap Style 2}
ann_df <- data.frame(phenotable$Sample_geo_accession,phenotable$sampleType)
ann_rownames <-unlist(ann_df[,1])
ann_df <- data.frame(ann_df[,-1])
rownames(ann_df) <- ann_rownames
colnames(ann_df) <-"sampleType"
ann_colors <- list(Group = c("KLF4" = "violet",
              "control" = "purple"))
pheatmap(deg_100_exp, 
      col = brewer.pal(10, 'RdYlGn'), # choose a colour scale for your data
      #color = colorRampPalette(rev(brewer.pal(n = 7, name =  "RdYlBu")))(100),
      cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
      clustering_distance_cols = 'euclidean',
      clustering_distance_rows = 'euclidean',
      clustering_method = 'ward.D',
      #annotation_row = gene_functions_df, # row (gene) annotations
      annotation_col = ann_df, # column (sample) annotations
      annotation_colors = ann_colors, # colours for your annotations
      annotation_names_row = F, 
      annotation_names_col = F,
      angle_col = "45" , # sample names at an angle
      legend_breaks = c(1.135840 , 3.967026, 6.798213), # legend customisation / I got the no. from >range(deg_100_exp)
      legend_labels = c("Low", "Medium", "High"), # legend customisation
      show_colnames = T, show_rownames = T, # displaying column and row names
      fontsize_row = 5,          # row label font size
      fontsize_col = 9,          # column label font size 
      main = "Top 100 DEGs Heatmap", # a title for our heatmap
      cutree_rows = 4, cutree_cols = 2)
```

### 9- (b) Perform PCA on DEGs

```{r PCA again}

pca_result <- prcomp(t(deg_100_exp), scale. = TRUE)
pca_data <-pca_result$x
pca_result$x
# Plot 2D PCA
autoplot(pca_result, data = phenotable, colour = 'sampleType',frame = T,label = T, label.size = 3,shape="sampleType")

# Plot 3D PCA
mycolors <- ifelse(phenotable$sampleType == "KLF4", "red", "lightgreen")
plot3d(pca_result$x[, 1:3], col = mycolors, size = 12, type = "s", main = "3D PCA Plot")
```

### 9- (c) Volcano plot

```{r Volcano plot}
res_df <- data.frame(degs_KLF4)
res_df$gene_symbol <- rownames(res_df)

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)))
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
res_df$diffexpressed <- "NA"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
res_df$diffexpressed[res_df$logFC > 1 & res_df$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_df$diffexpressed[res_df$logFC < -1 & res_df$P.Value < 0.05] <- "DOWN"
# Explore a bit
head(res_df[order(res_df$P.Value) & res_df$diffexpressed == 'DOWN', ])
# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
res_df$delabel <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$P.Value), "gene_symbol"], 10), res_df$gene_symbol, NA)


ggplot(res_df, aes(x = logFC, y = -log(P.Value), col = res_df$diffexpressed),label=delabel) +
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
coord_cartesian(ylim = c(0,20), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Regulation', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2))+ # to customise the breaks in the x axis
ggtitle('KLF4 overexpressed cells vs Control')+  # Plot title 
geom_text_repel(label=res_df$delabel,max.overlaps = Inf) # To show all labels 

```
