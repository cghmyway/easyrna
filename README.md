# **EASYRNA V1.5**

## Easyrna ( R Package of User-Friendly Bulk RNA Analysis)

<p align="center">
  <img src="https://github.com/user-attachments/assets/12d8dde7-a662-4519-9b99-eb0289e08847" alt="Example Image" width="50%" style="margin-top: 600mm; margin-bottom: 600mm;">
                                                          
</p>



## Global Object
The ***EasyRNA*** R package is a locally-based auxiliary package for RNA-seq data analysis, which includes functions for differential analysis, gene name conversion, gene expression conversion, and plotting. It lowers the barrier to entry for differential analysis.


## How to install
The development version can be installed using devtools:
```r
# install.packages("devtools")
devtools:install_github("cghmyway/easyrna")
```

## Related protocols

### 🚀1. Before any analysis, prepare the local database.
📢 Prepare **GTF**. Obtained from the [ENSEMBL](https://ftp.ensembl.org/pub/release-112/gtf/) website.
```r
database <- easyrna::refer("/path/to/your/gtf")
```
### 🌱2. Transformation of gene expression estimation methods
#### 🌱2.1 Convert read count to TPM
!!!The GeneReadCount files from STAR or featureCounts need to be integrated into a data frame with gene names as rows and sample names as columns, as shown below !!!


<p align="center">
  <img src="https://github.com/user-attachments/assets/7ad07a07-fb0a-42de-b8b5-4b6dde0695e5" alt="Example Image" width="50%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>

```r
head(ReadCount.data)
data <- geneLength( df = ReadCount.data, gtf = "/path/to/gtf.file")
data_TPM <- countToTpm(expmat = data, effLen_col = "est_len")
```

#### 🌱2.2  Convert readcount to FPKM

```r
head(ReadCount.data)
data <- geneLength( df = ReadCount.data, gtf = "/path/to/gtf.file")
data_FPKM <- easyrna::countToFpkm(expmat = data,effLen_col = "est_len" )
```

#### 🌱2.3  Convert readcount to CPM

```r
data_CPM <- easyrna::countToCPM(expmat = Readcount.data)
```

### 🚀3. Conversion between ENSEMBL and SYMBOL names

```r
data$gene_symbol <- findname(query = rownames(data), 
                           dataset = dataset,
                              type = "ENSEMBL") # type: "ENSEMBL";"SYMBOL"
```

### 🚀4. Analysis of differentially expressed genes using Deseq2

```r

database <- refer("/path/to/gtf") # make database
# desCompare return a List including 
res <- easyrna::desCompare(readcount = data,  # ReadCount file
                         group1.name = "A",   # Treatment group name, such as "OV"
                       group1.number = 3,     # Number of group1 sample
                         group2.name = "B",   # Control group name, such as "Control"
                       group2.number = 3,     # Number of group2 sample
                             dataset = database,  # database made by refer function
                                type = "ENSEMBL",  # Type of gene name in ReadCount file
                    log2FC.threshold = 1,
                    pvalue.threshold = 0.05,
                      padj.threshold = 0.05,
                               top.n = 10)  # The number of genes that need to be labeled for downstream plotting.
# save result
# write.csv(res[["DEG_Dataframe"]], "DEGs.csv")

# DEGs number
# table(res[["DEG_Dataframe"]]$Group)

```

##
### 🤩 5.Visualization 🤩

#### 🌱5.1 Volcano for bulk RNA-seq

```r

plot_volcano(deg.data = res[["DEG_Dataframe"]], log_scale="log10p")

```


<p align="center">
  <img src="https://github.com/user-attachments/assets/33374d52-528c-401b-8a5c-8ad5a8275c12" alt="Example Image" width="30%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>

#### 🌱5.2 PCA analysis for bulk RNA-seq

```r

pca_plot(
data = data_tpm, # A data frame of the TPM or FPKM values of genes, with row names as gene names and column names as sample names.
groups = groups,   #	Group names
group_positions = group_positions, 	# A list for the distribution of samples of each group in the data frame.
colors = colors, # Point colors
ellipse = TRUE) # Whether to draw a confidence ellipse. Default:FALSE
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/d2f0d7fe-c1af-43bd-9091-1ee2feb4c438" alt="Example Image" width="30%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>

#### 🌱5.3 Plot for correlation heatmap

```r
plot_correlation_heatmap(
data = data_tpm,  # A data frame of the TPM or FPKM values of genes, with row names as gene names and column names as sample names.
method = "spearman", # Correlation calculation method: "pearson", "spearman", "kendall"
sample_annotation = sample_annotation, #  Sample annotation data frame. The row names are sample names, and the column names are different traits or groups
annotation_colors = list(Group=c("A"="red3","B"="blue3"),Gender=c("Male"="purple", "Female"="yellow2")))
```
Sample annotation file:
```r
> sample_annotation
       Group Gender  Source
a1_TPM     A   Male Tissue1
a2_TPM     A Female Tissue2
a3_TPM     A Female Tissue1
b1_TPM     B   Male Tissue3
b2_TPM     B   Male Tissue2
b3_TPM     B Female Tissue1
```


<p align="center">
  <img src="https://github.com/user-attachments/assets/03ceec49-6751-47fb-8c21-781dc061fc81" alt="Example Image" width="50%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>

#### 🌱5.4 Plot of Gene expression between groups

```r
# Using boxplot
plot_gene_boxplot(data = data_tpm, groups = c("A","B"), group_positions = list(1:3,4:6), gene_of_interest = c("ENSGALG00000012172"), method = "t.test")

```
<p align="center">
  <img src="https://github.com/user-attachments/assets/6be56e82-ba26-4f2d-854b-ec700075348a" alt="Example Image" width="30%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>


```r
# Using violin
plot_gene_violin(data = data_tpm, groups = c("A","B"), group_positions = list(1:3,4:6), gene_of_interest = c("ENSGALG00000012172"), method = "t.test")

```

<p align="center">
  <img src="https://github.com/user-attachments/assets/a4d2175c-87bf-4774-8e18-b8809760892e" alt="Example Image" width="30%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>


#### 🌱5.5 GO/KEGG

```r
# GO
run_GO_analysis(des_res = des_res, topn = 15,  up_color = "red4",down_color = "green4",font_size = 3)

# des_res: The results of differentially expressed genes from desCompare function.
```
<p align="center">
  <img src="https://github.com/user-attachments/assets/04784a14-27ad-491b-977f-ad1b77a56af0" alt="Example Image" width="30%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>





```r
# KEGG
run_KEGG_analysis(des_res, pvalueCutoff = 0.05, qvalueCutoff = 0.05, font_size =2, point_size = 8)
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/f94fa911-86c2-42da-9fff-dc318e91ca52" alt="Example Image" width="30%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
</p>




###🚀 6. Visualization of alternative splicing results

📢Provides wrapper functions for rMATs result visualization

```r
event_files <- list(SE = "SE.MATS.JC.txt",RI = "RI.MATS.JC.txt",MXE = "MXE.MATS.JC.txt",A5SS = "A5SS.MATS.JC.txt",A3SS = "A3SS.MATS.JC.txt");
process_rMATS_results(event_files)

# param:
# event_files: A list showing the Alternative Splicing Event files
# pvalue_cutoff: p-value threshold for differential splicing events
# diff_cutoff: PSI threshold for differential splicing events
# orgdb: GO and KEGG enrichment analysis databases. Default: org.Hs.eg.db
# go_analysis: Whether to conduct GO analysis
# kegg_analysis: Whether to conduct  KEGG analysis
# save_results: Whether to save the result in csv file
# output_prefix: Output file prefix
```



