# **EASYRNA V1.0**

## Easyrna ( R Package of User-Friendly Bulk RNA Analysis)

<p align="center">
  <img src="https://github.com/user-attachments/assets/c2a36d3b-dce8-434d-b1b0-2abdbfbcd53d" alt="Example Image" width="50%" style="margin-top: 600mm; margin-bottom: 600mm;">
                                                          
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

#### 1. Before any analysis, prepare the local database.
prepare **GTF**. Obtained from the [ENSEMBL](https://ftp.ensembl.org/pub/release-112/gtf/) website.
```r
database <- easyrna::refer("/path/to/your/gtf")
```
#### 🌱2. Transformation of gene expression estimation methods
#### 🌱2.1 Convert read count to TPM
!!!The GeneReadCount files from STAR or featureCounts need to be integrated into a data frame with gene names as rows and sample names as columns, as shown below !!!


<p align="center">
  <img src="https://github.com/user-attachments/assets/581a826f-e3f0-4e38-a4e6-d7ace013c987" alt="Example Image" width="50%" style="margin-top: 60px; margin-bottom: 60px;">
                                                          
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

#### 3 Conversion between ENSEMBL and SYMBOL names

```r
data$gene_symbol <- findname(query = rownames(data), 
                           dataset = dataset,
                              type = "ENSEMBL") # type: "ENSEMBL";"SYMBOL"
```













