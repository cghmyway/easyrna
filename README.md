# **EASYRNA V1.0**

## easyrna (user-friendly Bulk RNA analysis R package)

## Global Object
**easyrna** The EasyRNA R package is a locally-based auxiliary package for RNA-seq data analysis, which includes functions for differential analysis, gene name conversion, gene expression conversion, and plotting. It lowers the barrier to entry for differential analysis.


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

#### 2.1 Convert read count to TPM


#### 2.2  Convert read count to FPKM


#### 2.3  Convert read count to CPM

#### 3. Conversion between ENSEMBL and SYMBOL names




