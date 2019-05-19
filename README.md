# scTHI
single cell Tumor Host Interaction 
### Introduction
The detection of tumor-host interaction is a major challange in cancer understanding. More recent advances in single cell next generation sequencing have provided new insights in unveiling tumor microenvironment (TME) cells comunications.
The scTHI package provides some useful functions to identify and visualize the cell types making up the TME, and a rank-based approach to discover interesting ligand-receptor interactions establishing among malignant and non tumor cells.
 

### Quick start
```{r,message=FALSE,warning=FALSE, eval=FALSE, echo=TRUE}
library(devtools)
install_github("miccec/yaGST")
install_github("miccec/scTHI")
library(scTHI)
```
### Input Data
As input, the scTHI package expects a matrix of count data (or normalized counts) obtained from single cell RNA-seq experiment, where rows are genes presented with Hugo Symbols and columns are cells. For a practical demonstration, let's use the scRNA-seq data available on the Broad Institute Single-Cell Portal (https://portals.broadinstitute.org/single_cell/study/single-cell-analysis-in-pediatric-midline-gliomas-withhistone-h3k27m-mutation) [1]. The dataset provides 3,321 scRNA-seq profiles from six primary Glioma (H3K27M-glioma), and includes both malignant cells and several tumor microenvironment cell types, as immune cells and oligodendrocytes. Our goal is to identify the significant ligand-receptor interactions that are established between cancer and immune cells present in the tumor microenvironment. For this reason, after downloading H3K27M-glioma data we preprocessed them, selecting only one sample of interest (i.e. BCH836), removing zero genes, trasforming data in TMP e log2, and appling quantile normalization. Processed data can be loaded as follows:

```{r,message=FALSE,warning=FALSE, eval=FALSE, echo=TRUE}
library(scTHI)
```
```{r,message=FALSE,warning=FALSE, eval=FALSE, echo=TRUE}
load(system.file("extdata", "H3K27_PatientBCH836.RData", package = "scTHI", mustWork = TRUE))
```
H3K27.meta includes the annotation of each cell, so as to identify the tumor cells and the immune cells in BCH836 patient:
```
## 
##          Filter     Immune cell       Malignant Oligodendrocyte 
##               2              53             438              34
```

```{r}
Malignant <- rownames(H3K27.meta)[H3K27.meta$Type == "Malignant"]
Immune <- rownames(H3K27.meta)[H3K27.meta$Type == "Immune cell"]
```
