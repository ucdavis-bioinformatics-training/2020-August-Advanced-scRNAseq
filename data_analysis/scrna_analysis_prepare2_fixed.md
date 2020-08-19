### Create a new RStudio project, PART 2

Using the same RStudtio project you've used so far.

In the R console run the following commands
```r
if (!any(rownames(installed.packages()) == "devtools")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("devtools")
}
library(devtools)

## install SeuratWrappers
devtools::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

BiocManager::install(c("Biobase", "SingleCellExperiment", "batchelor", "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SummarizedExperiment", "pcaMethods"))

## Possible you may have to install Rtools if the below fails
## http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/#1

## Mac users may also experience installation problems due to Xcode or gfortran.
## Instructions: https://cole-trapnell-lab.github.io/monocle3/docs/installation/
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)

devtools::install_github("velocyto-team/velocyto.R")
library(velocyto.R)

sessionInfo()
```

### Download the template Markdown workshop document for Monocle and open it.

In the R console run the following commands
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-August-Advanced-scRNAseq/master/data_analysis/adv_scrnaseq_monocle.Rmd", "monocle.Rmd")
```

Additionally, download the following data files for monocle analysis.

```r
download.file("https://github.com/ucdavis-bioinformatics-training/2020-August-Advanced-scRNAseq/raw/master/datasets/monocle3_expression_matrix.rds", "monocle3_expression_matrix.rds")
download.file("https://github.com/ucdavis-bioinformatics-training/2020-August-Advanced-scRNAseq/raw/master/datasets/monocle3_cell_metadata.rds", "monocle3_cell_metadata.rds")
download.file("https://github.com/ucdavis-bioinformatics-training/2020-August-Advanced-scRNAseq/raw/master/datasets/monocle3_gene_metadata.rds", "monocle3_gene_metadata.rds")
```

### Download the template Markdown workshop document for VDJ and open it.

In the R console run the following command
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-August-Advanced-scRNAseq/master/data_analysis/VDJ_Analysis.Rmd", "VDJ_Analysis.Rmd")
```

### Download the template Markdown workshop document for Velocity analysis and open it.

In the R console run the following command
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-August-Advanced-scRNAseq/master/data_analysis/Velocyto.Rmd", "Velocyto.Rmd")
```
