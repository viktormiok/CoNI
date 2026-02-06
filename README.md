<img src="https://github.com/viktormiok/viktormiok.me/blob/main/software/coni.png" align="right" height="200" width="200">

[![CRAN status](https://www.r-pkg.org/badges/version/ragt2ridges)](https://cran.r-project.org/package=ragt2ridges) ![](https://img.shields.io/badge/languages-R_and_Python-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/ragt2ridges) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/ragt2ridges)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/ragt2ridges) ![GitHub](https://img.shields.io/github/license/viktormiok/ragt2ridges)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/ragt2ridges) 

CoNI
========

## Correlation Guided Network Integration

### Short package description:
`CoNI` is a practical and easy-to-use R package for the unsupervised integration of numerical omics datasets. Our tool is based on partial correlations to identify putative confounding variables for a set of paired dependent variables. `CoNI` combines two omics datasets in an integrated, complex hypergraph-like network, represented as a weighted undirected graph, a bipartite graph, or a hypergraph structure. 

<img src="https://github.com/viktormiok/CoNI/blob/master/CONI_abstract.jpeg" align="top" height="480" width="1100">

### Installation
Before installing `CoNI`, a few dependencies are necessary:
```r
dependencies<-c("igraph", "doParallel", "cocor", "tidyverse", "foreach", "ggrepel", "gplots", "gridExtra", "plyr", "ppcor", "tidyr", "Hmisc")

`%notin%`<-Negate(`%in%`)
for(package in dependencies){
  if(package %notin% rownames(installed.packages())){
    install.packages(package,dependencies = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("genefilter")
```

CoNI requires Python 3. Ensure that Python 3 is installed and correctly added to your system PATH.

## References

Publications related to __`CoNI`__ include:

 - Monroy Kuhn, J.M., **Miok, V.**, Lutter, D. (2022),
   "[Correlation guided Network Integration (CoNI), an R package for integrating numerical omics data that allows multiform graph representations to study molecular interaction networks](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac042/6603108)". 
    *Bioinformatics Advances*, 2(1): 0-42.
   
Please cite the relevant publications if you use __`CoNI`__.
