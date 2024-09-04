<img src="https://github.com/viktormiok/viktormiok.me/blob/main/software/coni.png" align="right" height="200" width="200">

CoNI
========

## Correlation Guided Network Integration

### Short package description:
`CoNI` is a practical R package for the unsupervised integration of numerical omics datasets. Our tool is based on partial correlations to identify putative confounding variables for a set of paired dependent variables. `CoNI` combines two omics datasets in an integrated, complex hypergraph-like network, represented as a weighted undirected graph, a bipartite graph, or a hypergraph structure. 

<img src="https://github.com/viktormiok/CoNI/blob/master/CONI_abstract.jpeg" align="top" height="480" width="1100">

### Installation
Before installing `CoNI` a few dependencies are necessary:
```r
dependencies<-c("igraph", "doParallel", "cocor", "tidyverse", "foreach","ggrepel", "gplots", "gridExtra", "plyr", "ppcor", "tidyr", "Hmisc")

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

Python3 is also required to run CoNI. Make sure it is installed and it is in your path. 


