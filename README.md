<img src="https://github.com/viktormiok/viktormiok.me/blob/main/software/coni.png" align="right" height="200" width="200">

[![CRAN status](https://www.r-pkg.org/badges/version/ragt2ridges)](https://cran.r-project.org/package=ragt2ridges) ![](https://img.shields.io/badge/languages-R_and_Python-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/ragt2ridges) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/ragt2ridges)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/ragt2ridges) ![GitHub](https://img.shields.io/github/license/viktormiok/ragt2ridges)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/ragt2ridges) 

CoNI
========
- [Overview](#overview)
  * [Application](#application)
- [Installation](#installation)
- [Docker](#docker)
- [Data](#data)
- [Tutorials](#tutorials)
- [License](#license)
- [References](#references)
  
## Correlation Guided Network Integration

### Overview:
`CoNI` is an R package designed for the unsupervised integration of two numerical omics datasets. It leverages partial correlation analysis within a correlation-based framework to reveal how features from one dataset modulate associations between features in the second dataset, facilitating the exploration of potential biological interactions. By calculating partial correlations, CoNI can identify putative confounding variables in one omics layer (the “linker” dataset) that potentially explain or drive the observed correlations among paired features in another layer (the “vertex” dataset). In simpler terms, CoNI finds features in one data layer that could be causing or masking correlations in the other layer, helping researchers tease out hidden associations without any supervised guidance.

After identifying these cross-omics relationships, `CoNI` generates an interpretable network that integrates the two data layers. This integrated network can be represented in three complementary forms for analysis and visualization: as a weighted undirected graph (with features as nodes and confounding links as weighted edges), as a bipartite graph (where nodes are separated into two types corresponding to each omics layer and connections show their interactions), or as a hypergraph (in which a single confounding feature may connect multiple target features simultaneously). These versatile network representations make it easier to explore complex interactions across the two omics datasets, helping users to identify key multi-omics relationships, detect confounders, and generate new hypotheses about underlying biological interactions.

<img src="https://github.com/viktormiok/CoNI/blob/master/CONI_abstract.jpeg" align="top" height="480" width="1100">

#### Application
CoNI is designed to integrate two numerical omics datasets measured on the same samples (e.g., gene expression + metabolomics). Typical examples:

- Transcriptomics (RNA‑seq) + metabolomics
- Proteomics + metabolomics
- Microbiome abundance + metabolite profiles
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

## Docker

If your system configuration makes installing __`CoNI`__ natively difficult, a Docker container is an alternative way to get the __`CoNI`__ package running.

**Note:** Docker Machine has Memory and CPU limits on Mac OS X. To control it, please check instructions either for [CLI](https://stackoverflow.com/questions/32834082/how-to-increase-docker-machine-memory-mac/32834453#32834453) or for [Docker Desktop](https://docs.docker.com/docker-for-mac/#advanced).

For building a Docker image from the Dockerfile, download the [Dockerfile](https://github.com/viktormiok/CoNI/blob/main/Dockerf) (available in this repo) and run the following command to create it:
```
docker build -t CoNI.
```
This will create a __`CoNI`__ Docker image on your system (Please be patient, as the build may take approximately 30–50 minutes to complete).
You can then run it using the following command:
```
docker run -d -p 8787:8787 -e PASSWORD=pass --name tigaR -it CoNI
```
## License

__`CoNI`__ is distributed under the MIT license. Please read the license before using __`CoNI`__, distributed in the `LICENSE` file.

## References

Publications related to __`CoNI`__ include:

 - Monroy Kuhn, J.M., **Miok, V.**, Lutter, D. (2022),
   "[Correlation guided Network Integration (CoNI), an R package for integrating numerical omics data that allows multiform graph representations to study molecular interaction networks](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac042/6603108)". 
    *Bioinformatics Advances*, 2(1): 0-42.
   
Please cite the relevant publications if you use __`CoNI`__.
