# NicheDE
NicheDE is a method that detects context dependent changes in gene expression in spatial transcriptomic data.
In particular, given an index and a niche cell type, we define an (index,niche) niche gene as a gene as a gene that is up or down-regulated in 
the index cell type when in the presence of the niche cell type.\
NicheDE identifies (index,niche)+ genes by first characterizing a cell's neighborhood into a vector called the effective niche. 
Then, by regressing gene expression on the effective niche,testing if certain coefficients are equal to 0 is equivalent to 
determining if the gene is an (index,niche) niche gene. \
For more information about this method please check out [the manuscript](https://www.biorxiv.org/content/10.1101/2023.01.03.522646v1?rss=1)

# Installation
You can install NicheDE with the following code
```
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install.packages("devtools")
options(timeout=9999999)
devtools::install_github("Kmason23/NicheDE") # install
library(NicheDE)
```
# Vignette
Click the link below for a tutorial on using nicheDE\
[NicheDE applied to Visium data of liver metastasized colorectal carcinoma](https://github.com/Kmason23/NicheDE/blob/master/vignettes/Niche_DE_introduction.Rmd)

# citation
Kaishu Mason, Anuja Sathe, Paul Hess, Jiazhen Rong, Chi-Yun Wu, Emma Furth, Hanlee P. Ji, Nancy Zhang Niche differential gene expression analysis in spatial transcriptomics data identifies context-dependent cell-cell interactions (Biorxiv 2023:  https://doi.org/10.1101/2023.01.03.522646
