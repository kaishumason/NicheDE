# NicheDE
NicheDE is a method that detects context dependent changes in gene expression in spatial transcriptomic data.
In particular, nicheDE characterizes a cell's neighborhood into a vector called the effective niche. 
By regressing gene expression on the effective niche, nicheDE can determine if a gene is up or down regulated
in a specified cell type in the context of specific spatial niches.\
For more information about this method please check out [the manuscript](https://www.biorxiv.org/content/10.1101/2023.01.03.522646v1?rss=1)

# Installation
You can install NicheDE with the following code
```
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install.packages("devtools")
options(timeout=9999999)
devtools::install_github("Kmason23/NicheDE") # install
library(Clonalscope)
```
# Vignette
Click the link below for a tutorial on using nicheDE
[NicheDE applied to Visium data of liver metastasized colorectal carcinoma](https://www.biorxiv.org/content/10.1101/2023.01.03.522646v1?rss=1)
