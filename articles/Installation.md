# Installation

## Installation

You can install NicheDE with the following code

    Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
    install.packages("devtools")
    options(timeout=9999999)
    devtools::install_github("kaishumason/NicheDE") # install
    library(nicheDE)
