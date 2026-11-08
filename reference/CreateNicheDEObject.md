# CreateNicheDEObject

This function creates a niche-DE object

## Usage

``` r
CreateNicheDEObject(
  counts_mat,
  coordinate_mat,
  library_mat,
  deconv_mat,
  sigma,
  Int = T
)
```

## Arguments

- counts_mat:

  Counts matrix. Dimension should be cells/spots by genes

- coordinate_mat:

  Coordinate matrix. It will be scaled such that the median nearest
  neighbor distance is 100.

- library_mat:

  Matrix indicating average expression profile for each cell type in the
  sample

- deconv_mat:

  Deconvolution or cell type assignment matrix of data

- sigma:

  List of kernel bandwidths to use in calculating the effective niche

- Int:

  Boolean of if counts data supplied is integer. Default is true. When
  performing niche-DE, Negative binomial regression is performed if
  True. Linear regression with a gene specific variance is performed if
  False.

## Value

A niche-DE object
