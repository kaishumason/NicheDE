# CreateNicheDEObjectFromSeurat

This function creates a niche-DE object from a seurat object

## Usage

``` r
CreateNicheDEObjectFromSeurat(
  seurat_object,
  assay,
  library_mat,
  deconv_mat,
  sigma,
  Int = T
)
```

## Arguments

- seurat_object:

  A spatial seurat object.Coordinate matrix will be extracted via the
  seurat function 'GetTissueCoordinates'. The coordiantes will be scaled
  such that the median nearest neighbor distance is 100.

- assay:

  The assay from which to extract the counts matrix from. The counts
  matrix will be extracted from the counts slot.

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
