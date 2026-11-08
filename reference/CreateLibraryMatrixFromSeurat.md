# CreateLibraryMatrixFromSeurat

This function creates a expression profile matrix for single cell data
from a Seurat object

## Usage

``` r
CreateLibraryMatrixFromSeurat(seurat_object, assay)
```

## Arguments

- seurat_object:

  A seurat object

- assay:

  The assay from which we want to extract the counts matrix to calculate
  the average expression profile.

## Value

A library matrix with the average expression profile per cell type
