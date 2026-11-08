# CreateLibraryMatrix

This function creates a expression profile matrix for single cell data

## Usage

``` r
CreateLibraryMatrix(data, cell_type)
```

## Arguments

- data:

  Single cell rna-seq counts matrix. Dimension should be cells/spots by
  genes

- cell_type:

  Cell_type assignment matrix. First column is cell names and second
  column is cell type assignment.

## Value

A library matrix with the average expression profile per cell type
