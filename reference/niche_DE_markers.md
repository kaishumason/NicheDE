# Get Niche-DE marker genes(test)

This function returns genes marker genes in the index cell type when
near the first niche cell type realtive to the second one

## Usage

``` r
niche_DE_markers(object, index, niche1, niche2, alpha = 0.05)
```

## Arguments

- object:

  A niche-DE object

- index:

  The index cell type which we want to find marker genes for

- niche1:

  The niche cell type for the marker genes found

- niche2:

  The niche we wish to compare (index,niche1) patterns to

- alpha:

  The level at which to perform the Benjamini Hochberg correction.
  Default value is 0.05.

## Value

A vector of genes that are niche marker genes for the index cell type
near the niche1 cell type relative to the niche2 cell type
