# Getniche genes for the given index and niche cell types at the desired test.level

This function returns genes that show niche patterns at the desired
test.level

## Usage

``` r
get_niche_DE_genes(
  object,
  test.level,
  index,
  niche,
  positive = T,
  alpha = 0.05
)
```

## Arguments

- object:

  A niche-DE object

- test.level:

  At which test.level to return genes (gene level ("G"), cell type
  level("CT"), interaction level("I"))

- index:

  The index cell type

- niche:

  The niche cell type

- positive:

  Logica indicating whether to return upregulated/positive niche-DE
  genes(True) or downregulated/negative niche-DE genes(False). Default
  value is True.

- alpha:

  The level at which to perform the Benjamini Hochberg correction.
  Default value = 0.05

## Value

A vector of genes that are niche significant at the desired FDR,
test.level, index cell type, and niche cell type
