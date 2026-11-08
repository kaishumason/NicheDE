# Niche_DE

This function performs niche-DE but not in parallel. This would be
useful when the dataset is too large to run across multiple clusters as
currently the sparse daa matrix must be converted to a dense matrix. In
general, users should default to using the standard niche_DE function.

## Usage

``` r
niche_DE_no_parallel(
  object,
  C = 150,
  M = 10,
  gamma = 0.8,
  print = T,
  Int = T,
  batch = T,
  self_EN = F
)
```

## Arguments

- object:

  A niche-DE object

- C:

  Minimum total expression of a gene needed for the model to run.
  Default value is 150.

- M:

  Minimum number of spots containing the index cell type with the niche
  cell type in its effective niche for (index,niche) niche patterns to
  be investigated. Default value is 10.

- gamma:

  Percentile a gene needs to be with respect to expression in the index
  cell type in order for the model to investigate niche patterns for
  that gene in the index cell. Default value is 0.8 (80th percentile)

- print:

  Logical if function should print progress report (kernel, gene \#)

- Int:

  Logical if data is count data. If True, negative binomial regression
  is performed. Otherwise linear regression with a gene specific
  variance is applied. Default is True.

- batch:

  Logical if an indicator should be included for each batch in the
  regression. Default is True.

- Self_EN:

  Should niche interactions between the same cell types be considered.
  Default is False.

## Value

A niche-DE object with niche-DE analysis performed
