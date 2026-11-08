# Niche_DE

This function performs niche-DE

## Usage

``` r
niche_DE(
  object,
  num_cores,
  outfile,
  C = 150,
  M = 10,
  gamma = 0.8,
  print = T,
  Int = T,
  batch = T,
  self_EN = F,
  G = 1
)
```

## Arguments

- object:

  A niche-DE object

- num_cores:

  The number of cores to parallelize on

- outfile:

  The file path to write progress on

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

- G:

  Number of gigabytes each core should hold. If the counts matrix is
  bigger than G, it is split into multiple chunks such that the size of
  each chunk is less than G gigabytes.

- Self_EN:

  Should niche interactions between the same cell types be considered.
  Default is False.

## Value

A niche-DE object with niche-DE analysis performed
