# CalculateEffectiveNiche

This function calculates the effective niche of a niche-DE object

## Usage

``` r
CalculateEffectiveNicheLargeScale(
  object,
  batch_size = 1000,
  cutoff = 0.05,
  standardize = TRUE
)
```

## Arguments

- object:

  A niche-DE object

- batch_size:

  Number of cells to calculate the effective niche for at a time

- cutoff:

  Minimum kernel similarity. Similarities below this value get truncated
  to 0

- standardize:

  A logical indicator for if the effective niche should be z-normalized
  by column

## Value

A niche-DE object the effective niche calculated. The effective niche is
a list with each entry corresponding to a kernel bandwidth
