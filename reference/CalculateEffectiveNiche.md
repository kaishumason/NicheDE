# CalculateEffectiveNiche

This function calculates the effective niche of a niche-DE object

## Usage

``` r
CalculateEffectiveNiche(object, cutoff = 0.05)
```

## Arguments

- object:

  A niche-DE object

- cutoff:

  Minimum kernel similarity. Similarities below this value get truncated
  to 0

## Value

A niche-DE object the effective niche calculated. The effective niche is
a list with each entry corresponding to a kernel bandwidth
