# Tutorial

## Niche-DE Tutorial

This page contains tutorials for using the nicheDE package in R. Each
step in the niche-DE process is described in detail; from the initial
deconvolution step to downstream processes such as niche ligand-receptor
analysis. Click the links below to visit a specific step in the
pipeline.

### Installing Niche-DE

- [Installation](https://kmason23.github.io/NicheDE/articles/Installation.md):
  How to download Niche-DE using devtools.

### Deconvolution

- [Deconvolution](https://kmason23.github.io/NicheDE/articles/Deconvolution.md):
  A step by step guide for performing cell type deconvolution of your
  spatial transcriptomic data using RCTD.
- [Average Expression
  Profile](https://kmason23.github.io/NicheDE/articles/Average_Expression_Profile_Calculation.md):
  A guide on how to calculate the average expression profile matrix from
  your reference dataset.

### Creating Your Niche-DE Object

- Object Creation: A guide on how to create your own niche-DE object
  starting from a [spatial seurat
  object](https://kmason23.github.io/NicheDE/articles/Object_Creation_Seurat.md)
  or [raw data
  matrices](https://kmason23.github.io/NicheDE/articles/Object_Creation.md).

- [Effective Niche
  Calculation](https://kmason23.github.io/NicheDE/articles/Effective_Niche_Calculation.md):
  A guide on how to calculate the effective niche as well as suggested
  guidelines on how to choose the kernel bandwidth.

- [Merging and
  Filtering](https://kmason23.github.io/NicheDE/articles/Merging_and_Filtering.md):
  A guide on how to merge multiple niche-DE objects or filter out
  observations from a niche-DE object.

### Performing Niche-DE

- [Performing
  Niche-DE](https://kmason23.github.io/NicheDE/articles/Performing_Niche_DE.md):
  A guide on how to run niche-DE once you have your niche-DE object.

### Downstream Analysis

- [Getting Niche-DE
  Genes](https://kmason23.github.io/NicheDE/articles/Getting_Niche_DE_Genes.md):
  A guide on how to get niche-DE genes between a given index and niche
  cell type as well as how to perform pathway enrichment analysis given
  these niche-DE genes.
- [Marker
  Genes](https://kmason23.github.io/NicheDE/articles/Marker_Genes.md): A
  guide on performing marker gene analysis using niche-DE
- [Ligand Receptor
  Analysis](https://kmason23.github.io/NicheDE/articles/Ligand_Receptor_Analysis.md):
  A guide on performing ligand-receptor analysis using niche-DE.
