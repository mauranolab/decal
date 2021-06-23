# DECAL

**DECAL** (**D**ifferential **E**xpression analysis of **C**lonal
**A**lterations **L**ocal effects) provide you tools to conduct differential
expression analysis of single-cell perturbations to potential interacting
genes.
Similar to other expression analysis tools, it models gene expression using
a _Negative Binomial_ (or _Gamma-Poisson_) regression, modeling each gene
UMI count by the cell total count and the cell alteration status.

## Features

- `decal` is compatible with **tidyverse** analysis.
- Includes a suit of simulation functions to generate your own dataset based
  `decal` model and evaluate your statistical power.
- Package has few dependencies, requiring only `MASS`, `fastglm` and `Matrix`.
- Can evaluate specific clonal alteration and gene effect, instead of
  modeling all genes and all alterations. This allow us to quickly investigate
  a large number of interactions skipping unlikely effects.

## Instalation

To install `decal` package current version, open your R terminal and type:

```{r installation, eval = FALSE}
## Install remotes if not available with
## install.install.packages("remotes")
remotes::install_github("mauranolab/decal")
```

## Quick Start

`decal` is the package main function which estimates each gene dispersion
($\theta$) parameter and fit a _Negative Binomial_ regression to each gene
and alteration pair specified to evaluate the perturbation statistical
significance.
It requires three parameters a UMI count matrix, a list specifying your clones
cell (count matrix columns) composition, and a table indicating the potential
alteration and affected gene to evaluate.

Below we show a quick example of `decal` analysis using a simulated dataset
included with our package.

```{r quick_start}
library(decal)

data("sim_decal")
perturbations <- sim_decal$perturbations
count <- sim_decal$count
clone <- sim_decal$clone
```

Our dataset (`sim_decal`) is composed by the three required parameters:

- a UMI count matrix (`count`);
- a list dividing our cells by clonal population (`clone`);
- and a table (`perturbations`) indicating potential interactions to evaluate,
  specified by the perturbed `clone` and candidate `gene`.

With these inputs, we can run the full analysis using the `decal` function as
follow:

```{r quick_start_run}
res <- decal(perturbations, count, clone)
```

`decal` function conducts the whole differential analysis, by estimating each
gene dispersion and fitting a negative binomial regression to evaluate each
interaction statistical significance.
Finally, it updates `perturbations` table to include the following columns:

- `n0` and `n1`: number of non-perturbed and perturbed cells, respectively.
- `x0` and `x1`: average UMI count among perturbed and non-perturbed cells,
  respectively.
- `mu`: average UMI count among all cells.
- `xb`: expected average UMI count of perturbed cells.
- `theta`: estimated _negative binomial_ dispersion parameter.
- `z`: estimated perturbation z-score.
- `lfc`: _log2 fold-change_ of perturbed cells gene expression.
- `pvalue` and `p_adjusted`: perturbation _t-test_ significance values.

## Acknowledgement

- The negative binomial dispersion estimation process was extracted and modified
  from [`sctransform`](https://github.com/ChristophH/sctransform) package and
  described by _Hafemeister & Satija (2019)_