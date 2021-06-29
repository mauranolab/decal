# DECAL <img src="hex.png" align="right" width="150px"/>

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

## Installation

To install `decal` package current version, open your R terminal and type:

```{r installation, eval = FALSE}
## Install remotes if not available with
## install.install.packages("remotes")
remotes::install_github("mauranolab/decal")
```

## Using `decal`

The package main function (`decal`) runs the whole statistical analysis by
fitting a _Negative Binomial_ (or _Gamma-Poisson_) regression to each gene
and alteration pair specified and evaluate the perturbation statistical
significance for a single-cell perturbation experiment. You can run `decal`
as follow:

```{r}
decal(perturbations, count, clone)
```

The three parameters required are the following:

- a table specifying the population and gene to evaluate perturbation
  (`perturbations`).
- a UMI count matrix that each column correspond to a cell and each row
  a gene or feature (`count`).
- a list of cells specifying your clones composition (`clone`). This is
  represented by a R list composed by character (or integer) vectors named
  after the clone ID and indicating the cells belonging to this clone.

`decal` returns a deep copy of the `perturbations` table with the following
additional columns:

- `n0` and `n1`: number of non-perturbed and perturbed cells, respectively.
- `x0` and `x1`: average UMI count among perturbed and non-perturbed cells,
  respectively.
- `mu`: average UMI count among all cells.
- `xb`: expected average UMI count of perturbed cells.
- `theta`: estimated _negative binomial_ dispersion parameter.
- `z`: estimated perturbation z-score.
- `lfc`: _log2 fold-change_ of perturbed cells gene expression.
- `pvalue` and `p_adjusted`: perturbation _t-test_ significance values.


## Quick Start

Below we include a example of `decal` analysis on a simulated dataset included
with our package.

```{r quick_start}
library(decal)

data("sim_decal")
perturbations <- sim_decal$perturbations
count <- sim_decal$count
clone <- sim_decal$clone

decal(perturbations, count, clone)
```

Further details exploring this example are included in our _quick-start_
vignette.

## Acknowledgement

- The negative binomial dispersion estimation process was extracted and modified
  from [`sctransform`](https://github.com/ChristophH/sctransform) package and
  described by _Hafemeister & Satija (2019)_