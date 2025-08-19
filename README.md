# solanum-aa

This repository contains source code associated with the manuscript:

Muir CD, WS Lim, D Wang. Plasticity and adaptation to high light intensity amplify the advantage of amphistomatous leaves. 2025. *In review*.

## Author contributions

* [Chris Muir](https://cdmuir.netlify.app): Conceptualization, Methodology, Investigation, Visualization, Funding acquisition, Writing
* Wei Shen Lim: Methodology, Investigation, Writing – review & editing
* Dachuan Wang: Investigation, Writing – review & editing

## Contents

This repository has the following file folders:

- `data`: data files
- `figures`: figures generated from original artwork and *R* code
- `objects`: saved objects generated from *R* code
- `ms`: manuscript input (e.g. `ms.qmd` and `solanum-aa.bib`) and output (`ms.pdf`) files
- `r`: *R* scripts for all data processing and analysis

## Prerequisites:

To run code and render manuscript:

- [*R*](https://cran.r-project.org/) version >4.5.0 and [*RStudio*](https://www.posit.co/) (recommended)
- [LaTeX](https://www.latex-project.org/): you can install the full version or try [**tinytex**](https://yihui.org/tinytex/)
- [GNU Make](https://www.gnu.org/software/make/): In terminal, you can just type `make paper` to render the manuscript. You can also use it to re-run all scripts.

Before running scripts, you'll need to install the following *R* packages:

```
source("r/install-packages.R")
```

To fit **brms** model, set up [**cmdstanr**](https://mc-stan.org/cmdstanr/).

## Downloading data and code 

1. Download or clone this repository to your machine.

```
git clone git@github.com:cdmuir/solanum-aa.git
```

2. Open `solanum-aa.Rproj` in [RStudio](https://www.posit.co/)

## Rendering manuscript

### Software requirements

At minimum, you will need [R](https://cran.r-project.org/) version 4.5.0 or greater installed on your machine. Install additional packages by running `r/install-packages.R`.

### Rendering manuscript with pre-saved outout

Open `ms/ms.qmd` and knit using [RStudio](https://www.posit.co/).

You can also run the following code from the terminal:

```{terminal}
quarto render ms/ms.qmd
```
