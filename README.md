# Metacognitive Computations for Information Search -- Code

## Overview

Welcome! This repository contains the code to reproduce the results from the paper **"Metacognitive Computations for Information Search: Confidence in Control"**, authored by Lion Schulz, Steve Fleming and Peter Dayan. Don't hesitate to get in touch if you have any questions or comments (corresponding author: Lion Schulz).

## Code organisation

The code in this repository is organized as follows:

* Directly in the directory: Rmarkdowns reproducing the plots in the computations and results section of the paper, and the appendix.
  * `MetaSearch_01_Introduction.Rmd` produces the figures from the introduction
  * `MetaSearch_02_Results.Rmd` produces most theoretic results and the accompanying figures
  * `MetaSearch_03_Results_Data.Rmd` produces additional theoretic results and analyzes the data
  * `MetaSearch_04_Appendix.Rmd` produces the appendix results
* In `backend/`:
  * Several `global_[...].R` files which load the necessary libraries, custom functions and settings for plotting
    * `global_all.R` calls all of these and just serves as a shortcut for the Rmarkdowns
    * `global_functions.R` sources or specifies custom functions
    * `global_libraries.R` loads all the necessary R libraries for the project
    * `global_setttings.R` defines settings for the plotting in the main Rmarkdowns. This includes colour schemes, line sizes, etc. It also swithces on anti-aliasing in Windows.
  * `backend/data_analysis`: Functions for the data analysis.
  * `backend/models/`: Model code, divided by model class. Functions used by all models are in the `joint` folder.
  * `backend/plotting/`: Several `.R` files that specify code for specific figures in the Rmarkdowns.
  * `backend/hmetad` contains the code for fitting meta-d'. This is an adapted version from https://github.com/metacoglab/HMeta-d
* `output/` is the folder into which data and plots are saved. This contains  high resolution model runs.
* `data/` contains data from Schulz et al. (2020, PNAS) both in unprocessed form as well as with analyses applied. It also contains simulations.

The repository also includes files necessary for the `renv` environment management package (see more below). These include the `renv/` folder, and the `.Rprofile`, `renv.lock` files


## Models

The code for the two main models is in `backend/models/`. This code is written to work for the results in the paper. As such, it returns the relevant summary  statistic given the relevant parameters. Thus, we have two main files for the two classes of models:

* `model_Pd.R` for the postdecisional model
* `model_2o.R` for the vanilla second-order model
  * `model_2o_subjective.R` adapts the second-order model to allow for over- and underconfidence 

Both functions rely on `model_functions_joint.R`. The second-order model also uses `model_2o_functions.R` which contains a ported version of the original second-order confidence.


## Reproducing the figures

#### (1) Download Github directory

First, download the entire Github repository and R project and not just a single file. There are dependencies between folders and the environment within the project which is itself in a managed environment.

#### (2) Open the R project

This project uses the `renv` library to control the environment and to make the code more reproducible. To use this, you first need to launch the R project `MetaSearch_paper.Rproj` using
[RStudio](https://rstudio.com/products/rstudio/). Then, `renv` should automatically bootstrap itself, thereby downloading and installing the appropriate version of `renv` into the project library.

After this has completed, you then use `renv::restore()` to restore the libraries used in the project locally on your machine. You can find more information on `renv` [here](https://rstudio.github.io/renv/articles/renv.html) and on the collaboration functionality used [here](https://rstudio.github.io/renv/articles/collaborating.html). Note that `renv` takes care of the libraries (and their versions), but not of the specific version of R used. You will get a pointer upon calling `renv::restore()` on the R version used here (4.1.0.). For the best possible compatibility, we advise using this specific version of R, although others should also work.

#### (3) Running the Rmarkdowns

The main code for the figures in the paper is organized in [RMarkdowns](https://rmarkdown.rstudio.com/), where R code blocks ("chunks") are embedded within a file written in markdown. Output produced by these chunks is produced inline in the markdown.

