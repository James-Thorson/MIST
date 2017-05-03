Description
=============
MIST
* Is an R package for implementing a Multispecies Interaction Spatio-Temporal (MIST) model.
* Uses identical input/output and model-diagnostic tools as `VAST` and other packages at www.FishStats.org 

Background
=============
* This tool is designed to estimate pairwise interactions among species using spatially referenced data from multiple years (spatio-temporal data)
* It decomposes variation into (1) covariance among species in a given year, and (2) impacts of species density on per-capita productivity of every other species
* Spatial and spatiotemporal variation are approximated as Gaussian Markov random fields (Thorson Skaug Kristensen Shelton Ward Harms Banante 2014 Ecology), which imply that correlations in spatial variation decay as a function of distance.  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.260143.svg)](https://doi.org/10.5281/zenodo.260143)

Installation Instructions
=============
This function depends on R version >=3.1.1 and a variety of other tools.

First, install the "devtools" package from CRAN

    # Install and load devtools package
    install.packages("devtools")
    library("devtools")

Second, please install the following:
* TMB (Template Model Builder): https://github.com/kaskr/adcomp
* INLA (integrated nested Laplace approximations): http://www.r-inla.org/download

Note: at the moment, TMB and INLA can be installed using the commands 

    # devtools command to get TMB from GitHub
    install_github("kaskr/adcomp/TMB") 
    # source script to get INLA from the web
    source("http://www.math.ntnu.no/inla/givemeINLA.R")  
    
Next, please install the SpatialDFA from this GitHub repository using a function in the "devtools" package:

    # Install package
    install_github("james-thorson/MIST") 
    # Load package
    library(MIST)

Known installation/usage issues
=============
none

Example code
=============
Please see examples folder for single-species and multi-species examples of how to run the model.

Description of package
=============
### Please cite if using the software
* 2017. Thorson, J., Munch, S, and Swain, D.  Estimating partial regulation in community dynamics using spatio-temporal models. Ecology. 98(5): 1277-1289. DOI: 10.1002/ecy.1760. URL: http://onlinelibrary.wiley.com/doi/10.1002/ecy.1760/abstract



