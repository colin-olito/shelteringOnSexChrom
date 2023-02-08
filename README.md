# Consequences of partially recessive deleterious genetic variation for the evolution of inversions suppressing recombination between sex chromosomes.

## Overview

This is a GitHub repository for the development of a theoretical population genetics research project that is now published under the title "*Consequences of partially recessive deleterious genetic variation for the evolution of inversions suppressing recombination between sex chromosomes.*" (doi: 10.1111/evo.14496).
 Here you can find all of the necessary code to reproduce the simulations and main figures presented in the published paper and appendices. A version of record for this repository is archived on Zenodo. [![DOI](https://zenodo.org/badge/449208001.svg)](https://zenodo.org/badge/latestdoi/449208001).


## Abstract

The evolution of suppressed recombination between sex chromosomes is widely hypothesized to be driven by sexually antagonistic selection (SA), where tighter linkage between the sex-determining gene(s) and nearby SA loci is favoured when it couples male-beneficial alleles to the proto-Y chromosome, and female-beneficial alleles to the proto-X. Despite limited empirical evidence, the SA selection hypothesis overshadows several alternatives, including an incomplete but often-repeated "sheltering hypothesis" which suggests that expansion of the sex-linked region (SLR) reduces homozygous expression of partially recessive deleterious mutations at selected loci. Here, we use population genetic models to evaluate the consequences of deleterious mutational variation for the evolution of neutral chromosomal inversions expanding the SLR on proto-Y chromosomes. We find that SLR-expanding inversions face a race against time: lightly-loaded inversions are initially beneficial, but eventually become deleterious as they accumulate new mutations, and must fix before this window of opportunity closes. The outcome of this race is strongly influenced by inversion size, the mutation rate, and the dominance coefficient of deleterious mutations. Yet, small inversions have elevated fixation probabilities relative to neutral expectations for biologically plausible parameter values. Our results demonstrate that deleterious genetic variation can plausibly drive recombination suppression in small steps and would be most consistent with empirical patterns of small evolutionary strata or gradual recombination arrest.



## Citing information
*Please cite the paper as*:

Olito, C., B. Hansson, S. Ponnikas, J.K. Abbott. 2022. Consequences of partially recessive deleterious genetic variation for the evolution of inversions suppressing recombination between sex chromosomes. *Evolution* 76: 1320--1330. doi: 10.1111/evo.14496.

The published article, along with all supplementary material, is Open Access and [freely available through the publisher](https://doi.org/10.1111/evo.14496). You can also contact me directly if you would like a copy. 

A version of record for this repository at the time of acceptance is archived on Zenodo [![DOI](https://zenodo.org/badge/449208001.svg)](https://zenodo.org/badge/latestdoi/449208001).

##  Instructions

This repository provides all code necessary to (1) rerun the simulations and (2) produce figures as .pdf's. To do so, please follow these basic steps:

1. Clone the repo using the following: `git clone https://https://github.com/colin-olito/shelteringOnSexChrom`. Alternatively, on the project main page on GitHub, click on the green button `clone` or `download` and then click on `Download ZIP`.  
2. Check that you have a recent version of [`R`](https://www.r-project.org/) installed. 
3. Make sure that the working directory for your R session is the root directory of this repo (e.g., `shelteringOnSexChrom-master/`).
4. Run `./R/run-simulations-delMut.R` either interactively in R or in terminal.
5. Run `./R/run-simulations-PartialFullSib.R` either interactively in R or in terminal.
6. *Note*: We use CM fonts in the figures. To do this, be sure to correctly install the `R` font packages `extrafont` and `fontcm`. Alternatively, comment out L.4-5 in `./R/functions-figures`, and change the default font family to 'Arial' by swapping L.24 & L.25.
7. Run `makeFigs.R` (up to L.140), which will read the simulation output files and generate the main figures in the paper and supplementary material.  



## Repository structure and contents 

The directories/files in this repository needed to reproduce the results for this study are as follows:  

- **`R`**   
	- `functions-figures-inversion-delMut.R`  
	- `functions-figures.R`  
	- `run-simulations-3Loc-PartialFullSib.R`  
	- `run-simulations-PartialFullSib.R`  
	- `run-simulations-delMut.R`  
	- `simulations-3Loc-PartialFullSib.R`  
	- `simulations-PartialFullSib.R`  
	- `simulations-inversions-delMut.R`  
- **`data`***  
- **`figures`***  
- `makeFigs.R`  
- `LICENSE.txt`   

**Note:** * `Data` and `figures` directories will be created locally the first time `run-Simulations-delMut.R` is run (*if needed*).


### File & variable descriptions

Plotting function files
- `functions-figures-inversion-delMut.R`: plotting functions for figures specific to simulations for inversions with deleterious mutations.  
- `functions-figures.R`: general plotting functions.   

Simulation function files
- `simulations-3Loc-PartialFullSib.R`: workhorse simulation functions for 3-locus version of Charlesworth & Wall (1999) model.   
- `simulations-PartialFullSib.R`: workhorse simulation functions for 2-locus version of Charlesworth & Wall (1999) model.   
- `simulations-inversions-delMut.R`: workhorse simulation functions deterministic and W-F simulations of inversions under deleterious mutation pressure.   

Executables
- `run-simulations-3Loc-PartialFullSib.R`: executable functions for 3-locus version of Charlesworth & Wall (1999) model.   
- `run-simulations-PartialFullSib.R`: executable functions for 2-locus version of Charlesworth & Wall (1999) model.   
- `run-simulations-delMut.R`: executable functions for deterministic and W-F simulations of inversions under deleterious mutation pressure.   
- `makeFigs.R`: executable functions to create .pdf figures using simulation output files.

License    
- `LICENSE.txt`: MIT license for this repository.  


## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the inversionSize github [issues page](https://github.com/colin-olito/shelteringOnSexChrom/issues). If you would like to report a bug/issue, and do not have a github account (and don't want to get one), please send a brief email detailing the problem you encountered to colin.olito at biol dot lu dot se.



## Licence information

This repository is provided by the authors under the MIT License ([MIT](https://opensource.org/licenses/MIT)).