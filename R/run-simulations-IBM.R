################################################################
#  RUN INDIVIDUAL-BASED SIMULATIONS and CREATE OUTPUT DATA 
#  
#  R code for INDIVIDUAL-BASED forward simulations. 
#  Generates output data as .csv files saved to ./output/
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		


#####################
##  Dependencies
rm(list=ls())
source('R/functions-IBM-Simulations.R')


########################################
# Fixation Probability ~ Inv. Size
########################################
#' Note: This code runs individual based
#' simulations to create data that was
#' intended to reproduce essential parts
#' of Fig. 2, especially the results for 
#' small inversions. After spending a lot
#' of time, it has become obvious that the
#' necessary computation times to estimate 
#' the fixation probabilities for reasonable
#' parameter values is prohibitively long.
#' I have retained the code here in case it
#' is useful for future studies or other 
#' researchers. Advice to speed things up 
#' is most welcome.
#' 

#' IMPORTANT: 
#' In order to get reasonably 'equilibrium-like' behaviour
#' you need to use population sizes and dominance coefficients
#' that result in relatively strong selection relative to drift
#' FOR THE DELETERIOUS ALLELES AT EACH LOCUS. If not, the 
#' allele frequency distribution becomes too 'U' shaped, with 
#' frequencies of 0 for most loci, and a few with v.high frequencies.
#' 
#' After extensive fiddling, it seems a reasonable combination of 
#' parameter values is N = 10,000 and h = 0.1, with s = 0.01...
#' N = 5000 could maybe be an absolute lower limit....


numberOfCluster  <-  2*(detectCores() - 1)

# test using N = 1000
# makeDataPrFixIBMSimParallel(h = 0.25, s = 0.01, Us.factor = 2,
# 					nTot = 10^4, N = 10^3, invSize = c(0.05), burnin=100,
# 					nCluster = numberOfCluster, extraFileID = "_x0.05_freeR_parallel")


## N = 10000
# x = 0.05
makeDataPrFixIBMSimParallel(h = 0.1, s = 0.01, Us.factor = 5,
					nTot = 10^4, N = 10^4, invSize = 0.05,  burnin=10000,
					nCluster = numberOfCluster, extraFileID = "_x0.05_freeR_parallel")
# x = 0.1
makeDataPrFixIBMSimParallel(h = 0.1, s = 0.01, Us.factor = 5,
					nTot = 10^4, N = 10^4, invSize = 0.1,  burnin=10000,
					nCluster = numberOfCluster, extraFileID = "_x0.1_freeR_parallel")
# x = 0.2
makeDataPrFixIBMSimParallel(h = 0.1, s = 0.01, Us.factor = 5,
					nTot = 10^4, N = 10^4, invSize = 0.2,  burnin=10000,
					nCluster = numberOfCluster, extraFileID = "_x0.2_freeR_parallel")
