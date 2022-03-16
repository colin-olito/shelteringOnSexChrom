################################################################
#  RUN SIMULATIONS and CREATE OUTPUT DATA FOR PLOTTING
#  
#  R code for forward simulations. Generates output data
#  as .csv files saved to ./output/
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		


#####################
##  Dependencies
rm(list=ls())
source('R/simulations-inversions-delMut.R')

######################################
#' Create output directories if they
#' do not already exist
dataDirectoryExists  <-  dir.exists("./data")

if(!dataDirectoryExists) {
	dir.create("./data")
}

figuresDirectoryExists  <-  dir.exists("./figures")

if(!figuresDirectoryExists) {
	dir.create("./figures")
}

########################################
# Fixation Probability ~ Inv. Size
########################################
#' Note: these simulations create data
#' to produce Fig. 2 showing the fixation 
#' probability of different sized neutral 
#' inversions expanding the SLR on Y 
#' chromosomes under deleterious mutation 
#' pressure. Uses multilocus recursions

makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^2), Nfname="_N100_deterministic_q")

makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^3), Nfname="_N1k_deterministic_q")

makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_deterministic_q")

makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_deterministic_q")




makeDataPrFixInvSizeDetqI(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^3), Nfname="_N1k_det_qI")

makeDataPrFixInvSizeDetqI(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_det_qI")

makeDataPrFixInvSizeDetqI(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_det_qI")




#  More strongly recessive deleterious mutations
#  Interestingly, relation between P(fix) ~ x can change for lower levels of average
#  load (smaller U/s)
makeDataPrFixInvSize(h = 0.01, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^3), Nfname="_N1k_deterministic_q")

makeDataPrFixInvSize(h = 0.01, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_deterministic_q")


#  Less-recessive deleterious mutations (h = 0.25)
#  Interestingly, relation between P(fix) ~ x can change for lower levels of average
#  load (smaller U/s)
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^3), Nfname="_N1k_deterministic_q")

makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_deterministic_q")


# Autosomal inversions - Pr(fix | x)
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
						 nTot = 10^4, N.vals = c(10^3, 10^4))

makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
						 nTot = 10^4, N.vals = c(10^3, 10^4))

makeDataAutoPrFixInvSize(h = 0.5, s = 0.01, Us.factor.vals = c(2, 5, 10),
						 nTot = 10^4, N.vals = c(10^2, 10^3, 10^4))

