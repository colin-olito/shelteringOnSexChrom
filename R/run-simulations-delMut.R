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





# Estimate Time to fixation and Time inversion remains beneficial
makeDataTimeBen_Fix(h = 0.1, s = 0.01, Us.factor = 2, x = 0.2,
					nTot = 10^4, N.vals = c(500, 1000, 2000, 3000, 4000, 5000, 10000), Nfname = "_N10k")

makeDataTimeBen_Fix(h = 0.1, s = 0.01, Us.factor = 2, x = 0.5,
					nTot = 10^4, N.vals = c(500, 1000, 2000, 3000, 4000, 5000, 10000), Nfname = "_N10k")

makeDataTimeBen_Fix(h = 0.1, s = 0.01, Us.factor = 2, x = 0.8,
					nTot = 10^4, N.vals = c(500, 1000, 2000, 3000, 4000, 5000, 10000), Nfname = "_N10k")










# Hybrid sims


#####################
##  Dependencies
rm(list=ls())
source('R/simulations-inversions-delMut.R')

makeDataPrFixInvSizeHybridSim(h = 0.1, s = 0.01, Us.factor.vals = c(2),
					 			nTot = 10^4, N.vals = c(10^2), Nfname="_N100")

makeDataPrFixInvSizeHybridSim(h = 0.1, s = 0.01, Us.factor.vals = c(2, 10),
					 			nTot = 10^4, N.vals = c(10^3), Nfname="_N1k")

makeDataPrFixInvSizeHybridSim(h = 0.1, s = 0.01, Us.factor.vals = c(2, 10),
					 			nTot = 10^4, N.vals = c(4*10^3), Nfname="_N4k")

makeDataPrFixInvSizeHybridSim(h = 0.1, s = 0.01, Us.factor.vals = c(2, 10),
					 			nTot = 10^4, N.vals = c(10^4), Nfname="_N10k")






