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


#  Partially deleterious mutations (h = 0.25)

# N = 10k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us10_deterministic_q")

# N = 100k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us10_deterministic_q")

# N = 500k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us10_deterministic_q")

# N = 1mil, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us10_deterministic_q")






# h = 0.1
# N = 10k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us10_deterministic_q")

# N = 100k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us10_deterministic_q")

# N = 500k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us10_deterministic_q")

# N = 1mil, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(2),
					 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us2_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(5),
					 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us5_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, Us.factor.vals = c(10),
					 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us10_deterministic_q")





#  More strongly recessive deleterious mutations
#  Interestingly, relation between P(fix) ~ x can change for lower levels of average
#  load (smaller U/s)
makeDataPrFixInvSize(h = 0.01, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^3), Nfname="_N1k_deterministic_q")

makeDataPrFixInvSize(h = 0.01, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_deterministic_q")







# Autosomal inversions - Pr(fix | x)

# h = 0.25, N = 10k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us2")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us5")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_Us10")
# h = 0.25, N = 100k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us2")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us5")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_Us10")
# h = 0.25, N = 500k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
						 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us2")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
						 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us5")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
						 nTot = 10^4, N.vals = c(5*10^5), Nfname="_N500k_Us10")

# h = 0.25, N = 1 million, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us2")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(5),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us5")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(10),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_Us10")

