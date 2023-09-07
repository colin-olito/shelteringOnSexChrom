################################################################
#'  RUN SIMULATIONS and CREATE OUTPUT DATA FOR PLOTTING
#'  
#'  PrFix simulations for SLR-Expanding inversions, using 
#'  Recursions syled after Nei
#'		


#####################
##  Dependencies
rm(list=ls())
source('R/simulations-inversions-delMut-RECODE-Nei.R')

######################################
#' Create output directories if they
#' do not already exist
dataDirectoryExists  <-  dir.exists("./data/RECODE")

if(!dataDirectoryExists) {
	dir.create("./data/RECODE")
}

figuresDirectoryExists  <-  dir.exists("./figures/RECODE")

if(!figuresDirectoryExists) {
	dir.create("./figures/RECODE")
}


########################################
# Fixation Probability ~ Inv. Size
########################################
#' Note: these simulations use the EXPANDED
#' W-F model to generate the necessary data
#' to produce Fig. D1


#  Partially deleterious mutations (h = 0.25)
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U02_Nei_TEST")

# N = 10k, broken up by U value
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U02_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U05_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U1_Nei")

# N = 25k, broken up by U value
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U02_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U05_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U1_Nei")

# N = 100k, broken up by U value
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U02_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U05_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U1_Nei")

# N = 250k, broken up by U value
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U02_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U05_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U1_Nei")

# N = 1mil, broken up by U value
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U02_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U05_Nei")
makeDataPrFixInvSize_SLR_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U1_Nei")







##############################################
#' Timeseries of allele & inversion frequency
#' dynamics
#' 
#' Generates very large files
##############################################
rm(list=ls())
source('R/simulations-inversions-delMut-RECODE-Nei.R')

# makeData_SLR_Nei_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
# 							  nTot = 10^4, N = 10^4, Nfname = "_N10k")

# makeData_SLR_Nei_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
# 							  nTot = 10^4, N = 10^5, Nfname = "_N100k")

# makeData_SLR_Nei_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
# 							  nTot = 10^4, N = 5*10^5, Nfname = "_N500k")

# makeData_SLR_Nei_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
# 							  nTot = 10^4, N = 10^6, Nfname = "_N1mil")
