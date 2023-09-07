################################################################
#'  RUN SIMULATIONS and CREATE OUTPUT DATA FOR PLOTTING
#'  
#'  PrFix simulaitons for Autosomal inversions, using EXACT
#'  Recursions
#'		


#####################
##  Dependencies
rm(list=ls())
source('R/Auto-inversions-NeiMyVersion.R')

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
#' Note: these simulations create data
#' to produce Fig. 2 showing the fixation 
#' probability of different sized neutral 
#' inversions expanding the SLR on Y 
#' chromosomes under deleterious mutation 
#' pressure. Uses multilocus recursions




##  Partially recessive deleterious mutations (h = 0.25)

# N = 10k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U1_Nei")

# N = 25k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U1_Nei")

# N = 100k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U1_Nei")

# N = 250k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U1_Nei")

# N = 1mil, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.25, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U1_Nei")





##  Strongly recessive deleterious mutations (h = 0.25)

# N = 10k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^4, Nfname="_N10k_U1_Nei")

# N = 25k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 2.5*10^4, Nfname="_N25k_U1_Nei")

# N = 100k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^5, Nfname="_N100k_U1_Nei")

# N = 250k, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 2.5*10^5, Nfname="_N250k_U1_Nei")

# N = 1mil, broken up by U value
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.02),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U02_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.05),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U05_Nei")
makeDataExactAutoPrFixInvSize_Nei(h = 0.1, s = 0.01, U.vals = c(0.1),
						nTot = 10^4, N = 10^6, Nfname="_N1mil_U1_Nei")







##############################################
#' Timeseries of allele & inversion frequency
#' dynamics
#' 
#' Generate very larege files.
##############################################
rm(list=ls())
source('R/Auto-inversions-NeiMyVersion.R')

# makeData_AutoExact_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
# 							  nTot = 10^4, N = 10^5, Nfname = "_N100k")

# makeData_AutoExact_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
# 							  nTot = 10^4, N = 5*10^5, Nfname = "_N500k")

# makeData_AutoExact_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
# 							  nTot = 10^4, N = 10^6, Nfname = "_N1mil")
