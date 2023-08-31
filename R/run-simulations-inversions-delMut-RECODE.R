################################################################
#'  RUN SIMULATIONS and CREATE OUTPUT DATA FOR PLOTTING
#'  
#'  PrFix simulations after RECODE to troubleshoot odd results
#'		


#####################
##  Dependencies
rm(list=ls())
source('R/simulations-inversions-delMut-RECODE.R')

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


#  Partially deleterious mutations (h = 0.25)

# N = 10k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U02_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U05_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U1_deterministic_q")

# N = 100k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U02_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U05_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U1_deterministic_q")

# N = 500k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U02_deterministic_qHat")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U05_deterministic_qHat")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U1_deterministic_qHat")

# N = 1mil, broken up by Us.factor
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U02_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U05_deterministic_q")
makeDataPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U1_deterministic_q")






# h = 0.1
# N = 10k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U02_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U05_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U1_deterministic_q")

# N = 100k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U02_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U05_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U1_deterministic_q")

# N = 500k, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U02_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U05_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U1_deterministic_q")

# N = 1mil, broken up by Us.factor
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U02_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U05_deterministic_q")
makeDataPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U1_deterministic_q")




##########################################
#' Note: these simulations estimate fixation 
#' probabilities for different sized neutral 
#' inversions expanding the SLR on Y 
#' chromosomes under deleterious mutation 
#' pressure. Uses deterministic recursions
#' with assumption that deleterious mutations
#' on X and standard Y chromosomes remain
#' at equilibrium.
##########################################

## Partially recessive mutations (h = 0.25)

# N = 10k, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U02_test")

# N = 100k, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U02")
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U05")
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U1")

# N = 500k, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U02")
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U05")
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 5*10^5, Nfname="_N500k_U1")

# N = 1mil, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U02")
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U05")
makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U1")

## Strongly recessive mutations (h = 0.1)

# N = 10k, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.1, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U02_test")

# N = 100k, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.1, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U02")
makeDataPrFixInvSizeHaploid(h = 0.1, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U05")
makeDataPrFixInvSizeHaploid(h = 0.1, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U1")

# N = 1mil, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.1, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U02")
makeDataPrFixInvSizeHaploid(h = 0.1, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U05")
makeDataPrFixInvSizeHaploid(h = 0.1, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U1")


## Strongly recessive mutations (h = 0.01)

# N = 10k, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.01, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^4, Nfname="_N10k_U02_test")

# N = 100k, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.01, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U02")
makeDataPrFixInvSizeHaploid(h = 0.01, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U05")
makeDataPrFixInvSizeHaploid(h = 0.01, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^5, Nfname="_N100k_U1")

# N = 1mil, broken up by Us.factor
makeDataPrFixInvSizeHaploid(h = 0.01, s = 0.01, U.vals = c(0.02),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U02")
makeDataPrFixInvSizeHaploid(h = 0.01, s = 0.01, U.vals = c(0.05),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U05")
makeDataPrFixInvSizeHaploid(h = 0.01, s = 0.01, U.vals = c(0.1),
					 nTot = 10^4, N = 10^6, Nfname="_N1mil_U1")




############################################
#' Autosomal inversions - Pr(fix | x)
#' -- Recursions from Conallon & Olito 2021
############################################

# h = 0.25, N = 10k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_U02")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_U05")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_U1")
# h = 0.25, N = 100k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_U02")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_U05")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_U1")

# h = 0.25, N = 250k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(2.5*10^5), Nfname="_N250k_U02")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(2.5*10^5), Nfname="_N250k_U05")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(2.5*10^5), Nfname="_N250k_U1")

# h = 0.25, N = 1 million, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_U02")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_U05")
makeDataAutoPrFixInvSize(h = 0.25, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_U1")


# h = 0.1, N = 10k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_U02")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_U05")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_U1")

# h = 0.1, N = 100k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_U02")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_U05")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_U1")

# h = 0.1, N = 250k, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(2.5*10^5), Nfname="_N250k_U02")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(2.5*10^5), Nfname="_N250k_U05")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(2.5*10^5), Nfname="_N250k_U1")

# h = 0.1, N = 1 million, broken up by Us.factor
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.02),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_U02")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.05),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_U05")
makeDataAutoPrFixInvSize(h = 0.1, s = 0.01, U.vals = c(0.1),
						 nTot = 10^4, N.vals = c(10^6), Nfname="_N1mil_U1")





##############################################
#' Timeseries of allele & inversion frequency
#' dynamics
##############################################

makeData_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
					nTot = 10^4, N = 10^5, Nfname = "_N100k")

makeData_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
					nTot = 10^4, N = 5*10^5, Nfname = "_N500k")

makeData_WFDynamics(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
					nTot = 10^4, N = 10^6, Nfname = "_N1mil")
