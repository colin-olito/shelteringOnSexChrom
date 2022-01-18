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
source('R/simulations-PartialFullSib.R')


########################################
# Heterozygote advantage
########################################
#' Note: these simulations create data
#' to produce a plot similar to Fig. 1a
#' in Charlesworth & Wall (1999), but 
#' for our scenario where the selected
#' locus is in the PAR

s.vals  <-  c(0.1, 0.2, 0.1, 0.2)
t.vals  <-  c(0.1, 0.2, 0.2, 0.1)

makeDataHetAdvFig(generations = 10000, threshold=10^-8, 
				s.vals=s.vals, t.vals=t.vals, r = 1/2)

# What does linkage do?
makeDataHetAdvFig(generations = 10000, threshold=10^-8, 
				s.vals=s.vals, t.vals=t.vals, r = 0.1)

makeDataHetAdvFig(generations = 10000, threshold=10^-8, 
				s.vals=s.vals, t.vals=t.vals, r = 0.05)

########################################
# Sheltering
########################################

r.vals = c(0.001, 0.01, 0.1, 0.5)
h.vals = c(0, 0.1)
s.vals = c(0.01, 0.05)

makeDataShelteringFig(generations = 10000, threshold=10^-8, 
					  s.vals=s.vals, h.vals=h.vals, r.vals=r.vals, u = 10^-5) 


r.vals = c(0.001, 0.01, 0.1, 0.5)
h.vals = c(0, 0.1)
sf.vals = c(0.01, 0.05)
sm.vals = c(0.05, 0.01)

makeDataShelteringSexSpecific(generations = 10000, threshold=10^-8, 
					  sf.vals=sf.vals, sm.vals=sm.vals, h.vals=h.vals, r.vals=r.vals, u = 10^-5) 
