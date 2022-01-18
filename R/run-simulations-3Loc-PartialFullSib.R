################################################################
#  RUN SIMULATIONS and CREATE OUTPUT DATA FOR PLOTTING
#  
#  R code for forward 3 locus simulations. Generates output data
#  as .csv files saved to ./data/
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		


#####################
##  Dependencies
rm(list=ls())
source('R/simulations-3Loc-PartialFullSib.R')


########################################
# Heterozygote advantage sanity check
########################################
#' Note: these simulations create data
#' to reproduce a plot similar to Fig. 1a
#' in Charlesworth & Wall (1999). We keep
#' the same selection coefficients at the
#' A locus, assume the B locus is neutral,
#' and that all three loci, SDL, A , & B
#' recombine freely (q = r = 1/2)

s.vals  <-  rbind(c(0.1, 0.2, 0.1, 0.2),
				  c(0, 0, 0, 0))
t.vals  <-  rbind(c(0.1, 0.2, 0.2, 0.1),
				  c(0, 0, 0, 0))

make3LocusHetAdvSimData(generations = 1000, threshold=10^-8, 
						s.vals=s.vals, t.vals=t.vals, q=1/2, r=1/2)


# What does linkage do?
s.vals  <-  rbind(c(0.05, 0.1, 0.05, 0.1),
				  c(0.05, 0.1, 0.05, 0.1))
t.vals  <-  rbind(c(0.05, 0.1, 0.1, 0.05),
				  c(0.05, 0.1, 0.1, 0.05))

# make3LocusHetAdvSimData(generations = 1000, threshold=10^-8, 
#						s.vals=s.vals, t.vals=t.vals, q=1/2, r=1/2)

make3LocusHetAdvSimData(generations = 1000, threshold=10^-8, 
						s.vals=s.vals, t.vals=t.vals, q=0.05, r=1/2)

make3LocusHetAdvSimData(generations = 1000, threshold=10^-8, 
						s.vals=s.vals, t.vals=t.vals, q=0.01, r=1/2)

make3LocusHetAdvSimData(generations = 1000, threshold=10^-8, 
						s.vals=s.vals, t.vals=t.vals, q=0.05, r=0.01)

make3LocusHetAdvSimData(generations = 1000, threshold=10^-8, 
						s.vals=s.vals, t.vals=t.vals, q=0.01, r=0.01)




########################################
# Sheltering
########################################

q.vals = c(0.5, 0.1, 0.01, 0.001)
r.vals = c(0.5, 0.5, 0.5, 0.5)
hf.vals = rbind(c(0, 0.1),
				c(0, 0.1))
sf.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))
hm.vals = rbind(c(0, 0.1),
				c(0, 0.1))
sm.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))

makeDataShelteringFig(generations = 5000, threshold=10^-8, 
					  hf.vals=hf.vals, sf.vals=sf.vals, 
					  hm.vals=hm.vals, sm.vals=sm.vals,
					  q.vals=q.vals, r.vals=r.vals,
					  u=10^-5, v=10^-5) 

q.vals = c(0.5, 0.1, 0.01, 0.001)
r.vals = c(0.5, 0.5, 0.5, 0.5)
hf.vals = rbind(c(0.01, 0.25),
				c(0.01, 0.25))
sf.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))
hm.vals = rbind(c(0.01, 0.25),
				c(0.01, 0.25))
sm.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))

makeDataShelteringFig(generations = 5000, threshold=10^-8, 
					  hf.vals=hf.vals, sf.vals=sf.vals, 
					  hm.vals=hm.vals, sm.vals=sm.vals,
					  q.vals=q.vals, r.vals=r.vals,
					  u=10^-5, v=10^-5, filename="rExplore") 


q.vals = c(0.5, 0.001)
r.vals = c(0.5, 0.5)
hf.vals = rbind(c(0.1),
				c(0.1))
sf.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))
hm.vals = rbind(c(0.1),
				c(0.1))
sm.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))

makeDataShelteringFig(generations = 5000, threshold=10^-8, 
					  hf.vals=hf.vals, sf.vals=sf.vals, 
					  hm.vals=hm.vals, sm.vals=sm.vals,
					  q.vals=q.vals, r.vals=r.vals,
					  u=10^-5, v=10^-5, filename = "qExplore") 


q.vals = c(0.5)
r.vals = c(0.5)
hf.vals = rbind(c(0.1),
				c(0.1))
sf.vals = rbind(c(0.02),
				c(0.01))
hm.vals = rbind(c(0.1),
				c(0.1))
sm.vals = rbind(c(0.02),
				c(0.01))

makeDataShelteringFig(generations = 5000, threshold=10^-8, 
					  hf.vals=hf.vals, sf.vals=sf.vals, 
					  hm.vals=hm.vals, sm.vals=sm.vals,
					  q.vals=q.vals, r.vals=r.vals,
					  u=10^-5, v=10^-5, filename = "sAlarge_sBsmall")


q.vals = c(0.5)
r.vals = c(0.0)
hf.vals = rbind(c(0.25, 0.01),
				c(0.25, 0.01))
sf.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))
hm.vals = rbind(c(0.25, 0.01),
				c(0.25, 0.01))
sm.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))

makeDataShelteringFig(generations = 5000, threshold=10^-8, 
					  hf.vals=hf.vals, sf.vals=sf.vals, 
					  hm.vals=hm.vals, sm.vals=sm.vals,
					  q.vals=q.vals, r.vals=r.vals,
					  u=10^-5, v=10^-5, filename = "q0.5_r0.0_hsCompare")

q.vals = c(0.5, 0.5, 0.5, 0.5)
r.vals = c(0.5, 0.1, 0.01, 0.001)
hf.vals = rbind(c(0.25, 0.01),
				c(0.25, 0.01))
sf.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))
hm.vals = rbind(c(0.25, 0.01),
				c(0.25, 0.01))
sm.vals = rbind(c(0.01, 0.05),
				c(0.01, 0.05))

makeDataShelteringFig(generations = 5000, threshold=10^-8, 
					  hf.vals=hf.vals, sf.vals=sf.vals, 
					  hm.vals=hm.vals, sm.vals=sm.vals,
					  q.vals=q.vals, r.vals=r.vals,
					  u=10^-5, v=10^-5, filename = "q0.5_rhsCompare")


q.vals = c(0.5, 0.5, 0.001, 0.001)
r.vals = c(0.5, 0.001, 0.5, 0.001)
hf.vals = rbind(c(0.25, 0.01),
				c(0.25, 0.01))
sf.vals = rbind(c(0.01),
				c(0.01))
hm.vals = rbind(c(0.25, 0.01),
				c(0.25, 0.01))
sm.vals = rbind(c(0.01),
				c(0.01))

makeDataShelteringFig(generations = 5000, threshold=10^-8, 
					  hf.vals=hf.vals, sf.vals=sf.vals, 
					  hm.vals=hm.vals, sm.vals=sm.vals,
					  q.vals=q.vals, r.vals=r.vals,
					  u=10^-5, v=10^-5, filename = "qrhsCompare")