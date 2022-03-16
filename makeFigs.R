#'  Functions to generate figures for: 
#'    
#'  Title: 	A clarification of 'sheltering hypotheses' 
#' 			for the evolution of suppressed recombination 
#' 			between sex chromosomes 
#'
#'
#'  Author: Colin Olito
#'
#'
#'  NOTES: Run this file, either from terminal using Rscript,
#'		  or interactively in R. This should create all the 
#'		  figures needed to correctly compile the mansucript
#'		  LaTeX file.  
#'
rm(list=ls())
###############
# DEPENDENCIES
###############
source('./R/functions-figures.R')
source('./R/functions-figures-inversions-delMut.R')


########################
# Figures for the paper
########################




#  Fig.1 REVISED -- Illustration of deterministic dynamics
toPdf(deterministicDominanceIllustration(), 
			figPath(name='deterministicDomIllusFigCol.pdf'), width=10, height=10)
embed_fonts(figPath(name='deterministicDomIllusFigCol.pdf'))



#  Fig.2 REVISED -- Pr(fix) ~ x for Autosomal vs. SLR-expanding
toPdf(PrFixFig4Panel(), 
			figPath(name='PrFixFig4Panel_h0_25.pdf'), width=8, height=8)
embed_fonts(figPath(name='PrFixFig4Panel_h0_25.pdf'))


# Depricated figures
#  Fig.1 -- Illustration of deterministic dynamics
#toPdf(deterministicFig(), 
#			figPath(name='deterministicFig.pdf'), width=5, height=10)
#embed_fonts(figPath(name='deterministicFig.pdf'))

#toPdf(deterministicFigLessRecessive(), 
#			figPath(name='deterministicFigLessRec.pdf'), width=5, height=10)
#embed_fonts(figPath(name='deterministicFigLessRec.pdf'))

#toPdf(deterministicFigStrongRecessive(), 
#			figPath(name='deterministicFigStrongRec.pdf'), width=5, height=10)
#embed_fonts(figPath(name='deterministicFigStrongRec.pdf'))


#  Fig.2 -- Pr(fix) ~ x for Autosomal vs. SLR-expanding
#toPdf(PrFixFig(), 
#			figPath(name='PrFixFig_h0_25.pdf'), width=10, height=5)
#embed_fonts(figPath(name='PrFixFig_h0_25.pdf'))







########################
# Supplementary Figures
########################


#  Fig.S1 -- Deterministic Dynamics Overview Figure
toPdf(deterministicSuppFig(), 
			figPath(name='deterministicSuppFig_h0_1.pdf'), width=12, height=10)
embed_fonts(figPath(name='deterministicSuppFig_h0_1.pdf'))


#  Fig.S2 -- Deterministic Dynamics Overview Figure (h=0.01)
toPdf(deterministicSuppFig_h0.01(), 
			figPath(name='deterministicSuppFig_h0_01.pdf'), width=12, height=10)
embed_fonts(figPath(name='deterministicSuppFig_h0_01.pdf'))


#  Fig.S2 -- Deterministic Dynamics Overview Figure (h=0.25)
toPdf(deterministicSuppFig_h0.25(), 
			figPath(name='deterministicSuppFig_h0_25.pdf'), width=12, height=10)
embed_fonts(figPath(name='deterministicSuppFig_h0_25.pdf'))



#  Fig.S3 -- Fixation Probability for less recessive del. mut.
toPdf(PrFixFigRecessive(), 
			figPath(name='PrFixSuppFix_h0_1.pdf'), width=6, height=6)
embed_fonts(figPath(name='PrFixSuppFix_h0_1.pdf'))

#  Fig.S4 -- Fixation Probability for strongly recessive del. mut.
toPdf(PrFixFigStrongRecessive(), 
			figPath(name='PrFixSuppFix_h0_01.pdf'), width=6, height=6)
embed_fonts(figPath(name='PrFixSuppFix_h0_01.pdf'))



#  Fig.S5 -- Accumulation of del. mutations on Auto vs. SLR inversions (h=0.25)
toPdf(delMutAccumulationAutoVsSLR(), 
			figPath(name='delMutAccumulationAutoVsSLR.pdf'), width=12, height=5)
embed_fonts(figPath(name='delMutAccumulationAutoVsSLR.pdf'))






#' Selective advantage of inversion under heterozygote advantage
#' (similar to Fig. 1a of Charlesworth & Wall 1999)
toPdf(HetAdv3Locus(), 
			figPath(name='HetAdv3LocusCW1999.pdf'), width=7, height=7)
embed_fonts(figPath(name='HetAdv3LocusCW1999.pdf'))

toPdf(HetAdv3LocusQR(), 
			figPath(name='HetAdv3LocusQR.pdf'), width=10, height=10)
embed_fonts(figPath(name='HetAdv3LocusQR.pdf'))

toPdf(sheltering3LocusSuppFig(), 
			figPath(name='sheltering3Locus_SuppFig.pdf'), width=7, height=7)
embed_fonts(figPath(name='sheltering3Locus_SuppFig.pdf'))





######################
# Exploratory Figures
######################

toPdf(PrFixNeFig(), 
			figPath(name='PrFixNeFig.pdf'), width=5, height=10)
embed_fonts(figPath(name='PrFixNeFig.pdf'))


toPdf(PrFixAutoFig(), 
			figPath(name='PrFixAutoFig.pdf'), width=10, height=5)
embed_fonts(figPath(name='PrFixAutoFig.pdf'))


#  Comparison of t_ben vs. t_fix
toPdf(sim_Tben_Tfix_Fig(), 
			figPath(name='sim_Tben_Tfix_Fig.pdf'), width=5.5, height=5.5)
embed_fonts(figPath(name='sim_Tben_Tfix_Fig.pdf'))

########################
#' Two locus Model Figs

#' Selective advantage of inversion under heterozygote advantage
#' (similar to Fig. 1a of Charlesworth & Wall 1999)
toPdf(heterozygoteAdvantageFig(), 
			figPath(name='hetAdvFig.pdf'), width=7, height=7)
embed_fonts(figPath(name='hetAdvFig.pdf'))

#' What does linkage do to the heterozygote advantage model?
toPdf(heterozygoteAdvantageFig(df = "./data/hetAdvSimData_r0.1.csv"), 
			figPath(name='hetAdvFig_r0.1.pdf'), width=7, height=7)
embed_fonts(figPath(name='hetAdvFig_r0.1.pdf'))

toPdf(heterozygoteAdvantageFig(df = "./data/hetAdvSimData_r0.05.csv"), 
			figPath(name='hetAdvFig_r0.05.pdf'), width=7, height=7)
embed_fonts(figPath(name='hetAdvFig_r0.05.pdf'))

 
#' Selective advantage of inversion under Sheltering Scenario
toPdf(shelteringFig(), 
			figPath(name='shelteringFig.pdf'), width=7, height=7)
embed_fonts(figPath(name='shelteringFig.pdf'))

#' Selective advantage of inversion under Sheltering Scenario
#' when there are sex-specific selection coefficients
toPdf(shelteringSexSpecificFig(), 
			figPath(name='shelteringSexSpecificFig.pdf'), width=7, height=7)
embed_fonts(figPath(name='shelteringSexSpecificFig.pdf'))





########################
#' Three Locus Model Figs

#' Selective advantage of inversion under heterozygote advantage
#' (similar to Fig. 1a of Charlesworth & Wall 1999)
toPdf(HetAdv3Locus(), 
			figPath(name='HetAdv3LocusCW1999.pdf'), width=7, height=7)
embed_fonts(figPath(name='HetAdv3LocusCW1999.pdf'))

toPdf(HetAdv3LocusQR(), 
			figPath(name='HetAdv3LocusQR.pdf'), width=10, height=10)
embed_fonts(figPath(name='HetAdv3LocusQR.pdf'))

#' Selective advantage of inversion under Sheltering Scenario
toPdf(sheltering3LocusFig(), 
			figPath(name='sheltering3LocusFig.pdf'), width=7, height=7)
embed_fonts(figPath(name='sheltering3LocusFig.pdf'))

#' Selective advantage of inversion under Sheltering Scenario
toPdf(sheltering3LocusFig2(), 
			figPath(name='sheltering3LocusFig2.pdf'), width=7, height=7)
embed_fonts(figPath(name='sheltering3LocusFig2.pdf'))
