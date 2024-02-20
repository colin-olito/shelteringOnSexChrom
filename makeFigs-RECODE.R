#'  Makefigs file after RECODE to check for mistakes
#'
#'
###############
# DEPENDENCIES
###############
rm(list=ls())
source('./R/functions-figures.R')
source('./R/functions-figures-inversions-delMut-RECODE.R')

######################################
#' Create figures directories if they
#' do not already exist
figuresDirectoriesExist  <-  dir.exists("./figures")
figuresRecodeDirExist    <-  dir.exists("./figures/RECODE")

if(!figuresDirectoriesExist) {
	dir.create("./figures")
}
if(!figuresRecodeDirExist) {
	dir.create("./figures/RECODE")
}

########################
# Figures for the paper
########################

#  Fig.1 CORRECTION -- Illustration of deterministic dynamics
#' Deterministic dynamics for Y-linked inversions 
#' Illustrating effect of dominance where h = {0.25, 0.1, 0.01}
toPdf(deterministicDominanceIllustration(wHap = FALSE, x=0.3), 
			figPath(name='./RECODE/deterministicDomIllusFigCol_BCharlesworth_RECODE_R1.pdf'), width=12, height=12)
embed_fonts(figPath(name='./RECODE/deterministicDomIllusFigCol_BCharlesworth_RECODE_R1.pdf'))


#' Fig.2 CORRECTION -- Pr(fix) ~ x for Autosomal vs. SLR-expanding
#' 					using Eq. approximation for both Autosomal
#' 					and SLR-expanding inversions
#rm(list=ls())
#source('./R/functions-figures.R')
#source('./R/functions-figures-inversions-delMut-RECODE.R')
toPdf(PrFixFig4Panel_Approx(), 
			figPath(name='./RECODE/PrFixFig4Panel_Approx_h0_25_RECODE_Haploid.pdf'), width=8, height=8)
embed_fonts(figPath(name='./RECODE/PrFixFig4Panel_Approx_h0_25_RECODE_Haploid.pdf'))



########################
# Supplementary Figures 
########################

#' Supp. fig. -- Pr(fix) ~ x for Autosomal vs. SLR-expanding
#' 					using Exact recursions for both Autosomal & SLR-
#' 					expanding inversions
toPdf(PrFixFig4Panel_Exact(suppVersion=TRUE), 
			figPath(name='./RECODE/PrFixFig4Panel_h0_25_RECODE_Exact.pdf'), width=8, height=8)
embed_fonts(figPath(name='./RECODE/PrFixFig4Panel_h0_25_RECODE_Exact.pdf'))



#' Deterministic dynamics for Y-linked inversions 
#' Illustrating effect of mutation rate where U = {0.02, 0.05, 0.1}
toPdf(deterministicMutRateIllustration(wHap=FALSE), 
			figPath(name='./RECODE/deterministicMutRateIllusFigCol_BCharlesworth_RECODE_R1.pdf'), width=12, height=12)
embed_fonts(figPath(name='./RECODE/deterministicMutRateIllusFigCol_BCharlesworth_RECODE_R1.pdf'))



#'  Supplementary Figures -- Deterministic Dynamics Overview Figures
#' Partial recessivity (h = 0.25) 
toPdf(deterministicSuppFig(h=0.25, U = 0.02), 
			figPath(name='./RECODE/deterministicSuppFig_h0_25_U02_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/deterministicSuppFig_h0_25_U02_RECODE.pdf'))

toPdf(deterministicSuppFig(h=0.25, U = 0.05), 
			figPath(name='./RECODE/deterministicSuppFig_h0_25_U05_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/deterministicSuppFig_h0_25_U05_RECODE.pdf'))

toPdf(deterministicSuppFig(h=0.25, U = 0.1), 
			figPath(name='./RECODE/deterministicSuppFig_h0_25_U1_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/deterministicSuppFig_h0_25_U1_RECODE.pdf'))


#' Strongly recessive (h = 0.1) 
toPdf(deterministicSuppFig(h=0.1, U = 0.02), 
			figPath(name='./RECODE/deterministicSuppFig_h0_1_U02_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/deterministicSuppFig_h0_1_U02_RECODE.pdf'))

toPdf(deterministicSuppFig(h=0.1, U = 0.05), 
			figPath(name='./RECODE/deterministicSuppFig_h0_1_U05_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/deterministicSuppFig_h0_1_U05_RECODE.pdf'))

toPdf(deterministicSuppFig(h=0.1, U = 0.1), 
			figPath(name='./RECODE/deterministicSuppFig_h0_1_U1_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/deterministicSuppFig_h0_1_U1_RECODE.pdf'))


#######################
# Autosomal Inversions
rm(list=ls())
 source('./R/functions-figures-inversions-delMut-RECODE.R')
 source('./R/Auto-Inversions-NeiMyVersion.R')

toPdf(Nei_1967_recessiveMutations_AutoExact(), 
			figPath(name='./RECODE/AutoExact-Nei67-Fig1.pdf'), width=7, height=7)
embed_fonts(figPath(name='./RECODE/AutoExact-Nei67-Fig1.pdf'))

#' Deterministic dynamics of Autosomal Inversions using exact recursions
#' Illustrating effect of mutation rate where U = {0.02, 0.05, 0.1}
toPdf(deterministicMutRateIllustration_AutoExact(), 
			figPath(name='./RECODE/AutoExact-MutRateIllusFigCol_RECODE.pdf'), width=10, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-MutRateIllusFigCol_RECODE.pdf'))

#' Illustrating effect of dominance where h = {0.25, 0.1, 0.05}
toPdf(deterministicDominanceIllustration_AutoExact(), 
			figPath(name='./RECODE/AutoExact-DominanceIllusFigCol_RECODE.pdf'), width=10, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-DominanceIllusFigCol_RECODE.pdf'))





#'  Supplementary Figures -- Deterministic Dynamics Overview Figures
#' Partial recessivity (h = 0.25) 

toPdf(deterministicSuppFig_AutoExact(h=0.25, U = 0.02), 
			figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_25_U02_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_25_U02_RECODE.pdf'))

toPdf(deterministicSuppFig_AutoExact(h=0.25, U = 0.05), 
			figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_25_U05_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_25_U05_RECODE.pdf'))

toPdf(deterministicSuppFig_AutoExact(h=0.25, U = 0.1), 
			figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_25_U1_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_25_U1_RECODE.pdf'))


#' Strongly recessive (h = 0.1) 

toPdf(deterministicSuppFig_AutoExact(h=0.1, U = 0.02), 
			figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_1_U02_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_1_U02_RECODE.pdf'))

toPdf(deterministicSuppFig_AutoExact(h=0.1, U = 0.05), 
			figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_1_U05_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_1_U05_RECODE.pdf'))

toPdf(deterministicSuppFig_AutoExact(h=0.1, U = 0.1), 
			figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_1_U1_RECODE.pdf'), width=12, height=10)
embed_fonts(figPath(name='./RECODE/AutoExact-deterministicSuppFig_h0_1_U1_RECODE.pdf'))




## WF simulation figure for h = 0.1
rm(list=ls())
source('./R/functions-figures-inversions-delMut-RECODE.R')

#PrFixFig4Panel_Exact(suppVersion=FALSE)
#toPdf(PrFixFig4Panel_SuppFig(), 
#			figPath(name='./RECODE/PrFixFig4Panel_SuppFig_h0_1_RECODE.pdf'), width=8, height=8)
#embed_fonts(figPath(name='./RECODE/PrFixFig4Panel_SuppFig_h0_1_RECODE.pdf'))

toPdf(PrFixFig2Panel_SuppFig(), 
			figPath(name='./RECODE/PrFixFig2Panel_SuppFig_h0_1_RECODE.pdf'), width=8, height=5)
embed_fonts(figPath(name='./RECODE/PrFixFig2Panel_SuppFig_h0_1_RECODE.pdf'))





## Figures examining frequency dynamics during WF simulations.
rm(list=ls())
source('./R/functions-figures.R')
source('./R/functions-figures-inversions-delMut-RECODE.R')

# Overview of all inversion sizes
toPdf(SLR_WF_DynamicsSuppFig(), 
			figPath(name='./RECODE/SLR_Nei_WF_Dynamics_SuppFig.pdf'), width=50, height=14)
embed_fonts(figPath(name='./RECODE/SLR_Nei_WF_Dynamics_SuppFig.pdf'))



rm(list=ls())
 source('./R/functions-figures-inversions-delMut-RECODE.R')
 source('./R/Auto-Inversions-NeiMyVersion.R')


# Overview of all inversion sizes
toPdf(AutoExact_WF_DynamicsFigSuppFig(), 
			figPath(name='./RECODE/AutoExact_WF_Dynamics_SuppFig.pdf'), width=50, height=14)
embed_fonts(figPath(name='./RECODE/AutoExact_WF_Dynamics_SuppFig.pdf'))


toPdf(AutoExact_WF_DynamicsFig(invSize=0.8), 
			figPath(name='./RECODE/AutoExact_WF_Dynamics_x_0.8.pdf'), width=7, height=14)
embed_fonts(figPath(name='./RECODE/AutoExact_WF_Dynamics_x_0.8.pdf'))





toPdf(perspectiveInvRelFit(), 
			figPath(name='./RECODE/perspectiveInvRelFit.pdf'), width=7, height=6)
embed_fonts(figPath(name='./RECODE/perspectiveInvRelFit.pdf'))
