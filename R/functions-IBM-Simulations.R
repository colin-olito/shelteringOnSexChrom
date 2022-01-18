# Individual-Based Simulation Code: 
# Inversions expanding SLR w/ recessive deleterious mutations

rm(list=ls())


# install.packages("foreach")
library(foreach)
#install.packages("doParallel")
library(doParallel)
#install.packages("doSNOW")
library(doSNOW)


mutateReplace  <-  function(x, p, ...) {
	x[x == 0]     <-  rbinom(n=length(x[x == 0]), size=1, prob=p)
	x
}

relFitness  <-  function(mat1, mat2, W) {
	genotype  <-  mat1 + mat2 + 1
	W.i       <-  as.vector(apply(X=genotype, MARGIN=1, function(x){ prod(W[x])} ))
	W.i/mean(W.i)
}

recombine.sampOvules  <-  function(v1, v2, x, n) {
	if(runif(1) < x) {
		breakpoint        <-  floor(n*runif(1, min=1/n, max=1-(1/n)))
		temp              <-  v1[breakpoint:n]
		v1[breakpoint:n]  <-  v2[breakpoint:n]
		v2[breakpoint:n]  <-  temp
	}
	get(ifelse(runif(1) > 1/2, "v1", "v2"))
}

freeRecombine.sampOvules  <-  function(v1, v2, n) {
	shuffle  <-  rbinom(n=n, size=1, prob=1/2)
	c1       <-  rep(0,times=n)
	c2       <-  rep(0,times=n)
	c1[shuffle==0]  <-  v1[shuffle==0]
	c1[shuffle==1]  <-  v2[shuffle==1]
	c2[shuffle==0]  <-  v2[shuffle==0]
	c2[shuffle==1]  <-  v1[shuffle==1]
	get(ifelse(runif(1) > 1/2, "c1", "c2"))
}

recombine.sampSperm  <-  function(v1, v2, x, n, inversion, sons) {
	if(inversion) {
		get(ifelse(sons, "v2", "v1"))
	}
	else {
		if(runif(1) < x) {
			breakpoint        <-  floor(n*runif(1, min=1/n, max=1-(1/n)))
			temp              <-  v1[breakpoint:n]
			v1[breakpoint:n]  <-  v2[breakpoint:n]
			v2[breakpoint:n]  <-  temp
		}
		get(ifelse(sons, "v2", "v1"))
	}
}

#a=c(1,1,1,1,1,1,1,1,1,1,1,1)
#b=c(2,2,2,2,2,2,2,2,2,2,2,2)
#freeRecombine.sampSperm(v1=a, v2=b, n=length(a), inversion=FALSE, sons=FALSE)

#a=matrix(c(1,1,1,1,1,1,1,1,1,1,1,1), nrow=2)
#b=matrix(c(2,2,2,2,2,2,2,2,2,2,2,2), nrow=2)
#t(apply(1:2, function(X) freeRecombine.sampSperm(v1=a[X,], v2=b[X,], n=ncol(a), inversion=FALSE, sons=FALSE)))


freeRecombine.sampSperm  <-  function(v1, v2, n, inversion, sons) {
	if(inversion) {
		return(get(ifelse(sons, "v2", "v1")))
	}
	if(!inversion) {
		shuffle  <-  rbinom(n=n, size=1, prob=1/2)
		c1       <-  rep(0,times=n)
		c2       <-  rep(0,times=n)
		c1[shuffle==0]  <-  v1[shuffle==0]
		c1[shuffle==1]  <-  v2[shuffle==1]
		c2[shuffle==0]  <-  v2[shuffle==0]
		c2[shuffle==1]  <-  v1[shuffle==1]
		return(get(ifelse(sons, "c2", "c1")))
	}
}


freeRecombine.sampOvules.FL  <-  function(m1, m2, n, moms) {
	f1.mat  <-  matrix(NA,nrow=nrow(m2), ncol=ncol(m2))
	for(i in 1:nrow(m2)) {
		shuffle  <-  rbinom(n=n, size=1, prob=1/2)
		c1       <-  rep(0,times=n)
		c2       <-  rep(0,times=n)
		c1[shuffle==0]  <-  m1[moms[i],][shuffle==0]
		c1[shuffle==1]  <-  m2[moms[i],][shuffle==1]
		c2[shuffle==0]  <-  m2[moms[i],][shuffle==0]
		c2[shuffle==1]  <-  m1[moms[i],][shuffle==1]
		f1.mat[i,]      <-  get(ifelse(runif(1) > 1/2, "c2", "c1"))
	}
	return(f1.mat)
}

freeRecombine.sampSperm.FL  <-  function(m1, m2, n, dads, inverted, sons) {
	
	f1.mat  <-  matrix(NA,nrow=nrow(m2), ncol=ncol(m2))
	for(i in 1:nrow(m2)) {
		if(inverted[i]) {
			f1.mat[i,]  <-  get(ifelse(sons, "m2", "m1"))[dads[i],]
		}
		if(!inverted[i]) {
			shuffle  <-  rbinom(n=n, size=1, prob=1/2)
			c1       <-  rep(0,times=n)
			c2       <-  rep(0,times=n)
			c1[shuffle==0]  <-  m1[dads[i],][shuffle==0]
			c1[shuffle==1]  <-  m2[dads[i],][shuffle==1]
			c2[shuffle==0]  <-  m2[dads[i],][shuffle==0]
			c2[shuffle==1]  <-  m1[dads[i],][shuffle==1]
			f1.mat[i,]  <-  get(ifelse(sons, "c2", "c1"))
		}
	}
	return(f1.mat)
}

#testY.f1 <-  freeRecombine.sampSperm.FL(m1=Xm, m2=Y, dads=dads[1:(N/2)], n=n, inverted=dadInverted, sons=TRUE)

#testY.f1[invertedY[4],]
#sum(testY.f1[invertedY[4],])

makeDataPrFixIBMSimParallel  <-  function(h = 0.1, s = 0.01, Us.factor = 2,
										  nTot = 10^3, N = 10^3, invSize = 0.05, 
										  nCluster = 4, extraFileID = "",
										  burnin=100) {

	# Set number of cores to parallelize over
	numCluster  <-  2*(detectCores() - 1)
	if(nCluster < numCluster && !is.na(nCluster) ) {
		numCluster  <-  nCluster
	}
#	registerDoParallel(numCluster)
	cl  <-  makeCluster(numCluster, type="SOCK")
	registerDoSNOW(cl)

	# Fitness effects for indiv. loci
	# W.i  <-  c(w_AA, w_Aa, w_aa)
	W.i  <-  c(1, 1 - h*s, 1 - s)

	#number of simulations 
	sims  = 100*N/2

	# Mutation rate and equilibrium frequencies prior to inversion
	U     <-  Us.factor*s
	u     <-  U/nTot
	qHat  <-  (U/(nTot*h*s))
			
	# number of loci captured by inversion
	n    <-  nTot*invSize

	# Initiate chromosome matrices (rows = individuals, cols = loci spanned by inversion)
	# and populate with mutations
	Xf.mat  <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # Maternally-inherited X's in females
	Xf.pat  <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # Paternally-inherited X's in females
	Xm      <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # X chromosomes in Males
	Y       <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # Y chromosomes
par(mfrow=c(2,2))
hist(colSums(rbind(Xf.mat, Xf.pat, Xm, Y)), xlim=c(0,(10*qHat*N)))
hist(colSums(Y), xlim=c(0,(10*qHat*N/2)))
	
mean(colSums(rbind(Xf.mat, Xf.pat, Xm, Y))/(2*N)	)
	# burnin to equilbrium prior to inversion
	for(t in 1:burnin) {
		cat('\r', paste(sprintf("burnin generation %s / ", t), burnin, " is complete"))
		# Relative fitness & sampling probs. 
		Wf.i          <-  relFitness(mat1=Xf.mat, mat2=Xf.pat, W=W.i)
		Wm.i          <-  relFitness(mat1=Xm,     mat2=Y,      W=W.i)
		# Sample parent gametes
		moms          <-  sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wf.i/(N/2))
		dads          <-  sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wm.i/(N/2))
		# Meiosis and fertilization
		# for sons
#		Xm.f1         <-  t(sapply(1:nrow(Xm), function(X) freeRecombine.sampOvules(v1=Xf.mat[moms[X],], v2=Xf.pat[moms[X],], n=n), simplify=TRUE))
#		Y.f1          <-  t(sapply(1:nrow(Y), function(X) freeRecombine.sampSperm(v1=Xm[dads[X],], v2=Y[dads[X],], n=n, inversion=dadInverted[X], sons=TRUE)))
		Xm.f1         <-  freeRecombine.sampOvules.FL(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[1:(N/2)])
		Y.f1          <-  freeRecombine.sampSperm.FL(m1=Xm, m2=Y, n=n, dads=dads[1:(N/2)], inverted=rep(FALSE, times=(N/2)), sons=TRUE)
		# for daughters
#		Xf.mat.f1     <-  t(sapply(1:nrow(Xf.mat), function(X) freeRecombine.sampOvules(v1=Xf.mat[moms[(N/2)+X],], v2=Xf.pat[moms[(N/2)+X],], n=n), simplify=TRUE ))
#		Xf.pat.f1     <-  t(sapply(1:nrow(Xf.pat), function(X) freeRecombine.sampSperm(v1=Xm[dads[(N/2)+X],], v2=Y[dads[(N/2)+X],], n=n, inversion=dadInverted[(N/2)+X], sons=FALSE), simplify=TRUE ))
		Xf.mat.f1     <-  freeRecombine.sampOvules.FL(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[((N/2)+1):N])
		Xf.pat.f1     <-  freeRecombine.sampSperm.FL(m1=Xm, m2=Y, n=n, dads=dads[((N/2)+1):N], inverted=rep(FALSE, times=(N/2)), sons=FALSE)
		# Mutation
		Xf.mat.f1[Xf.mat.f1 == 0]  <-  rbinom(n=length(Xf.mat.f1[Xf.mat.f1 == 0]), size=1, prob=u)
		Xf.pat.f1[Xf.pat.f1 == 0]  <-  rbinom(n=length(Xf.pat.f1[Xf.pat.f1 == 0]), size=1, prob=u)
		Xm.f1[Xm.f1 == 0]          <-  rbinom(n=length(Xm.f1[Xm.f1 == 0]),         size=1, prob=u)
		Y.f1[Y.f1 == 0]            <-  rbinom(n=length(Y.f1[Y.f1 == 0]),           size=1, prob=u)
		# F1's become adults in next generation
		Xm      <-  Xm.f1
		Y       <-  Y.f1
		Xf.mat  <-  Xf.mat.f1
		Xf.pat  <-  Xf.pat.f1
	}
hist(colSums(rbind(Xf.mat, Xf.pat, Xm, Y)), col=2)#, add=T)
hist(colSums(Y), col=2)#, add=T)

	# assign initial state chromosome matrices
	Xm.init      <-  Xm
	Y.init       <-  Y
	Xf.mat.init  <-  Xf.mat
	Xf.pat.init  <-  Xf.pat

# hist(rowSums(rbind(Xm.init, Y.init,Xf.mat.init, Xf.pat.init)))
# hist(colSums(rbind(Xm.init, Y.init,Xf.mat.init, Xf.pat.init)), breaks=40)
# table(colSums(rbind(Xm.init, Y.init, Xf.mat.init, Xf.pat.init)))
# 1/sum(rowSums(Y) == 0)

	# Export list for %dopar%
	ex.vec  <-  c("mutateReplace", "relFitness", 
				  "freeRecombine.sampOvules", "freeRecombine.sampSperm",
				  "freeRecombine.sampOvules.FL", "freeRecombine.sampSperm.FL")

	# Progress
	progress <- function(n) cat('\r', paste(sprintf("simulation %d /", n), sims, " is complete"))
	opts <- list(progress = progress)

	# Loop over replicate simulations
	fix  <-  foreach(i=icount(sims), .combine= '+', 
					 .options.snow=opts, .export = ex.vec) %dopar% {

		# shuffle alleles at each locus among individuals to give
		# unique initial equilibrium state
		Xf.mat  <-  apply(Xf.mat.init, MARGIN=2, function(x) x[sample(seq_along(x))])
		Xf.pat  <-  apply(Xf.pat.init, MARGIN=2, function(x) x[sample(seq_along(x))])
		Xm      <-  apply(Xm.init,     MARGIN=2, function(x) x[sample(seq_along(x))])
		Y       <-  apply(Y.init,      MARGIN=2, function(x) x[sample(seq_along(x))])

		# Single-copy inversion mutation occurs on randomly chosen Y chromosome 
		invertedYs  <-  sample(c(1:(N/2)), size=1)
		YI.t        <-  length(invertedYs)/(N/2)
print(sum(Y[invertedYs,]))

		# Run forward simulation
		while(YI.t*(1 - YI.t) > 0) {

			# Relative fitness & sampling probs. 
			Wf.i         <-  relFitness(mat1=Xf.mat, mat2=Xf.pat, W=W.i)
			Wm.i         <-  relFitness(mat1=Xm,     mat2=Y,      W=W.i)
#hist(Wm.i)
#abline(v=1, col=1, lwd=3)
#abline(v=Wm.i[invertedYs], col=2, lwd=3)
#print(mean(Wm.i[invertedYs]))
sum(Y[invertedYs,])
rbind(Y[invertedYs,],Xm[invertedYs,])
hist(rowSums(Y))
abline(v=mean(rowSums(Y)))
qHat*n
			# Sample N mothers & fathers with probabilities weighted by their relative fitnesses
			moms         <-  sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wf.i/(N/2))
			dads         <-  sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wm.i/(N/2))
			# Which dads carry the inversion, which sons will inherit it, calculate inversion frequency among sons
			dadInverted  <-  dads %in% invertedYs
			invertedYs   <-  which(dadInverted)[which(dadInverted) <= (N/2)]
			YI.t         <-  length(invertedYs)/(N/2)
			# Meiosis and fertilization
			# for sons
#			Xm.f1         <-  t(sapply(1:nrow(Xm), function(X) freeRecombine.sampOvules(v1=Xf.mat[moms[X],], v2=Xf.pat[moms[X],], n=n), simplify=TRUE))
#			Y.f1          <-  t(sapply(1:nrow(Y), function(X) freeRecombine.sampSperm(v1=Xm[dads[X],], v2=Y[dads[X],], n=n, inversion=dadInverted[X], sons=TRUE)))
			Xm.f1        <-  freeRecombine.sampOvules.FL(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[1:(N/2)])
			Y.f1         <-  freeRecombine.sampSperm.FL(m1=Xm, m2=Y, n=n, dads=dads[1:(N/2)], inverted=dadInverted[1:(N/2)], sons=TRUE)
#print(length(invertedYs))
			# for daughters
#			Xf.mat.f1     <-  t(sapply(1:nrow(Xf.mat), function(X) freeRecombine.sampOvules(v1=Xf.mat[moms[(N/2)+X],], v2=Xf.pat[moms[(N/2)+X],], n=n), simplify=TRUE ))
#			Xf.pat.f1     <-  t(sapply(1:nrow(Xf.pat), function(X) freeRecombine.sampSperm(v1=Xm[dads[(N/2)+X],], v2=Y[dads[(N/2)+X],], n=n, inversion=dadInverted[(N/2)+X], sons=FALSE), simplify=TRUE ))
			Xf.mat.f1     <-  freeRecombine.sampOvules.FL(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[((N/2)+1):N])
			Xf.pat.f1     <-  freeRecombine.sampSperm.FL(m1=Xm, m2=Y, n=n, dads=dads[((N/2)+1):N], inverted=dadInverted[((N/2)+1):N], sons=FALSE)
			# Mutation
			Xf.mat.f1[Xf.mat.f1 == 0]  <-  rbinom(n=length(Xf.mat.f1[Xf.mat.f1 == 0]), size=1, prob=u)
			Xf.pat.f1[Xf.pat.f1 == 0]  <-  rbinom(n=length(Xf.pat.f1[Xf.pat.f1 == 0]), size=1, prob=u)
			Xm.f1[Xm.f1 == 0]          <-  rbinom(n=length(Xm.f1[Xm.f1 == 0]),         size=1, prob=u)
			Y.f1[Y.f1 == 0]            <-  rbinom(n=length(Y.f1[Y.f1 == 0]),           size=1, prob=u)
			# F1's become adults in next generation
			Xm      <-  Xm.f1
			Y       <-  Y.f1
			Xf.mat  <-  Xf.mat.f1
			Xf.pat  <-  Xf.pat.f1
		}
		YI.t
	}
	# Calculate proportion of fixations
	PrFix  <-  fix/sims

	# Close cluster
	stopCluster(cl)

	# Export results as a dataframe
	filename  <-  paste("./data/PrFixIBM_h", h, "_s", s, "_N", N, "_Ufac", Us.factor, "_", extraFileID, ".csv", sep="")
	d  <-  data.frame(
										"h"        =  rep(h, times=length(PrFix)),
										"s"        =  rep(s, times=length(PrFix)),
										"N"        =  rep(N, times=length(PrFix)),
										"U"        =  rep(U, times=length(PrFix)),
										"Ufac"     =  rep(Us.factor, times=length(PrFix)),
										"x"        =  invSize,
										"PrFix"    =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}





PoC_invasionFitness  <-  function(h = 0.1, s = 0.01, Us.factor = 5,
								  nTot = 10^4, N = 4*10^3, invSize = 0.05, 
								  nCluster = 4, extraFileID = "",
								  burnin=1000) {


	# Fitness effects for indiv. loci
	# W.i  <-  c(w_AA, w_Aa, w_aa)
	W.i  <-  c(1, 1 - h*s, 1 - s)

	#number of simulations 
	nreps  = 10^3

	# Mutation rate and equilibrium frequencies prior to inversion
	U        <-  Us.factor*s
	u        <-  U/nTot
	qHat     <-  (U/(nTot*h*s))
	rBar.eq  <-  U*invSize/(h*s)

	# number of loci captured by inversion
	n    <-  nTot*invSize

	# Initiate chromosome matrices (rows = individuals, cols = loci spanned by inversion)
	# and populate with mutations
	Xf.mat  <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # Maternally-inherited X's in females
	Xf.pat  <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # Paternally-inherited X's in females
	Xm      <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # X chromosomes in Males
	Y       <-  matrix(rbinom((n*N/2), size=1, prob=qHat), nrow=N/2, ncol=n) # Y chromosomes

	PrR.init   <-  as.vector(table(rowSums(Y)))/(nrow(Y))
	qBar.init  <-  mean(colSums(rbind(Xf.mat, Xf.pat, Xm, Y))/(2*N))
	rBar.init  <-  mean(rowSums(Y))
#par(mfrow=c(2,2))
#hist(colSums(rbind(Xf.mat, Xf.pat, Xm, Y)), xlim=c(0,(10*qHat*N)))
#mean(colSums(rbind(Xf.mat, Xf.pat, Xm, Y))/(2*N))
#mean(rowSums(Y))
#hist(rowSums(Y))
	# burnin to equilbrium prior to inversion
	for(t in 1:burnin) {
		cat('\r', paste(sprintf("burnin generation %s / ", t), burnin, " is complete"))
		# Relative fitness & sampling probs. 
		Wf.i          <-  relFitness(mat1=Xf.mat, mat2=Xf.pat, W=W.i)
		Wm.i          <-  relFitness(mat1=Xm,     mat2=Y,      W=W.i)
		# Sample parent gametes
		moms          <-  sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wf.i/(N/2))
		dads          <-  sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wm.i/(N/2))
		# Meiosis and fertilization
		# for sons
		Xm.f1         <-  freeRecombine.sampOvules.FL(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[1:(N/2)])
		Y.f1          <-  freeRecombine.sampSperm.FL(m1=Xm, m2=Y, n=n, dads=dads[1:(N/2)], inverted=rep(FALSE, times=(N/2)), sons=TRUE)
		# for daughters
		Xf.mat.f1     <-  freeRecombine.sampOvules.FL(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[((N/2)+1):N])
		Xf.pat.f1     <-  freeRecombine.sampSperm.FL(m1=Xm, m2=Y, n=n, dads=dads[((N/2)+1):N], inverted=rep(FALSE, times=(N/2)), sons=FALSE)
		# Mutation
		Xf.mat.f1[Xf.mat.f1 == 0]  <-  rbinom(n=length(Xf.mat.f1[Xf.mat.f1 == 0]), size=1, prob=u)
		Xf.pat.f1[Xf.pat.f1 == 0]  <-  rbinom(n=length(Xf.pat.f1[Xf.pat.f1 == 0]), size=1, prob=u)
		Xm.f1[Xm.f1 == 0]          <-  rbinom(n=length(Xm.f1[Xm.f1 == 0]),         size=1, prob=u)
		Y.f1[Y.f1 == 0]            <-  rbinom(n=length(Y.f1[Y.f1 == 0]),           size=1, prob=u)
		# F1's become adults in next generation
		Xm      <-  Xm.f1
		Y       <-  Y.f1
		Xf.mat  <-  Xf.mat.f1
		Xf.pat  <-  Xf.pat.f1
	}

#hist(colSums(rbind(Xf.mat, Xf.pat, Xm, Y)), col=2)#, add=T)
#hist(colSums(rbind(Xf.mat, Xf.pat, Xm, Y)), col=2, xlim=c(0,10), breaks=100)	
#mean(colSums(rbind(Xf.mat, Xf.pat, Xm, Y))/(2*N)	)
#qHat
# colSums(Y)
#hist(rowSums(Y), col=2)#, add=T)
	rBar.calc  <-  mean(rowSums(Y))
	qBar.calc  <-  mean(colSums(rbind(Xf.mat, Xf.pat, Xm, Y))/(2*N))
	nDel.eq    <-  colSums(rbind(Xf.mat, Xf.pat, Xm, Y))

	# assign initial state chromosome matrices
	Xm.init      <-  Xm
	Y.init       <-  Y
	Xf.mat.init  <-  Xf.mat
	Xf.pat.init  <-  Xf.pat

	# Pr(r | x)
	PrR.pois  <-  dpois(c(0:max(rowSums(Y))), lambda=mean(rowSums(Y)))	
	PrR.calc  <-  as.vector(table(rowSums(Y)))/(nrow(Y))
# plot(PrR.pois, type='l',  lwd=2)
# lines(PrR.calc, col=2, lwd=2)
# sum(PrR.pois)

output  <-  c(NA, NA)
	# loop over replicate instances of new inversion
	for(i in 1:nreps){
		cat('\r', paste(sprintf("replicate %s / ", i), nreps, " is complete"))
	
		# shuffle alleles at each locus among individuals to give
		# unique initial equilibrium state
		Xf.mat  <-  apply(Xf.mat.init, MARGIN=2, function(x) x[sample(seq_along(x))])
		Xf.pat  <-  apply(Xf.pat.init, MARGIN=2, function(x) x[sample(seq_along(x))])
		Xm      <-  apply(Xm.init,     MARGIN=2, function(x) x[sample(seq_along(x))])
		Y       <-  apply(Y.init,      MARGIN=2, function(x) x[sample(seq_along(x))])
	
		# Single-copy inversion mutation occurs on randomly chosen Y chromosome 
		invertedYs  <-  sample(c(1:(N/2)), size=1)
		r  <-  sum(Y[invertedYs,])

		# Calculate relative fitness in 1st generation
		Wf.i         <-  relFitness(mat1=Xf.mat, mat2=Xf.pat, W=W.i)
		Wm.i         <-  relFitness(mat1=Xm,     mat2=Y,      W=W.i)

		wI.initial   <-  Wm.i[invertedYs]
	
	output  <-  rbind(output, c(r, wI.initial))
	}
	colnames(output)  <-  c("r", "wI")
	output  <-  data.frame(output[-1,])
	res  <-  list("Y.eq"       =  Y.init,
				  "out"        =  output,
				  "qHat"       =  qHat,
				  "qBar.init"  =  qBar.init,
				  "qBar.calc"  =  qBar.calc,
				  "rBar.eq"    =  rBar.eq,
				  "rBar.init"  =  rBar.init,
				  "rBar.calc"  =  rBar.calc,
				  "PrR.init"   =  PrR.init,
				  "PrR.pois"   =  PrR.pois,
				  "PrR.calc"   =  PrR.calc)
	return(res)
}


# output  <-  PoC_invasionFitness(h = 0.25, s = 0.01, Us.factor = 2,
# 								nTot = 10^4, N = 4*10^3, invSize = 0.1, 
# 								burnin=1000)
# 
# dat  <-  output$out
# 
# 
# wI.bars  <-  data.frame(aggregate(dat$wI, list(dat$r), mean))
# wI.bars  <-  data.frame(cbind(wI.bars, aggregate(dat$wI, list(dat$r), sd)[,2]))
# colnames(wI.bars)  <-  c("r", "wI.bar", "wI.sd")
# 
# par(mfrow=c(2,2))
# plot(wI.bars$wI.bar ~ r, data=wI.bars, ylim=c(min(wI.bar)*0.99, max(wI.bar)*1.01))
# arrows(x0=wI.bars$r, y0=wI.bars$wI.bar-wI.bars$wI.sd, x1=wI.bars$r, y1=wI.bars$wI.bar+wI.bars$wI.sd, code=3, angle=90, length=0.1)
# abline(h=1, lwd=2, col=1)
# 
# plot(output$PrR.pois, pch=21, bg='grey70', cex=1.5, ylim=c(0,max(output$PrR.pois,output$PrR.calc)*1.05))
# points(output$PrR.calc, pch=21, bg='tomato')
# 
# 
# optimisticPrFix  <-  2*(wI.bars[,2] - 1)
# optimisticPrFix[optimisticPrFix < 0]  <-  0
# optimisticPrFix  <-  optimisticPrFix*c(output$PrR.calc,rep(0.00001, times=length(optimisticPrFix)-length(output$PrR.calc)))
# plot(optimisticPrFix, ylim=c(min(optimisticPrFix), max(optimisticPrFix, (2/N))*1.05))
# abline(h=c((2/N), 0), lwd=2, col=c(2,1))
# 
# sum(optimisticPrFix)
# 
# output