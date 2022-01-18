# Simulation code, partial full-sib mating model

#' Here we present exact recursions for the frequency
#' of all possible mating types for the genetic system
#' described in the main text (1 sex determining locus, 
#' one selected locus), and code to simulate the selective
#' advantage of a rare Y-linked inversion under different
#' selection scenarios.
#' 
#' The life cycle proceeds as follows: 
#'  fertilization --> mutation --> selection --> meiosis --> fertilization


# euclidian distance between two vectors (for checking equilibrium)
eucDist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
} 


##################################
##################################
#' Genotypes
#' 
#' Females
#' x1:  XA/XA
#' x2:  XA/Xa
#' x3:  Xa/Xa
#' 
#' Males
#' y1:   XA/YA
#' y2c:  XA/Ya
#' y2cI: XA/YaI
#' y2t:  Xa/YA
#' y3:   Xa/Ya
#' y3I:  Xa/YaI

#' Mating combination frequencies
#' 
#' z1:  x1 * y1    --  XA/XA * XA/YA
#' z2:  x1 * y2c   --  XA/XA * XA/Ya
#' z3:  x1 * y2cI  --  XA/XA * XA/YaI
#' z4:  x1 * y2t   --  XA/XA * Xa/YA
#' z5:  x1 * y3    --  XA/XA * Xa/Ya
#' z6:  x1 * y3I   --  XA/XA * Xa/YaI
#' z7:  x2 * y1    --  XA/Xa * XA/YA
#' z8:  x2 * y2c   --  XA/Xa * XA/Ya
#' z9:  x2 * y2cI  --  XA/Xa * XA/YaI
#' z10: x2 * y2t   --  XA/Xa * Xa/YA
#' z11: x2 * y3    --  XA/Xa * Xa/Ya
#' z12: x2 * y3I   --  XA/Xa * Xa/YaI
#' z13: x3 * y1    --  Xa/Xa * XA/YA
#' z14: x3 * y2c   --  Xa/Xa * XA/Ya
#' z15: x3 * y2cI  --  Xa/Xa * XA/YaI
#' z16: x3 * y2t   --  Xa/Xa * Xa/YA
#' z17: x3 * y3    --  Xa/Xa * Xa/Ya
#' z18: x3 * y3I   --  Xa/Xa * Xa/YaI
#' 
#' Vector of mating combination frequencies
#c(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18)


#' Define mutation matrix
u.matrix  <-  function(u) {
	as.matrix(
		rbind(
			c((1-u)^2, ((1-u)*u)/2, 0, ((1-u)*u)/2, 0, 0, u*(1-u), (u^2)/2, 0, (u^2)/2, 0, 0, 0, 0, 0, 0, 0, 0),
			c(0, (1-u)*(1-(u/2)), 0, 0, ((1-u)*u)/2, 0, 0, u*(1-(u/2)), 0, 0, (u^2)/2, 0, 0, 0, 0, 0, 0, 0),
			c(0, 0, (1-u)*(1-(u/2)), 0, 0, (1-u)*u/2, 0, 0, u*(1-(u/2)), 0, 0, (u^2)/2, 0, 0, 0, 0, 0, 0),
			c(0, 0, 0, (1-u)*(1-(u/2)), (1-u)*u/2, 0, 0, 0, 0, u*(1-(u/2)), (u^2)/2, 0, 0, 0, 0, 0, 0, 0),
			c(0, 0, 0, 0, (1-u), 0, 0, 0, 0, 0, u, 0, 0, 0, 0, 0, 0, 0),
			c(0, 0, 0, 0, 0, (1-u), 0, 0, 0, 0, 0, u, 0, 0, 0, 0, 0, 0),
			c(0, 0, 0, 0, 0, 0, (1-(u/2))*(1-u), (1-(u/2))*u/2, 0, (1-(u/2))*u/2, 0, 0, (u/2)*(1-u), (u/2)^2, 0, (u/2)^2, 0, 0),
			c(0, 0, 0, 0, 0, 0, 0, (1-(u/2))^2, 0, 0, (1-(u/2))*u/2, 0, 0, (u/2)*(1-(u/2)), 0, 0, (u/2)^2, 0),
			c(0, 0, 0, 0, 0, 0, 0, 0, (1-(u/2))^2, 0, 0, (1-(u/2))*u/2, 0, 0, u/2*(1-(u/2)), 0, 0, (u/2)^2),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, (1-(u/2))^2, (1-(u/2))*u/2, 0, 0, 0, 0, u/2*(1-(u/2)), (u/2)^2, 0),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-(u/2)), 0, 0, 0, 0, 0, u/2, 0),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-(u/2)), 0, 0, 0, 0, 0, u/2),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-u, u/2, 0, u/2, 0, 0),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-(u/2), 0, 0, u/2, 0),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-(u/2), 0, 0, u/2),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-(u/2), u/2, 0),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
			c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
		)
	)
}

#' Define female offspring frequency array (meiosis & random mating)
off.X.array  <-  function(r) {
	as.matrix(
		rbind(
			c(1, 0, 0),
			c((1-r), r, 0),
			c(1, 0, 0),
			c(r, (1-r), 0),
			c(0, 1, 0),
			c(0, 1, 0),
			c(1/2, 1/2, 0),
			c((1-r)/2, 1/2, r/2),
			c(1/2, 1/2, 0),
			c(r/2, 1/2, (1-r)/2),
			c(0, 1/2, 1/2),
			c(0, 1/2, 1/2),
			c(0, 1, 0),
			c(0, (1-r), r),
			c(0, 1, 0),
			c(0, r, (1-r)),
			c(0, 0, 1),
			c(0, 0, 1)
			)
		)
}


#' Define male offspring frequency array (meiosis & random mating)
off.Y.array  <-  function(r) {
	as.matrix(
		rbind(
			  c(1, 0, 0, 0, 0, 0),
			  c(r, (1-r), 0, 0, 0, 0),
			  c(0, 0, 1, 0, 0, 0),
			  c((1-r), r, 0, 0, 0, 0),
			  c(0, 1, 0, 0, 0, 0),
			  c(0, 0, 1, 0, 0, 0),
			  c(1/2, 0, 0, 1/2, 0, 0),
			  c(r/2, (1-r)/2, 0, r/2, (1-r)/2, 0),
			  c(0, 0, 1/2, 0, 0, 1/2),
			  c((1-r)/2, r/2, 0, (1-r)/2, r/2, 0),
			  c(0, 1/2, 0, 0, 1/2, 0),
			  c(0, 0, 1/2, 0, 0, 1/2),
			  c(0, 0, 0, 1, 0, 0),
			  c(0, 0, 0, r, (1-r), 0),
			  c(0, 0, 0, 0, 0, 1),
			  c(0, 0, 0, (1-r), r, 0),
			  c(0, 0, 0, 0, 1, 0),
			  c(0, 0, 0, 0, 0, 1)
			)
		)
}

#' Define sib-mating array (frequency of mating types resulting from 
#' full-sib matings)
#' 
#' illustration of how sib.Z.array was calculated:
#' test  <-  matrix(0, ncol=18,nrow=18)
#' r  <-  0.01
#' off.X  <-  off.X.array(r)
#' off.Y  <-  off.Y.array(r)
#' for(i in 1:nrow(off.X)) {
#' 	for(j in 1:ncol(off.X)) {
#' 		for(k in 1:ncol(off.Y)) {
#' 			test[i,(((j-1)*ncol(off.Y))+k)]  <-  off.X[i,j] * off.Y[i,k]
#' 		}
#' 	}
#' }
#' t(test) == sib.Z.array(r)

sib.Z.array  <-  function(r) {
	as.matrix(
		rbind(
			  c(1, (1-r)*r, 0, r*(1-r), 0, 0, 1/4, ((1-r)*r)/4, 0, (r*(1-r))/4, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, (1-r)^2, 0, r^2, 0, 0, 0, (1-r)^2/4, 0, r^2/4, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 1, 0, 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 1/4, ((1-r)*r)/4, 0, (r*(1-r))/4, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, (1-r)^2/4, 0, r^2/4, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, r^2, 0, (1-r)^2, 0, 0, 1/4, r/4, 0, (1-r)/4, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, r*(1-r), 0, (1-r)*r , 1 , 0, 0, (1-r)/4, 0, r/4, 1/4, 0, 0, 0 , 0, 0, 0 , 0),
			  c(0, 0, 0, 0, 0, 1, 0, 0, 1/4, 0, 0, 1/4, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 0, 0 , 0, 0, 1/4, r/4, 0, (1-r)/4, 0, 0, 1 ,  (1-r)*r, 0, r*(1-r) , 0 , 0),
			  c(0, 0, 0, 0, 0, 0, 0, (1-r)/4, 0, r/4, 1/4, 0, 0, (1-r)^2, 0, r^2, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, 0, 1/4, 0, 0, 1/4, 0, 0, 1, 0, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, r^2/4, 0, (1-r)^2/4, 0, 0, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, (r*(1-r))/4, 0, ((1-r)*r)/4, 1/4, 0, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, r^2/4, 0, (1-r)^2/4, 0, 0, 0, r^2, 0, (1-r)^2, 0, 0),
			  c(0, 0, 0, 0, 0, 0, 0, (r*(1-r))/4, 0, ((1-r)*r)/4, 1/4, 0, 0, r*(1-r), 0, (1-r)*r, 1, 0),
			  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 0, 1)
			)
		)
}




#' Fitness Vectors
#' Wf  <-  c(wx1, wx2, wx3)
#' Wm  <-  c(wy1, wy2c, wy2cI, wy2t, wy3, wy3I)

#' Define vector of fitness expressions for each mating combination
Wvec  <-  function(Wf, Wm) {
	Wvec  <-  c()
	for(i in 1:length(Wf)) {
		for(j in 1:length(Wm)) {
			Wvec  <-  c(Wvec,Wf[i]*Wm[j])
		}
	}
	return(Wvec)
}


# z.0  <-  c(1/12,1/12,0, 1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0)
#' Recursion giving frequency of each mating type in the next generation
#' 
Z.pr  <-  function(z, u=10^-5, Wf, Wm, r=1/2, beta=0.15, ...) {

	# Step through the life-cycle
	# [1] Mutation 
	z.m  <-  z %*% u.matrix(u)

	# [2] [3] Selection, Meiosis, Mating,
	# Outcrossing
	# female offspring genotypic frequencies
#	X  <-  Wf * z.m %*% off.X.array(r)
	X  <-  z.m %*% off.X.array(r) * Wf
	X  <-  X/norm(t(X), type="1")
	# male offspring genotypic frequencies
#	Y  <-  Wm * z.m %*% off.Y.array(r)
	Y  <-  z.m %*% off.Y.array(r) * Wm
	Y  <-  Y/norm(t(Y), type="1")
	# Outcross mating-type frequencies
	z.pr.out  <-  c()
	for(i in 1:length(X)) {
		for(j in 1:length(Y)) {
			z.pr.out  <-  c(z.pr.out, X[i]*Y[j])
		}
	}
	# Sib-mating
	z.pr.sib  <-  Wvec(Wf=Wf, Wm=Wm) * z.m %*% t(sib.Z.array(r))
#	z.pr.sib  <-  Wvec(Wf=Wf, Wm=Wm) * z.m %*% sib.Z.array(r)
	z.pr.sib  <-  z.pr.sib/norm(t(z.pr.sib), type="1")
	
	# [4] Overall mating-type frequencies in next generation
	z.pr  <-  (1 - beta)*z.pr.out + beta*z.pr.sib

	return(as.vector(z.pr))
}


#' Pared-down forward simulation of recursion equation to equilibrium
findEq  <-  function(threshold=10^-8, 
					 z.init, u=10^-5, Wf, Wm, r=1/2, beta=0.15, ...) {

	# Define delta values to check eq.
	delta  <-  1

	# Forward simulation
	t  <-  1
	z  <-  z.init
	while(delta >= threshold) {
		z.next      <-  Z.pr(z=z, u=u, Wf=Wf, Wm=Wm, r=r, beta=beta)
		delta       <-  eucDist(z.next, z)
		z           <-  z.next
		t           <-  t + 1
	}
	return(z)
}


#' Forward simulation of recursion equation to equilibrium
fwdRecSim  <-  function(generations = 1000, threshold=10^-8, 
						z.init, u=10^-5, Wf, Wm, r=1/2, beta=0.15, ...) {

	# Matrix to store frequencies
	outMat  <-  matrix(0, ncol=length(z.init), nrow=generations)
	# Define delta values to check eq.
	delta  <-  1

	# Forward simulation
	t  <-  1
	z  <-  z.init
	while(delta >= threshold && t <= generations) {
		z.next      <-  Z.pr(z=z, u=u, Wf=Wf, Wm=Wm, r=r, beta=beta)
		delta       <-  eucDist(z.next, z)
		z           <-  z.next
		outMat[t,]  <-  z
		t           <-  t + 1
	}

	# trim output Matrix
	outMat  <-  outMat[rowSums(outMat) != 0,]

	# Calculate frequencies of interest
	Xfreqs  <-  cbind((rowSums(outMat[,1:6])),
					  (rowSums(outMat[,7:12])),
					  (rowSums(outMat[,13:18]))
					  )
	Yfreqs  <-  cbind((rowSums(outMat[,c(1,7,13)])),
					  (rowSums(outMat[,c(2,8,14)])),
					  (rowSums(outMat[,c(3,9,15)])),
					  (rowSums(outMat[,c(4,10,16)])),
					  (rowSums(outMat[,c(5,11,17)])),
					  (rowSums(outMat[,c(6,12,18)]))
					  )
	aFreqX   <-  (Xfreqs[,2]/2) + Xfreqs[,3]
	aFreqY   <-  ((Yfreqs[,2] + Yfreqs[,3] + Yfreqs[,4])/2) + Yfreqs[,5] + Yfreqs[,6]
	aFreq    <-  (aFreqX + aFreqY)/2
	invFreq  <-  Yfreqs[,3] + Yfreqs[,6]

	results  <-  list(
					  "outMat"   =  outMat,
					  "Xfreqs"   =  Xfreqs,
					  "Yfreqs"   =  Yfreqs,
					  "aFreqX"   =  aFreqX,
					  "aFreqY"   =  aFreqY,
					  "aFreq"    =  aFreq,
					  "invFreq"  =  invFreq,
					  "z.eq"     =  outMat[nrow(outMat),]
					  )
	return(results)
}



#' Fast forward simulation of recursion equation to equilibrium
#' ends simulation when inversion reaches frequency of 0.01
fwdRecSimFast  <-  function(generations = 1000, invCutoff=0.01,
						z.init, u=10^-5, Wf, Wm, r=1/2, beta=0.15, ...) {

	# Matrix to store frequencies
	outMat  <-  matrix(0, ncol=length(z.init), nrow=generations)
	# Define delta values to check eq.
	invFreqTest  <-  0.001

	# Forward simulation
	t  <-  1
	z  <-  z.init
	while(invFreqTest <= invCutoff && invFreqTest > 0.0001 && t <= generations) {
		z.next      <-  Z.pr(z=z, u=u, Wf=Wf, Wm=Wm, r=r, beta=beta)
		invFreqTest <-  sum(z.next[c(3,6,9,12,15,18)])
		z           <-  z.next
		outMat[t,]  <-  z
		t           <-  t + 1
	}

	# trim output Matrix
	outMat  <-  outMat[rowSums(outMat) != 0,]

	# Calculate frequencies of interest
	Xfreqs  <-  cbind((rowSums(outMat[,1:6])),
					  (rowSums(outMat[,7:12])),
					  (rowSums(outMat[,13:18]))
					  )
	Yfreqs  <-  cbind((rowSums(outMat[,c(1,7,13)])),
					  (rowSums(outMat[,c(2,8,14)])),
					  (rowSums(outMat[,c(3,9,15)])),
					  (rowSums(outMat[,c(4,10,16)])),
					  (rowSums(outMat[,c(5,11,17)])),
					  (rowSums(outMat[,c(6,12,18)]))
					  )
	aFreqX   <-  (Xfreqs[,2]/2) + Xfreqs[,3]
	aFreqY   <-  ((Yfreqs[,2] + Yfreqs[,3] + Yfreqs[,4])/2) + Yfreqs[,5] + Yfreqs[,6]
	aFreq    <-  (aFreqX + aFreqY)/2
	invFreq  <-  Yfreqs[,3] + Yfreqs[,6]

	results  <-  list(
					  "outMat"   =  outMat,
					  "Xfreqs"   =  Xfreqs,
					  "Yfreqs"   =  Yfreqs,
					  "aFreqX"   =  aFreqX,
					  "aFreqY"   =  aFreqY,
					  "aFreq"    =  aFreq,
					  "invFreq"  =  invFreq,
					  "z.eq"     =  outMat[nrow(outMat),]
					  )
	return(results)
}




#' Loop over selection & Sib-mating rates to produce 
#' similar plot as Charlesworth & Wall (1999) for
#' heterozygote advantage
makeDataHetAdvFig  <-  function(generations = 10000, threshold=10^-8, 
								s.vals, t.vals, r = 1/2) {
	# data storage
	invSelAdv  <-  c()

	# parameters
	betas  <-  c(0:9)/10
	u      <-  0

	# Define arbitrary initial
	# frequencies in absence of inversion
	z.init  <-  c(1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0)

	# progress counter
	count  <-  0
	total  <-  length(betas)*length(s.vals)
	cat('\r', paste("Progress:", round(100*(count/total)), "% complete"))

	# loop over parameters
	for(i in 1:length(betas)) {

		for(j in 1:length(s.vals)) {

			# Define fitness expressions for heterozygote Adv.
			sf  <- s.vals[j]
			sm  <- s.vals[j]
			tf  <- t.vals[j]
			tm  <- t.vals[j]
			Wf  <-  c(1 - sf, 1, 1 - tf)
			Wm  <-  c(1-sm, 1, 1, 1, 1-tm, 1-tm)

			# Find eq. z values in absence of inversion	
			z.eq  <-  findEq(threshold=10^-8, z.init=z.init, u=u, Wf=Wf, Wm=Wm, r=r, beta=betas[i])

			# Introduce new inversion at 0.001 freq. Following Charlesworth & Wall,
			# introduce onto Xa/Ya male, at equal frequencies.
			z.intro     <-  z.eq
			for (m in c(6,12,18)) {
				z.intro[(m-1)]  <-  z.intro[(m-1)] - (0.001/4)
				z.intro[m]      <-  z.intro[m] + (0.001/4)
			}
#			z.intro[8]  <-  z.intro[8] - 0.001
#			z.intro[9]  <-  z.intro[9] + 0.001
			res  <-  fwdRecSimFast(generations = 10000, invCutoff=0.01, 
								z.init=z.intro, u=u, Wf=Wf, Wm=Wm, r=r, beta=betas[i])

			# Calculate asymp. rate of increase in log() frequency of inversion
			rareLogInv  <-  log(res$invFreq)
			diff  <-  c()
			for(k in 10:length(rareLogInv)) {
				diff  <-  c(diff, rareLogInv[k] - rareLogInv[k-1])
			}
			invSelAdv  <-  c(invSelAdv, mean(diff))
			# update progress counter
			count  <-  count + 1
			cat('\r', paste("Progress:", round(100*(count/total)), "% complete"))

		}
	}
	beta  <-  rep(betas, each=length(s.vals))
	s     <-  rep(s.vals, times=length(betas))
	t     <-  rep(t.vals, times=length(betas))

	results  <-  as.data.frame(cbind(beta, s, t, invSelAdv))

	# Export data file for plotting
	filename  <-  paste("./data/hetAdvSimData_r", r,".csv", sep="")
	write.csv(results, file=filename)
}




#' Loop over selection & Sib-mating rates to produce 
#' similar plot as Charlesworth & Wall (1999) for
#' heterozygote advantage
makeDataShelteringFig  <-  function(generations = 10000, threshold=10^-8, 
								    s.vals, h.vals, r.vals, u = 10^-5) {
	# data storage
	invSelAdv  <-  c()

	# parameters
	betas  <-  c(0:9)/10

	# Define arbitrary initial
	# frequencies in absence of inversion
	z.init  <-  c(1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0)

	# Progress counter
	count  <-  0
	total  <-  length(betas)*length(r.vals)*length(h.vals)*length(s.vals)
	cat('\r', paste("Progress:", round(100*(count/total)), "% complete"))

	# loop over parameters
	for(i in 1:length(betas)) {
		for(j in 1:length(r.vals)) {
			for(k in 1:length(h.vals)) {
				for(l in 1:length(s.vals)) {

					# Define fitness expressions for heterozygote Adv.
					hf  <- h.vals[k]
					hm  <- h.vals[k]
					sf  <- s.vals[l]
					sm  <- s.vals[l]
					Wf  <-  c(1, 1 - hf*sf, 1 - sf)
					Wm  <-  c(1, 1 - hm*sm, 1 - hm*sm, 1 - hm*sm, 1 - sm, 1 - sm)

					# Find eq. z values in absence of inversion	
					z.eq  <-  findEq(threshold=10^-8, z.init=z.init, u=u, Wf=Wf, Wm=Wm, r=r.vals[j], beta=betas[i])
	
					# Introduce new inversion on a randomly chosen Ya chromosome at 0.001 freq.
					z.intro            <-  z.eq
					newMut             <-  c(2,5,8,11,14,17)[rmultinom(n=1, size=1, prob=z.eq[c(2,5,8,11,14,17)])]
					z.intro[newMut]    <-  z.intro[newMut] - 0.001
					z.intro[newMut+1]  <-  z.intro[newMut+1] + 0.001
					# Run simulation
					res  <-  fwdRecSimFast(generations = 10000, invCutoff=0.01, 
											z.init=z.intro, u=u, Wf=Wf, Wm=Wm, r=r.vals[j], beta=betas[i])
			
					# Calculate asymp. rate of increase in log() frequency of inversion
					rareLogInv  <-  log(res$invFreq)
					diff  <-  c()
					for(t in 2:length(rareLogInv)) {
						diff  <-  c(diff, rareLogInv[t] - rareLogInv[t-1])
					}
					invSelAdv  <-  c(invSelAdv, mean(diff))

					# progress counter
					count  <-  count + 1
					cat('\r', paste("Progress:", round(100*(count/total)), "% complete"))
				}
			}
		}
	}

	# Parameter vectors for plotting
	beta  <-  rep(betas, each=length(r.vals)*length(h.vals)*length(s.vals))
	r     <-  rep(rep(r.vals, each=length(h.vals)*length(s.vals)), times=length(betas))
	h     <-  rep(rep(h.vals, each=length(s.vals)), times=length(betas)*length(r.vals))
	s     <-  rep(s.vals, times=length(betas)*length(r.vals)*length(h.vals))

	# Make Dataframe
	results  <-  as.data.frame(cbind(beta, r, h, s, invSelAdv))

	# Export data file for plotting
	filename  <-  paste("./data/shelteringData.csv", sep="")
	write.csv(results, file=filename)
}






#' Loop over selection & Sib-mating rates to produce 
#' similar plot as Charlesworth & Wall (1999) for
#' heterozygote advantage
makeDataShelteringSexSpecific  <-  function(generations = 10000, threshold=10^-8, 
								    		sf.vals, sm.vals, h.vals, r.vals, u = 10^-5) {
	# data storage
	invSelAdv  <-  c()

	# parameters
	betas  <-  c(0:9)/10

	# Define arbitrary initial
	# frequencies in absence of inversion
	z.init  <-  c(1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0,1/12,1/12,0)

	# Progress counter
	count  <-  0
	total  <-  length(betas)*length(r.vals)*length(h.vals)*length(sf.vals)
	cat('\r', paste("Progress:", round(100*(count/total)), "% complete"))

	# loop over parameters
	for(i in 1:length(betas)) {
		for(j in 1:length(r.vals)) {
			for(k in 1:length(h.vals)) {
				for(l in 1:length(sf.vals)) {

					# Define fitness expressions for heterozygote Adv.
					hf  <- h.vals[k]
					hm  <- h.vals[k]
					sf  <- sf.vals[l]
					sm  <- sm.vals[l]
					Wf  <-  c(1, 1 - hf*sf, 1 - sf)
					Wm  <-  c(1, 1 - hm*sm, 1 - hm*sm, 1 - hm*sm, 1 - sm, 1 - sm)

					# Find eq. z values in absence of inversion	
					z.eq  <-  findEq(threshold=10^-8, z.init=z.init, u=u, Wf=Wf, Wm=Wm, r=r.vals[j], beta=betas[i])
	
					# Introduce new inversion on a randomly chosen Ya chromosome at 0.001 freq.
					z.intro            <-  z.eq
					newMut             <-  c(2,5,8,11,14,17)[rmultinom(n=1, size=1, prob=z.eq[c(2,5,8,11,14,17)])]
					z.intro[newMut]    <-  z.intro[newMut] - 0.001
					z.intro[newMut+1]  <-  z.intro[newMut+1] + 0.001
					# Run simulation
					res  <-  fwdRecSimFast(generations = 10000, invCutoff=0.01, 
											z.init=z.intro, u=u, Wf=Wf, Wm=Wm, r=r.vals[j], beta=betas[i])
			
					# Calculate asymp. rate of increase in log() frequency of inversion
					rareLogInv  <-  log(res$invFreq)
					diff  <-  c()
					for(t in 2:length(rareLogInv)) {
						diff  <-  c(diff, rareLogInv[t] - rareLogInv[t-1])
					}
					invSelAdv  <-  c(invSelAdv, mean(diff))

					# progress counter
					count  <-  count + 1
					cat('\r', paste("Progress:", round(100*(count/total)), "% complete"))
				}
			}
		}
	}

	# Parameter vectors for plotting
	beta  <-  rep(betas, each=length(r.vals)*length(h.vals)*length(sf.vals))
	r     <-  rep(rep(r.vals, each=length(h.vals)*length(sf.vals)), times=length(betas))
	h     <-  rep(rep(h.vals, each=length(sf.vals)), times=length(betas)*length(r.vals))
	sf    <-  rep(sf.vals, times=length(betas)*length(r.vals)*length(h.vals))
	sm    <-  rep(sm.vals, times=length(betas)*length(r.vals)*length(h.vals))

	# Make Dataframe
	results  <-  as.data.frame(cbind(beta, r, h, sf, sm, invSelAdv))

	# Export data file for plotting
	filename  <-  paste("./data/shelteringDataSexSpecific.csv", sep="")
	write.csv(results, file=filename)
}


