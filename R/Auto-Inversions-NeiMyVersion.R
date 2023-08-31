#' Autosomal Recursions based on Eqs(16--23) of Nei, Kojima, Schaffer (1967).

qstar  <-  function(q.0, U, x, s, h) {
     (q.0*exp((U*x)/(1 - exp(-s*h)))) / (1 - q.0*(1 - exp((U*x)/(1 - exp(-h*s)))))
 }

I.genSol  <-  function(q.0, U, x, s, h, t) {
     (exp(((1 - exp(-h*s*t))*U*x) / (1 - exp(-h*s)))*q.0)/(1 - q.0 + exp(((1 - exp(-h*s*t))*U*x) / (1 - exp(-h*s)))*q.0)
 }


#' Nei et al. (1967) use x_t for inversion frequency, but since I've
#' already been using x for inversion size, I'm using I.
I.Prime  <-  function(qN, qI, qD, I, h, s, n, r, ...) {
	QI   <-  qI/I
	PI   <-  1 - QI
	QN   <-  qN/(1 - I)
	PN   <-  1 - QN
	QND  <-  qD/(1 - I)
	PND  <-  1 - QND
	if(I == 0) {
		QI  <-  0
		PI  <-  0
	}
	if(I == 1) {
		QN   <-  0
		PN   <-  0
		QND  <-  0
		PND  <-  0
	}
	( (I^2)*(1 - s)^r*(1 - s*(2*h*QI*PI + QI^2))^(n-r) +   I*(1 - I)*(1 - s*(  h*PND + QND))^r*(1 - s*(h*(QN*PI + PN*QI) + QN*QI))^(n-r) ) / 
	( (I^2)*(1 - s)^r*(1 - s*(2*h*QI*PI + QI^2))^(n-r) + 2*I*(1 - I)*(1 - s*(2*h*PND + QND))^r*(1 - s*(h*(QN*PI + PN*QI) + QN*QI))^(n-r) + ((1 - I)^2)*(1 - s*(2*h*QND*PND + QND^2))^r*(1 - s*(2*h*QN*PN + QN^2))^(n-r) )
}

qI.Prime.W  <-  function(qN, pN, qI, pI, I, Iprime, u, h, s, ...) {
	qNm  <-  qN + pN*u
	pNm  <-  pN*(1 - u)
	qIm  <-  qI + pI*u
	pIm  <-  pI*(1 - u)

	(qIm*(1 - s*(h*(pIm + pNm) + qIm + qNm) ) / 
			(I - qIm*s*(h*(pIm + pNm) + qIm + qNm) )) * Iprime
}

pI.Prime.W  <-  function(Iprime, qIprime.W) {
	Iprime - qIprime.W
}

qN.Prime.W  <-  function(qN, pN, qI, pI, I, Iprime, u, h, s, ...) {
	qNm  <-  qN + pN*u
	pNm  <-  pN*(1 - u)
	qIm  <-  qI + pI*u
	pIm  <-  pI*(1 - u)

	(qNm*(1 - s*(h*(pNm + pIm) + qNm + qIm) ) / 
			(1 - I - qNm*s*(h*(pNm + pIm) + qNm + qIm) )) * (1 - Iprime)
}

pN.Prime.W  <-  function(Iprime, qNprime.W) {
	1 - Iprime - qNprime.W
}

qN.Prime.D  <-  function(qN.D, pN.D, I, Iprime, u, h, s, ...) {
	qNm.D  <-  qN.D + pN.D*u
	pNm.D  <-  pN.D*(1 - u)

	(qNm.D*(1 - s*(h*pNm.D + I + qNm.D) ) / 
			(1 - I - qNm.D*s*(h*pNm.D + I + qNm.D) )) * (1 - Iprime)
}

pN.Prime.D  <-  function(Iprime, qNprime.D) {
	1 - Iprime - qNprime.D
}


autoInvRelFit  <-  function(qN, qI, qD, I, h, s, n, r, ...) {
	QI   <-  qI/I
	PI   <-  1 - QI
	QN   <-  qN/(1 - I)
	PN   <-  1 - QN
	QND  <-  qD/(1 - I)
	PND  <-  1 - QND

	I.prime  <-  ( (I^2)*(1 - s)^r*(1 - s*(2*h*QI*PI + QI^2))^(n-r) +   I*(1 - I)*(1 - s*(  h*PND + QND))^r*(1 - s*(h*(QN*PI + PN*QI) + QN*QI))^(n-r) ) / 
	( (I^2)*(1 - s)^r*(1 - s*(2*h*QI*PI + QI^2))^(n-r) + 2*I*(1 - I)*(1 - s*(2*h*PND + QND))^r*(1 - s*(h*(QN*PI + PN*QI) + QN*QI))^(n-r) + ((1 - I)^2)*(1 - s*(2*h*QND*PND + QND^2))^r*(1 - s*(2*h*QN*PN + QN^2))^(n-r) )

	I.prime/I
}

#' Function to calculate euclidean distance between two vectors
#' Convenient for checking if equilibrium has been reached
eucDist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
} 

#' Function to find equilibrium frequencies in the absence of inverions
findEqFullRec  <-  function(u, h, s, n, r, qHat_init, ...) {

	# Empty Frequency Vectors
	I.t     <- c()
	qt.I.W  <- c()
	pt.I.W  <- c()
	qt.W    <- c()
	pt.W    <- c()
	qt.D    <- c()
	pt.D    <- c()

	I.0  <-  10^-10
	qHat  <-  qHat_init

	# First generation (all loci at equilibrium frequencies when inversion arises)
	I.t[1]     <-  I.Prime(qN=qHat*(1 - I.0), qI=0, qD=qHat*(1 - I.0), I=I.0, h=h, s=s, n=n, r=r)
	qt.I.W[1]  <-  qI.Prime.W(qN=qHat*(1 - I.0), pN=(1 - I.0 - qHat*(1 - I.0)), qI=0, pI=(I.0 - qHat*I.0), I=I.0, Iprime=I.0, u=u, h=h, s=s)
	pt.I.W[1]  <-  pI.Prime.W(Iprime=I.0, qIprime.W=qt.I.W[1])
	qt.W[1]    <-  qN.Prime.W(qN=qHat*(1 - I.0), pN=(1 - I.0 - qHat*(1 - I.0)), qI=0, pI=(I.0 - qHat*I.0), I=I.0, Iprime=I.0, u=u, h=h, s=s)
	pt.W[1]    <-  pN.Prime.W(Iprime=I.0, qNprime.W=qt.W[1])
	qt.D[1]    <-  qN.Prime.D(qN.D=qHat*(1 - I.0), pN.D=(1 - I.0 - qHat*(1 - I.0)), I=I.0, Iprime=I.0, u=u, h=h, s=s)
	pt.D[1]    <-  pN.Prime.D(Iprime=I.0, qNprime.D=qt.D[1])

	# Subsequent generations
	i=2
	while(i < 10^5) {
		I.t[i]     <-  I.Prime(qN=qt.W[i-1], qI=qt.I.W[i-1], qD=qt.D[i-1], I=I.0, h=h, s=s, n=n, r=r)
		qt.I.W[i]  <-  qI.Prime.W(qN=qt.W[i-1], pN=pt.W[i-1], qI=qt.I.W[i-1], pI=pt.I.W[i-1], I=I.0, Iprime=I.0, u=u, h=h, s=s)
		pt.I.W[i]  <-  pI.Prime.W(Iprime=I.0, qIprime.W=qt.I.W[i])
		qt.W[i]    <-  qN.Prime.W(qN=qt.W[i-1], pN=pt.W[i-1], qI=qt.I.W[i-1], pI=pt.I.W[i-1], I=I.0, Iprime=I.0, u=u, h=h, s=s)
		pt.W[i]    <-  pN.Prime.W(Iprime=I.0, qNprime.W=qt.W[i])
		qt.D[i]    <-  qN.Prime.D(qN.D=qt.D[i-1], pN.D=pt.D[i-1], I=I.0, Iprime=I.0, u=u, h=h, s=s)
		pt.D[i]    <-  pN.Prime.D(Iprime=I.0, qNprime.D=qt.D[i])
		diff       <-  eucDist(qt.W[i], qt.W[i-1])
		i  <-  i + 1
	}
	# Return results as df
	results  <-  data.frame(
							"qt.W"    =  qt.W,
							"pt.W"    =  pt.W,
							"qt.D"    =  qt.D,
							"pt.D"    =  pt.D,
							"qt.I.W"  =  qt.I.W,
							"pt.I.W"  =  pt.I.W,
							"I.t"     =  I.t
							)
	return(results[nrow(results),])
}


#' Function to generate deterministic frequency dynamics 
#' for AUTOSOMAL inversions
#' Used to generate Fig.C1, S1, S3, S10-S15
makeAutoDeterministicFigSimData_Nei  <-  function(U = 0.02, r = 0, x = 0.2, h = 0.25, s = 0.01, I.0 = (2/5000), nTot=10^4, generations, ...) {

	# Parameters
	u     <-  U/nTot
	if(h == 0) {
		qHat  <-  sqrt(u/s)  # for h = 0		
	} else{qHat  <-  (U/(nTot*h*s)) } # for h > 0
	n     <-  nTot*x

	# find equilibrium frequencies prior to inversion mutation
	init_eqs  <-  findEqFullRec(u=u, h=h, s=s, n=n, r=r, qHat_init=qHat)
	qHat  <-  init_eqs$qt.W

	# Empty Frequency Vectors
	I.t     <- c()
	qt.I.W  <- c()
	pt.I.W  <- c()
	qt.W    <- c()
	pt.W    <- c()
	qt.D    <- c()
	pt.D    <- c()

	# First generation (all loci at equilibrium frequencies when inversion arises)
	I.t[1]     <-  I.Prime(qN=qHat*(1 - I.0), qI=0, qD=qHat*(1 - I.0), I=I.0, h=h, s=s, n=n, r=r)
	qt.I.W[1]  <-  qI.Prime.W(qN=qHat*(1 - I.0), pN=(1 - I.0 - qHat*(1 - I.0)), qI=0, pI=(I.0 - qHat*I.0), I=I.0, Iprime=I.t[1], u=u, h=h, s=s)
	pt.I.W[1]  <-  pI.Prime.W(Iprime=I.t[1], qIprime.W=qt.I.W[1])
	qt.W[1]    <-  qN.Prime.W(qN=qHat*(1 - I.0), pN=(1 - I.0 - qHat*(1 - I.0)), qI=0, pI=(I.0 - qHat*I.0), I=I.0, Iprime=I.t[1], u=u, h=h, s=s)
	pt.W[1]    <-  pN.Prime.W(Iprime=I.t[1], qNprime.W=qt.W[1])
	qt.D[1]    <-  qN.Prime.D(qN.D=qHat*(1 - I.0), pN.D=(1 - I.0 - qHat*(1 - I.0)), I=I.0, Iprime=I.t[1], u=u, h=h, s=s)
	pt.D[1]    <-  pN.Prime.D(Iprime=I.t[1], qNprime.D=qt.D[1])

	# Subsequent generations
	i=2
	while(i < generations+1) {
		I.t[i]     <-  I.Prime(qN=qt.W[i-1], qI=qt.I.W[i-1], qD=qt.D[i-1], I=I.t[i-1], h=h, s=s, n=n, r=r)
		qt.I.W[i]  <-  qI.Prime.W(qN=qt.W[i-1], pN=pt.W[i-1], qI=qt.I.W[i-1], pI=pt.I.W[i-1], I=I.t[i-1], Iprime=I.t[i], u=u, h=h, s=s)
		pt.I.W[i]  <-  pI.Prime.W(Iprime=I.t[i], qIprime.W=qt.I.W[i])
		qt.W[i]    <-  qN.Prime.W(qN=qt.W[i-1], pN=pt.W[i-1], qI=qt.I.W[i-1], pI=pt.I.W[i-1], I=I.t[i-1], Iprime=I.t[i], u=u, h=h, s=s)
		pt.W[i]    <-  pN.Prime.W(Iprime=I.t[i], qNprime.W=qt.W[i])
		qt.D[i]    <-  qN.Prime.D(qN.D=qt.D[i-1], pN.D=pt.D[i-1], I=I.t[i-1], Iprime=I.t[i], u=u, h=h, s=s)
		pt.D[i]    <-  pN.Prime.D(Iprime=I.t[i], qNprime.D=qt.D[i])
		i  <-  i + 1
	}
	# Return results as df
	results  <-  data.frame(
							"qt.W"    =  qt.W,
							"pt.W"    =  pt.W,
							"qt.D"    =  qt.D,
							"pt.D"    =  pt.D,
							"qt.I.W"  =  qt.I.W,
							"pt.I.W"  =  pt.I.W,
							"I.t"     =  I.t
							)
	return(results)
}



##########################################################################
#' Wright-Fisher Simulations to estimate Inversion fixation probabilities


# Function to estimate Pr(fix | x) for different AUTOSOMAL inversion sizes (x)
makeDataExactAutoPrFixInvSize_Nei  <-  function(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
											nTot = 10^4, N = 10^4, Nfname = "") {

	# Containers
	PrFix      <-  c()
	rFixedInv  <-  c()

	# Empty Frequency Vectors
	I.t     <- c()
	qt.I.W  <- c()
	pt.I.W  <- c()
	qt.W    <- c()
	pt.W    <- c()
	qt.D    <- c()
	pt.D    <- c()

	# Containers for # del mutations on
	# fixed inversions 
	rFixedInvTab  <-  data.frame()
	rInvSize      <-  c()
	rNs           <-  c()
	rUs           <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	I.0   <-  1/(2*N)

	#number of simulations 
	sims  = 200*N


	# Loop over Us factor values
	for(j in 1:length(U.vals)) {
		# mutation rate & initial qHat
		U     <-  U.vals[j]
		u     <-  U/nTot
		qHat  <-  (U/(nTot*h*s))
		init_eqs  <-  findEqFullRec(u=u, h=h, s=s, n=1, r=0, qHat_init=qHat)
		qHat  <-  init_eqs$qt.W

		# Loop over inversion size
		for(k in 1:length(invSize)) {
			
			# Draw random values for # loci captured by inversion
			ns   <-  rpois(sims, lambda=nTot*invSize[k])
			# Draw random value for # del. mutations captured by inversion (r) given x
			rs  <-  c()
			for(m in 1:length(ns)) {
				rs[m]  <-  sum(rbinom(ns[m], size=1, prob=qHat))	
			}

			# counter for fixations
			fix   = 0
			rFixedInv  <-  c()

			# Loop over replicate simulations
			for(l in 1:sims){

				# take randomly drawn n and r values
				n  <-  ns[l]
				r  <-  rs[l]

				# Assign frequencies
				# Note implicit assumption of equal initial
				# frequencies in XOv, XSp, Y chromosomes
				I.t[1]     <-  I.Prime(qN=qHat*(1 - I.0), qI=0, qD=qHat*(1 - I.0), I=I.0, h=h, s=s, n=n, r=r)
				qt.I.W[1]  <-  qI.Prime.W(qN=qHat*(1 - I.0), pN=(1 - I.0 - qHat*(1 - I.0)), qI=0, pI=(I.0 - qHat*I.0), I=I.0, Iprime=I.t[1], u=u, h=h, s=s)
				pt.I.W[1]  <-  pI.Prime.W(Iprime=I.t[1], qIprime.W=qt.I.W[1])
				qt.W[1]    <-  qN.Prime.W(qN=qHat*(1 - I.0), pN=(1 - I.0 - qHat*(1 - I.0)), qI=0, pI=(I.0 - qHat*I.0), I=I.0, Iprime=I.t[1], u=u, h=h, s=s)
				pt.W[1]    <-  pN.Prime.W(Iprime=I.t[1], qNprime.W=qt.W[1])
				qt.D[1]    <-  qN.Prime.D(qN.D=qHat*(1 - I.0), pN.D=(1 - I.0 - qHat*(1 - I.0)), I=I.0, Iprime=I.t[1], u=u, h=h, s=s)
				pt.D[1]    <-  pN.Prime.D(Iprime=I.t[1], qNprime.D=qt.D[1])


				# Run forward simulation
				while(I.t[1]*(1 - I.t[1]) > 0) {
					#expected frequency after selection
					I.t[2]     <-  I.Prime(qN=qt.W[1], qI=qt.I.W[1], qD=qt.D[1], I=I.t[1], h=h, s=s, n=n, r=r)
					qt.I.W[2]  <-  qI.Prime.W(qN=qt.W[1], pN=pt.W[1], qI=qt.I.W[1], pI=pt.I.W[1], I=I.t[1], Iprime=I.t[2], u=u, h=h, s=s)
					pt.I.W[2]  <-  pI.Prime.W(Iprime=I.t[2], qIprime.W=qt.I.W[2])
					qt.W[2]    <-  qN.Prime.W(qN=qt.W[1], pN=pt.W[1], qI=qt.I.W[1], pI=pt.I.W[1], I=I.t[1], Iprime=I.t[2], u=u, h=h, s=s)
					pt.W[2]    <-  pN.Prime.W(Iprime=I.t[2], qNprime.W=qt.W[2])
					qt.D[2]    <-  qN.Prime.D(qN.D=qt.D[1], pN.D=pt.D[1], I=I.t[1], Iprime=I.t[2], u=u, h=h, s=s)
					pt.D[2]    <-  pN.Prime.D(Iprime=I.t[2], qNprime.D=qt.D[2])

					#binomial sampling
					I.t[1]     <-  rbinom(1, 2*N, I.t[2]) / (2*N)

					# Frame shift
					qt.I.W[1]  <-  qt.I.W[2]
					pt.I.W[1]  <-  pt.I.W[2]
					qt.W[1]    <-  qt.W[2]
					pt.W[1]    <-  pt.W[2]
					qt.D[1]    <-  qt.D[2]
					pt.D[1]    <-  pt.D[2]
				}
				if(I.t[1] == 1) {
					fix        <-  fix + 1
					rFixedInv  <-  c(rFixedInv, rs[l])
				}
			}

			PrFix         <-  c(PrFix, (fix/sims))
			rTab          <-  as.data.frame(table(rFixedInv))
			rTab$Freq     <-  rTab$Freq/sims
			rFixedInvTab  <-  rbind(rFixedInvTab, rTab)
			rNs           <-  c(rNs, rep(N, times = nrow(rTab)))
			rUs           <-  c(rUs, rep(U.vals[j], times = nrow(rTab)))
			rInvSize      <-  c(rInvSize, rep(invSize[k], times = nrow(rTab)))
								cat('\r', paste("U: ", j, "/", length(U.vals), 
								", x: ", k, "/", length(invSize), " complete", sep=""))
		}
	}

	# Index variables
	Ns        <-  rep(N, times=(length(U.vals)*length(invSize)))
	Us        <-  rep(U.vals, each=(length(invSize)))
	invSizes  <-  rep(invSize, times=length(U.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./data/RECODE/PrFixFig_AutoExact_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

	filename  <-  paste("./data/RECODE/PrFixFig_AutoExact_rFixedInv_h", h, "_s", s, Nfname, ".csv", sep="")
	r.d       <-  as.data.frame(cbind(rNs, rUs, rInvSize, rFixedInvTab))
	write.csv(r.d, file=filename, row.names=FALSE)

}








####################################################################
####################################################################
#' Function to generate timeseries of allele frequencies during W-F
#' Simulations to look at what is going on with feedback between 
#' inversion frequency & del. allele frequencies on each chromosome
#' class
makeData_AutoExact_WFDynamics  <-  function(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
											nTot = 10^4, N = 10^5, Nfname = "") {

	# Empty Vectors for 
	# Allele/Inversion timeseries
	I.tseries      <-  c()
	q.I.W.tseries  <-  c()
	p.I.W.tseries  <-  c()
	q.W.tseries    <-  c()
	p.W.tseries    <-  c()
	q.D.tseries    <-  c()
	p.D.tseries    <-  c()
	fixCounter     <-  c()
	rTracker       <-  c()
	invSizeVar     <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	I.0   <-  1/(2*N)

	# maxSims
	maxSims  <-  10*N

		# mutation rate & initial qHat
		u     <-  U/nTot
		qHat  <-  (U/(nTot*h*s))
		# Find equilibrium prior to inversion 
		init_eqs  <-  findEqFullRec(u=u, h=h, s=s, n=1, r=0, qHat_init=qHat)
		qHat  <-  init_eqs$qt.W


		# Loop over inversion size
		for(k in 1:length(invSize)) {
			
			# Counter for fixations
			fix = 0
			rep = 0
			# Loop over replicate simulations
			while(fix < nFix & rep < maxSims){

				# Draw random values for # loci captured by inversion
				n   <-  rpois(1, lambda=nTot*invSize[k])
				# Draw random value for # del. mutations captured by inversion (r) given x
				r  <-  sum(rbinom(n, size=1, prob=qHat))	

				# Empty Vectors for 
				# Allele/Inversion Frequencies
				I.t     <- c()
				qt.I.W  <- c()
				pt.I.W  <- c()
				qt.W    <- c()
				pt.W    <- c()
				qt.D    <- c()
				pt.D    <- c()

				# Assign frequencies
				# Note implicit assumption of equal initial
				# frequencies in XOv, XSp, Y chromosomes
				I.t[1]     <-  I.0 
				qt.I.W[1]  <-  0 
				pt.I.W[1]  <-  (I.0 - qHat*I.0) 
				qt.W[1]    <-  qHat*(1 - I.0) 
				pt.W[1]    <-  (1 - I.0 - qHat*(1 - I.0)) 
				qt.D[1]    <-  qHat*(1 - I.0) 
				pt.D[1]    <-  (1 - I.0 - qHat*(1 - I.0)) 
				

				# Run forward simulation
				i = 2
				while(I.t[i-1]*(1 - I.t[i-1]) > 0) {
					#expected frequency after selection
					I.t.next    <-  I.Prime(qN=qt.W[i-1], qI=qt.I.W[i-1], qD=qt.D[i-1], I=I.t[i-1], h=h, s=s, n=n, r=r)
					#binomial sampling
					I.t[i]     <-  rbinom(1, 2*N, I.t.next) / (2*N)
					#expected frequency after selection
					qt.I.W[i]   <-  qI.Prime.W(qN=qt.W[i-1], pN=pt.W[i-1], qI=qt.I.W[i-1], pI=pt.I.W[i-1], I=I.t[i-1], Iprime=I.t[i], u=u, h=h, s=s)
					pt.I.W[i]   <-  pI.Prime.W(Iprime=I.t[i], qIprime.W=qt.I.W[i-1])
					qt.W[i]     <-  qN.Prime.W(qN=qt.W[i-1], pN=pt.W[i-1], qI=qt.I.W[i-1], pI=pt.I.W[i-1], I=I.t[i-1], Iprime=I.t[i], u=u, h=h, s=s)
					pt.W[i]     <-  pN.Prime.W(Iprime=I.t[i], qNprime.W=qt.W[i-1])
					qt.D[i]     <-  qN.Prime.D(qN.D=qt.D[i-1], pN.D=pt.D[i-1], I=I.t[i-1], Iprime=I.t[i], u=u, h=h, s=s)
					pt.D[i]     <-  pN.Prime.D(Iprime=I.t[i], qNprime.D=qt.D[i-1])

					i  <-  i + 1
				}

				if(I.t[i-1] == 1) {
					fix        <-  fix + 1 
	
					# Concatenate timeseries
					I.tseries      <-  c(I.tseries, I.t)
					q.I.W.tseries  <-  c(q.I.W.tseries, qt.I.W)
					p.I.W.tseries  <-  c(p.I.W.tseries, pt.I.W)
					q.W.tseries    <-  c(q.W.tseries, qt.W)
					p.W.tseries    <-  c(p.W.tseries, pt.W)
					q.D.tseries    <-  c(q.D.tseries, qt.D)
					p.D.tseries    <-  c(p.D.tseries, pt.D)
					fixCounter     <-  c(fixCounter, rep(fix, times=length(I.t)))
					rTracker       <-  c(rTracker, rep(r, times=length(I.t)))
					invSizeVar     <-  c(invSizeVar, rep(invSize[k], times = length(I.t)))
				}

				rep  <-  rep+1
			}

			# Print Progress
			cat('\r', paste("x: ", k, "/", length(invSize), " complete", sep=""))

		}

	# Export Results Dataframe
	filename  <-  paste("./data/RECODE/AutoExact-WF-Dynamics_h", h, "_s", s, "_U", U, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"       =  rep(h, times=length(I.tseries)),
										"s"       =  rep(s, times=length(I.tseries)),
										"N"       =  rep(N, times=length(I.tseries)),
										"U"       =  rep(U, times=length(I.tseries)),
										"x"       =  invSizeVar,
										"I.t"     =  I.tseries,
										"qt.I.W"  =  q.I.W.tseries,
										"pt.I.W"  =  p.I.W.tseries,
										"qt.W"    =  q.W.tseries,
										"pt.W"    =  p.W.tseries,
										"qt.D"    =  q.D.tseries,
										"pt.D"    =  p.D.tseries,
										"fixRep"  =  fixCounter,
										"r"       =  rTracker
										)
	write.csv(d, file=filename, row.names=FALSE)

}
