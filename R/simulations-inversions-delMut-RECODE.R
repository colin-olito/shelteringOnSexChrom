# Simulation code: 
# Inversions expanding SLR w/ recessive deleterious mutations

#' Due to some strange inflation of fixation probabilities in 
#' my Wright-Fisher simulations at higher population sizes, 
#' I've decided to completely recode my recursions and W-F
#' simulation functions to make sure that I haven't just 
#' missed another small typo somewhere. These functions will
#' independently produce results that should be comparagble to
#' the deterministic and W-F results from the earlier version.



#' Exact Recursions.
#' This model requires a system of 8 recursions:
#' q.Prime.Xf.D, q.Prime.Xm.D, q.Prime.Y.D,
#' q.Prime.Xf.W, q.Prime.Xm.W, q.Prime.Y.W, q.Prime.YI.W, and
#' Y.I.Prime


# Recursions for loci where inversion captures a WT allele
q.Prime.Xf.W  <-  function(q.Xf.W, q.Xm.W, u, h, s) {
	(2*(1 - h*s)*u + q.Xm.W*(1 - u - s*(h + u*(2 - 3*h))) + q.Xf.W*(1 - u - s*(h + 2*q.Xm.W*(1 - h) + u*(2 - 3*h - 4*(1 - h)*q.Xm.W)))) / 
		(2*(1 - s*(q.Xm.W*u + q.Xf.W*(q.Xm.W + u*(1 - 2*q.Xm.W))) - h*s*(q.Xf.W + q.Xm.W*(1 - 2*q.Xf.W) + u*(2 - 3*q.Xm.W - q.Xf.W*(3 - 4*q.Xm.W)))))
}

q.Prime.Xm.W  <-  function(q.Xf.W, q.Y.W, YI, q.I.W, u, h, s) {
	(q.Xf.W*(1 + YI*(1 - 2*q.I.W*s) - h*s*(1 + YI*(1 - 2*q.I.W))) + u*(2 - 2*h*s - q.Xf.W*(1 + 2*s - 3*h*s) - YI*(q.Xf.W*(1 - h*s) + 2*(1 - h)*q.I.W*s*(1 - 2*q.Xf.W))) + q.Y.W*(1 - YI)*(1 - u - s*(h + q.Xf.W*(2 - 2*h) + u*(2 - 3*h - 4*(1 - h)*q.Xf.W)))) / 
		(2 - s*q.Xf.W*(2*q.Y.W*(1 - YI) + 2*q.I.W*YI) - 2*h*s*(q.Xf.W + q.Y.W*(1 - 2*q.Xf.W) + YI*(q.I.W - q.Y.W)*(1 - 2*q.Xf.W)) - 2*s*u*(2*h + q.Y.W*(1 - 3*h) + q.Xf.W*(1 - 3*h - q.Y.W*(2 - 4*h)) + YI*(q.I.W - q.Y.W)*(1 - 2*q.Xf.W - h*(3 - 4*q.Xf.W))))
}

q.Prime.Y.W  <-  function(q.Y.W, q.Xf.W, u, h, s) {
	(q.Y.W*(1 - u - s*(h + 2*q.Xf.W*(1 - h) + u*(2 - 3*h - 4*(1 - h)*q.Xf.W))) + q.Xf.W*(1 - u - s*(h + u*(2 - 3*h))) + 2*(1 - h*s)*u) / 
		(2*(1 - s*q.Xf.W*u - q.Y.W*s*(q.Xf.W + u*(1 - 2*q.Xf.W)) - h*s*(q.Y.W + q.Xf.W*(1 - 2*q.Y.W) + u*(2 - 3*q.Xf.W - q.Y.W*(3 - 4*q.Xf.W)))))
}

q.Prime.YI.W  <-  function(q.I.W, q.Xf.W, u, h, s) {
	(q.I.W*(1 - u - s*(q.Xf.W + u*(1 - 2*q.Xf.W) + h*(1 - q.Xf.W)*(1 - 2*u))) + u*(1 - s*(q.Xf.W + h*(1 - q.Xf.W)))) / 
		(1 - s*q.Xf.W*u - q.I.W*s*(q.Xf.W + u*(1 - 2*q.Xf.W)) - h*s*(q.I.W + q.Xf.W*(1 - 2*q.I.W) + u*(2 - 3*q.Xf.W - q.I.W*(3 - 4*q.Xf.W))))
}

# Recursions for loci where inversion captures a Del. allele
q.Prime.Xf.D  <-  function(q.Xf.D, q.Xm.D, u, h, s) {
	(q.Xf.D*(1 - u - s*(h + 2*q.Xm.D*(1 - h) + u*(2 - 3*h - 4*(1 - h)*q.Xm.D))) + q.Xm.D*(1 - u - s*(2*u + h*(1 - 3*u))) + 2*(1 - h*s)*u) / 
		(2*(1 - s*q.Xm.D*u - s*q.Xf.D*(q.Xm.D + u*(1 - 2*q.Xm.D)) - h*s*(q.Xf.D + q.Xm.D*(1 - 2*q.Xf.D) + u*(2 - 3*q.Xm.D - q.Xf.D*(3 - 4*q.Xm.D)))))
}

q.Prime.Xm.D  <-  function(q.Xf.D, q.Y.D, YI, u, h, s) {
	(q.Xf.D*(1 - h*s*(1 - YI) + (1 - 2*s)*YI) + (2 - 2*h*s - (1 + 2*s - 3*h*s)*q.Xf.D - (2*(1 - h)*s*(1 - 2*q.Xf.D) + (1 - h*s)*q.Xf.D)*YI)*u + q.Y.D*(1 - YI)*(1 - u - s*(h + (2 - 2*h)*q.Xf.D + (2 - 3*h - 4*(1 - h)*q.Xf.D)*u))) / 
		(2 - s*q.Xf.D*(2*q.Y.D*(1 - YI) + 2*YI) - 2*h*s*(q.Y.D*(1 - 2*q.Xf.D) + q.Xf.D + (1 - q.Y.D)*(1 - 2*q.Xf.D)*YI) - 2*s*(2*h + (1 - 3*h)*q.Y.D + (1 - 3*h - (2 - 4*h)*q.Y.D)*q.Xf.D + (1 - q.Y.D)*(1 - h*(3 - 4*q.Xf.D) - 2*q.Xf.D)*YI)*u)
}

q.Prime.Y.D  <-  function(q.Xf.D, q.Y.D, u, h, s) {
	(2*(1 - h*s)*u + q.Xf.D*(1 - u - s*(h + u*(2 - 3*h))) + q.Y.D*(1 - u - s*(h + q.Xf.D*(2 - 2*h) + (2 - 3*h - 4*(1 - h)*q.Xf.D)*u))) / 
		(2*(1 - s*q.Xf.D*u - q.Y.D*s*(u + q.Xf.D*(1 - 2*u)) - h*s*(q.Y.D + q.Xf.D*(1 - 2*q.Y.D) + (2 - 3*q.Xf.D - q.Y.D*(3 - 4*q.Xf.D))*u)))
}

# Recursion for Inversion frequency (Eqs. 1 & 2 from paper)
wBarY  <-  function(n, r, s, h, YI, q.I.W, q.Y.W, q.Y.D, q.Xf.W, q.Xf.D) {
	p.Xf.D  <-  (1 - q.Xf.D)
	p.Xf.W  <-  (1 - q.Xf.W)
	p.I.W   <-  (1 - q.I.W)
	p.Y.D   <-  (1 - q.Y.D)
	p.Y.W   <-  (1 - q.Y.W)
		   YI*( ( 1 - s*(h*p.Xf.D + q.Xf.D) )^(r) * ( 1 - s*(h*(p.I.W*q.Xf.W + q.I.W*p.Xf.W) + q.I.W*q.Xf.W))^(n-r)) +
	 (1 - YI)*( ( 1 - s*(h*(p.Xf.D*q.Y.D + q.Xf.D*p.Y.D) + q.Xf.D*q.Y.D) )^(r) * ( 1 - s*(h*(p.Xf.W*q.Y.W + q.Xf.W*p.Y.W) + q.Xf.W*q.Y.W) )^(n - r) )
}

YI.prime  <-  function(n, r, s, h, YI, q.I.W, q.Y.W, q.Y.D, q.Xf.W, q.Xf.D) {
	p.Xf.D  <-  (1 - q.Xf.D)
	p.Xf.W  <-  (1 - q.Xf.W)
	p.I.W   <-  (1 - q.I.W)
	YI*( ( 1 - s*(h*p.Xf.D + q.Xf.D) )^(r) * ( 1 - s*(h*(p.I.W*q.Xf.W + q.I.W*p.Xf.W) + q.I.W*q.Xf.W))^(n-r)) / wBarY(n=n, r=r, s=s, h=h, YI=YI, q.I.W=q.I.W, q.Y.W=q.Y.W, q.Y.D=q.Y.D, q.Xf.W=q.Xf.W, q.Xf.D=q.Xf.D)
}

invRelFit  <-  function(n, r, s, h, q.I.W, q.Y.W, q.Y.D, q.Xf.W, q.Xf.D){
	p.I.W  <-  1 - q.I.W
	p.Y.D  <-  1 - q.Y.D
	p.Y.W  <-  1 - q.Y.W
		  (( 1 - s*(h*(1 - q.Xf.D) + q.Xf.D) )^r * (1 - s*(h*(p.I.W*q.Xf.W + q.I.W*(1 - q.Xf.W)) + q.I.W*q.Xf.W))^(n-r)) /
				( (1 - s*(h*((1 - q.Xf.D)*q.Y.D + q.Xf.D*p.Y.D) + q.Xf.D*q.Y.D))^r * (1 - s*(h*((1 - q.Xf.W)*q.Y.W + q.Xf.W*p.Y.W) + q.Xf.W*q.Y.W))^(n - r) )
}

invRelFit2  <-  function(n, r, s, h, YI, q.I.W, q.Y.W, q.Y.D, q.Xf.W, q.Xf.D) {
	p.Xf.D  <-  (1 - q.Xf.D)
	p.Xf.W  <-  (1 - q.Xf.W)
	p.I.W   <-  (1 - q.I.W)
	YI.prime  <-  YI*( ( 1 - s*(h*p.Xf.D + q.Xf.D) )^(r) * ( 1 - s*(h*(p.I.W*q.Xf.W + q.I.W*p.Xf.W) + q.I.W*q.Xf.W))^(n-r)) / wBarY(n=n, r=r, s=s, h=h, YI=YI, q.I.W=q.I.W, q.Y.W=q.Y.W, q.Y.D=q.Y.D, q.Xf.W=q.Xf.W, q.Xf.D=q.Xf.D)
	YI.prime/YI
}


#################################
#################################

#' Function to calculate euclidean distance between two vectors
#' Convenient for checking if equilibrium has been reached
eucDist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
} 

#' Function to find equilibrium frequencies in the absence of inverions
findEq  <- function(r, h, s, n, u, qHat, qHatDel) {

	# Empty Frequency Vectors
	YI.t     <-  c()
	qt.Xf.W  <-  c()
	qt.Xm.W  <-  c()
	qt.YI.W  <-  c()
	qt.Y.W   <-  c()
	qt.Xf.D  <-  c()
	qt.Xm.D  <-  c()
	qt.Y.D   <-  c()

	# Find equilibrium prior to inversion 
	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]     <-  round(YI.prime(n=n, r=r, s=s, h=h, YI=0, q.I.W=0, q.Y.W=qHat, q.Y.D=qHatDel, q.Xf.W=qHat, q.Xf.D=qHatDel), digits=16)
	qt.Xf.W[1]  <-  round(q.Prime.Xf.W(u=u, s=s, h=h, q.Xf.W=qHat, q.Xm.W=qHat), digits=16)
	qt.Xm.W[1]  <-  round(q.Prime.Xm.W(u=u, s=s, h=h, YI=0, q.I.W=0, q.Y.W=qHat, q.Xf.W=qHat), digits=16)
	qt.YI.W[1]  <-  round(q.Prime.YI.W(u=u, s=s, h=h, q.I.W=0, q.Xf.W=qHat), digits=16)
	qt.Y.W[1]   <-  round(q.Prime.Y.W(u=u, s=s, h=h, q.Y.W=qHat, q.Xf.W=qHat), digits=16)
	qt.Xf.D[1]  <-  round(q.Prime.Xf.D(u=u, s=s, h=h, q.Xf.D=qHatDel, q.Xm.D=qHatDel), digits=16)
	qt.Xm.D[1]  <-  round(q.Prime.Xm.D(u=u, s=s, h=h, YI=0, q.Y.D=qHatDel, q.Xf.D=qHatDel), digits=16)
	qt.Y.D[1]   <-  round(q.Prime.Y.D(u=u, s=s, h=h, q.Xf.D=qHatDel, q.Y.D=qHatDel), digits=16)

	# Subsequent generations
	i=2
	diff  <-  1
	while(diff > 1e-12) {
		YI.t[i]     <-  round(YI.prime(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.I.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		qt.Xf.W[i]  <-  round(q.Prime.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[i-1], q.Xm.W=qt.Xm.W[i-1]), digits=16)
		qt.Xm.W[i]  <-  round(q.Prime.Xm.W(u=u, s=s, h=h, YI=YI.t[i-1], q.I.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=16)
		qt.Xf.D[i]  <-  round(q.Prime.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Xm.D=qt.Xm.D[i-1]), digits=16)
		qt.Xm.D[i]  <-  round(q.Prime.Xm.D(u=u, s=s, h=h, YI=YI.t[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		qt.YI.W[i]  <-  round(q.Prime.YI.W(u=u, s=s, h=h, q.I.W=qt.YI.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=16)
		qt.Y.W[i]   <-  round(q.Prime.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=16)
		qt.Y.D[i]   <-  round(q.Prime.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Y.D=qt.Y.D[i-1]), digits=16)
		diff  <-  eucDist(c(qt.Xf.W[i-1],qt.Xm.W[i-1],qt.Xf.D[i-1],qt.Xm.D[i-1],qt.Y.W[i-1],qt.Y.D[i-1]),
						  c(qt.Xf.W[i],qt.Xm.W[i],qt.Xf.D[i],qt.Xm.D[i],qt.Y.W[i],qt.Y.D[i]))
		i  <-  i + 1
	}
	res  <-  list(
				  "q.Xf.W.init" = qt.Xf.W[i-1],
				  "q.Xm.W.init" = qt.Xm.W[i-1],
				  "q.Xf.D.init" = qt.Xf.D[i-1],
				  "q.Xm.D.init" = qt.Xm.D[i-1],
				  "q.Y.W.init"  = qt.Y.W[i-1],
				  "q.Y.D.init"  = qt.Y.D[i-1]
				  )
	return(res)
}




#' Function to generate deterministic frequency dynamics 
#' for SLR-expanding inversions
makeDeterministicFigSimData  <-  function(U = 0.02, r = 0, x = 0.2, h = 0.25, s = 0.01, generations = 10^4, ...) {
	# Parameters
	nTot  <-  10^4
	u     <-  U/nTot
	qHat  <-  (U/(nTot*h*s))
	if(r == 0) {
		qHatDel  <-  0
	} else {qHatDel  <-  qHat}
	n     <-  nTot*x
	N     <-  5000
	YI.0  <-  2/N

	# Empty Frequency Vectors
	YI.t       <- c()
	qt.Xf.W    <- c()
	qt.Xm.W    <- c()
	qt.YI.W    <- c()
	qt.Y.W     <- c()
	qt.Xf.D    <- c()
	qt.Xm.D    <- c()
	qt.Y.D     <- c()
	relW.YI.t  <- c()
	Wbar.Y.t   <- c()

	# Find equilibrium prior to inversion 
	eqs  <-  findEq(r=r, h=h, s=s, n=n, u=u, qHat=qHat, qHatDel=qHatDel)

	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]       <-  round(YI.prime(n=n, r=r, s=s, h=h, YI=YI.0, q.I.W=0, q.Y.W=eqs$q.Y.W.init, q.Y.D=eqs$q.Y.D.init, q.Xf.W=eqs$q.Xf.W.init, q.Xf.D=eqs$q.Xf.D.init), digits=12)
	qt.Xf.W[1]    <-  round(q.Prime.Xf.W(u=u, s=s, h=h, q.Xf.W=eqs$q.Xf.W.init, q.Xm.W=eqs$q.Xm.W.init), digits=12)
	qt.Xm.W[1]    <-  round(q.Prime.Xm.W(u=u, s=s, h=h, YI=YI.0, q.I.W=0, q.Y.W=eqs$q.Y.W.init, q.Xf.W=eqs$q.Xf.W.init), digits=12)
	qt.YI.W[1]    <-  round(q.Prime.YI.W(u=u, s=s, h=h, q.I.W=0, q.Xf.W=eqs$q.Xf.W.init), digits=12)
	qt.Y.W[1]     <-  round(q.Prime.Y.W(u=u, s=s, h=h, q.Y.W=eqs$q.Y.W.init, q.Xf.W=eqs$q.Xf.W.init), digits=12)
	qt.Xf.D[1]    <-  round(q.Prime.Xf.D(u=u, s=s, h=h, q.Xf.D=eqs$q.Xf.D.init, q.Xm.D=eqs$q.Xm.D.init), digits=12)
	qt.Xm.D[1]    <-  round(q.Prime.Xm.D(u=u, s=s, h=h, YI=YI.0, q.Y.D=eqs$q.Y.D.init, q.Xf.D=eqs$q.Xf.D.init), digits=12)
	qt.Y.D[1]     <-  round(q.Prime.Y.D(u=u, s=s, h=h, q.Xf.D=eqs$q.Xf.D.init, q.Y.D=eqs$q.Y.D.init), digits=12)
	relW.YI.t[1]  <-  round(invRelFit(n=n, r=r, s=s, h=h, q.I.W=0, q.Y.W=eqs$q.Y.W.init, q.Y.D=eqs$q.Y.D.init, q.Xf.W=eqs$q.Xf.W.init, q.Xf.D=eqs$q.Xf.D.init), digits=12)
	Wbar.Y.t[1]   <-  round(wBarY(n=n, r=r, s=s, h=h, YI=YI.0, q.I.W=0, q.Y.W=eqs$q.Y.W.init, q.Y.D=eqs$q.Y.D.init, q.Xf.W=eqs$q.Xf.W.init, q.Xf.D=eqs$q.Xf.D.init), digits=12)


	# Subsequent generations
	i=2
	while(i < generations+1) {
		YI.t[i]       <-  round(YI.prime(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.I.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=12)
		qt.Xf.W[i]    <-  round(q.Prime.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[i-1], q.Xm.W=qt.Xm.W[i-1]), digits=12)
		qt.Xm.W[i]    <-  round(q.Prime.Xm.W(u=u, s=s, h=h, YI=YI.t[i-1], q.I.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=12)
		qt.Xf.D[i]    <-  round(q.Prime.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Xm.D=qt.Xm.D[i-1]), digits=12)
		qt.Xm.D[i]    <-  round(q.Prime.Xm.D(u=u, s=s, h=h, YI=YI.t[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=12)
		qt.YI.W[i]    <-  round(q.Prime.YI.W(u=u, s=s, h=h, q.I.W=qt.YI.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=12)
		qt.Y.W[i]     <-  round(q.Prime.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=12)
		qt.Y.D[i]     <-  round(q.Prime.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Y.D=qt.Y.D[i-1]), digits=12)
		relW.YI.t[i]  <-  round(invRelFit(n=n, r=r, s=s, h=h, q.I.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=12)
		Wbar.Y.t[i]   <-  round(wBarY(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.I.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=12)
		i  <-  i + 1
	}
	# Tack on initial freq.
	YI.t       <-  c(YI.0,YI.t)
	qt.Xf.W    <-  c(eqs$q.Xf.W.init, qt.Xf.W)
	qt.Xm.W    <-  c(eqs$q.Xm.W.init, qt.Xm.W)
	qt.Xf.D    <-  c(eqs$q.Xf.D.init, qt.Xf.D)
	qt.Xm.D    <-  c(eqs$q.Xm.D.init, qt.Xm.D)
	qt.YI.W    <-  c(0, qt.YI.W)
	qt.Y.W     <-  c(eqs$q.Y.W.init, qt.Y.W)
	qt.Y.D     <-  c(eqs$q.Y.D.init, qt.Y.D)
	relW.YI.t  <-  c(relW.YI.t[1], relW.YI.t)
	Wbar.Y.t   <-  c(Wbar.Y.t[1], Wbar.Y.t)
	# Return results as df
	results  <-  data.frame(
							"YI.t"        =  YI.t,
							"qt.Xf.W"     =  qt.Xf.W,
							"qt.Xm.W"     =  qt.Xm.W,
							"qt.Xf.D"     =  qt.Xf.D,
							"qt.Xm.D"     =  qt.Xm.D,
							"qt.YI.W"     =  qt.YI.W,
							"qt.Y.W"      =  qt.Y.W,
							"qt.Y.D"      =  qt.Y.D,
							"rel.w.YI.t"  =  relW.YI.t,
							"Wbar.Y.t"    =  Wbar.Y.t
		 					)
	return(results)
}


################################
################################
# Wright-Fisher Simulations

#' Function to estimate Pr(fix | x) for different inversion sizes (x)
#' EXPANDED WF simulation model: uses full set of deterministic
#' recursions to model changes in deleterious allele frequencies.
makeDataPrFixInvSize  <-  function(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
									nTot = 10^4, N = 10^4, Nfname = "") {

	# Containers
	PrFix      <-  c()
	rFixedInv  <-  c()
	YI.t       <-  c()
	qt.Xf.W    <-  c()
	qt.Xm.W    <-  c()
	qI.wt.t    <-  c()
	qY.wt.t    <-  c()
	qt.Xf.D    <-  c()
	qt.Xm.D    <-  c()
	qY.del.t   <-  c()

	rFixedInvTab  <-  data.frame()
	rInvSize      <-  c()
	rNs           <-  c()
	rUs           <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	YI.0   <-  2/N

	#number of simulations 
	sims  = 200*N/2


	# Loop over Us factor values
	for(j in 1:length(U.vals)) {
		# mutation rate & initial qHat
		U     <-  U.vals[j]
		u     <-  U/nTot
		qHat  <-  (U/(nTot*h*s))

		# Find equilibrium prior to inversion 
		eqs  <-  findEq(r=0, h=h, s=s, n=nTot, u=u, qHat=qHat, qHatDel=qHat)
		qHat  <-  eqs$q.Xf.W.init

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
				YI.t[1]     <-  YI.0
				qt.Xf.W[1]  <-  qHat
				qt.Xm.W[1]  <-  qHat
				qt.I.W[1]   <-  0
				qt.Y.W[1]   <-  qHat
				qt.Xf.D[1]  <-  qHat
				qt.Xm.D[1]  <-  qHat
				qt.Y.D[1]   <-  qHat

				# Run forward simulation
				while(YI.t[1]*(1 - YI.t[1]) > 0) {
					#expected frequency after selection
					qt.Xf.W[2]  <-  round(q.Prime.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[1], q.Xm.W=qt.Xm.W[1]), digits=12)
					qt.Xf.D[2]  <-  round(q.Prime.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[1], q.Xm.D=qt.Xm.D[1]), digits=12)
					qt.Xm.W[2]  <-  round(q.Prime.Xm.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[1], q.Y.W=qt.Y.W[1], YI=YI.t[1], q.I.W=qt.I.W[1]), digits=12)
					qt.Xm.D[2]  <-  round(q.Prime.Xm.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[1], q.Y.D=qt.Y.D[1], YI=YI.t[1]), digits=12)
					qt.Y.W[2]   <-  round(q.Prime.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[1], q.Xf.W=qt.Xf.W[1]), digits=12)
					qt.Y.D[2]   <-  round(q.Prime.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[1], q.Y.D=qt.Y.D[1]), digits=12)
					qt.I.W[2]   <-  round(q.Prime.YI.W(u=u, s=s, h=h, q.I.W=qt.I.W[1], q.Xf.W=qt.Xf.W[1]), digits=12)
					YI.t[2]     <-  round(YI.prime(n=n, r=r, s=s, h=h, YI=YI.t[1], q.I.W=qt.I.W[1], q.Y.W=qt.Y.W[1], q.Y.D=qt.Y.D[1], q.Xf.W=qt.Xf.W[1], q.Xf.D=qt.Xf.D[1]), digits=12)

					#binomial sampling
					YI.t[1]     <-  rbinom(1, (N/2), YI.t[2]) / (N/2)

					# Frame shift
					qt.Xf.W[1]  <-  qt.Xf.W[2]
					qt.Xf.D[1]  <-  qt.Xf.D[2]
					qt.Xm.W[1]  <-  qt.Xm.W[2]
					qt.Xm.D[1]  <-  qt.Xm.D[2]
					qt.Y.W[1]   <-  qt.Y.W[2]
					qt.Y.D[1]   <-  qt.Y.D[2]
					qt.I.W[1]   <-  qt.I.W[2]
				}
				if(YI.t[1] == 1) {
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
	filename  <-  paste("./data/RECODE/PrFixFig_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

	filename  <-  paste("./data/RECODE/PrFixFig_rFixedInv_h", h, "_s", s, Nfname, ".csv", sep="")
	r.d       <-  as.data.frame(cbind(rNs, rUs, rInvSize, rFixedInvTab))
	write.csv(r.d, file=filename, row.names=FALSE)

}





######################################
#' Equilibrium Approximation WF model

#' Deterministic fitness expressions for
#' equilibrium approximation WF model 
w.YI.X  <-  function(n, r, t, u, s, h, qHat) {
	((qHat*(1 - s) + (1 - qHat)*(1 - s*h))^r) * (1 - u*(2 - exp(-s*h*t)))^(n - r)
}
w.Y.X  <-  function(n, u) {
	(1 - 2*u)^(n)
}


####################################################################
#' Function to estimate Pr(fix | x) for different inversion sizes (x)
#' Using equilibrium approximation WF simulations
makeDataPrFixInvSizeHaploid  <-  function(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
																					nTot = 10^4, N = 10^4, Nfname = "") {

	# Containers
	PrFix      <-  c()
	rFixedInv  <-  c()
	YI.t       <-  c()

	rFixedInvTab  <-  data.frame()
	rInvSize      <-  c()
	rNs           <-  c()
	rUs           <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	YI.0   <-  2/N

	#number of simulations 
	sims  = 200*N/2


	# Loop over Us factor values
	for(j in 1:length(U.vals)) {
		# mutation rate & initial qHat
		U     <-  U.vals[j]
		u     <-  U/nTot
		qHat  <-  u/(h*s)

		# Loop over inversion size
		for(k in 1:length(invSize)) {


			# Draw random values for # loci captured by inversion
			ns   <-  rpois(sims, lambda=nTot*invSize[k])

			# Draw random value for # del. mutations captured by inversion (r) given x
			rs  <-  c()
			for(m in 1:length(ns)) {
				rs[m]  <-  sum(rbinom(n=ns[m], size=1, prob=qHat))
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
				YI.t  <-  YI.0
				t     <-  1

				# Run forward simulation
				while(YI.t*(1 - YI.t) > 0) {
					# Fitness expressions
							w_YI.X  <-  w.YI.X(n=n, r=r, t=t, u=u, s=s, h=h, qHat=qHat)
							w_Y.X   <-  w.Y.X(n=n, u=u)
							w_avg   <-  YI.t*w_YI.X + (1 - YI.t)*w_Y.X
          		YI.sel  <-  YI.t*w_YI.X/w_avg

					#binomial sampling
					YI.t     <-  rbinom(1, (N/2), YI.sel) / (N/2)
					# time counter
					t  <-  t + 1
				}
				if(YI.t == 1) {
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
	filename  <-  paste("./data/RECODE/PrFixFig_Haploid_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

	filename  <-  paste("./data/RECODE/PrFixFig_Haploid_rFixedInv_h", h, "_s", s, Nfname, ".csv", sep="")
	r.d       <-  as.data.frame(cbind(rNs, rUs, rInvSize, rFixedInvTab))
	write.csv(r.d, file=filename, row.names=FALSE)

}





##################################
#' Make deterministic time-series
#' using haploid model
makeDeterministicSimHaploidData  <-  function(U = 0.02, r = 0, x = 0.2, 
																								h = 0.25, s = 0.01, generations = 10^4, ...) {
	# Parameters
	nTot  <-  10^4
	u     <-  U/nTot
	qHat  <-  (U/(nTot*h*s))
	n     <-  nTot*x
	N     <-  5000
	YI.0  <-  2/N

	# Empty Frequency Vectors
	YI.t        <- rep(NA, times=generations+1)
	rel.w.YI.t  <- rep(NA, times=generations)
	w_avg       <- rep(NA, times=generations)


	# Run deterministic forward simulation
	t  <-  1
	for(i in 1:generations) {
		if(i==1) {
			YI.t[i]  <-  YI.0
		}
		# Fitness expressions
		w_YI.X      <-  w.YI.X(n=n, r=r, t=t, u=u, s=s, h=h, qHat=qHat)
		w_Y.X       <-  w.Y.X(n=n, u=u)
		rel.w.YI.t[i]  <-  w_YI.X/w_Y.X
		w_avg[i]       <-  YI.t[i]*w_YI.X + (1 - YI.t[i])*w_Y.X
		YI.t[i+1]      <-  YI.t[i]*w_YI.X/w_avg[i]

		# time counter
		t  <-  t + 1
	}
	YI.t  <-  YI.t[1:generations]
	# Return results as df
	results  <-  data.frame(
							"YI.t"        =  YI.t,
							"rel.w.YI.t"  =  rel.w.YI.t,
							"Wbar.Y.t"    =  w_avg
		 					)
	return(results)
}




########################################
##  AUTOSOMAL W-F Simulations using
#' equilibrium approximation model 
#' (after Connallon & Olito 2019)

#' Deterministic time-dependent fitness expressions
w.II  <-  function(n, nd, u, sdHom, h, Ud, x, t) {
	((1 - sdHom)^nd)*(1 - 2*u*(1 - exp(-sdHom*h*t)))^(n - nd)
}
w.IS  <-  function(n, nd, u, sdHom, h, Ud, x, t) {
	((u/(sdHom*h))*(1 - sdHom) + (1 - (u/(sdHom*h)))*(1 - sdHom*h))^nd * (1 - u*(2 - exp(-sdHom*h*t)))^(n - nd)
}
w.SS  <-  function(n, nd, u, sdHom, h, Ud, x, t) {
	(1 - 2*u)^n
}

####################################################################
#' Function to estimate Pr(fix | x) for different inversion sizes (x)
#' Using equilibrium approximation AUTOSOMAL WF simulations
makeDataAutoPrFixInvSize  <-  function(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
																	 nTot = 10^4, N.vals = c(10^3, 10^4), Nfname = "") {

	# Containers
	PrFix   <-  c()
	PrFix0  <-  c()
	q.t     <-  c()
	
	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N     <-  N.vals[i]
		q.0   <-  1/(2*N)

		#number of simulations 
		sims  = 100*N

		# Loop over Us factor values
		for(j in 1:length(U.vals)) {
			U     <-  U.vals[j]
			u     <-  U/nTot
			qHat  <-  u/(h*s)

				# Loop over inversion size
				for(k in 1:length(invSize)) {
					fix   = 0
					fix0  = 0

						# Draw random values for # loci captured by inversion
						ns   <-  rpois(sims, lambda=nTot*invSize[k])
						# Draw random value for # del. mutations captured by inversion (r) given x
						rs  <-  c()
						for(m in 1:length(ns)) {
							rs[m]  <-  sum(rbinom(ns[m], size=1, prob=qHat))	
						}

					# Loop over replicate simulations
					for(l in 1:sims){

						# Take randomly drawn value for # loci captured by inversion
						n   <-  ns[l]

						# Take randomly drawn value for # del. mutations captured
						# by inversion (r) given x
						r  <-  rs[l]

						# Assign frequencies
						# Note implicit assumption of equal
						# frequencies at all loci
						q.t  <-  q.0
						p.t  <-  1 - q.t

						# Run forward simulation
						t  <-  1
						while(q.t*p.t > 0) {

							#expected frequency after selection
							w_II    <-  w.II(n = n, nd=r, sdHom=s, h=h, u=u, Ud=U, x=invSize[k], t=t)
							w_IS    <-  w.IS(n = n, nd=r, sdHom=s, h=h, u=u, Ud=U, x=invSize[k], t=t)
							w_SS    <-  w.SS(n = n, nd=r, sdHom=s, h=h, u=u, Ud=U, x=invSize[k], t=t)
							w.avg   <-  p.t^2*w_SS + 2*p.t*q.t*w_IS + q.t^2*w_II
							SS.sel  <-  p.t^2*w_SS/w.avg
							IS.sel  <-  2*p.t*q.t*w_IS/w.avg
							II.sel  <-  q.t^2*w_II/w.avg
							
							# multinomial sampling
							drift  <-  rmultinom(1, N, c(SS.sel, IS.sel, II.sel)) / N          
          		p.t    <-  drift[1] + drift[2]/2
          		q.t    <-  drift[2]/2 + drift[3]

							# time counter
							t  <-  t + 1
						}
						if(q.t == 1){
							fix  <-  fix + 1
							if(r  == 0){
								fix0  <-  fix0 + 1
							}
						}
					}
				PrFix  <-  c(PrFix, (fix/sims))
				PrFix0  <-  c(PrFix0, (fix0/sims))
				cat('\r', paste("N: ", i, "/", length(N.vals), 
											", U: ", j, "/", length(U.vals), 
											", x: ", round(100*(k/length(invSize))), "% complete", sep=""))
				}

		}

	}

	# Index variables
	Ns        <-  rep(N.vals, each=(length(U.vals)*length(invSize)))
	Us        <-  rep((U.vals), each=(length(invSize)), times=length(N.vals))
	invSizes  <-  rep(invSize, times=length(N.vals)*length(U.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./data/RECODE/PrFixAutoFig_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"x"      =  invSizes,
										"PrFix"  =  PrFix,
										"PrFix0"  =  PrFix0
										)
	write.csv(d, file=filename, row.names=FALSE)

}












####################################################################
####################################################################
#' Function to generate timeseries of allele frequencies during W-F
#' Simulations to look at what is going on with feedback between 
#' inversion frequency & del. allele frequencies on each chromosome
#' class. Used to explore what was going on in expanded WF simulations
makeData_WFDynamics  <-  function(h = 0.25, s = 0.01, U = 0.02, nFix = 3,
																	nTot = 10^4, N = 10^5, Nfname = "") {

	# Empty Vectors for 
	# Allele/Inversion timeseries
	YI.tseries      <-  c()
	q.Xf.W.tseries  <-  c()
	q.Xm.W.tseries  <-  c()
	qI.wt.tseries   <-  c()
	qY.wt.tseries   <-  c()
	q.Xf.D.tseries  <-  c()
	q.Xm.D.tseries  <-  c()
	qY.del.tseries  <-  c()
	fixCounter      <-  c()
	rTracker        <-  c()
	invSizeVar      <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	YI.0   <-  2/N

	# maxSims
	maxSims  <-  10*N

		# mutation rate & initial qHat
		u     <-  U/nTot
		qHat  <-  (U/(nTot*h*s))
		# Find equilibrium prior to inversion 
		eqs  <-  findEq(r=0, h=h, s=s, n=nTot, u=u, qHat=qHat, qHatDel=qHat)
		qHat  <-  eqs$q.Xf.W.init


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
				YI.t     <-  c()
				qt.Xf.W  <-  c()
				qt.Xm.W  <-  c()
				qt.I.W   <-  c()
				qt.Y.W   <-  c()
				qt.Xf.D  <-  c()
				qt.Xm.D  <-  c()
				qt.Y.D   <-  c()

				# Assign frequencies
				# Note implicit assumption of equal initial
				# frequencies in XOv, XSp, Y chromosomes
				YI.t[1]     <-  YI.0
				qt.Xf.W[1]  <-  qHat
				qt.Xm.W[1]  <-  qHat
				qt.I.W[1]   <-  0
				qt.Y.W[1]   <-  qHat
				qt.Xf.D[1]  <-  qHat
				qt.Xm.D[1]  <-  qHat
				qt.Y.D[1]   <-  qHat

				# Run forward simulation
				i = 2
					while(YI.t[i-1]*(1 - YI.t[i-1]) > 0) {
					#expected frequency after selection
					qt.Xf.W[i]  <-  round(q.Prime.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[i-1], q.Xm.W=qt.Xm.W[i-1]), digits=12)
					qt.Xf.D[i]  <-  round(q.Prime.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Xm.D=qt.Xm.D[i-1]), digits=12)
					qt.Xm.W[i]  <-  round(q.Prime.Xm.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[i-1], q.Y.W=qt.Y.W[i-1], YI=YI.t[i-1], q.I.W=qt.I.W[i-1]), digits=12)
					qt.Xm.D[i]  <-  round(q.Prime.Xm.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Y.D=qt.Y.D[i-1], YI=YI.t[i-1]), digits=12)
					qt.Y.W[i]   <-  round(q.Prime.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=12)
					qt.Y.D[i]   <-  round(q.Prime.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Y.D=qt.Y.D[i-1]), digits=12)
					qt.I.W[i]   <-  round(q.Prime.YI.W(u=u, s=s, h=h, q.I.W=qt.I.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=12)
					YI.sel      <-  round(YI.prime(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.I.W=qt.I.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=12)
					#binomial sampling
					YI.t[i]     <-  rbinom(1, (N/2), YI.sel) / (N/2)
					i  <-  i + 1
				}
				if(YI.t[i-1] == 1) {
					fix        <-  fix + 1 
	
					# Concatenate timeseries
					YI.tseries      <-  c(YI.tseries, YI.t)
					q.Xf.W.tseries  <-  c(q.Xf.W.tseries, qt.Xf.W)
					q.Xm.W.tseries  <-  c(q.Xm.W.tseries, qt.Xm.W)
					qI.wt.tseries   <-  c(qI.wt.tseries, qt.I.W)
					qY.wt.tseries   <-  c(qY.wt.tseries, qt.Y.W)
					q.Xf.D.tseries  <-  c(q.Xf.D.tseries, qt.Xf.D)
					q.Xm.D.tseries  <-  c(q.Xm.D.tseries, qt.Xm.D)
					qY.del.tseries  <-  c(qY.del.tseries, qt.Y.D)
					fixCounter      <-  c(fixCounter, rep(fix, times=length(YI.t)))
					rTracker        <-  c(rTracker, rep(r, times=length(YI.t)))
					invSizeVar  <-  c(invSizeVar, rep(invSize[k], times = length(YI.t)))
				}

				rep  <-  rep+1
			}

			# Print Progress
			cat('\r', paste("x: ", k, "/", length(invSize), " complete", sep=""))

		}

	# Export Results Dataframe
	filename  <-  paste("./data/RECODE/SLR-WF-Dynamics_h", h, "_s", s, "_U", U, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"        =  rep(h, times=length(YI.tseries)),
										"s"        =  rep(s, times=length(YI.tseries)),
										"N"        =  rep(N, times=length(YI.tseries)),
										"U"        =  rep(U, times=length(YI.tseries)),
										"x"        =  invSizeVar,
										"YI.t"     =   YI.tseries,
										"qt.Xf.W"  =   q.Xf.W.tseries,
										"qt.Xm.W"  =   q.Xm.W.tseries,
										"qt.I.W"   =   qI.wt.tseries,
										"qt.Y.W"   =   qY.wt.tseries,
										"qt.Xf.D"  =   q.Xf.D.tseries,
										"qt.Xm.D"  =   q.Xm.D.tseries,
										"qt.Y.D"   =   qY.del.tseries,
										"fixRep"   =   fixCounter,
										"r"        =   rTracker
										)
	write.csv(d, file=filename, row.names=FALSE)

}
