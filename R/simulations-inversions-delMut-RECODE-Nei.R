# Simulation code: 
# Inversions expanding SLR w/ recessive deleterious mutations

#' Due to some strange inflation of fixation probabilities in 
#' my Wright-Fisher simulations at higher population sizes, 
#' I've decided to completely recode my recursions and W-F
#' simulation functions to make sure that I haven't just 
#' missed another small typo somewhere. These functions will
#' independently produce results that should be comparable to
#' the deterministic and W-F results from the earlier version.

#' MODIFIED VERSION OF ORIGINAL RECURSIONS FOLLOWING THE LOGIC
#' OF NEI ET AL. (1967), I.E., NOT BACK-TRANSFORMING FROM 
#' ABSOLUTE TO RELATIVE FREQUENCIES UNTIL THE NEXT GENERATION  
#' 

#' Exact Recursions.
#' This model requires a system of 8 recursions:
#' q.Prime2.Xf.D, q.Prime2.Xm.D, q.Prime2.Y.D,
#' q.Prime2.Xf.W, q.Prime2.Xm.W, q.Prime2.Y.W, q.Prime2.YI.W, and
#' Y.I.Prime


# Recursions for loci where inversion captures a WT allele
q.Prime2.Xf.W  <-  function(q.Xf.W, q.Xm.W, u, h, s) {
	(2*(-1 + h*s)*u + q.Xm.W*(-1 + u + s*(h + 2*u - 3*h*u)) + q.Xf.W*(-1 + u + s*(h + 2*q.Xm.W - 2*h*q.Xm.W + (2 - 3*h + 4*(-1 + h)*q.Xm.W)*u))) / 
		(2*(-1 + q.Xm.W*s*u + q.Xf.W*s*(q.Xm.W + u - 2*q.Xm.W*u) + h*s*(q.Xf.W + q.Xm.W - 2*q.Xf.W*q.Xm.W + (2 - 3*q.Xm.W + q.Xf.W*(-3 + 4*q.Xm.W))*u)))
}

q.Prime2.Xm.W  <-  function(q.Xf.W, q.Y.W, YI, q.YI.W, p.YI.W, u, h, s) {
	(2*(-1 + (h + q.YI.W - h*q.YI.W)*s)*u + q.Y.W*(-1 + u + s*(h + 2*u - 3*h*u)) + q.Xf.W*(-1 - p.YI.W*(-1 + h*s)*(-1 + u) + u + s*(h + 2*q.Y.W - 2*h*q.Y.W + (2 - 3*h + 4*(-1 + h)*q.Y.W)*u) + q.YI.W*(-1 + u + s*(2 - h - 4*u + 3*h*u)))) / 
		(2*(-1 + q.Y.W*s*u + q.YI.W*s*(q.Xf.W + u - 2*q.Xf.W*u) + q.Xf.W*s*(q.Y.W + u - 2*q.Y.W*u) + h*s*(q.YI.W + q.Xf.W - 2*q.YI.W*q.Xf.W + q.Y.W - 2*q.Xf.W*q.Y.W + (2 - 3*q.Xf.W + q.YI.W*(-3 + 4*q.Xf.W) - 3*q.Y.W + 4*q.Xf.W*q.Y.W)*u)))
}

q.Prime2.Y.W  <-  function(q.Xf.W, q.Y.W, p.Y.W, YI, YIprime, u, h, s) { #, q.YI.W, p.YI.W
	p.Xf.W  <-  (1 - q.Xf.W)

	# mutation
	q.Xf.W  <-  q.Xf.W + p.Xf.W*u
	p.Xf.W  <-  p.Xf.W*(1 - u)
	q.Y.W   <-  q.Y.W + p.Y.W*u
	p.Y.W   <-  p.Y.W*(1 - u)

	(q.Y.W*(1 - s*(h*p.Xf.W + q.Xf.W) ) / 
			(1 - YI - q.Y.W*s*(h*p.Xf.W + q.Xf.W) ) )* (1 - YIprime)
}

p.Prime2.Y.W  <-  function(YIprime, qYprime.W) {
		1 - YIprime - qYprime.W
}

q.Prime2.YI.W  <-  function(q.YI.W, p.YI.W, q.Xf.W, YI, YIprime, u, h, s) { #, q.Y.W
	p.Xf.W  <-  (1 - q.Xf.W)

	# mutation
	q.Xf.W  <-  q.Xf.W + p.Xf.W*u
	p.Xf.W  <-  p.Xf.W*(1 - u)
	q.YI.W  <-  q.YI.W + p.YI.W*u
	p.YI.W  <-  p.YI.W*(1 - u)

	(q.YI.W*(1 - s*(h*p.Xf.W + q.Xf.W) ) / 
			(YI - q.YI.W*s*(h*p.Xf.W + q.Xf.W) ) ) * YIprime
}

pI.Prime2.W  <-  function(YIprime, qYIprime.W) {
	YIprime - qYIprime.W
}


# Recursions for loci where inversion captures a Del. allele

q.Prime2.Xf.D  <-  function(q.Xf.D, q.Xm.D, u, h, s) {
	(2*(-1 + h*s)*u + q.Xm.D*(-1 + u + s*(h + 2*u - 3*h*u)) +q.Xf.D*(-1 + u + s*(h + 2*q.Xm.D - 2*h*q.Xm.D + (2 - 3*h + 4*(-1 + h)*q.Xm.D)*u))) / 
		(2*(-1 + q.Xm.D*s*u + q.Xf.D*s*(q.Xm.D + u - 2*q.Xm.D*u) + h*s*(q.Xf.D + q.Xm.D - 2*q.Xf.D*q.Xm.D + (2 - 3*q.Xm.D + q.Xf.D*(-3 + 4*q.Xm.D))*u)))
}

q.Prime2.Xm.D  <-  function(q.Xf.D, q.Y.D, YI, u, h, s) {
	(2*(-1 + s*(h + YI*(1 - h)))*u + q.Xf.D*(-1 - YI + 2*s*(q.Y.D + YI) - h*s*(-1 + 2*q.Y.D + YI) + u + YI*u + s*(2 - 4*q.Y.D - 4*YI + h*(-3 + 4*q.Y.D + 3*YI))*u) + q.Y.D*(-1 + u + s*(h + 2*u - 3*h*u))) / 
		(2*(-1 + s*(q.Y.D + YI)*u + q.Xf.D*s*(q.Y.D + YI + u - 2*(q.Y.D + YI)*u) + h*s*(q.Y.D + YI + 2*u - 3*(q.Y.D + YI)*u + q.Xf.D*(1 - 2*q.Y.D - 2*YI - 3*u + 4*(q.Y.D + YI)*u))))
}

q.Prime2.Y.D  <-  function(q.Xf.D, q.Y.D, YI, YIprime, u, h, s) {
	p.Xf.D  <-  (1 - q.Xf.D)
	p.Y.D   <-  1 - YI - q.Y.D
	# mutation
	q.Xf.D  <-  q.Xf.D + p.Xf.D*u
	p.Xf.D  <-  p.Xf.D*(1 - u)
	q.Y.D   <-  q.Y.D + p.Y.D*u
	p.Y.D   <-  p.Y.D*(1 - u)

	(q.Y.D*(1 - s*(h*p.Xf.D + q.Xf.D) ) / 
			(1 - YI - q.Y.D*s*(h*p.Xf.D + q.Xf.D) ) )* (1 - YIprime)
}




# Recursion for Inversion frequency*(Eqs. 1 & 2 from paper)
wBarY2  <-  function(n, r, s, h, YI, q.YI.W, q.Y.W, q.Y.D, q.Xf.W, q.Xf.D) {

	p.Xf.W  <-  (1 - q.Xf.W)
	p.Xf.D  <-  (1 - q.Xf.D)

	Q.Y.W   <-  q.Y.W/(1 - YI)
	Q.YI.W   <-  q.YI.W/YI
	Q.Y.D   <-  q.Y.D/(1 - YI)

	P.Y.W   <-  (1 - Q.Y.W)
	P.YI.W   <-  (1 - Q.YI.W)
	P.Y.D   <-  (1 - Q.Y.D)

	if(YI == 0) {
		Q.YI.W  <-  0
		P.YI.W  <-  0
	}
	if(YI == 1) {
		Q.Y.W  <-  0
		P.Y.W  <-  0
		Q.Y.D  <-  0
		P.Y.D  <-  0
	}

		   YI*( ( 1 - s*(h*p.Xf.D + q.Xf.D) )^(r) * ( 1 - s*(h*(P.YI.W*q.Xf.W + Q.YI.W*p.Xf.W) + Q.YI.W*q.Xf.W))^(n-r)) +
	 (1 - YI)*( ( 1 - s*(h*(p.Xf.D*Q.Y.D + q.Xf.D*P.Y.D) + q.Xf.D*Q.Y.D) )^(r) * ( 1 - s*(h*(p.Xf.W*Q.Y.W + q.Xf.W*P.Y.W) + q.Xf.W*Q.Y.W) )^(n - r) )
}

YI.prime2  <-  function(n, r, s, h, YI, q.YI.W, q.Y.W, q.Y.D, q.Xf.W, q.Xf.D) {
	p.Xf.W  <-  (1 - q.Xf.W)
	p.Xf.D  <-  (1 - q.Xf.D)

	Q.Y.W   <-  q.Y.W/(1 - YI)
	Q.YI.W   <-  q.YI.W/YI
	Q.Y.D   <-  q.Y.D/(1 - YI)

	P.Y.W   <-  (1 - Q.Y.W)
	P.YI.W   <-  (1 - Q.YI.W)
	P.Y.D   <-  (1 - Q.Y.D)

	if(YI == 0) {
		Q.YI.W  <-  0
		P.YI.W  <-  0
	}
	if(YI == 1) {
		Q.Y.W  <-  0
		P.Y.W  <-  0
		Q.Y.D  <-  0
		P.Y.D  <-  0
	}
	
	(YI*( ( 1 - s*(h*p.Xf.D + q.Xf.D) )^(r) * ( 1 - s*(h*(P.YI.W*q.Xf.W + Q.YI.W*p.Xf.W) + Q.YI.W*q.Xf.W))^(n-r))) / 
		(YI*( ( 1 - s*(h*p.Xf.D + q.Xf.D) )^(r) * ( 1 - s*(h*(P.YI.W*q.Xf.W + Q.YI.W*p.Xf.W) + Q.YI.W*q.Xf.W))^(n-r)) +
   (1 - YI)*( ( 1 - s*(h*(p.Xf.D*Q.Y.D + q.Xf.D*P.Y.D) + q.Xf.D*Q.Y.D) )^(r) * ( 1 - s*(h*(p.Xf.W*Q.Y.W + q.Xf.W*P.Y.W) + q.Xf.W*Q.Y.W) )^(n - r) ) )
}

#' Inversion relative fitness
invRelFit2  <-  function(n, r, s, h, q.YI.W, q.Y.W, q.Y.D, q.Xf.W, q.Xf.D, YI){
	p.Xf.W  <-  (1 - q.Xf.W)
	p.Xf.D  <-  (1 - q.Xf.D)

	Q.Y.W   <-  q.Y.W/(1 - YI)
	Q.YI.W   <-  q.YI.W/YI
	Q.Y.D   <-  q.Y.D/(1 - YI)

	P.Y.W   <-  (1 - Q.Y.W)
	P.YI.W   <-  (1 - Q.YI.W)
	P.Y.D   <-  (1 - Q.Y.D)

	if(YI == 0) {
		Q.YI.W  <-  0
		P.YI.W  <-  0
	}
	if(YI == 1) {
		Q.Y.W  <-  0
		P.Y.W  <-  0
		Q.Y.D  <-  0
		P.Y.D  <-  0
	}

		   ( ( 1 - s*(h*p.Xf.D + q.Xf.D) )^(r) * ( 1 - s*(h*(P.YI.W*q.Xf.W + Q.YI.W*p.Xf.W) + Q.YI.W*q.Xf.W))^(n-r)) /
	 ( ( 1 - s*(h*(p.Xf.D*Q.Y.D + q.Xf.D*P.Y.D) + q.Xf.D*Q.Y.D) )^(r) * ( 1 - s*(h*(p.Xf.W*Q.Y.W + q.Xf.W*P.Y.W) + q.Xf.W*Q.Y.W) )^(n - r) )
}




#################################
#################################

#' Function to calculate euclidean distance between two vectors
#' Convenient for checking if equilibrium has been reached
eucDist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
} 

#' Function to find equilibrium frequencies in the absence of inverions
findEqFullRec2  <- function(r, h, s, n, u, qHat_init) {

	# Empty Frequency Vectors
	YI.t     <-  c()
	qt.Xf.W  <-  c()
	qt.Xm.W  <-  c()
	qt.YI.W  <-  c()
	pt.YI.W  <-  c()
	qt.Y.W   <-  c()
	pt.Y.W   <-  c()
	qt.Xf.D  <-  c()
	qt.Xm.D  <-  c()
	qt.Y.D   <-  c()

	# Find equilibrium prior to inversion 
	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=0, q.YI.W=0, q.Y.W=qHat_init, q.Y.D=qHat_init, q.Xf.W=qHat_init, q.Xf.D=qHat_init), digits=16)
	qt.Xf.W[1]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qHat_init, q.Xm.W=qHat_init), digits=16)
	qt.Xm.W[1]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=0, q.YI.W=0, p.YI.W = 0, q.Y.W=qHat_init, q.Xf.W=qHat_init), digits=16)
	qt.YI.W[1]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=0, p.YI.W=0, q.Xf.W=qHat_init, YI=0, YIprime=0), digits=16) #, q.Y.W=qHat_init
	pt.YI.W[1]  <-  round(pI.Prime2.W(YIprime=0, qYIprime.W=0), digits=16)						
	qt.Y.W[1]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qHat_init, p.Y.W=(1- qHat_init), q.Xf.W=qHat_init, YI=0, YIprime=0), digits=16) #, p.YI.W=0, q.YI.W=0
	pt.Y.W[1]   <-  round(p.Prime2.Y.W(YIprime=0, qYprime.W=qt.Y.W[1]), digits=16)						
	qt.Xf.D[1]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qHat_init, q.Xm.D=qHat_init), digits=16)
	qt.Xm.D[1]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=0, q.Y.D=qHat_init, q.Xf.D=qHat_init), digits=16)
	qt.Y.D[1]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qHat_init, q.Y.D=qHat_init, YI=0, YIprime=0), digits=16)

	# Subsequent generations
	diff  <-  1
	i=2
	while(diff > 1e-12) {
		YI.t[i]     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=0, q.YI.W=0, q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		qt.Xf.W[i]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[i-1], q.Xm.W=qt.Xm.W[i-1]), digits=16)
		qt.Xm.W[i]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=0, q.YI.W=0, p.YI.W = 0, q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=16)
		qt.YI.W[i]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=0, p.YI.W=0, q.Xf.W=qt.Xf.W[i-1], YI=0, YIprime=0), digits=16) #, q.Y.W=qt.Y.W[i-1]
		pt.YI.W[i]  <-  round(pI.Prime2.W(YIprime=0, qYIprime.W=0), digits=16)						
		qt.Y.W[i]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[i-1], p.Y.W=(1-YI.t[i-1]-qt.Y.W[i-1]), q.Xf.W=qt.Xf.W[i-1], YI=0, YIprime=0), digits=16) #, p.YI.W=0, q.YI.W=0
		pt.Y.W[i]   <-  round(p.Prime2.Y.W(YIprime=0, qYprime.W=qt.Y.W[i]), digits=16)						
		qt.Xf.D[i]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Xm.D=qt.Xm.D[i-1]), digits=16)
		qt.Xm.D[i]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=0, q.Y.D=qt.Y.D[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		qt.Y.D[i]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Y.D=qt.Y.D[i-1], YI=0, YIprime=0), digits=16)

		diff        <-  eucDist(qt.Xf.W[i], qt.Xf.W[i-1])
		i  <-  i + 1
	}
	res  <-  list(
				  "q.Xf.W.init" = qt.Xf.W[i-1],
				  "q.Xm.W.init" = qt.Xm.W[i-1],
				  "q.Xf.D.init" = qt.Xf.D[i-1],
				  "q.Xm.D.init" = qt.Xm.D[i-1],
				  "q.Y.W.init"  = qt.Y.W[i-1],
				  "p.Y.W.init"  = pt.Y.W[i-1],
				  "q.Y.D.init"  = qt.Y.D[i-1]
				  )
	return(res)
}



#' Function to generate deterministic frequency dynamics 
#' for SLR-expanding inversions
#' Used to generate Fig.1
makeDeterministicFigSimData_SLR_Nei  <-  function(U = 0.02, r = 0, x = 0.2, h = 0.25, s = 0.01, generations = 10^4, ...) {
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
	YI.t     <-  c()
	qt.Xf.W  <-  c()
	qt.Xm.W  <-  c()
	qt.YI.W  <-  c()
	pt.YI.W  <-  c()
	qt.Y.W   <-  c()
	pt.Y.W   <-  c()
	qt.Xf.D  <-  c()
	qt.Xm.D  <-  c()
	qt.Y.D   <-  c()
	relW.YI.t  <-  c()
	Wbar.Y.t   <-  c()

	# find equilibrium frequencies prior to inversion mutation
	eqs   <-  findEqFullRec2(u=u, h=h, s=s, n=n, r=r, qHat_init=qHat)
	qHat  <-  eqs$q.Xf.W.init

	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=YI.0, q.YI.W=0, q.Y.W=qHat*(1-YI.0), q.Y.D=qHat*(1 - YI.0), q.Xf.W=qHat, q.Xf.D=qHat), digits=16)
	qt.Xf.W[1]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qHat, q.Xm.W=qHat), digits=16)
	qt.Xm.W[1]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=YI.0, q.YI.W=0, p.YI.W=YI.0, q.Y.W=qHat*(1-YI.0), q.Xf.W=qHat), digits=16)
	qt.YI.W[1]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=0, p.YI.W=YI.0, q.Xf.W=qHat, YI=YI.0, YIprime=YI.t[1]), digits=16)
	pt.YI.W[1]  <-  round(pI.Prime2.W(YIprime=YI.0, qYIprime.W=qt.YI.W[1]), digits=16)						
	qt.Y.W[1]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qHat*(1-YI.0), p.Y.W=(1-YI.0-qHat), q.Xf.W=qHat, YI=YI.0, YIprime=YI.t[1]), digits=16) #, p.YI.W=0, q.YI.W=0
	pt.Y.W[1]   <-  round(p.Prime2.Y.W(YIprime=YI.t[1], qYprime.W=qt.Y.W[1]), digits=16)						
	qt.Xf.D[1]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qHat, q.Xm.D=qHat), digits=16)
	qt.Xm.D[1]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=YI.0, q.Y.D=qHat*(1 - YI.0), q.Xf.D=qHat), digits=16)
	qt.Y.D[1]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qHat, q.Y.D= qHat*(1 - YI.0), YI=YI.0, YIprime=YI.t[1]), digits=16)
	relW.YI.t[1]  <-  round(invRelFit2(n=n, r=r, s=s, h=h, YI=YI.0, q.YI.W=0, q.Y.W=eqs$q.Y.W.init*(1 - YI.0), q.Y.D=eqs$q.Y.D.init*(1 - YI.0), q.Xf.W=eqs$q.Xf.W.init, q.Xf.D=eqs$q.Xf.D.init), digits=12)
	Wbar.Y.t[1]   <-  round(wBarY2(n=n, r=r, s=s, h=h, YI=YI.0, q.YI.W=0, q.Y.W=eqs$q.Y.W.init*(1 - YI.0), q.Y.D=eqs$q.Y.D.init*(1 - YI.0), q.Xf.W=eqs$q.Xf.W.init, q.Xf.D=eqs$q.Xf.D.init), digits=12)


	# Subsequent generations
	i=2
	while(i < generations+1) {
		YI.t[i]     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.YI.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		qt.Xf.W[i]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[i-1], q.Xm.W=qt.Xm.W[i-1]), digits=16)
		qt.Xm.W[i]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=YI.t[i-1], q.YI.W=qt.YI.W[i-1], p.YI.W=pt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=16)
		qt.YI.W[i]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=qt.YI.W[i-1], p.YI.W=pt.YI.W[i-1], q.Xf.W=qt.Xf.W[i-1], YI=YI.t[i-1], YIprime=YI.t[i]), digits=16) #, q.Y.W=qt.Y.W[i-1]
		pt.YI.W[i]  <-  round(pI.Prime2.W(YIprime=YI.t[i], qYIprime.W=qt.YI.W[i]), digits=16)						
		qt.Y.W[i]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[i-1], p.Y.W=pt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1], YI=YI.t[i-1], YIprime=YI.t[i]), digits=16) #, p.YI.W=0, q.YI.W=0
		pt.Y.W[i]   <-  round(p.Prime2.Y.W(YIprime=YI.t[i], qYprime.W=qt.Y.W[i]), digits=16)						
		qt.Xf.D[i]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Xm.D=qt.Xm.D[i-1]), digits=16)
		qt.Xm.D[i]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=YI.t[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		qt.Y.D[i]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Y.D=qt.Y.D[i-1], YI=YI.t[i-1], YIprime=YI.t[i]), digits=16)
		relW.YI.t[i]  <-  round(invRelFit2(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.YI.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		Wbar.Y.t[i]   <-  round(wBarY2(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.YI.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
		i  <-  i + 1
	}

	# Tack on initial freq.
	YI.t       <-  c(YI.0,YI.t)
	qt.Xf.W    <-  c(eqs$q.Xf.W.init, qt.Xf.W)
	qt.Xm.W    <-  c(eqs$q.Xm.W.init, qt.Xm.W)
	qt.Xf.D    <-  c(eqs$q.Xf.D.init, qt.Xf.D)
	qt.Xm.D    <-  c(eqs$q.Xm.D.init, qt.Xm.D)
	qt.YI.W    <-  c(0, qt.YI.W)
	pt.YI.W    <-  c(0, pt.YI.W)
	qt.Y.W     <-  c(eqs$q.Y.W.init, qt.Y.W)
	pt.Y.W     <-  c(eqs$p.Y.W.init, pt.Y.W)
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
							"pt.YI.W"     =  pt.YI.W,
							"qt.Y.W"      =  qt.Y.W,
							"pt.Y.W"      =  pt.Y.W,
							"qt.Y.D"      =  qt.Y.D,
							"rel.w.YI.t"  =  relW.YI.t,
							"Wbar.Y.t"    =  Wbar.Y.t
		 					)
	return(results)
}






##########################################################################
#' Wright-Fisher Simulations to estimate Inversion fixation probabilities


#' Function to estimate Pr(fix | x) for different inversion sizes (x)
#' EXPANDED WF simulation model: uses full set of deterministic
#' recursions to model changes in deleterious allele frequencies.
makeDataPrFixInvSize_SLR_Nei  <-  function(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
											nTot = 10^4, N = 10^4, Nfname = "") {

	# Containers
	PrFix      <-  c()
	rFixedInv  <-  c()

	# Empty Frequency Vectors
	YI.t     <-  c()
	qt.Xf.W  <-  c()
	qt.Xm.W  <-  c()
	qt.YI.W  <-  c()
	pt.YI.W  <-  c()
	qt.Y.W   <-  c()
	pt.Y.W   <-  c()
	qt.Xf.D  <-  c()
	qt.Xm.D  <-  c()
	qt.Y.D   <-  c()

	# Containers for # del mutations on
	# fixed inversions 
	rFixedInvTab  <-  data.frame()
	rInvSize      <-  c()
	rNs           <-  c()
	rUs           <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	YI.0   <-  2/N

	#number of simulations 
	sims  = 200*N


	# Loop over Us factor values
	for(j in 1:length(U.vals)) {
		# mutation rate & initial qHat
		U     <-  U.vals[j]
		u     <-  U/nTot
		qHat  <-  (U/(nTot*h*s))
		init_eqs  <-  findEqFullRec2(u=u, h=h, s=s, n=1, r=0, qHat_init=qHat)
		qHat  <-  init_eqs$q.Xf.W.init

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
				YI.t[1]     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=YI.0, q.YI.W=0, q.Y.W=qHat*(1-YI.0), q.Y.D=qHat*(1 - YI.0), q.Xf.W=qHat, q.Xf.D=qHat), digits=16)
				qt.Xf.W[1]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qHat, q.Xm.W=qHat), digits=16)
				qt.Xm.W[1]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=YI.0, q.YI.W=0, p.YI.W=YI.0, q.Y.W=qHat*(1-YI.0), q.Xf.W=qHat), digits=16)
				qt.YI.W[1]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=0, p.YI.W=YI.0, q.Xf.W=qHat, YI=YI.0, YIprime=YI.t[1]), digits=16)
				pt.YI.W[1]  <-  round(pI.Prime2.W(YIprime=YI.0, qYIprime.W=qt.YI.W[1]), digits=16)						
				qt.Y.W[1]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qHat*(1-YI.0), p.Y.W=(1-YI.0-qHat), q.Xf.W=qHat, YI=YI.0, YIprime=YI.t[1]), digits=16) #, p.YI.W=0, q.YI.W=0
				pt.Y.W[1]   <-  round(p.Prime2.Y.W(YIprime=YI.t[1], qYprime.W=qt.Y.W[1]), digits=16)						
				qt.Xf.D[1]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qHat, q.Xm.D=qHat), digits=16)
				qt.Xm.D[1]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=YI.0, q.Y.D=qHat*(1 - YI.0), q.Xf.D=qHat), digits=16)
				qt.Y.D[1]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qHat, q.Y.D= qHat*(1 - YI.0), YI=YI.0, YIprime=YI.t[1]), digits=16)

				# Run forward simulation
				while(YI.t[1]*(1 - YI.t[1]) > 0) {
					#expected frequency after selection
					YI.t[2]     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=YI.t[1], q.YI.W=qt.YI.W[1], q.Y.W=qt.Y.W[1], q.Y.D=qt.Y.D[1], q.Xf.W=qt.Xf.W[1], q.Xf.D=qt.Xf.D[1]), digits=16)
					qt.Xf.W[2]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[1], q.Xm.W=qt.Xm.W[1]), digits=16)
					qt.Xm.W[2]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=YI.t[1], q.YI.W=qt.YI.W[1], p.YI.W=pt.YI.W[1], q.Y.W=qt.Y.W[1], q.Xf.W=qt.Xf.W[1]), digits=16)
					qt.YI.W[2]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=qt.YI.W[1], p.YI.W=pt.YI.W[1], q.Xf.W=qt.Xf.W[1], YI=YI.t[1], YIprime=YI.t[2]), digits=16) #, q.Y.W=qt.Y.W[1]
					pt.YI.W[2]  <-  round(pI.Prime2.W(YIprime=YI.t[2], qYIprime.W=qt.YI.W[2]), digits=16)						
					qt.Y.W[2]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[1], p.Y.W=pt.Y.W[1], q.Xf.W=qt.Xf.W[1], YI=YI.t[1], YIprime=YI.t[2]), digits=16) #, p.YI.W=0, q.YI.W=0
					pt.Y.W[2]   <-  round(p.Prime2.Y.W(YIprime=YI.t[2], qYprime.W=qt.Y.W[2]), digits=16)						
					qt.Xf.D[2]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[1], q.Xm.D=qt.Xm.D[1]), digits=16)
					qt.Xm.D[2]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=YI.t[1], q.Y.D=qt.Y.D[1], q.Xf.D=qt.Xf.D[1]), digits=16)
					qt.Y.D[2]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[1], q.Y.D=qt.Y.D[1], YI=YI.t[1], YIprime=YI.t[2]), digits=16)

					#binomial sampling
					YI.t[1]     <-  rbinom(1, (N/2), YI.t[2]) / (N/2)
					
					# Frame Shift
					qt.Xf.W[1]  <-  qt.Xf.W[2]
					qt.Xm.W[1]  <-  qt.Xm.W[2]
					qt.YI.W[1]  <-  qt.YI.W[2]
					pt.YI.W[1]  <-  pt.YI.W[2]
					qt.Y.W[1]   <-  qt.Y.W[2] 
					pt.Y.W[1]   <-  pt.Y.W[2] 
					qt.Xf.D[1]  <-  qt.Xf.D[2]
					qt.Xm.D[1]  <-  qt.Xm.D[2]
					qt.Y.D[1]   <-  qt.Y.D[2] 
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
	filename  <-  paste("./data/RECODE/PrFixFig_SLR-Nei_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

	filename  <-  paste("./data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h", h, "_s", s, Nfname, ".csv", sep="")
	r.d       <-  as.data.frame(cbind(rNs, rUs, rInvSize, rFixedInvTab))
	write.csv(r.d, file=filename, row.names=FALSE)

}









####################################################################
####################################################################
#' Function to generate timeseries of allele frequencies during W-F
#' Simulations to look at what is going on with feedback between 
#' inversion frequency & del. allele frequencies on each chromosome
#' class
makeData_SLR_Nei_WFDynamics  <-  function(h = 0.25, s = 0.01, U = 0.02, nFix = 1,
											nTot = 10^4, N = 10^5, Nfname = "") {

	# Empty Vectors for 
	# Allele/Inversion timeseries
	YI.tseries      <-  c()
	q.Xf.W.tseries  <-  c()
	q.Xm.W.tseries  <-  c()
	q.YI.W.tseries  <-  c()
	p.YI.W.tseries  <-  c()
	q.Y.W.tseries   <-  c()
	p.Y.W.tseries   <-  c()
	q.Xf.D.tseries  <-  c()
	q.Xm.D.tseries  <-  c()
	q.Y.D.tseries   <-  c()
	invRelFit.tseries  <-  c()
	fixCounter      <-  c()
	rTracker        <-  c()
	invSizeVar      <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	YI.0   <-  1/(2*N)

	# maxSims
	maxSims  <-  10*N

		# mutation rate & initial qHat
		u     <-  U/nTot
		qHat  <-  (U/(nTot*h*s))
		# Find equilibrium prior to inversion 
		init_eqs  <-  findEqFullRec2(u=u, h=h, s=s, n=1, r=0, qHat_init=qHat)
		qHat  <-  init_eqs$q.Xf.W.init

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
				qt.YI.W  <-  c()
				pt.YI.W  <-  c()
				qt.Y.W   <-  c()
				pt.Y.W   <-  c()
				qt.Xf.D  <-  c()
				qt.Xm.D  <-  c()
				qt.Y.D   <-  c()
				invRelFit.t  <-  c()

				# Assign frequencies
				# Note implicit assumption of equal initial
				# frequencies in XOv, XSp, Y chromosomes
				YI.t[1]     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=YI.0, q.YI.W=0, q.Y.W=qHat*(1-YI.0), q.Y.D=qHat*(1 - YI.0), q.Xf.W=qHat, q.Xf.D=qHat), digits=16)
				qt.Xf.W[1]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qHat, q.Xm.W=qHat), digits=16)
				qt.Xm.W[1]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=YI.0, q.YI.W=0, p.YI.W=YI.0, q.Y.W=qHat*(1-YI.0), q.Xf.W=qHat), digits=16)
				qt.YI.W[1]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=0, p.YI.W=YI.0, q.Xf.W=qHat, YI=YI.0, YIprime=YI.t[1]), digits=16)
				pt.YI.W[1]  <-  round(pI.Prime2.W(YIprime=YI.0, qYIprime.W=qt.YI.W[1]), digits=16)						
				qt.Y.W[1]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qHat*(1-YI.0), p.Y.W=(1-YI.0-qHat), q.Xf.W=qHat, YI=YI.0, YIprime=YI.t[1]), digits=16) #, p.YI.W=0, q.YI.W=0
				pt.Y.W[1]   <-  round(p.Prime2.Y.W(YIprime=YI.t[1], qYprime.W=qt.Y.W[1]), digits=16)						
				qt.Xf.D[1]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qHat, q.Xm.D=qHat), digits=16)
				qt.Xm.D[1]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=YI.0, q.Y.D=qHat*(1 - YI.0), q.Xf.D=qHat), digits=16)
				qt.Y.D[1]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qHat, q.Y.D= qHat*(1 - YI.0), YI=YI.0, YIprime=YI.t[1]), digits=16)
				invRelFit.t[1]  <-  round(invRelFit2(n=n, r=r, s=s, h=h, YI=YI.0, q.YI.W=0, q.Y.W=qHat*(1-YI.0), q.Y.D=qHat*(1 - YI.0), q.Xf.W=qHat, q.Xf.D=qHat), digits=16)				

				# Run forward simulation
				i = 2
				while(YI.t[i-1]*(1 - YI.t[i-1]) > 0) {
					# expected inversion frequency after selection
					YI.next     <-  round(YI.prime2(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.YI.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
					# binomial of inversion sampling
					YI.t[i]     <-  rbinom(1, (N/2), YI.next) / (N/2)
					# Calculate relative fitness
					invRelFit.t[i]  <-  round(invRelFit2(n=n, r=r, s=s, h=h, YI=YI.t[i-1], q.YI.W=qt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.W=qt.Xf.W[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
					# expected frequency of del. mutations in next generation
					qt.Xf.W[i]  <-  round(q.Prime2.Xf.W(u=u, s=s, h=h, q.Xf.W=qt.Xf.W[i-1], q.Xm.W=qt.Xm.W[i-1]), digits=16)
					qt.Xm.W[i]  <-  round(q.Prime2.Xm.W(u=u, s=s, h=h, YI=YI.t[i-1], q.YI.W=qt.YI.W[i-1], p.YI.W=pt.YI.W[i-1], q.Y.W=qt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1]), digits=16)
					qt.YI.W[i]  <-  round(q.Prime2.YI.W(u=u, s=s, h=h, q.YI.W=qt.YI.W[i-1], p.YI.W=pt.YI.W[i-1], q.Xf.W=qt.Xf.W[i-1], YI=YI.t[i-1], YIprime=YI.t[i]), digits=16)
					pt.YI.W[i]  <-  round(pI.Prime2.W(YIprime=YI.t[i], qYIprime.W=qt.YI.W[i]), digits=16)						
					qt.Y.W[i]   <-  round(q.Prime2.Y.W(u=u, s=s, h=h, q.Y.W=qt.Y.W[i-1], p.Y.W=pt.Y.W[i-1], q.Xf.W=qt.Xf.W[i-1], YI=YI.t[i-1], YIprime=YI.t[i]), digits=16)
					pt.Y.W[i]   <-  round(p.Prime2.Y.W(YIprime=YI.t[i], qYprime.W=qt.Y.W[i]), digits=16)						
					qt.Xf.D[i]  <-  round(q.Prime2.Xf.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Xm.D=qt.Xm.D[i-1]), digits=16)
					qt.Xm.D[i]  <-  round(q.Prime2.Xm.D(u=u, s=s, h=h, YI=YI.t[i-1], q.Y.D=qt.Y.D[i-1], q.Xf.D=qt.Xf.D[i-1]), digits=16)
					qt.Y.D[i]   <-  round(q.Prime2.Y.D(u=u, s=s, h=h, q.Xf.D=qt.Xf.D[i-1], q.Y.D=qt.Y.D[i-1], YI=YI.t[i-1], YIprime=YI.t[i]), digits=16)

					i  <-  i + 1
				}
				if(YI.t[i-1] == 1) {
					fix        <-  fix + 1 
	
					# Concatenate timeseries
					YI.tseries      <-  c(YI.tseries, YI.t)
					invRelFit.tseries  <-  c(invRelFit.tseries, invRelFit.t)
					q.Xf.W.tseries  <-  c(q.Xf.W.tseries, qt.Xf.W)
					q.Xm.W.tseries  <-  c(q.Xm.W.tseries, qt.Xm.W)
					q.YI.W.tseries  <-  c(q.YI.W.tseries, qt.YI.W)
					p.YI.W.tseries  <-  c(p.YI.W.tseries, pt.YI.W)
					q.Y.W.tseries   <-  c(q.Y.W.tseries, qt.Y.W)
					p.Y.W.tseries   <-  c(p.Y.W.tseries, pt.Y.W)
					q.Xf.D.tseries  <-  c(q.Xf.D.tseries, qt.Xf.D)
					q.Xm.D.tseries  <-  c(q.Xm.D.tseries, qt.Xm.D)
					q.Y.D.tseries   <-  c(q.Y.D.tseries, qt.Y.D)
					fixCounter      <-  c(fixCounter, rep(fix, times=length(YI.t)))
					rTracker        <-  c(rTracker, rep(r, times=length(YI.t)))
					invSizeVar      <-  c(invSizeVar, rep(invSize[k], times = length(YI.t)))
				}

				rep  <-  rep+1
			}

			# Print Progress
			cat('\r', paste("x: ", k, "/", length(invSize), " complete", sep=""))

		}

	# Export Results Dataframe
	filename  <-  paste("./data/RECODE/SLR-Nei-WF-Dynamics_h", h, "_s", s, "_U", U, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"       =  rep(h, times=length(YI.tseries)),
										"s"       =  rep(s, times=length(YI.tseries)),
										"N"       =  rep(N, times=length(YI.tseries)),
										"U"       =  rep(U, times=length(YI.tseries)),
										"x"       =  invSizeVar,
										"YI.t"    = YI.tseries,
										"invRelFit" = invRelFit.tseries,
										"qt.Xf.W" = q.Xf.W.tseries,
										"qt.Xm.W" = q.Xm.W.tseries,
										"qt.YI.W" = q.YI.W.tseries,
										"pt.YI.W" = p.YI.W.tseries,
										"qt.Y.W"  = q.Y.W.tseries,
										"pt.Y.W"  = p.Y.W.tseries,
										"qt.Xf.D" = q.Xf.D.tseries,
										"qt.Xm.D" = q.Xm.D.tseries,
										"qt.Y.D"  = q.Y.D.tseries,
										"fixRep"  =  fixCounter,
										"r"       =  rTracker
										)
	write.csv(d, file=filename, row.names=FALSE)

}
