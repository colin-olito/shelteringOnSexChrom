# Simulation code: 
# Inversions expanding SLR w/ recessive deleterious mutations

#' Here we present exact recursions for the deterministic 
#' change in deleterious allele as well as inversion 
#' frequencies. Functions are provided that are necessary
#' To produce Figs. 1 & 2 in the main text.
#' 


##  EXACT RECURSIONS
wBarY  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del) {
	pt.I.wt    <-  1 - qt.I.wt
	pt.Y.del   <-  1 - qt.Y.del
	pt.Y.wt    <-  1 - qt.Y.wt
		  YI.t*(( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )^r) * (1 - s*(h*(pt.I.wt*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))^(n-r) +
	(1 - YI.t)*( (1 - s*(h*((1 - XaOv.t.del)*qt.Y.del + XaOv.t.del*pt.Y.del) + XaOv.t.del*qt.Y.del)^r) * ((1 - s*(h*((1 - XaOv.t.wt)*qt.Y.wt + XaOv.t.wt*pt.Y.wt) + XaOv.t.wt*qt.Y.wt))^(n - r)) )
}

YI.prime  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del, ...) {
	(YI.t*(( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )^r) * (1 - s*(h*((1 - qt.I.wt)*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))^(n-r)) / wBarY(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qt.I.wt, qt.Y.wt=qt.Y.wt, qt.Y.del=qt.Y.del, XaOv.t.wt=XaOv.t.wt, XaOv.t.del=XaOv.t.del)
}

qI.wt.prime  <-  function(n, r, u, s, h, YI.t, qt.I.wt, XaOv.t.wt, XaOv.t.del) {
		(qt.I.wt*(1 - u - s*(XaOv.t.wt + u*(1 - 2*XaOv.t.wt) + h*(1 - XaOv.t.wt)*(1 - 2*u))) + u*(1 - s*(XaOv.t.wt + h*(1 - XaOv.t.wt)))) / 
			(1 - s*XaOv.t.wt*u - qt.I.wt*s*(XaOv.t.wt + u*(1 - 2*XaOv.t.wt)) - h*s*(qt.I.wt + XaOv.t.wt*(1 - 2*qt.I.wt) + u*(2 - 3*XaOv.t.wt - qt.I.wt*(3 - 4*XaOv.t.wt))))
}

qY.wt.prime  <-  function(u, s, h, qt.Y.wt, XaOv.t.wt) {
	(qt.Y.wt*(1 - u - s*(h + 2*XaOv.t.wt*(1 - h) + u*(2 - 3*h - 4*XaOv.t.wt*(1 - h)))) + XaOv.t.wt*(1 - u - s*(h + u*(2 - 3*h))) + 2*u*(1 - h*s)) / 
		(2*(1 - s*XaOv.t.wt*u - qt.Y.wt*s*(XaOv.t.wt + u*(1 - 2*XaOv.t.wt)) - h*s*(qt.Y.wt + XaOv.t.wt*(1 - 2*qt.Y.wt) + u*(2 - 3*XaOv.t.wt - qt.Y.wt*(3 - 4*XaOv.t.wt)))))
}

qY.del.prime  <-  function(u, s, h, qt.Y.del, XaOv.t.del) {
	(qt.Y.del*(1 - u - s*(h + 2*XaOv.t.del*(1 - h ) + u*(2 - 3*h - 4*XaOv.t.del*(1 - h)))) + XaOv.t.del*(1 - u - s*(h + u*(2 - 3*h))) + 2*u*(1 - h*s)) / 
		(2*(1 - s*XaOv.t.del*u - qt.Y.del*s*(XaOv.t.del + u*(1 - 2*XaOv.t.del)) - h*s*(qt.Y.del + XaOv.t.del*(1 - 2*qt.Y.del) + u*(2 - 3*XaOv.t.del - qt.Y.del*(3 - 4*XaOv.t.del)))))
}

XaOv.wt.prime  <-  function(u, s, h, XaOv.t.wt, XaSp.t.wt) {
	(2*u*(1 - h*s) + XaSp.t.wt*(1 - u - s*(h + u*(2 - 3*h))) + XaOv.t.wt*(1 - u - s*(h + 2*XaSp.t.wt*(1 - h) + u*(2 - 3*h - 4*XaSp.t.wt*(1 - h))))) / 
		(2*(1 - s*(XaSp.t.wt*u + XaOv.t.wt*(XaSp.t.wt + u*(1 - 2*XaSp.t.wt))) - h*s*(XaOv.t.wt + XaSp.t.wt*(1 - 2*XaOv.t.wt) + u*(2 - 3*XaSp.t.wt - XaOv.t.wt*(3 - 4*XaSp.t.wt)))))
}

XaOv.del.prime  <-  function(u, s, h, XaOv.t.del, XaSp.t.del) {
	(XaOv.t.del*(1 - u - s*(h + 2*XaSp.t.del*(1 - h) + u*(2 - 3*h - 4*XaSp.t.del*(1 - h)))) + XaSp.t.del*(1 - u - s*(2*u + h*(1 - 3*u))) + 2*u*(1 - h*s)) / 
		(2*(1 - s*XaSp.t.del*u - s*XaOv.t.del*(XaSp.t.del + u*(1 - 2*XaSp.t.del)) - h*s*(XaOv.t.del + XaSp.t.del*(1 - 2*XaOv.t.del) + u*(2 - 3*XaSp.t.del - XaOv.t.del*(3 - 4*XaSp.t.del)))))
}

XaSp.wt.prime  <-  function(u, s, h, YI.t, qt.I.wt, qt.Y.wt, XaOv.t.wt) {
(XaOv.t.wt*(1 + YI.t*(1 - 2*qt.I.wt*s) - h*s*(1 + YI.t*(1 - 2*qt.I.wt))) + u*(2 - 2*h*s - XaOv.t.wt*(1 + 2*s - 3*h*s) - YI.t*(XaOv.t.wt*(1 - h*s) + 2*(1 - h)*qt.I.wt*s*(1 - 2*XaOv.t.wt))) + qt.Y.wt*(1 - YI.t)*(1 - u - s*(h + XaOv.t.wt*(2 - 2*h) + u*(2 - 3*h - 4*XaOv.t.wt*(1 - h))))) / 
	(2 - s*XaOv.t.wt*(2*qt.Y.wt*(1 - YI.t) + 2*qt.I.wt*YI.t) - 2*h*s*(XaOv.t.wt + qt.Y.wt*(1 - 2*XaOv.t.wt) + YI.t*(qt.I.wt - qt.Y.wt)*(1 - 2*XaOv.t.wt)) - 2*s*u*(2*h + qt.Y.wt*(1 - 3*h) + XaOv.t.wt*(1 - 3*h - qt.Y.wt*(2 - 4*h)) + YI.t*(qt.I.wt - qt.Y.wt)*(1 - 2*XaOv.t.wt - h*(3 - 4*XaOv.t.wt))))
}

XaSp.del.prime  <-  function(u, s, h, YI.t, qt.Y.del, XaOv.t.del) {
(XaOv.t.del*(1 - h*s + qt.Y.del*YI.t*(1 - s*(2 - h))) + u*(2 - XaOv.t.del*(1 + qt.Y.del*YI.t) - h*s*(2 - 3*XaOv.t.del)*(1 - qt.Y.del*YI.t) - 2*s*(XaOv.t.del + qt.Y.del*YI.t*(1 - 2*XaOv.t.del))) + qt.Y.del*(1 - YI.t)*(1 - u - s*(h + 2*XaOv.t.del*(1 - h) + u*(2 - 3*h - 4*XaOv.t.del*(1 - h))))) / 
	(2 - 2*qt.Y.del*s*XaOv.t.del*(1 - YI.t) - 2*s*(qt.Y.del*XaOv.t.del*YI.t - h*(qt.Y.del + XaOv.t.del*(1 - 2*qt.Y.del) + YI.t*(qt.Y.del - qt.Y.del)*(1 - 2*XaOv.t.del))) - 2*s*u*(qt.Y.del + h*(2 - 3*qt.Y.del) + XaOv.t.del*(1 - 3*h - 2*qt.Y.del*(1 - 2*h)) + YI.t*(qt.Y.del - qt.Y.del)*(1 - 2*XaOv.t.del - h*(3 - 4*XaOv.t.del))))
}

eucDist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
} 

#################
#' Function to generate fig.1 data
findEq  <- function(r, h, s, n, u, qHat, qHatDel) {

	# Empty Frequency Vectors
	YI.t        <- c()
	XaOv.wt.t   <- c()
	XaSp.wt.t   <- c()
	qI.wt.t     <- c()
	qY.wt.t     <- c()
	XaOv.del.t  <- c()
	XaSp.del.t  <- c()
	qY.del.t    <- c()

	# Find equilibrium prior to inversion 
	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=0, qt.I.wt=0, qt.Y.wt=qHat, qt.Y.del=qHatDel, XaOv.t.wt=qHat, XaOv.t.del=qHatDel), digits=9)
	XaOv.wt.t[1]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=qHat, XaSp.t.wt=qHat), digits=9)
	XaSp.wt.t[1]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=0, qt.I.wt=0, qt.Y.wt=qHat, XaOv.t.wt=qHat), digits=9)
	qI.wt.t[1]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=0, qt.I.wt=0, XaOv.t.wt=qHat, XaOv.t.del=qHatDel), digits=9)
	qY.wt.t[1]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qHat, XaOv.t.wt=qHat), digits=9)
	XaOv.del.t[1]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=qHatDel, XaSp.t.del=qHatDel), digits=9)
	XaSp.del.t[1]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=0, qt.Y.del=qHatDel, XaOv.t.del=qHatDel), digits=9)
	qY.del.t[1]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qHatDel, XaOv.t.del=qHatDel), digits=9)

	# Subsequent generations
	i=2
	diff  <-  1
	while(diff > 1e-8) {
		YI.t[i]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		XaOv.wt.t[i]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t[i-1], XaSp.t.wt=XaSp.wt.t[i-1]), digits=9)
		XaSp.wt.t[i]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		XaOv.del.t[i]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t[i-1], XaSp.t.del=XaSp.del.t[i-1]), digits=9)
		XaSp.del.t[i]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qI.wt.t[i]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qY.wt.t[i]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		qY.del.t[i]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		diff  <-  eucDist(c(XaOv.wt.t[i-1],XaSp.wt.t[i-1],XaOv.del.t[i-1],XaSp.del.t[i-1],qY.wt.t[i-1],qY.del.t[i-1]),
						  c(XaOv.wt.t[i],XaSp.wt.t[i],XaOv.del.t[i],XaSp.del.t[i],qY.wt.t[i],qY.del.t[i]))
		i  <-  i + 1
	}
	res  <-  list(
				  "XaOv.wt.init" = XaOv.wt.t[i-1],
				  "XaSp.wt.init" = XaSp.wt.t[i-1],
				  "XaOv.del.init" = XaOv.del.t[i-1],
				  "XaSp.del.init" = XaSp.del.t[i-1],
				  "qI.wt.init" = qI.wt.t[i-1],
				  "qY.wt.init" = qY.wt.t[i-1],
				  "qY.del.init" = qY.del.t[i-1]
				  )
	return(res)
}

makeDeterministicFigSimData  <-  function(r = 0, x = 0.2, h = 0.1, s = 0.01, Ufactor = 2, generations = 10^4, ...) {
	# Parameters
	U     <-  Ufactor*s
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
	YI.t        <- c()
	XaOv.wt.t   <- c()
	XaSp.wt.t   <- c()
	qI.wt.t     <- c()
	qY.wt.t     <- c()
	XaOv.del.t  <- c()
	XaSp.del.t  <- c()
	qY.del.t    <- c()
	wbarYI.t    <- c()

	# Find equilibrium prior to inversion 
	eqs  <-  findEq(r=r, h=h, s=s, n=n, u=u, qHat=qHat, qHatDel=qHatDel)
	
	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=eqs$qY.wt.init, qt.Y.del=eqs$qY.del.init, XaOv.t.wt=eqs$XaOv.wt.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	XaOv.wt.t[1]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=eqs$XaOv.wt.init, XaSp.t.wt=eqs$XaSp.wt.init), digits=9)
	XaSp.wt.t[1]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=eqs$qY.wt.init, XaOv.t.wt=eqs$XaOv.wt.init), digits=9)
	qI.wt.t[1]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, XaOv.t.wt=eqs$XaOv.wt.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	qY.wt.t[1]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=eqs$qY.wt.init, XaOv.t.wt=eqs$XaOv.wt.init), digits=9)
	XaOv.del.t[1]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=eqs$XaOv.del.init, XaSp.t.del=eqs$XaSp.del.init), digits=9)
	XaSp.del.t[1]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.0, qt.Y.del=eqs$qY.del.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	qY.del.t[1]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=eqs$qY.del.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	wbarYI.t[1]    <-  round((YI.t[1]/YI.0), digits=9)
	# Subsequent generations
	i=2
	while(i < generations+1) {
		YI.t[i]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		XaOv.wt.t[i]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t[i-1], XaSp.t.wt=XaSp.wt.t[i-1]), digits=9)
		XaSp.wt.t[i]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		XaOv.del.t[i]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t[i-1], XaSp.t.del=XaSp.del.t[i-1]), digits=9)
		XaSp.del.t[i]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qI.wt.t[i]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qY.wt.t[i]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		qY.del.t[i]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		wbarYI.t[i]    <-  round((YI.t[i]/ YI.t[i-1]), digits=9)
		i  <-  i + 1
	}
	# Tack on initial freq.
	YI.t        <- c(YI.0,YI.t)
	XaOv.wt.t   <-  c(eqs$XaOv.wt.init, XaOv.wt.t)
	XaSp.wt.t   <-  c(eqs$XaSp.wt.init, XaSp.wt.t)
	XaOv.del.t  <-  c(eqs$XaOv.del.init, XaOv.del.t)
	XaSp.del.t  <-  c(eqs$XaSp.del.init, XaSp.del.t)
	qI.wt.t     <-  c(0, qI.wt.t)
	qY.wt.t     <-  c(eqs$qY.wt.init, qY.wt.t)
	qY.del.t    <-  c(eqs$qY.del.init, qY.del.t)
	wbarYI.t    <-  c(NA, wbarYI.t)

	# Return results as df
	results  <-  data.frame(
							"YI.t"        =  YI.t,
							"XaOv.wt.t"   =  XaOv.wt.t,
							"XaSp.wt.t"   =  XaSp.wt.t,
							"XaOv.del.t"  =  XaOv.del.t,
							"XaSp.del.t"  =  XaSp.del.t,
							"qI.wt.t"     =  qI.wt.t,
							"qY.wt.t"     =  qY.wt.t,
							"qY.del.t"    =  qY.del.t,
							"wbar.YI.t"   =  wbarYI.t
							)
	return(results)
}



########################################
## Multilocus Wright-Fisher Simulations
## Using exact recursions derived in Mathematica using Xf, Xm, Y
## inluding sampling variance for XaOv, XaSp, etc.
## 
wBarY.multi  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del) {
	pt.I.wt    <-  1 - qt.I.wt
	pt.Y.del   <-  1 - qt.Y.del
	pt.Y.wt    <-  1 - qt.Y.wt
		  YI.t*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )) * prod(1 - s*(h*(pt.I.wt*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt)) +
	(1 - YI.t)*( prod(1 - s*(h*((1 - XaOv.t.del)*qt.Y.del + XaOv.t.del*pt.Y.del) + XaOv.t.del*qt.Y.del)) * prod((1 - s*(h*((1 - XaOv.t.wt)*qt.Y.wt + XaOv.t.wt*pt.Y.wt) + XaOv.t.wt*qt.Y.wt))) )
}

YI.multi.prime  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del, ...) {
	(YI.t*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del))) * prod(1 - s*(h*((1 - qt.I.wt)*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))) / wBarY.multi(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qt.I.wt, qt.Y.wt=qt.Y.wt, qt.Y.del=qt.Y.del, XaOv.t.wt=XaOv.t.wt, XaOv.t.del=XaOv.t.del)
}

# Function to estimate Pr(fix | x) for different inversion sizes (x)
# Looping over different Population size and U/s ratios
makeDataPrFixInvSize  <-  function(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
																	 nTot = 10^4, N.vals = c(10^3, 10^4), Nfname = "") {

	# Containers
	PrFix       <-  c()
	YI.t        <-  c()
	XaOv.wt.t   <-  c()
	XaSp.wt.t   <-  c()
	qI.wt.t     <-  c()
	qY.wt.t     <-  c()
	XaOv.del.t  <-  c()
	XaSp.del.t  <-  c()
	qY.del.t    <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N  <-  N.vals[i]
		YI.0   <-  2/N

		#number of simulations 
		sims  = 100*N/2

		# Loop over Us factor values
		for(j in 1:length(Us.factor.vals)) {
			U     <-  Us.factor.vals[j]*s
			u     <-  U/nTot
			qHat  <-  (U/(nTot*h*s))

			# Loop over inversion size
			for(k in 1:length(invSize)) {
				# # loci captured by inversion
				n   <-  nTot*invSize[k]

				# counter for fixations
				fix   = 0

				# Loop over replicate simulations
				for(l in 1:sims){


					# Draw random value for # del. mutations captured by inversion (r) given x
					r   <-  rpois(1, lambda=(U*invSize[k]/(s*h)))

					# Assign frequencies
					# Note implicit assumption of equal initial
					# frequencies in XOv, XSp, Y chromosomes
					YI.t        <- YI.0
					XaOv.wt.t   <- qHat
					XaSp.wt.t   <- qHat
					qI.wt.t     <- 0
					qY.wt.t     <- qHat
					XaOv.del.t  <- qHat
					XaSp.del.t  <- qHat
					qY.del.t    <- qHat
# t=1
# plot(NA, ylim=c(0,1.001), xlim=c(0,N))
# abline(h=1)
# text(x=(N/2), y=0.9, labels=paste(l, sep=''))
					# Run forward simulation
					while(YI.t*(1 - YI.t) > 0) {
						#expected frequency after selection
						XaOv.wt.t   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
						XaOv.del.t  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
						XaSp.wt.t   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						XaSp.del.t  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qY.wt.t     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						qY.del.t    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qI.wt.t     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
						YI.sel      <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
#points(YI.sel/YI.t ~ t, col=2)

						#binomial sampling
						YI.t        <-  rbinom(1, (N/2), YI.sel)/(N/2)

#points(YI.t ~ t)
#points(qI.wt.t ~ t, col=2)
#t=t+1
					}
					fix  <-  fix + YI.t
				}
			PrFix  <-  c(PrFix, (fix/sims))
			cat('\r', paste("N: ", i, "/", length(N.vals), 
										", U: ", j, "/", length(Us.factor.vals), 
										", x: ", round(100*(k/length(invSize))), "% complete", sep=""))
			}

		}

	}

	# Index variables
	Ns        <-  rep(N.vals, each=(length(Us.factor.vals)*length(invSize)))
	Us        <-  rep((Us.factor.vals*s), each=(length(invSize)), times=length(N.vals))
	Ufac      <-  rep((Us.factor.vals), each=(length(invSize)), times=length(N.vals))
	invSizes  <-  rep(invSize, times=length(N.vals)*length(Us.factor.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./data/PrFixFig_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"Ufac"   =  Ufac,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}



########################################
##  RECURSIONS W/ RANDOM GAMETE SAMPLING FOR X CHROMOSOMES
# Function to estimate Pr(fix | x) for different inversion sizes (x)
# Looping over different Population size and U/s ratios
makeDataPrFixInvSizeDetqI  <-  function(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
																	 nTot = 10^4, N.vals = c(10^3, 10^4), Nfname = "") {

	# Containers
	PrFix       <-  c()
	YI.t        <-  c()
	XaOv.wt.t   <-  c()
	XaSp.wt.t   <-  c()
	qI.wt.t     <-  c()
	qY.wt.t     <-  c()
	XaOv.del.t  <-  c()
	XaSp.del.t  <-  c()
	qY.del.t    <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N  <-  N.vals[i]
		YI.0   <-  2/N

		#number of simulations 
		sims  = 100*N/2

		# Loop over Us factor values
		for(j in 1:length(Us.factor.vals)) {
			U     <-  Us.factor.vals[j]*s
			u     <-  U/nTot
			qHat  <-  (U/(nTot*h*s))

			# Loop over inversion size
			for(k in 1:length(invSize)) {
				# # loci captured by inversion
				n   <-  nTot*invSize[k]

				# counter for fixations
				fix   = 0

				# Loop over replicate simulations
				for(l in 1:sims){


					# Draw random value for # del. mutations captured by inversion (r) given x
					r   <-  rpois(1, lambda=(U*invSize[k]/(s*h)))

					# Assign frequencies
					# Note implicit assumption of equal initial
					# frequencies in XOv, XSp, Y chromosomes
					YI.t        <- YI.0
					XaOv.wt.t   <- qHat
					XaSp.wt.t   <- qHat
					qI.wt.t     <- 0
					qY.wt.t     <- qHat
					XaOv.del.t  <- qHat
					XaSp.del.t  <- qHat
					qY.del.t    <- qHat
# t=1
# plot(NA, ylim=c(0,1.001), xlim=c(0,N))
# abline(h=1)
# text(x=(N/2), y=0.9, labels=paste(l, sep=''))
					# Run forward simulation
					while(YI.t*(1 - YI.t) > 0) {
						#expected frequency after selection
						#expected frequency after selection
						XaOv.wt.sel   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
						XaOv.del.sel  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
						XaSp.wt.sel   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						XaSp.del.sel  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qY.wt.sel     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						qY.del.sel    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qI.wt.t       <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
						YI.sel        <-  round(YI.multi.prime(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
#points(YI.sel/YI.t ~ t, col=2)

						#binomial sampling
						XaOv.wt.t   <-  rbinom((n-r), N, prob=XaOv.wt.sel)/N
						XaOv.del.t  <-  rbinom(r, N, XaOv.del.sel)/N
						YI.t        <-  rbinom(1, (N/2), YI.sel)/(N/2)
						qY.wt.t     <-  rbinom((n-r), (round((N/2)*(1 - YI.t))), qY.wt.sel)/(round((N/2)*(1 - YI.t)))
						qY.del.t    <-  rbinom(r, (round((N/2)*(1 - YI.t))), qY.del.sel)/(round((N/2)*(1 - YI.t)))
						XaSp.wt.t   <-  rbinom((n-r), (N/2), XaSp.wt.sel)/(N/2)
						XaSp.del.t  <-  rbinom(r, (N/2), XaSp.del.sel)/(N/2)

#points(YI.t ~ t)
#points(qI.wt.t ~ t, col=2)
#t=t+1
					}
					fix  <-  fix + YI.t
				}
			PrFix  <-  c(PrFix, (fix/sims))
			cat('\r', paste("N: ", i, "/", length(N.vals), 
										", U: ", j, "/", length(Us.factor.vals), 
										", x: ", round(100*(k/length(invSize))), "% complete", sep=""))
			}

		}

	}

	# Index variables
	Ns        <-  rep(N.vals, each=(length(Us.factor.vals)*length(invSize)))
	Us        <-  rep((Us.factor.vals*s), each=(length(invSize)), times=length(N.vals))
	Ufac      <-  rep((Us.factor.vals), each=(length(invSize)), times=length(N.vals))
	invSizes  <-  rep(invSize, times=length(N.vals)*length(Us.factor.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./data/PrFixFigDetqI_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"Ufac"   =  Ufac,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}



########################################
##  AUTOSOMAL RECURSIONS for comparison

w.II  <-  function(n, nd, u, sdHom, h, Ud, x, t) {
	((1 - sdHom)^nd)*(1 - 2*u*(1 - exp(-sdHom*h*t)))^(n - nd)
}
w.IS  <-  function(n, nd, u, sdHom, h, Ud, x, t) {
	((u/(sdHom*h))*(1 - sdHom) + (1 - (u/(sdHom*h)))*(1 - sdHom*h))^nd * (1 - u*(2 - exp(-sdHom*h*t)))^(n - nd)
}
w.SS  <-  function(n, nd, u, sdHom, h, Ud, x, t) {
	(1 - 2*u)^n
}

makeDataAutoPrFixInvSize  <-  function(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
																	 nTot = 10^4, N.vals = c(10^3, 10^4)) {

	# Containers
	PrFix  <-  c()
	q.t    <-  c()
	
	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N     <-  N.vals[i]
		q.0   <-  1/(2*N)

		#number of simulations 
		sims  = 500*N

		# Loop over Us factor values
		for(j in 1:length(Us.factor.vals)) {
			U     <-  Us.factor.vals[j]*s
			u     <-  U/nTot
			qHat  <-  (u/(h*s))

				# Loop over inversion size
				for(k in 1:length(invSize)) {
					fix   = 0

					# Loop over replicate simulations
					for(l in 1:sims){

						## initial frequencies
						# Draw random value for # loci captured by inversion
						n   <-  rpois(1,lambda=nTot*invSize[k])

						# Draw random value for # del. mutations captured by inversion (r) given x
						r  <-  sum(rbinom(n, size=1, prob=qHat))

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
							drift  <-  rmultinom(1, N, c(SS.sel, IS.sel, II.sel))/N          
          		p.t    <-  drift[1] + drift[2]/2
          		q.t    <-  drift[2]/2 + drift[3]

							# time counter
							t  <-  t + 1
						}
						if(q.t == 1){
						fix  <-  fix + 1
						}
					}
				PrFix  <-  c(PrFix, (fix/sims))
				cat('\r', paste("N: ", i, "/", length(N.vals), 
											", U: ", j, "/", length(Us.factor.vals), 
											", x: ", round(100*(k/length(invSize))), "% complete", sep=""))
				}

		}

	}

	# Index variables
	Ns        <-  rep(N.vals, each=(length(Us.factor.vals)*length(invSize)))
	Us        <-  rep((Us.factor.vals*s), each=(length(invSize)), times=length(N.vals))
	Ufac      <-  rep((Us.factor.vals), each=(length(invSize)), times=length(N.vals))
	invSizes  <-  rep(invSize, times=length(N.vals)*length(Us.factor.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./data/PrFixAutoFig_h", h, "_s", s, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"Ufac"   =  Ufac,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}





rollAve  <-  function(x, winSize) {

	rollMean  <-  c()
	for(i in 1:(length(x) - winSize)) {
		sub  <-  x[i:(i-1+winSize)]
		rollMean[i]  <-  mean(sub)
	}
	midpoints  <-  seq_along(rollMean) + (winSize/2)
	cbind(rollMean, midpoints)
}





# Function to estimate Pr(fix | x) for different inversion sizes (x)
# Looping over different Population size and U/s ratios
makeDataTimeBen_Fix  <-  function(h = 0.1, s = 0.01, Us.factor = 2, x = 0.2,
																	 nTot = 10^4, N.vals = c(10^3, 10^4), reps=30, Nfname = "_") {

	# Containers
	tBen  <-  c()
	tFix  <-  c()

	# inversion size
	invSize  <-  x

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N  <-  N.vals[i]
		YI.0   <-  2/N

		#number of simulations 
#		sims  = 100*N/2

		# Loop over Us factor values
			U     <-  Us.factor*s
			u     <-  U/nTot
			qHat  <-  (U/(nTot*h*s))


				fix   = 0
				tBen.N  <-  c()
				tFix.N  <-  c()
				
				# Loop over replicate simulations
				while(length(tBen.N) < reps | length(tFix.N) < reps){

					## initial frequencies
					# Draw random value for # loci captured by inversion
					n   <-  rpois(1,lambda=nTot*invSize)

					# Draw random initial frequencies for del. alleles at each locus
					qi  <-  rbinom(n, size = 2*N, prob=qHat)/(2*N)

					# Draw random value for # del. mutations captured by inversion (r) 
					# conditioned on the inversion being initially beneficial
					r  <-  n*qHat + 10
					while(r > n*qHat) {
						ri  <-  rbinom(n, size=1, prob=qi)
						r   <-  sum(ri)
					}
					# Assign frequencies
					# Note implicit assumption of equal initial
					# frequencies in XOv, XSp, Y chromosomes
					YI.t        <- YI.0
					XaOv.wt.t   <- qi[ri == 0]
					XaSp.wt.t   <- qi[ri == 0]
					qI.wt.t     <- 0
					qY.wt.t     <- qi[ri == 0]
					XaOv.del.t  <- qi[ri == 1]
					XaSp.del.t  <- qi[ri == 1]
					qY.del.t    <- qi[ri == 1]

					# Run forward simulation
					t  <-  1
					wbarYI.t  <-  c()
					while(YI.t*(1 - YI.t) > 0) {

						#expected frequency after selection
						XaOv.wt.sel   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
						XaOv.del.sel  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
						XaSp.wt.sel   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						XaSp.del.sel  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qY.wt.sel     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						qY.del.sel    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qI.wt.sel     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
						YI.sel        <-  round(YI.multi.prime(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
						wbarYI.t[t]   <-  (YI.sel/ YI.t)
						if(t == 1 && wbarYI.t[1] < 1) {
							break
						}

						#binomial sampling
						XaOv.wt.t   <-  rbinom((n-r), N, prob=XaOv.wt.sel)/N
						XaOv.del.t  <-  rbinom(r, N, XaOv.del.sel)/N
						YI.t        <-  rbinom(1, (N/2), YI.sel)/(N/2)
						qY.wt.t     <-  rbinom((n-r), (round((N/2)*(1 - YI.t))), qY.wt.sel)/(round((N/2)*(1 - YI.t)))
						qY.del.t    <-  rbinom(r, (round((N/2)*(1 - YI.t))), qY.del.sel)/(round((N/2)*(1 - YI.t)))
						qI.wt.t     <-  rbinom((n-r), (round((N/2)*YI.t)), qI.wt.sel)/(round((N/2)*YI.t))
						XaSp.wt.t   <-  rbinom((n-r), (N/2), XaSp.wt.sel)/(N/2)
						XaSp.del.t  <-  rbinom(r, (N/2), XaSp.del.sel)/(N/2)

						t  <- t + 1
					}
					if(t > 200) {
						rM  <-  rollAve(x=wbarYI.t, winSize=50)
						if(any(rM[,1] < 1)) {
						plot(wbarYI.t)
						lines(rM[,1] ~ rM[,2], lwd=2, col=4)
						abline(h=1)
						tBen.N  <-  c(tBen.N, rM[,2][rM[,1] < 1][1])
						}
					}
					fix  <-  fix + YI.t
					if(YI.t == 1) {
						tFix.N  <-  c(tFix.N, t-1)
					}
					cat('\r', paste("N: ", i, "/", length(N.vals), ", ",
											"tBen: ", length(tBen.N), "/", reps, ", ", 
										  "tFix: ", length(tFix.N), "/", reps, sep=""))
				}
				tBen[i]  <-  mean(tBen.N)
				tFix[i]  <-  mean(tFix.N)
			}
	
	# Export Results Dataframe
	filename  <-  paste("./data/timeBen_Fix_h", h, "_s", s, "_U", U, "_x", x, Nfname, ".csv", sep="")
#	filename  <-  paste("./data/PrFixFig_h", h, "_s", s, "_5k", ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(N.vals)),
										"s"      =  rep(s, times=length(N.vals)),
										"U"      =  rep(U, times=length(N.vals)),
										"x"      =  rep(x, times=length(N.vals)),
										"N"      =  N.vals,
										"tBen"  =  tBen,
										"tFix"  =  tFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}








####################################
####################################
## Hybrid simulations
####################################
####################################


# wBarY.multi  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del) {
# 	pt.I.wt    <-  1 - qt.I.wt
# 	pt.Y.del   <-  1 - qt.Y.del
# 	pt.Y.wt    <-  1 - qt.Y.wt
# 		  YI.t*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )) * prod(1 - s*(h*(pt.I.wt*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt)) +
# 	(1 - YI.t)*( prod(1 - s*(h*((1 - XaOv.t.del)*qt.Y.del + XaOv.t.del*pt.Y.del) + XaOv.t.del*qt.Y.del)) * prod((1 - s*(h*((1 - XaOv.t.wt)*qt.Y.wt + XaOv.t.wt*pt.Y.wt) + XaOv.t.wt*qt.Y.wt))) )
# }

#YI.multi.prime  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del, ...) {
#	(YI.t*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del))) * prod(1 - s*(h*((1 - qt.I.wt)*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))) / wBarY.multi(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qt.I.wt, qt.Y.wt=qt.Y.wt, qt.Y.del=qt.Y.del, XaOv.t.wt=XaOv.t.wt, XaOv.t.del=XaOv.t.del)
#}


mutateReplace  <-  function(x, p=u, ...) {
	n.wt.alleles  <-  length(x[x == 0])
	x[x == 0]     <-  rbinom(n.wt.alleles, size=1, prob=p)
	x
}


indYI.mating.mutation  <-  function(N, n, r, u, YI.Mat.t, YI.sel, ...) {
	# Realized number of inverted chromosomes in next generation
	nYI.next     <-  rbinom(1, (N/2), YI.sel)
	if(nYI.next == 0) {
		list(
				 "YI.t"         =  0,
				 "YI.Mat.next"  =  NA,
				 "qI.wt.next"   =  NA
				 )
	}
	else{
		# Empty matrix of inverted chromosomes
		YI.Mat.next  <-  matrix(0, nrow=nYI.next, ncol=(n-r))
		# random gamete sampling
		YI.Mat.next[1:nYI.next,]  <-  YI.Mat.t[sample(c(1:nrow(YI.Mat.t)), size=nYI.next, replace=TRUE),]
		# mutation
		YI.Mat.next  <-  apply(X=YI.Mat.next, MARGIN=2, function(x) mutateReplace(x, p=u))
		# calculate del. allele frequencies
		if(nYI.next == 1) {
			qI.wt.next   <-  YI.Mat.next
			YI.Mat.next  <-  matrix(YI.Mat.next, nrow=1)
		} else {
			qI.wt.next   <-  colSums(YI.Mat.next)/nYI.next
		}
		# return list of results
		list(
				 "YI.t"         =  nYI.next/(N/2),
				 "YI.Mat.next"  =  YI.Mat.next,
				 "qI.wt.next"   =  qI.wt.next
				 )
	}
}


# N=1000
# nTot = 10000
# n=1000
# U=0.2
# u=U/nTot
# r=1
# YI.t  <-  matrix(0, nrow=1, ncol=(n-r))
# YI.sel=0.004
# test  <-  indYI.mating.mutation(N=N, n=n, r=r, u=u, YI.Mat.t=YI.t, YI.sel=YI.sel) 
# nrow(test$YI.Mat.next)
# any(test$qI.wt.next > 0)
# test$YI.t

# Function to estimate Pr(fix | x) for different inversion sizes (x)
# Looping over different Population size and U/s ratios
makeDataPrFixInvSizeHybridSim  <-  function(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
																	 nTot = 10^4, N.vals = c(10^3, 10^4), Nfname = "") {

	# Containers
	PrFix       <-  c()
	YI.t        <-  c()
	XaOv.wt.t   <-  c()
	XaSp.wt.t   <-  c()
	qI.wt.t     <-  c()
	qY.wt.t     <-  c()
	XaOv.del.t  <-  c()
	XaSp.del.t  <-  c()
	qY.del.t    <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N  <-  N.vals[i]
		YI.0   <-  2/N

		#number of simulations 
		sims  = 100*N/2

		# Loop over Us factor values
		for(j in 1:length(Us.factor.vals)) {
			U     <-  Us.factor.vals[j]*s
			u     <-  U/nTot
			qHat  <-  (U/(nTot*h*s))

			# Loop over inversion size
			for(k in 1:length(invSize)) {
				fix   = 0

				# Loop over replicate simulations
				for(l in 1:sims){

					## initial frequencies
					# Draw random value for # loci captured by inversion
					n   <-  rpois(1,lambda=nTot*invSize[k])

					# Draw random initial frequencies for del. alleles at each locus
					qi  <-  rbinom(n, size = 2*N, prob=qHat)/(2*N)

					# Draw random value for # del. mutations captured by inversion (r) given x
					ri  <-  rbinom(n, size=1, prob=qi)
					r   <-  sum(ri)

					# Assign frequencies
					# Note implicit assumption of equal initial
					# frequencies in XOv, XSp, Y chromosomes
					YI.t        <- YI.0
					YI.Mat.t    <- matrix(0, nrow=1, ncol=(n-r))
					XaOv.wt.t   <- qi[ri == 0]
					XaSp.wt.t   <- qi[ri == 0]
					qI.wt.t     <- 0
					qY.wt.t     <- qi[ri == 0]
					XaOv.del.t  <- qi[ri == 1]
					XaSp.del.t  <- qi[ri == 1]
					qY.del.t    <- qi[ri == 1]

					# Run forward simulation
					while(YI.t*(1 - YI.t) > 0) {
						#expected frequency after selection
						XaOv.wt.sel   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
						XaOv.del.sel  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
						XaSp.wt.sel   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						XaSp.del.sel  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qY.wt.sel     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						qY.del.sel    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qI.wt.sel     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
						YI.sel        <-  round(YI.multi.prime(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)

						# ind. simulation of inverted Y's
						simYI  <-  indYI.mating.mutation(N=N, n=n, r=r, u=u, YI.Mat.t=YI.Mat.t, YI.sel=YI.sel) 
						YI.t        <-  simYI$YI.t
						qI.wt.t     <-  simYI$qI.wt.next
						YI.Mat.t    <-  simYI$YI.Mat.next
						# Binomial sampling of other loci x chromosome classes
						XaOv.wt.t   <-  rbinom((n-r), N, prob=XaOv.wt.sel)/N
						XaOv.del.t  <-  rbinom(r, N, XaOv.del.sel)/N
						qY.wt.t     <-  rbinom((n-r), (round((N/2)*(1 - YI.t))), qY.wt.sel)/(round((N/2)*(1 - YI.t)))
						qY.del.t    <-  rbinom(r, (round((N/2)*(1 - YI.t))), qY.del.sel)/(round((N/2)*(1 - YI.t)))
						XaSp.wt.t   <-  rbinom((n-r), (N/2), XaSp.wt.sel)/(N/2)
						XaSp.del.t  <-  rbinom(r, (N/2), XaSp.del.sel)/(N/2)
					}
					fix  <-  fix + YI.t
				}
			PrFix  <-  c(PrFix, (fix/sims))
			cat('\r', paste("N: ", i, "/", length(N.vals), 
										", U: ", j, "/", length(Us.factor.vals), 
										", x: ", round(100*(k/length(invSize))), "% complete", sep=""))
			}

		}

	}

	# Index variables
	Ns        <-  rep(N.vals, each=(length(Us.factor.vals)*length(invSize)))
	Us        <-  rep((Us.factor.vals*s), each=(length(invSize)), times=length(N.vals))
	Ufac      <-  rep((Us.factor.vals), each=(length(invSize)), times=length(N.vals))
	invSizes  <-  rep(invSize, times=length(N.vals)*length(Us.factor.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./data/PrFixHybridSimFig_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"Ufac"   =  Ufac,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}





###########################################
#

#YI.Mat.t    <- matrix(0, nrow=2, ncol=(n-r))

mutateReplace  <-  function(x, p=u, ...) {
	n.wt.alleles  <-  length(x[x == 0])
	x[x == 0]     <-  rbinom(n.wt.alleles, size=1, prob=p)
	x
}

#YI.Mat.mut  <-  mutateReplace(x=YI.Mat.t, p=u)

indYI.mating.mutation  <-  function(N, n, r, u, YI.Mat.t, YI.sel, ...) {
	# Realized number of inverted chromosomes in next generation
	nYI.next     <-  rbinom(1, (N/2), YI.sel)
	if(nYI.next == 0) {
		list(
				 "YI.t"         =  0,
				 "YI.Mat.next"  =  NA,
				 "qI.wt.next"   =  NA
				 )
	}
	else{
		# Empty matrix of inverted chromosomes
		YI.Mat.next  <-  matrix(0, nrow=nYI.next, ncol=(n-r))
		# random gamete sampling
		YI.Mat.next[1:nYI.next,]  <-  YI.Mat.t[sample(c(1:nrow(YI.Mat.t)), size=nYI.next, replace=TRUE),]
		# mutation
		YI.Mat.next  <-  apply(X=YI.Mat.next, MARGIN=2, function(x) mutateReplace(x, p=u))
		# calculate del. allele frequencies
		if(nYI.next == 1) {
			qI.wt.next   <-  YI.Mat.next
			YI.Mat.next  <-  matrix(YI.Mat.next, nrow=1)
		} else {
			qI.wt.next   <-  colSums(YI.Mat.next)/nYI.next
		}
		# return list of results
		list(
				 "YI.t"         =  nYI.next/(N/2),
				 "YI.Mat.next"  =  YI.Mat.next,
				 "qI.wt.next"   =  qI.wt.next
				 )
	}
}

#simYI  <-  indYI.mating.mutation(N=N, n=n, r=r, u=u, YI.Mat.t=YI.Mat.t, YI.sel=YI.sel) 





# Function to estimate Pr(fix | x) for different inversion sizes (x)
# Looping over different Population size and U/s ratios
makeTESTDataPrFixInvSize  <-  function(h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
																	 nTot = 10^4, N.vals = c(10^3, 10^4), Nfname = "") {

	# Containers
	PrFix       <-  c()
	YI.t        <-  c()
	XaOv.wt.t   <-  c()
	XaSp.wt.t   <-  c()
	qI.wt.t     <-  c()
	qY.wt.t     <-  c()
	XaOv.del.t  <-  c()
	XaSp.del.t  <-  c()
	qY.del.t    <-  c()

	# fitness expresions
	W  <-  c(1, 1 - (h*s), 1 - s)

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N  <-  N.vals[i]
		YI.0   <-  2/N

		#number of simulations 
		sims  = 100*N/2

		# Loop over Us factor values
		for(j in 1:length(Us.factor.vals)) {
			U     <-  Us.factor.vals[j]*s
			u     <-  U/nTot
			qHat  <-  (U/(nTot*h*s))
	
			# Loop over inversion size
			for(k in 1:length(invSize)) {
				# Draw random value for # loci captured by inversion
				n   <-  nTot*invSize[k]

				# counter for fixations
				fix   = 0

				# Loop over replicate simulations
				for(l in 1:sims){

					## initial frequencies
					# Draw random initial frequencies for del. alleles at each locus
					qi  <-  rbinom(n, size = 2*N, prob=qHat)/(2*N)

					# Draw random value for # del. mutations captured by inversion (r) given x
					ri  <-  rbinom(n, size=1, prob=qi)
					r   <-  sum(ri)

					# Assign frequencies
					# Note implicit assumption of equal initial
					# frequencies in XOv, XSp, Y chromosomes
					YI.t        <- YI.0
					XaOv.wt.t   <- qi[ri == 0]
					XaSp.wt.t   <- qi[ri == 0]
					qI.wt.t     <- 0
					qY.wt.t     <- qi[ri == 0]
					XaOv.del.t  <- qi[ri == 1]
					XaSp.del.t  <- qi[ri == 1]
					qY.del.t    <- qi[ri == 1]

plot(NA, ylim=c(0,1), xlim=c(0,N))
text(x=(N/2), y=0.9, labels=paste(l, sep=''))
points(YI.t ~ t)
points(qI.wt.t[1] ~ t, col=2)					# Run forward simulation
t=1
					while(YI.t*(1 - YI.t) > 0) {

						## Step through recursions:
						# 1) expected frequency after mutation & selection
						XaOv.wt.mut   <-  XaOv.wt.prime(u=u, s=0, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t)
						XaOv.del.mut  <-  XaOv.del.prime(u=u, s=0, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t)
						XaSp.wt.mut   <-  XaSp.wt.prime(u=u, s=0, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t)
						XaSp.del.mut  <-  XaSp.del.prime(u=u, s=0, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t)
						qY.wt.mut     <-  qY.wt.prime(u=u, s=0, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t)
						qY.del.mut    <-  qY.del.prime(u=u, s=0, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t)
						qI.wt.mut     <-  qI.wt.prime(n=n, r=r, u=u, s=0, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t)
						# 2) adult genotypic frequencies
						Fii.f.wt    <-  cbind((1 - XaOv.wt.mut)*(1 - XaSp.wt.mut), ((1 - XaOv.wt.mut)*XaSp.wt.mut + XaOv.wt.mut*(1 - XaSp.wt.mut)), XaOv.wt.mut*XaSp.wt.mut)
#						colnames(Fii.f.wt)  <-  c("AA", "Aa", "aa")
						Fii.f.del   <-  cbind((1 - XaOv.del.mut)*(1 - XaSp.del.mut), ((1 - XaOv.del.mut)*XaSp.del.mut + XaOv.del.mut*(1 - XaSp.del.mut)), XaOv.del.mut*XaSp.del.mut)
#						colnames(Fii.f.del)  <-  c("AA", "Aa", "aa")
						Fii.Y.wt    <-  cbind((1 - XaOv.wt.mut)*(1 - qY.wt.mut), ((1 - XaOv.wt.mut)*qY.wt.mut + XaOv.wt.mut*(1 - qY.wt.mut)), XaOv.wt.mut*qY.wt.mut)
#						colnames(Fii.Y.wt)  <-  c("AA", "Aa", "aa")
						Fii.Y.del   <-  cbind((1 - XaOv.del.mut)*(1 - qY.del.mut), ((1 - XaOv.del.mut)*qY.del.mut + XaOv.del.mut*(1 - qY.del.mut)), XaOv.del.mut*qY.del.mut)
#						colnames(Fii.Y.del)  <-  c("AA", "Aa", "aa")
						Fii.YI.wt   <-  cbind(cbind((1 - XaOv.wt.mut)*(1 - qI.wt.mut), ((1 - XaOv.wt.mut)*qI.wt.mut + XaOv.wt.mut*(1 - qI.wt.mut)), XaOv.wt.mut*qI.wt.mut))
#						colnames(Fii.YI.wt)  <-  c("AA", "Aa", "aa")
						Fii.YI.del  <-  cbind(cbind((1 - XaOv.del.mut)*(1 - 1), (1 - XaOv.del.mut)*1 , XaOv.del.mut*1))
#						colnames(Fii.YI.del)  <-  c("AA", "Aa", "aa")
						# 3) expected genotypic frequencies after selection
						E.YI.sel     <-  YI.multi.prime(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.mut, qt.Y.wt=qY.wt.mut, qt.Y.del=qY.del.mut, XaOv.t.wt=XaOv.wt.mut, XaOv.t.del=XaOv.del.mut)
						E.Fii.f.wt   <-  (t(apply(X=Fii.f.wt, MARGIN=1, function(x)  (c(x[1]*W[1], x[2]*W[2], x[3]*W[3])/ sum(c(x[1]*W[1], x[2]*W[2], x[3]*W[3]))) )))
						E.Fii.f.del  <-  (t(apply(X=Fii.f.del, MARGIN=1, function(x) (c(x[1]*W[1], x[2]*W[2], x[3]*W[3])/ sum(c(x[1]*W[1], x[2]*W[2], x[3]*W[3]))) )))
						E.Fii.Y.wt   <-  (t(apply(X=Fii.Y.wt, MARGIN=1, function(x)  (c(x[1]*W[1], x[2]*W[2], x[3]*W[3])/ sum(c(x[1]*W[1], x[2]*W[2], x[3]*W[3]))) )))
						E.Fii.Y.del  <-  (t(apply(X=Fii.Y.del, MARGIN=1, function(x) (c(x[1]*W[1], x[2]*W[2], x[3]*W[3])/ sum(c(x[1]*W[1], x[2]*W[2], x[3]*W[3]))) )))
						# 4) realized frequencies after multinomial sampling
						YI.t         <-  rbinom(1, (N/2), E.YI.sel)/(N/2)
						R.Fii.f.wt   <-  t(apply(X=E.Fii.f.wt, MARGIN=1, function(x) rmultinom(n=1, size=(N/2), p=c(x[1], x[2], x[3]))/(N/2)))
						R.Fii.Y.wt   <-  t(apply(X=E.Fii.Y.wt, MARGIN=1, function(x) rmultinom(n=1, size=(N/2)*(1-YI.t), p=c(x[1], x[2], x[3]))/((N/2)*(1-YI.t))))
						if(r > 0) {
							R.Fii.f.del  <-  t(apply(X=E.Fii.f.del, MARGIN=1, function(x) rmultinom(n=1, size=(N/2), p=c(x[1], x[2], x[3]))/(N/2)))
							R.Fii.Y.del  <-  t(apply(X=E.Fii.Y.del, MARGIN=1, function(x) rmultinom(n=1, size=(N/2)*(1-YI.t), p=c(x[1], x[2], x[3]))/((N/2)*(1-YI.t))))
						}
						# 5) Resulting allele frequencies among gametes
						XaOv.wt.t    <-  (R.Fii.f.wt[,3] + R.Fii.f.wt[,2])/2
						XaSp.wt.t    <-  ((R.Fii.Y.wt[,3] + R.Fii.Y.wt[,2])/2)*(1 - YI.t) + ((Fii.YI.wt[,3] + Fii.YI.wt[,2])/2)*YI.t
						qY.wt.t      <-  (R.Fii.Y.wt[,3] + R.Fii.Y.wt[,2])/2
						if(r > 0) {
							XaOv.del.t   <-  (R.Fii.f.del[,3] + R.Fii.f.del[,2])/2
							XaSp.del.t   <-  ((R.Fii.Y.del[,3] + R.Fii.Y.del[,2])/2)*(1 - YI.t) + ((Fii.YI.del[,3] + Fii.YI.del[,2])/2)*YI.t
							qY.del.t     <-  (R.Fii.Y.del[,3] + R.Fii.Y.del[,2])/2
						}
						# deterministic change in allele frequency at wt loci on inversion
						qI.wt.t      <-  qI.wt.prime(n=n, r=r, u=0, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.mut, XaOv.t.wt=XaOv.wt.mut, XaOv.t.del=XaOv.del.mut)
points(YI.t ~ t)
points(qI.wt.t[1] ~ t, col=2)
t=t + 1
					}
					fix  <-  fix + YI.t
				}
			PrFix  <-  c(PrFix, (fix/sims))
			cat('\r', paste("N: ", i, "/", length(N.vals), 
										", U: ", j, "/", length(Us.factor.vals), 
										", x: ", round(100*(k/length(invSize))), "% complete", sep=""))
			}

		}

	}

	# Index variables
	Ns        <-  rep(N.vals, each=(length(Us.factor.vals)*length(invSize)))
	Us        <-  rep((Us.factor.vals*s), each=(length(invSize)), times=length(N.vals))
	Ufac      <-  rep((Us.factor.vals), each=(length(invSize)), times=length(N.vals))
	invSizes  <-  rep(invSize, times=length(N.vals)*length(Us.factor.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./data/PrFixFigTEST_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"Ufac"   =  Ufac,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}

