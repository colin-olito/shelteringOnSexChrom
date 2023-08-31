
#' Recursions based on Nei, Kojima, Schaffer (1967).
#' Eqs(16--23)

rm(list=ls())
#' They use x_t for inversion frequency, but since I've
#' already been using x for inversion size, I'm using I.
I.Prime  <-  function(qN, qI, qD, I, h, s, n, r, ...) {
	QN  <-  qN/(1-I)
	QI  <-  qI/I
	( (I^2)*(1 - s*QI^2)^(n) +   I*(1 - I)*(1 - s*QI*QN)^(n) ) / 
	( (I^2)*(1 - s*QI^2)^(n) + 2*I*(1 - I)*(1 - s*QI*QN)^(n) + ((1 - I)^2)*(1 - s*QN^2)^(n))
}

qI.Prime.W  <-  function(qN, qI, I, Iprime, u, h, s, ...) {
	pI   <-  I - qI
	w_i  <-  1 - s*(qI + qN)^2
	( ( qI*(1 - s*(qN + qI)) + u*pI*w_i) / 
		(I - s*qI*(qN + qI)) ) * Iprime
}

pI.Prime.W  <-  function(Iprime, qIprime.W) {
	Iprime - qIprime.W
}

qN.Prime.W  <-  function(qN, qI, I, Iprime, u, h, s, ...) {
	pN   <-  1 - I - qN
	w_i  <-  1 - s*(qI + qN)^2
	( ( qN*(1 - s*(qN + qI)) + u*pN*w_i ) / 
		(1 - I - s*qN*(qN + qI)) ) * (1 - Iprime)
}

pN.Prime.W  <-  function(Iprime, qNprime.W) {
	1 - Iprime - qNprime.W
}

#' For now, can ignore what's going on with qN.D. I've removed
#' qD from the recursion I.Prime(), so it doesn't affect the 
#' dynamics. Right now, I'm just trying to reproduce the 
#' Nei et al. results. 
qN.Prime.D  <-  function(qN.D, I, Iprime, u, h, s, ...) {
	pN.D  <-  1 - I - qN.D
	w_i   <-  1 - s*(qN.D + I)^2
	( (qN.D*(1 - s*(qN.D + I)) + pN.D*u*w_i) / 
		(1 - I - qN.D*s*(qN.D + I)^2) ) * (1 - Iprime)
}

pN.Prime.D  <-  function(Iprime, qNprime.D) {
	1 - Iprime - qNprime.D
}


makeAutoDeterministicFigSimData_Nei  <-  function(U = 0.02, r = 0, x = 0.2, h = 0.25, s = 0.01, generations = 5*10^3, ...) {

#browser()
	# Parameters
	nTot  <-  10^4
	u     <-  U/nTot
#	qHat  <-  (U/(nTot*h*s)) # for h > 0
	qHat  <-  sqrt(u/s)  # for h = 0
	n     <-  nTot*x
	N     <-  10000
	I.0   <-  1/(2*N)

	# Empty Frequency Vectors
	I.t     <- c()
	qt.I.W  <- c()
	pt.I.W  <- c()
	qt.W    <- c()
	pt.W    <- c()
	qt.D    <- c()
	pt.D    <- c()

	# First generation (all loci at equilibrium frequencies when inversion arises)
	I.t[1]     <-  round(I.Prime(qN=qHat*(1-I.0), qI=0, qD=qHat*(1-I.0), I=I.0, h=h, s=s, n=n, r=r), digits=9)
	qt.I.W[1]  <-  round(qI.Prime.W(qN=qHat*I.0, qI=0, I=I.0, Iprime=I.t[1], u=u, h=h, s=s), digits=9)
	pt.I.W[1]  <-  round(pI.Prime.W(Iprime=I.t[1], qIprime.W=qt.I.W[1]), digits=9)
	qt.W[1]    <-  round(qN.Prime.W(qN=qHat*(1-I.0), qI=0, I=I.0, Iprime=I.t[1], u=u, h=h, s=s), digits=9)
	pt.W[1]    <-  round(pN.Prime.W(Iprime=I.t[1], qNprime.W=qt.W[1]), digits=9)
	qt.D[1]    <-  round(qN.Prime.D(qN.D=qHat*(1-I.0), I=I.0, Iprime=I.t[1], u=u, h=h, s=s), digits=9)
	pt.D[1]    <-  round(pN.Prime.D(Iprime=I.t[1], qNprime.D=qt.D[1]), digits=9)

	# Subsequent generations
	i=2
	while(i < generations+1) {
		I.t[i]     <-  round(I.Prime(qN=qt.W[i-1], qD=qt.D[i-1], qI=qt.I.W[i-1], I=I.t[i-1], h=h, s=s, n=n, r=r), digits=9)
		qt.I.W[i]  <-  round(qI.Prime.W(qN=qt.W[i-1], qI=qt.I.W[i-1], I=I.t[i], Iprime=I.t[i], u=u, h=h, s=s), digits=9)
		pt.I.W[i]  <-  round(pI.Prime.W(Iprime=I.t[i-1], qIprime.W=qt.I.W[i-1]), digits=9)
		qt.W[i]    <-  round(qN.Prime.W(qN=qt.W[i-1], qI=qt.I.W[i-1], I=I.t[i], Iprime=I.t[i], u=u, h=h, s=s), digits=9)
		pt.W[i]    <-  round(pN.Prime.W(Iprime=I.t[i-1], qNprime.W=qt.W[i-1]), digits=9)
		qt.D[i]    <-  round(qN.Prime.D(qN.D=qt.D[i-1], I=I.t[i], Iprime=I.t[i], u=u, h=h, s=s), digits=9)
		pt.D[i]    <-  round(pN.Prime.D(Iprime=I.t[i-1], qNprime.D=qt.D[i-1]), digits=9)

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


test  <-  makeAutoDeterministicFigSimData_Nei(U = 0.02, r = 0, x = 0.2, h = 0, s = 0.01, generations = 5*10^3)
head(test)

test$qt.W + test$pt.W + test$qt.I.W + test$pt.I.W

N=10000
U=0.02
s=0.01
nTot  <-  10^4
u     <-  U/nTot
#qHat  <-  (U/(nTot*h*s)) 
qHat  <-  sqrt(u/s) 
#qstr  <-  qstar(q.0=(1/(2*N)), U=U, x=x, s=s, h=h)
#qHat  <-  (U/nTot)/(s*h)
#qIgensol  <-  I.genSol(q.0=(1/(2*N)), U=U, x=x, s=s, h=h, t=c(1:10000)) 
par(mfrow=c(2,2))
plot(I.t ~ seq_along(I.t), type='l', ylim=c(0,1), data=test)
lines(qIgensol ~ seq_along(qIgensol), type='l', col=2, lty=2)
abline(h=qstr, lty=2)
plot(qt.W/(1 - I.t) ~ seq_along(qt.W), type='l', data=test)
abline(h=qHat, lty=2)
plot(qt.D/(1 - I.t) ~ seq_along(qt.D), type='l', data=test)
abline(h=qHat, lty=2)
plot(qt.I.W/I.t ~ seq_along(qt.I.W), type='l', data=test)
abline(h=qHat, lty=2)

#plot(test$qt.W * (1 - test$I.t), log='y')
#plot(test$qt.I.W * test$I.t, log='y') 