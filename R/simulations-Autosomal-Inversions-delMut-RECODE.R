#' Simulation code: 
#' Autosomal Inversions w/ recessive deleterious mutations,
#' taking feeback between inversion frequency and frequency
#' of deleterious mutations on standard arrangement 
#' chromosomes
#' 
#' These functions will produce results that should be 
#' comparable to the deterministic and W-F results from
#' the Autosomal model in my earlier results.

rm(list=ls())
#' Exact Recursions.
#' This model requires a system of 4 recursions:
#' q.Prime.W   -- freq. del. mut. on Standard Arrangement loci where inversion captures WT allele
#' q.Prime.D   -- freq. del. mut. on Standard Arrangement loci where inversion captures del. mut
#' qI.Prime.W  -- freq. del. mut. at loci within Inversion where a WT allele was initially captured
#' AI.Prime    -- freq. of Inverted chromosomes.

q.Prime.W  <-  function(q, qI, AI, u, h, s) {
#	( (q*(-1 + u) - u)*(1 + h*(1 + (-1 + AI)*q - AI*qI)*s*(-1 + u) - s*(((-1 + AI)*q - AI*qI)*(-1 + u) + u))) / 
#		(-1 - s*(q*(-1 + u) - u)*(((-1 + AI)*q - AI*qI)*(-1 + u) + u) + h*s*(-1 + u)*(2*(-1 + AI)*(q^2)*(-1 + u) - 2*u + AI*qI*(-1 + 2*u) + q*(-2 + AI + 2*AI*qI - 2*(-2 + AI + AI*qI)*u)))
((1 - AI)*((1 - h*s)*(1 - q - (1 - q)*u)*(q + (1 - q)*u) + (1 - s)*(q + (1 - q)*u)^2) + AI*((1 - h*s)*(q + (1 - q)*u)*(1 - qI - (1 - qI)*u) + (1 - s)*(q + (1 - q)*u)*(qI + (1 - qI)*u))) / 
	((1 - AI)*((1 - q - (1 - q)*u)^2 + 2*(1 - h*s)*(1 - q - (1 - q)*u)*(q + (1 - q)*u) + (1 - s)*(q + (1 - q)*u)^2) + AI*((1 - q - (1 - q)*u)*(1 - qI - (1 - qI)*u) + (1 - s)*(q + (1 - q)*u)*(qI + (1 - qI)*u) + (1 - h*s)*((q + (1 - q)*u)*(1 - qI - (1 - qI)*u) + (1 - q - (1 - q)*u)*(qI + (1 - qI)*u))))
}

q.Prime.D  <-  function(qD, AI, u, h, s) {
#	((qD*(-1 + u) - u)*(1 + AI*(-1 + h)*(-1 + qD)*s*(-1 + u) - s*(qD + h*(-1 + qD)*(-1 + u) + u - qD*u))) / 
#    	(-1 + AI*(-1 + qD)*s*(-1 + u)*(h + qD - 2*h*qD + (-1 + 2*h)*(-1 + qD)*u) - s*(qD*(-1 + u) - u)*(qD + 2*h*(-1 + qD)*(-1 + u) + u - qD*u))
(AI*(1 - s)*(qD + (1 - qD)*u) + (1 - AI)*((1 - h*s)*(1 - qD - (1 - qD)*u)*(qD + (1 - qD)*u) + (1 - s)*(qD + (1 - qD)*u)^2)) / 
	(AI*((1 - h*s)*(1 - qD - (1 - qD)*u) + (1 - s)*(qD + (1 - qD)*u)) + (1 - AI)*((1 - qD - (1 - qD)*u)^2 + 2*(1 - h*s)*(1 - qD - (1 - qD)*u)*(qD + (1 - qD)*u) + (1 - s)*(qD + (1 - qD)*u)^2))
}

q.Prime.I.W  <-  function(q, qI, AI, u, h, s) {
#	-(((qI*(-1 + u) - u)*(-1 + s*(((-1 + AI)*q - AI*qI)*(-1 + u) - h*(1 + (-1 + AI)*q - AI*qI)*(-1 + u) + u))) / 
#			(-1 - s*(qI*(-1 + u) - u)*(((-1 + AI)*q - AI*qI)*(-1 + u) + u) + h*s*(-1 + u)*((-1 + AI)*q*(1 + 2*qI*(-1 + u) - 2*u) - 2*AI*(qI^2)*(-1 + u) - 2*u + (1 + AI)*qI*(-1 + 2*u))))
((1 - AI)*((1 - h*s)*(1 - q - (1 - q)*u)*(qI + (1 - qI)*u) + (1 - s)*(q + (1 - q)*u)*(qI + (1 - qI)*u)) + AI*((1 - h*s)*(1 - qI - (1 - qI)*u)*(qI + (1 - qI)*u) + (1 - s)*(qI + (1 - qI)*u)^2)) / 
	(AI*((1 - qI - (1 - qI)*u)^2 + 2*(1 - h*s)*(1 - qI - (1 - qI)*u)*(qI + (1 - qI)*u) + (1 - s)*(qI + (1 - qI)*u)^2) + (1 - AI)*((1 - q - (1 - q)*u)*(1 - qI - (1 - qI)*u) + (1 - s)*(q + (1 - q)*u)*(qI + (1 - qI)*u) + (1 - h*s)*((q + (1 - q)*u)*(1 - qI - (1 - qI)*u) + (1 - q - (1 - q)*u)*(qI + (1 - qI)*u))))
}

WBar  <-  function(q, qD, qI, AI, u, h, s, n, r) {
#	A   <-  1 - AI
#	p   <-  1 - q
#	pD  <-  1 - qD
#	pI  <-  1 - qI
	((AI^2)*(1 - s) + 2*A*AI*(qD*(1 - s) + pD*(1 - h*s)) + (A^2)*(pD^2 + (qD^2)*(1 - s) + 2*pD*qD*(1 - h*s)))^r *
		((A^2)*(p^2 + (q^2)*(1 - s) + 2*p*q*(1 - h*s)) + (AI^2)*(pI^2 + (qI^2)*(1 - s) + 2*pI*qI*(1 - h*s)) + 2*A*AI*(p*pI + q*qI*(1 - s) + (q*pI + p*qI)*(1 - h*s)))^(n - r)
}

AI.Prime  <-  function(q, qD, qI, AI, u, h, s, n, r, ...) {

	A   <-  1 - AI
	p   <-  1 - q
	pD  <-  1 - qD
	pI  <-  1 - qI
	A   <-  1 - AI
#	q   <-  q/A
#	qD  <-  qD/A
#	qI  <-  qI/AI
#	p   <-  (1 - AI - q)
#	pD  <-  (1 - AI - qD)/A
#	pI  <-  (AI - qI)

	W.SS   <-  (1 - s*(2*h*qD*pD + qD^2))^r * (1 - s*(2*h*q*p + q^2))^(n-r)
	W.IS   <-  (1 - s*(2*h*pD + qD))^r * (1 - s*(h*(q*pI + p*qI) + q*qI))^(n-r)
	W.II   <-  (1 - s)^r * (1 - s*(2*h*qI*pI + qI^2))^(n-r)
	W.avg  <-  (A^2)*W.SS + 2*A*AI*W.IS + (AI^2)*W.II
	
	SS.sel  <-  (A^2)*W.SS/W.avg
	IS.sel  <-  2*A*AI*W.IS/W.avg
	II.sel  <-  (AI^2)*W.II/W.avg

	freqs   <-  c(SS.sel, IS.sel, II.sel)
	AI.pr   <-  freqs[2]/2 + freqs[3]
	return(AI.pr)
}




makeAutoDeterministicFigSimData  <-  function(U = 0.02, r = 0, x = 0.2, h = 0.25, s = 0.01, generations = 5*10^3, ...) {

#browser()
	# Parameters
	nTot  <-  10^4
	u     <-  U/nTot
	qHat  <-  0.001987262#(U/(nTot*h*s))
	n     <-  nTot*x
	N     <-  10000
	AI.0  <-  1/(2*N)

	# Empty Frequency Vectors
	AI.t       <- c()
	qt.W       <- c()
	qt.D       <- c()
	qt.I.W     <- c()

	# First generation (all loci at equilibrium frequencies when inversion arises)
	qt.W[1]       <-  round(q.Prime.W(q=qHat, qI=0, AI=AI.0, u=u, h=h, s=s), digits=9)
	qt.D[1]       <-  round(q.Prime.D(qD=qHat, AI=AI.0, u=u, h=h, s=s), digits=9)
	qt.I.W[1]     <-  round(q.Prime.I.W(q=qHat, qI=0, AI=AI.0, u=u, h=h, s=s), digits=9)
	AI.t[1]       <-  round(AI.Prime(q=qHat, qD=qHat, qI=0, AI=AI.0, u=u, h=h, s=s, n=n, r=r), digits=9)


	# Subsequent generations
	i=2
	while(i < generations+1) {
		qt.W[i]       <-  round(q.Prime.W(q=qt.W[i-1], qI=qt.I.W[i-1], AI=AI.t[i-1], u=u, h=h, s=s), digits=9)
		qt.D[i]       <-  round(q.Prime.D(qD=qt.D[i-1], AI=AI.t[i-1], u=u, h=h, s=s), digits=9)
		qt.I.W[i]     <-  round(q.Prime.I.W(q=qt.W[i-1], qI=qt.I.W[i-1], AI=AI.t[i-1], u=u, h=h, s=s), digits=9)
		AI.t[i]       <-  round(AI.Prime(q=qt.W[i-1], qD=qt.D[i-1], qI=qt.I.W[i-1], AI=AI.t[i-1], u=u, h=h, s=s, n=n, r=r), digits=9)
		i  <-  i + 1
	}
	# Return results as df
	results  <-  data.frame(
							"qt.W"       =  qt.W,
							"qt.D"       =  qt.D,
							"qt.I.W"     =  qt.I.W,
							"AI.t"       =  AI.t
		 					)
	return(results)
}


qstar  <-  function(q.0, U, x, s, h) {
     (q.0*exp((U*x)/(1 - exp(-s*h)))) / (1 - q.0*(1 - exp((U*x)/(1 - exp(-h*s)))))
 }

I.genSol  <-  function(q.0, U, x, s, h, t) {
     (exp(((1 - exp(-h*s*t))*U*x) / (1 - exp(-h*s)))*q.0)/(1 - q.0 + exp(((1 - exp(-h*s*t))*U*x) / (1 - exp(-h*s)))*q.0)
 }


test  <-  makeAutoDeterministicFigSimData(U = 0.02, r = 0, x = 0.2, h = 0.25, s = 0.01, generations = 2*10^4)
N=10000
U=0.02
h=0.25
qstr  <-  qstar(q.0=(1/(2*N)), U=U, x=x, s=s, h=h)
qHat  <-  (U/nTot)/(s*h)
qIgensol  <-  I.genSol(q.0=(1/(2*N)), U=U, x=x, s=s, h=h, t=c(1:10000)) 
par(mfrow=c(2,2))
plot(AI.t ~ seq_along(AI.t), type='l', ylim=c(0,qstr), data=test)
lines(qIgensol ~ seq_along(AI.t), type='l', col=2, lty=2, data=test)
abline(h=qstr, lty=2)
plot(qt.W ~ seq_along(qt.W), type='l', data=test)
abline(h=qHat, lty=2)
plot(qt.D ~ seq_along(qt.D), type='l', data=test)
abline(h=qHat, lty=2)
plot(qt.I.W ~ seq_along(qt.I.W), type='l', data=test)
abline(h=qHat, lty=2)







