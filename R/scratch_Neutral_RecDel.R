# Approx. Deterministic trajectory 

YInext  <-  function(YIt, U, s, n, r, t, x) {
	qHat  <-  sqrt(U/(n*s))
	qt    <-  qHat*(1 - exp(-(sqrt(s*U)/n)*t))
	YIt*exp(-(n*x - r)*s*qt*qHat - r*s*qHat) / ((1 - YIt)*exp(-U*x) + YIt*exp(-(n*x - r)*s*qt*qHat - r*s*qHat))
}

# Parameters
U    <-  0.1
s    <-  0.1
n     <-  1000
r     <-  3
qHat  <-  sqrt(U/n*s)
n*qHat
x     <-  0.2
N     <-  1000
t     <-  c(1:19999)
YI0   <-  2/N

YIt  <- c()
for(i in 1:length(t)) {
	if(i == 1) {
		YIt[i]  <-  YInext(YIt=YI0, U=U, s=s, n=n, r=r, t=t[i], x=x)
	}else{YIt[i]  <-  YInext(YIt = YIt[i-1], U=U, s=s, n=n, r=r, t=t[i], x=x) }
}
if(YIt[i] == 1) {
	YIt[i]  <-  YIt[i] - 0.00001
}
YIt  <- c(YI0,YIt)
plot(YIt ~ seq_along(YIt), ylim=c(0,1), type='l', lwd=2)


rbernoulli <- function(n,p) runif(n) < p


n   <-  2000000
mu  <-  10^-5
h   <-  0.1
qHat <-  c()
mut  <-  c()
for(i in 1:n) {
	s  <-  rnorm(n=1,mean=0.01, s=0.001)
	qHat[i]  <-  (mu/h*s)
	mut[i]   <-  rbernoulli(n=1,p=qHat[i])
}

sum(mut)
sum(qHat)
n*(mu/h*0.01)
sum(rpois(n=2000000,lambda=(mu/h*0.01)))

# Deterministic trajectory 

# Parameters
U    <-  0.1
s    <-  0.01
n     <-  1000
r     <-  10
qHat  <-  sqrt(U/n*s)
x     <-  0.1
N     <-  1000
YI0   <-  2/N

#number of simulations 
sims  = 100*N/2
fix   = 0
for(i in 1:sims){
	YI  <-  YI0
	t   <-  0
	while(YI*(1 - YI) > 0) {

		#expected frequency after selection          
		qt    <-  qHat*(1 - exp(-(sqrt(s*U)/n)*t))
        w.avg = ((1 - YI)*exp(-U*x) + YI*exp(-(n*x - r)*s*qt*qHat - r*s*qHat))
        F.sel = YI + YI*(1 - YI)*(exp((r*U - n*r*s*sqrt(U/(n*s)) - exp(-t*sqrt((s*U)/n))*U*(r - n*x))/n)-1)/w.avg



        #binomial sampling
        YI = rbinom(1, (N/2), F.sel)/(N/2)
        t = t + 1
    }
    if(YI == 1){
    	fix  <-  fix + 1
    }
    cat('\r', paste("Progress:", round(100*(i/sims)), "% complete"))
}    

fix/sims



# Parameters used in biorXiv paper
# u  <-  10^-8
# n  <-  20*10^6
# n*u = 0.2
# x  <-  0.1
#n*x*qHat
# n*x = 2*10^6




# loop over r to get a profile of Pr(fix)

# Parameters
U    <-  0.1
s    <-  0.01
n     <-  1000
r     <-  c(0:29)
qHat  <-  sqrt(U/n*s)
x     <-  0.1
N     <-  1000
YI0   <-  2/N

PrFix  <-  c()
#number of simulations 
sims  = 100*N/2

for(j in 1:length(r)) {
fix   = 0
for(i in 1:sims){
	YI  <-  YI0
	t   <-  0
	while(YI*(1 - YI) > 0) {

		#expected frequency after selection          
		qt    <-  qHat*(1 - exp(-(sqrt(s*U)/n)*t))
        w.avg = ((1 - YI)*exp(-U*x) + YI*exp(-(n*x - r[j])*s*qt*qHat - r[j]*s*qHat))
        F.sel = YI + YI*(1 - YI)*(exp((r[j]*U - n*r[j]*s*sqrt(U/(n*s)) - exp(-t*sqrt((s*U)/n))*U*(r[j] - n*x))/n)-1)/w.avg



        #binomial sampling
        YI = rbinom(1, (N/2), F.sel)/(N/2)
        t = t + 1
    }
    if(YI == 1){
    	fix  <-  fix + 1
    }
}    

PrFix[j]  <-  fix/sims
cat('\r', paste("Progress:", round(100*(j/length(r))), "% complete"))
}

par(mfrow=c(1,2))
plot(PrFix ~ r)
PrR  <-  ppois(r,lambda=n*U*x)
plot(PrFix*PrR ~ r, type='l')

NrR  <-  rpois((100*N/2),lambda=n*U*x)
range(NrR)
sum(NrR*PrFix[NrR+1])/sims




#######################

# Parameters
h    <-  0.1
s    <-  0.01
U    <-  10*s
nTot  <-  10000
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
N     <-  1000
YI.0   <-  2/N



#number of simulations 
sims  = 100*N/2

for(j in 1:length(r)) {
fix   = 0
for(i in 1:sims){
	YI  <-  YI0
	t   <-  0
	while(YI*(1 - YI) > 0) {

		#expected frequency after selection          
		qt    <-  qHat*(1 - exp(-(sqrt(s*U)/n)*t))
        w.avg = ((1 - YI)*exp(-U*x) + YI*exp(-(n*x - r[j])*s*qt*qHat - r[j]*s*qHat))
        F.sel = YI + YI*(1 - YI)*(exp((r[j]*U - n*r[j]*s*sqrt(U/(n*s)) - exp(-t*sqrt((s*U)/n))*U*(r[j] - n*x))/n)-1)/w.avg



        #binomial sampling
        YI = rbinom(1, (N/2), F.sel)/(N/2)
        t = t + 1
    }
    if(YI == 1){
    	fix  <-  fix + 1
    }
}    

PrFix[j]  <-  fix/sims
cat('\r', paste("Progress:", round(100*(j/length(r))), "% complete"))
}

par(mfrow=c(1,2))
plot(PrFix ~ r)
PrR  <-  ppois(r,lambda=n*U*x)
plot(PrFix*PrR ~ r, type='l')