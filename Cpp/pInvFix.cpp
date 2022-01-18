#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

// example compiling commands
// g++ -std=c++11 trajectory.cpp `pkg-config --libs gsl` command for compiling
// clang++ -std=c++11 -stdlib=libc++ example.cpp -o example_program

// random number generator
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

// Rfunction
//	mutateReplace  <-  function(x, p, ...) {
//		x[x == 0]     <-  rbinom(n=length(x[x == 0]), size=1, prob=p)
//		x
//	}
int mutation (int x, double p) 
{
	for(int i=0; i != sizeof( x ); i++) {
		if(x[i]==0){
			x[i] = unsigned int gsl_ran_binomial(const gsl_rng *r, double p, unsigned int 1);
		}
	}
	return x;
}

// Rfunction
//	relFitness  <-  function(mat1, mat2, W) {
//		genotype  <-  mat1 + mat2 + 1
//		W.i       <-  as.vector(apply(X=genotype, MARGIN=1, function(x){ prod(W[x])} ))
//		W.i/mean(W.i)
//	}
double relFitness(int mat1, int mat2, double W) 
{
	int genotype = mat1 + mat2 + 1;
	double W_i[sizeof(mat1)];
	for (int i=0; i != sizeof( W_i ); i++) 
	{
		auto indFit = accumulate(begin(genotype[i]), end(genotype[i]), 1, multiplies<>());
		W_i[i] = indFit;
	}
	W_i = W_i/gsl_stats_mean(W_i, 1, sizeof(W_i));
	return W_i
}

// Rfunction
//	freeRecombine.sampOvules.FL  <-  function(m1, m2, n, moms) {
//		f1.mat  <-  matrix(NA,nrow=nrow(m2), ncol=ncol(m2))
//		for(i in 1:nrow(m2)) {
//			shuffle  <-  rbinom(n=n, size=1, prob=1/2)
//			c1       <-  rep(0,times=n)
//			c2       <-  rep(0,times=n)
//			c1[shuffle==0]  <-  m1[moms[i],][shuffle==0]
//			c1[shuffle==1]  <-  m2[moms[i],][shuffle==1]
//			c2[shuffle==0]  <-  m2[moms[i],][shuffle==0]
//			c2[shuffle==1]  <-  m1[moms[i],][shuffle==1]
//			f1.mat[i,]      <-  get(ifelse(runif(1) > 1/2, "c2", "c1"))
//		}
//		return(f1.mat)
//	}
int freeRecombine_sampOvules(int m1, int m2, int n, int moms)
{
	// Random distributions
	double gsl_ran_flat(const gsl_rng * r, double a, double b);      // a = lower bound, b= upper bound
	unsigned int gsl_ran_binomial(const gsl_rng *r, double p, unsigned int k) // p = prob. success, k = n_trials

	int n_ind = sizeof(m2)
	int f1_mat[int n_ind][int n];
	for (int i=0; i != n; i++) 
	{
		int shuffle[n] =  int gsl_ran_binomial(r, 0.5, n)
		int std::vector<int> c1(n, 0);
		int std::vector<int> c2(n, 0);
		for (auto j = c1.begin(); j <= c1.end(); ++j) 
		{
			switch (shuffle[j]) 
			{								// if (shuffle[j] == 0) {
				case 0:						// c1[j] = m1[moms[i]][j];
					c1[j] = m1[moms[i]][j];	// c2[j] = m2[moms[i]][j];
					c2[j] = m2[moms[i]][j];	// }
				case 1:						// else {
					c1[j] = m2[moms[i]][j];	// c1[j] = m2[moms[i]][j];
					c2[j] = m1[moms[i]][j];	// c2[j] = m1[moms[i]][j];
			}								// }
		}
		double randGamete = gsl_ran_flat(r, 0.0, 1.0);
		if (randGamete < 0.5) 
		{
			f1_mat[i] = c1;
		}
		else 
		{
			f1_mat[i] = c2;
		}
	}
}




// Rfunction
//	freeRecombine.sampSperm.FL  <-  function(m1, m2, n, dads, inverted, sons) {	
//		f1.mat  <-  matrix(NA,nrow=nrow(m2), ncol=ncol(m2))
//		for(i in 1:nrow(m2)) {
//			if(inverted[i]) {
//				f1.mat[i,]  <-  get(ifelse(sons, "m2", "m1"))[dads[i],]
//			}
//			if(!inverted[i]) {
//				shuffle  <-  rbinom(n=n, size=1, prob=1/2)
//				c1       <-  rep(0,times=n)
//				c2       <-  rep(0,times=n)
//				c1[shuffle==0]  <-  m1[dads[i],][shuffle==0]
//				c1[shuffle==1]  <-  m2[dads[i],][shuffle==1]
//				c2[shuffle==0]  <-  m2[dads[i],][shuffle==0]
//				c2[shuffle==1]  <-  m1[dads[i],][shuffle==1]
//				f1.mat[i,]  <-  get(ifelse(sons, "c2", "c1"))
//			}
//		}
//		return(f1.mat)
//	}
int freeRecombine_sampSperm(int m1, int m2, int n, int dads, int inverted, int sons)
{
	// Random distributions
	unsigned int gsl_ran_binomial(const gsl_rng *r, double p, unsigned int k) // p = prob. success, k = n_trials

	int n_ind = sizeof(m2)
	int f1_mat[int n_ind][int n];
	for (int i=0; i != n; i++) 
	{
		if (inverted[i]) 
		{
			if (sons[i]==1) 
			{
				f1_mat[i] = m2[i];
			}
			else 
			{
				f1_mat[i] = m1[i];
			}
		}
		else
		{
			int shuffle[n] =  int gsl_ran_binomial(r, 0.5, n)
			int std::vector<int> c1(n, 0);
			int std::vector<int> c2(n, 0);
			for (auto j = c1.begin(); j <= c1.end(); ++j) 
			{
				switch (shuffle[j]) 
				{								// if (shuffle[j] == 0) {
					case 0:						// c1[j] = m1[moms[i]][j];
						c1[j] = m1[dads[i]][j];	// c2[j] = m2[moms[i]][j];
						c2[j] = m2[dads[i]][j];	// }
					case 1:						// else {
						c1[j] = m2[dads[i]][j];	// c1[j] = m2[moms[i]][j];
						c2[j] = m1[dads[i]][j];	// c2[j] = m1[moms[i]][j];
				}								// }
			}
			if (sons[i]) 
			{
				f1_mat[i] = c1;
			}
			else 
			{
				f1_mat[i] = c2;
			}
		}
}


/* 
Using main() for burnin 
*/

int main(double h = 0.1, double s = 0.01, int Us_factor = 2,
		 int nTot = 10^3, int N = 10^3, double invSize = 0.05, 
		 int burnin=100)
{
	// Random distributions
	unsigned int gsl_ran_binomial(const gsl_rng *r, double p, unsigned int k) // p = prob. success, k = n_trials

	unsigned int gsl_ran_discrete_t *gsl_ran_discrete_preproc(size_t K, const double *P)


	// Fitness effects for indiv. loci
	// W.i  <-  c(w_AA, w_Aa, w_aa)
	double W_i  <-  {1, 1 - h*s, 1 - s};

	//number of simulations ** only relevant for full simulation code
	//	sims  = 100*N/2

	// Mutation rate and equilibrium frequencies prior to inversion
	double U = Us_factor*s;
	double u = U/nTot;
	double qHat = U/(nTot*h*s);
			
	// number of loci captured by inversion
	int n = nTot*invSize;

	// initiate vector of moms & dads
	int indSeq[N/2];
	for (int i=0 , i<(N/2), i++)
	{
		indSeq[i] = i+1;
	}

 	// Initiate chromosome matrices (rows = individuals, cols = loci spanned by inversion)
	// and populate with mutations
	int Xf_mat[N/2][n] = gsl_ran_binomial(r, qHat, (n*N/2)); // Maternally-inherited X's in females
	int Xf_pat[N/2][n] = gsl_ran_binomial(r, qHat, (n*N/2)); // Paternally-inherited X's in females
	int Xm [N/2][n]    = gsl_ran_binomial(r, qHat, (n*N/2)); // X chromosomes in Males
	int Y [N/2][n]     = gsl_ran_binomial(r, qHat, (n*N/2)); // Y chromosomes

	// initiate variables used in generation loop
	double Wf_i
	double Wm_i

	// burnin to equilbrium prior to inversion
	for (t, t<burnin, t++) {
		// print progres to screen
		printf("burnin generation %s / %d is complete \r", t, burnin)
		// Relative fitness & sampling probs. 
		double Wf_i = relFitness(Xf_mat, Xf_pat, W_i)
		double Wm_i = relFitness(Xm, Y, W_i)
>> CONTINUE TRANSLATION FROM HERE <<
		// Sample parent gametes
		int moms = sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wf_i/(N/2))
		int dads = sample(c(1:(N/2)), size=N, replace=TRUE, prob=Wm_i/(N/2))
		// Meiosis and fertilization
		// for sons
		Xm.f1         <-  freeRecombine_sampOvules(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[1:(N/2)])
		Y.f1          <-  freeRecombine_sampSperm(m1=Xm, m2=Y, n=n, dads=dads[1:(N/2)], inverted=rep(FALSE, times=(N/2)), sons=TRUE)
		// for daughters
		Xf.mat.f1     <-  freeRecombine_sampOvules(m1=Xf.mat, m2=Xf.pat, n=n, moms=moms[((N/2)+1):N])
		Xf.pat.f1     <-  freeRecombine_sampSperm(m1=Xm, m2=Y, n=n, dads=dads[((N/2)+1):N], inverted=rep(FALSE, times=(N/2)), sons=FALSE)
		// Mutation
		Xf.mat.f1[Xf.mat.f1 == 0]  <-  rbinom(n=length(Xf.mat.f1[Xf.mat.f1 == 0]), size=1, prob=u)
		Xf.pat.f1[Xf.pat.f1 == 0]  <-  rbinom(n=length(Xf.pat.f1[Xf.pat.f1 == 0]), size=1, prob=u)
		Xm.f1[Xm.f1 == 0]          <-  rbinom(n=length(Xm.f1[Xm.f1 == 0]),         size=1, prob=u)
		Y.f1[Y.f1 == 0]            <-  rbinom(n=length(Y.f1[Y.f1 == 0]),           size=1, prob=u)
		// F1's become adults in next generation
		Xm      <-  Xm.f1
		Y       <-  Y.f1
		Xf.mat  <-  Xf.mat.f1
		Xf.pat  <-  Xf.pat.f1
	}

	// assign initial state chromosome matrices
	Xm.init      <-  Xm
	Y.init       <-  Y
	Xf.mat.init  <-  Xf.mat
	Xf.pat.init  <-  Xf.pat
}


hist(colSums(rbind(Xf.mat, Xf.pat, Xm, Y)), col=2)#, add=T)
hist(colSums(Y), col=2)#, add=T)
