/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004 Anthony Brockwell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * -------------------------------------------------------------------
 *
 * This random number generation module began
 * (around 1995) as a simple collection of routines generating
 * different types of random variables by applying various transformations
 * to uniform random number generators.
 *
 * Sections of this code will probably be replaced with
 * functions from the GSL library as they become available
 * in the future.  (One might expect that given the effort
 * put into the GSL routines by the developers, those routines
 * will be more efficient than the routines in this library.)
 *
 * Development history:
 *
 * sometime 2000: added the Gamma random number generator.
 *            see "Squeeze methods for generating Gamma variates",
 *            JASA 1980, Schmeiser and Lal.
 * Mar. 2001: added a Poisson random number generator, see
 *            Algorithm 8 of Kemp & Kemp, Poisson Random Variate
 *            Generation, Applied Statistics, JRSS 1991.
 * Apr. 2001: added a Beta random number generator -
 *            (see Gelman, Carlin, Stern, Rubin, p. 481.)
 * Aug. 2001: geometric: See Gelman et. al. Appendix,
 *            also fixed bug with Poisson generator when $\lambda < 0.5$.
 * Jan. 2002: started adding density functions
 * Apr. 2002: added multivariate normal and conditional MVN routines
 * Jun. 2002: removed reliance on LAPACK++ since it had bugs
 *            in finding eigenvalues (critical for MVN routines!)
 *            still uses LAPACK library
 * May  2003: incorporated into mslib library and replaced
 *            R250 generator with SPRNG, which can be used for
 *            parallel generation of random numbers, see
 *            http://sprng.cs.fsu.edu/Version2.0/install.html
 * Jul  2004: replaced SPRNG with GSL (no parallel support, but
 *            eliminates dependence on MPI library) 
*/

#include "abrand.h"
#define SIMPLE_SPRNG
#define USE_MPI
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#ifdef USE_LAPACK
extern "C" {
  void dsyev_(char *job, char *uplo, int *n, double *a,
	      int *lda, double *evals, double *work, int *lwork,
	      int *output);
}
#endif
const double zero=1e-10;


namespace mslib {

gsl_rng *the_generator;

void init_generator(long int sd)
{
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;  // get default type

  the_generator = gsl_rng_alloc(T);
  gsl_rng_set(the_generator, sd);
}

// uniform random number generator, returns a double
// in the interval [0,1)
// most other generation routines are based on this
double unif_random()
{
  return gsl_rng_uniform(the_generator);
}


// Code for sampling of a binomial random variable
// based on Alg. in Statistical Computing, Kennedy & Gentle,
// Dekker, p. 222

int bino_random(int n, double p)  
{
  if ((n<0) || (p<0) || (p>1))
    {
      cout << "Bad Binomial Parameters! (" << n << "," << p << ")" << endl;
      return 0;
    }

  int i, a, count=0, m=n, k=0;
  double q=p, y=0.0, h=1.0, s, g, z;

  while (m>40)
    {
      if (m%2==0)
	{
	  --m;
	  if (unif_random()<=q)
	    ++k;
	}
      a = (m+1)/2;
      s = beta_random(a, a);
      g = h*s;
      z = y+g;
      if (z<=p)
	{
	  y = z;	  h -= g;
	  q = (p-z)/h;    k += a;
	}
      else
	{
	  h = g;
	  q = (p-y)/h;
	}
      m = a-1;
    }
  
  for (i=0; i<m; ++i)
    count += (unif_random() <= q);
  return count + k;
}

// Acceptable code for sampling from a t-distribution
double t_random(int df)
{
  double xx=0.0,yy;
  xx = chi2_random(df)/df;
  yy = norm_random()/sqrt(xx);
  return yy;
}

// Code for sampling from an exponential distribution
double exp_random(double lambda)
{
  double u=unif_random();
  return (-log(u)/lambda);
}

// Standard Box-Muller transformation for generating
// a realization of a normal random variable
double norm_random()
  {
    double pi=M_PI;
    double x1,x2,w;
    static double y1,y2;
    
    static int toggle=0;
    toggle = !toggle;
    
    if (toggle) {
      do {
        x1 = 2.0 * unif_random() - 1.0;
        x2 = 2.0 * unif_random() - 1.0;
        w = x1 * x1 + x2 * x2;
      } while ( w >= 1.0 );
      w = sqrt( (-2.0 * log( w ) ) / w );
      y1 = x1 * w;
      y2 = x2 * w;
    }

    if (toggle)
      return y1;
    else
      return y2;
  }

// normal random variable with specified mean and variance
double norm_random(double mu, double sig2)
{
  return (norm_random()*sqrt(sig2)+mu);
}

double norm_density(double x, double mu, double sig2)
{
  double tx = 1.0/sqrt(2*M_PI*sig2);
  tx *= exp(-0.5*(x-mu)*(x-mu)/sig2);
  return tx;
}

double lognorm_density(double x, double mu, double sig2)
{ 
  return ( -log(sqrt(2*M_PI*sig2)) - 0.5*(x-mu)*(x-mu)/sig2 ); 
}

// Here is a new attempt to write a function which
// generates realizations of Gamma random variables.
// This one is an implementation of Algorithm G2PE


double gammadensity(double x, double alpha)		// beta assumed = 1
{
  double tx;
  tx = exp(-x) * pow(x,alpha-1) / gsl_sf_gamma(alpha);
  return tx;
}

double gammafof(double x, double x3)
{
  return exp(x3*log(x/x3)+x3-x);
}

double gam_random(double alpha, double beta)
// use the method of rejection:
// we need a comparison function (see above)
{
  bool done=false, restart;
  double x,u,v;
  double x1,x2,x3,x4,x5,D,lambdal,f1,f2,f4,p1,p2,p3,lambdar;

  if (alpha<1)
    {
      cout << "Gamma Badness: Alpha= " << alpha << "< 1!"  << endl;
      return 0;
    }
  if (fabs(alpha-1.0)<1e-6)
    return exp_random(beta)*alpha;

  // initialization: Steps 1 and 2
  x3 = alpha-1;
  D = sqrt(x3);
  lambdal = 1;
  x1 = f2 = x2 = 0;
  if (D<x3)
    {
      x2 = x3-D;
      lambdal = 1-x3/x2;
      x1 = x2+1/lambdal;
      f2 = gammafof(x2,x3);
    }
  x4 = x3 + D;
  lambdar = 1-x3/x4;
  x5 = x4 + 1/lambdar;
  f4 = gammafof(x4,x3);
  p1 = x4-x2;
  p2 = p1 - f2/lambdal;
  p3 = p2 + f4/lambdar;

  int iteration=0;

  // iteration to get Gamma
  while (!done) {
    ++iteration;
    if (iteration==1000)
      {
	cout << "Gamma badness with " << alpha << "," << beta << "!" << endl;
	return 0;
      }
    restart = false;
    u = unif_random()*p3;   v = unif_random();
    if (u>p1)
      {
        if (u>p2)
          { // Step 5
            u=(u-p2)/(p3-p2);
            x=x4-log(u)/lambdar;
            v=v*f4*u;
            if (v<=f4*(x5-x)/(x5-x4))
              done = true;
          }
        else
          { // Step 4
            u=(u-p1)/(p2-p1);
            x=x2-log(u)/lambdal;
            if (x<0)
              restart=true;
            v = v*f2*u;
            if (v<=f2*(x-x1)/(x2-x1))
              done = true;
          }
      }
    else
      {
        x = x2 + u;
        if ((x>x3) && (v <= f4+(x4-x)*(1-f4)/(x4-x3)))
          done = true;
      }
    // Step 6.
    if ((!restart) && (!done))
      {
        if (log(v)<=x3*log(x/x3)+x3-x)
          done=true;
      }
  }

  return (x/beta);
}

// using the Gamma generator, we can get Chi^2 random variables
double chi2_random(double dof)
{
  double a=dof/2.0, b=0.5;
  return gam_random(a,b);
}

int specialprandom(double lambda)
{
  // For values of lambda < 0.5, the algorithm below doesn't work.
  // It is pretty simple to simulate, though.
  double u = unif_random(), pj = exp(-lambda), cumsum = pj;
  int found=0, j=0;
  while (!found)
    {
      if (u <= cumsum)
	found = 1;
      else
	{
	  ++j;
	  pj *= lambda / j;
	  cumsum += pj;
	}
    }
  return j;
}

int pois_random(double lambda)
{
  if (lambda < 0.5)
    return specialprandom(lambda);

  int i;
  static double lastlambda=-1.0;
  static double prr, frr, prlam, frlam;
  double r = (int) (lambda+0.5), tx;
  double alpha = lambda - r, g = 1.0/sqrt(2*M_PI*r);

  // 1. compute the various constants if lambda has changed
  if (fabs(lambda-lastlambda)>1e-10) {
    prr = g*(1-1/(12*r+0.5+293.0/(720*r)));
    tx = 12.0*r + 138134432.0/105880005.0;
    tx = 30557.0/4508.0 / tx;
    tx += 12*r + 15.0/14.0;
    tx = 23.0/15.0 / tx;
    frr = 0.5 + 2*g/3*(1-tx);
    prlam = prr*(r+2*alpha/3 - alpha*alpha/4 - alpha*alpha/(18*r))
      /(r+2*alpha/3 + alpha*alpha/4 - alpha*alpha/(18*r));
    frlam = frr - alpha*prr*(r+alpha/2-alpha*alpha/60-alpha*alpha/(20*r))
      /(r+alpha/2+3*alpha*alpha/20-alpha*alpha/(20*r));
    lastlambda = lambda;
  }

  double u = unif_random(), p;
  int x;

  // Squeeze step
  int downup = 0;
  if (u <= 0.5)
    downup = -1;
  else
    if (u >= 0.5+7*g/6)
      downup = 1;
    else
      if (u > frlam)
	downup = 1;
      else
	downup = -1;

  if (downup==-1) // downward search
    {
      if (u < prlam)
	return (int)(r+0.5);
      p = prlam;
      for (i=0 ; i<r ; ++i)
	{
	  u -= p;
	  p = (r-i)*p/lambda;
	  if (u < p)
	    return (int)(r-i-1+0.5);
	}
      return 0; // failsafe?
    }
  else // upward search
    {
      u = 1-u;
      p = prlam;
      for (i=(int)(r+1+0.5) ; i<999999999 ; ++i)
	{
	  p = p*lambda/i;
	  if (u<p)
	    return (int)(i+0.5);
	  u -= p;
	}
      return 999999999;  // failsafe?
    }
}

double beta_random(double a, double b)
{
  if ((a<1) || (b<1))
    {
      cout << "Beta Badness with " << a << ", " << b << endl;
      return 0;
    }
   
  // first get a couple of chi^2 r.v.s
  double xa = chi2_random(2*a), xb = chi2_random(2*b);
  return xa/(xa+xb);
}

double logbeta_density(double x, double a, double b)
{
  double tx;
  tx = gsl_sf_lngamma(a+b) - gsl_sf_lngamma(a) - gsl_sf_lngamma(b)
    + (a-1)*log(x) + (b-1)*log(1.0-x);
  return tx;
}

double cauchy_random(double location)
{
  double u = unif_random()*M_PI;
  return tan(u) + location;
}

int geo_random(double p)
{
  double beta = p/(1.0-p);
  double ll = exp_random(beta);
  return pois_random(ll)+1;
}
   
//-----------------------------------------------------
// Here is something pretty handy too:
// approximations to the cdf and inverse cdf of a
// normal random variable.
//-----------------------------------------------------
// These functions use an approximation given
// by Abramowitz and Stegun (1964), Handbook of Mathematical Functions.

double norm_cdf(double x)
{
  static double oneonroot2pi = 1.0/(sqrt(2.0*M_PI));
  static double b1=0.31938153,b2=-.35656378,b3=1.7814779;
  static double b4=-1.8212560,b5=1.33027443,p=.2316419;
  double x1=fabs(x);
  double t=1.0/(1+p*x1),phi;
  phi = oneonroot2pi * exp(-x1*x1/2) * (((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t);
  if (x>0.0) phi=1-phi;
  return(phi);
}

double norm_cdf(double x, double mu, double sig2)
{
  return norm_cdf((x-mu)/sqrt(sig2));
}

double norm_invcdf(double x)  {
  static double c0=2.515517,c1=0.802853,c2=0.010328;
  static double d1=1.432788,d2=0.189269,d3=0.001308;
  double p1,t,phiinv;
  if (x>0.5) p1=1-x;  else p1=x;
  t = sqrt(-log(p1*p1));
  phiinv = t-((c2*t+c1)*t+c0)/(((d3*t+d2)*t+d1)*t+1);
  if (x<0.5) phiinv=-phiinv;
  return(phiinv);
}


double exp_density(double x, double lambda)
{
  return lambda*exp(-lambda*x);
}

double logexp_density(double x, double lambda)
{
  return log(lambda) - lambda*x;
}

double exp_cdf(double x, double lambda)
{
  if (x>0)
    return 1.0-exp(-lambda*x);
  else
    return 0.0;
}

double exp_invcdf(double x, double lambda)
{
  return (-log(1-x)/lambda);
}

double logmvn_density(Vector& x, Vector *mn, Matrix *cov, 
		     Matrix *siginv, double covdet)
{
  double tx, det;
  int dim = x.nrows();
  
  // now compute density
  Vector tv = x - (*mn), tv2 = tv;
  tx = log(2*M_PI)*(-0.5*dim);
  if (siginv==NULL) {
    if (cov->singular(&det))
      {
	cout << "Singular Matrix in LOGMVNDENSITY!" << endl;
      }
    else
      tx -= 0.5*log(det);  
    cov->symm_solve_system(&tv2);
    tx -= 0.5*tv.dot(tv2);
  } else {
    tx -= 0.5*log(covdet);
    tx -= 0.5*tv.dot((*siginv)*tv2);
  }
  return tx;
}

Matrix GramSchmidtBasis(Matrix& h)
  // returns basis with $h$ at the bottom,
  // rows (except for $h$) are normalized
{
  int dim=h.ncols(), basisdim=h.nrows(), i, j, kk;
  Vector *xi = new Vector[dim+basisdim],
    *ui = new Vector[dim+basisdim],
    *yi = new Vector[dim+basisdim];
  Matrix V(dim,dim);
  double tx;

  // step 1: Find V, m0, E
  // We get V by Gram-Schmidt orthogonalization

  // start with dim+basisdim x's, identity Matrix along with h Matrix
  for (i=0 ; i<basisdim ; ++i)
    xi[i]=h.extract_row(i);
  for (i=0 ; i<dim ; ++i)
    {
      xi[i+basisdim] = xi[0];  xi[i+basisdim].zeroes();
      xi[i+basisdim][i] = 1.0;
    }
  
  // then apply Gram-Schmidt procedure to construct
  // an orthonormal basis, including h
  yi[0]=xi[0];   ui[0]=yi[0]/yi[0].norm();
  for (i=1 ; i<dim+basisdim ; ++i)
    {
      yi[i] = xi[i];
      for (j=0 ; j<i ; ++j)
	yi[i] = yi[i] - ui[j]*(xi[i].dot(ui[j]));
      ui[i] = yi[i];
      tx = ui[i].norm();
      if (tx>zero)
	ui[i] /= tx;
    }

  // now we can fill in V with the non-zero uis
  for (i=0,j=0 ; i<dim+basisdim ; ++i)
    if (ui[i].norm()>zero)
      {
	// copy into V
	for (kk=0 ; kk<dim ; ++kk)
	  V[dim-1-j][kk] = ui[i][kk];
	++j;
      }

  for (i=0 ; i<basisdim ; ++i)
    for (kk=0 ; kk<dim ; ++kk)
      V[dim-1-i][kk] = h[basisdim-1-i][kk];

  delete[] xi;
  delete[] ui;
  delete[] yi;
  return V;
}

Matrix cov_msqrt(Matrix& sig2, int reduce, int& error)  // returns square root of cov. Matrix
{
  int i,j,d = sig2.nrows();
  Matrix temp(d,d),*temp2;
  static double scratch[256];
  static int scratchsize=256;
  error = 0;  // no error

  if (d==1)
    {
      if (sig2[0][0]<0)
	error = 1;
      else
        temp[0][0]=sqrt(sig2[0][0]);
      return temp;
    }

  int output;
  Vector evals(d);
  Matrix evecs(d,d);
#ifdef USE_LAPACK
  evecs = sig2;
  dsyev_("V","U",&d,&evecs[0][0],&d,&evals[0],scratch,&scratchsize,
	 &output);
  evecs=evecs.transpose();
#else
  complex *cevals = sig2.eigenvalues();
  for (i=0 ; i<d ; ++i)
    evals[i] = cevals[i].real();  // they are only real anyway
  delete[] cevals;
  evecs = sig2.symm_eigenvectors(&evals[0]);
#endif

  if (!reduce)
    for (i=0 ; i<d ; ++i)
      for (j=0 ; j<d ; ++j)
	if (evals[j]>zero)
	  temp[i][j] = evecs[i][j]*sqrt(evals[j]);
	else if (fabs(evals[j])<=zero)
	  temp[i][j] = 0;
        else
	  error = 1;
  else
    {
      int realj;
      // same procedure but eliminate zero columns
      for (j=0,realj=0 ; j<d ; ++j)
	if (fabs(evals[j])>zero)
	  {
	    for (i=0 ; i<d ; ++i)
	      temp[i][realj] = evecs[i][j]*sqrt(evals[j]);
	    ++realj;
	  }
      if (realj==0)
	{
	  realj=1;
	  for (j=0 ; j<d ; ++j)
	    temp[0][j]=0;
	}
      temp2 = new Matrix(d,realj);
      for (j=0 ; j<realj ; ++j)
	for (i=0 ; i<d ; ++i)
	  (*temp2)[i][j] = temp[i][j];
      temp = *temp2;
      delete temp2;
    }

  return temp;
}

Vector mvn_with_sqrtvar(Vector& mean, Matrix& sqrtsig2)
{
  int i,d = sqrtsig2.ncols();
  Vector temp(d);

  for (i=0 ; i<d ; ++i)
    temp[i] = norm_random(0,1);
  temp = sqrtsig2*temp + mean;
  return temp;
}

Vector mvn_random(Vector& mean, Matrix& sig2, int *error)   
  // Matrix argument is cov. instead of sqrt. of cov.
{
  int i,j,d = mean.nrows(),localerror;
  Matrix sqrtsig2; 
  sqrtsig2 = cov_msqrt(sig2,0,localerror);
  if (error!=NULL)
    *error = localerror;
  return mvn_with_sqrtvar(mean,sqrtsig2);
}

double condnic(double mean, double sd, double k)
  // returns normal given it is greater than k
{
  double kpos = (k-mean)/sd, var=sd*sd;
  double tx, logu;
  if (kpos<1)
    {
      // repeat draw
      tx = norm_random(mean, var);
      for ( ; tx < k ; )
	tx = norm_random(mean, var);
    }
  else
    {
      // accept/reject with exponential envelope
      double delta = sd;
      double k2=k+delta;
      double logf1 = lognorm_density(k, mean, var),
	logf2 = lognorm_density(k2, mean, var);
      double lambda = (logf1-logf2)/delta;
      double logalpha = log(lambda) - logf1;
      int ok = 0;
      
      while (!ok)
	{
	  tx = exp_random(lambda)+k;
	  logu = logalpha + lognorm_density(tx, mean, var)
	    - log(lambda) + lambda*(tx-k);
	  ok=(unif_random()<=exp(logu));
	}
    }
  return tx;
}


void covfixup(Matrix& cov)
{
  int i,j,n=cov.nrows();
  for (i=0 ; i<n ; ++i)
    cov[i][i] = fabs(cov[i][i]);
  for (i=0 ; i<n ; ++i)
    for (j=0 ; j<n ; ++j)
      if (fabs(cov[i][j])<zero)
	cov[i][j] = 0.0;
}

Vector cond_mvn_random(Vector& mean, Matrix& var, 
		 Matrix& Hs, Vector& rels, Vector& ks,
		 Matrix *suppliedP, Matrix *suppliedPinv,
		 int fillinP)
  // This is a general purpose routine which draws from
  // a multivariate normal (mean,var) X, given H rel. k
  // "rel" is a vector of integers specifying a relationship
  // to hold for the corresponding elements of HX and k,
  // -1 is <, 0 is =, 1 is >

  // If suppliedP, suppliedPinv are NOT NULL, then
  //   1. if fillinP is 1, they are computed, used, and filled in
  //   2. if fillinP is 0, they are simply used
  // Note that they only depend on H and the dimension,
  // and supplying them saves solving a system of eqns as
  // well as a call to GramSchmidtBasis.
{
  int i,j,jj,iter,neqcon,ncon=Hs.nrows(),nineqcon,satisfied,dim=mean.nrows();
  Matrix P,H;
  Vector tv,rel,k;
  double tx,tvar;

  k = ks;   rel = rels;   H = Hs;

  // 0. Preprocess so that equality constraints are at the bottom
  int negate;
  Matrix H2(ncon,dim);
  Vector rel2(ncon), k2(ncon);
  for (i=0,j=0 ; i<ncon ; ++i)
    {
      if ((int)(floor(rel[i]+0.5))!=0)
	{
	  if (floor(rel[i]+0.5)==-1)
	    negate = -1;
	  else
	    negate = 1;
	  for (jj=0 ; jj<dim ; ++jj)
	    H2[j][jj] = negate*H[i][jj];
	  rel2[j] = 1;  // they are all translated to >
	  k2[j] = negate*k[i];
	  ++j;
	}
    }
  for (i=0 ; i<ncon ; ++i)
    {
      if ((int)(floor(rel[i]+0.5))==0)
	{
	  for (jj=0 ; jj<dim ; ++jj)
	    H2[j][jj] = H[i][jj];
	  rel2[j] = rel[i];
	  k2[j] = k[i];
	  ++j;
	}
    }
  H = H2;
  k = k2;
  rel = rel2;

  // 1. determine conditional normal based only on
  //    equality constraints

  if (suppliedP==NULL)
    P = GramSchmidtBasis(H);
  else
    if (fillinP)
      {
	P = GramSchmidtBasis(H);
	*suppliedP = P;
	*suppliedPinv = P.inverse();
      }
    else
      P = *suppliedP;
  
  Vector nmean = P*mean, nmean2;
  Matrix nvar = P*var*P.transpose(), nvar2;
  Vector Y(dim);

  Matrix sig11, sig12, sig22, s22inv;
  Vector mu1, mu2;

  for (neqcon=0 ; neqcon<ncon && rel[ncon-1-neqcon]==0 ; ++neqcon);  
  nineqcon = ncon - neqcon;
  if (neqcon>0)
    {
      // fill in last few bits of Y before proceeding
      for (i=0 ; i<neqcon ; ++i)
	Y[i+dim-neqcon] = k[i+ncon-neqcon];
      // replace nmean and nvar by conditional things
      Vector x2(neqcon);
      sig11 = nvar.submatrix(0,0,dim-neqcon,dim-neqcon);
      sig12 = nvar.submatrix(0,dim-neqcon,dim-neqcon,dim);
      sig22 = nvar.submatrix(dim-neqcon, dim-neqcon, dim, dim);
      mu1.resize(dim-neqcon);
      mu2.resize(neqcon);
      for (i=0 ; i<dim-neqcon ; ++i)
	mu1[i] = nmean[i];
      for (i=0 ; i<neqcon ; ++i) {
	mu2[i] = nmean[i+dim-neqcon];
	x2[i] = k[i+nineqcon];
      }
      s22inv = sig22.inverse();
      nmean = mu1 + sig12*s22inv*(x2-mu2);
      nvar = sig11 - sig12*s22inv*sig12.transpose();
    }

  // do the remaining part by repetition
  int newdim=nmean.nrows(),error;
  satisfied = 0;
  Vector X(newdim),Xsub,modk;
  Matrix sv, Q(newdim,newdim);
  iter = 0;

  if (newdim>0)
    {
      if (nineqcon>0) 
	{
	  // peel off one constraint
	  tx = nmean[newdim-1];
	  tvar = nvar[newdim-1][newdim-1];
	  sig11 = nvar.submatrix(0,0,newdim-1,newdim-1);
	  sig12 = nvar.submatrix(0,newdim-1,newdim-1,newdim);
	  sig22 = nvar.submatrix(newdim-1,newdim-1,newdim,newdim);
	  s22inv = sig22.inverse();
	  if (newdim>1)
	    nvar2 =  sig11 - sig12*s22inv*sig12.transpose();
	  else
	    nvar2 = sig22;
	  sv = cov_msqrt(nvar2, 1, error);
	  mu1.resize(newdim-1);
	  for (i=0 ; i<newdim-1 ; ++i)
	    mu1[i] = nmean[i];
	  mu2.resize(1);
	  mu2[0] = nmean[newdim-1];	
	}
      iter = 0;
      while (!satisfied)
	{
	  ++iter;
	  if (iter==10000)
	    cout << "Warning: Iteration exceeding 10000." << endl;
	  if (nineqcon>0)
	    {
	      // further reduce the problem
	      X[newdim-1] = condnic(tx,sqrt(tvar),k[nineqcon-1]);
	      // and recompute nmean, nvar and newdim
	      if (newdim>1)
		nmean2 = mu1 + sig12*(1.0/sig22[0][0])*(X[newdim-1]-mu2[0]);
	    }

	  // now we just need nmean, nvar, newdim
	  if (newdim>1)
	    {
	      if (nineqcon==0)
		{
		  sv = cov_msqrt(nvar, 1, error);
		  Xsub = mvn_with_sqrtvar(nmean, sv);
		  for (i=0 ; i<Xsub.nrows() ; ++i)
		    X[i] = Xsub[i];
		}
	      else
		{
		  Xsub = mvn_with_sqrtvar(nmean2, sv);
		  for (i=0 ; i<Xsub.nrows() ; ++i)
		    X[i] = Xsub[i];
		}
	    }

	  // and now we check constraints
	  satisfied = 1;
	  for (i=0 ; (i<nineqcon) && (satisfied) ; ++i)
	    if (X[i+dim-ncon]<k[i])
	      satisfied = 0;
	}
    }

  // copy the result into the first part of Y
  for (i=0 ; i<newdim ; ++i)
    Y[i] = X[i];

  // back-transform
  if (suppliedPinv == NULL)
    P.solve_system(&Y);
  else
    Y = (*suppliedPinv)*Y;
  return Y;
}

Vector diri_random(Vector& nums)
{
  int i,dim=nums.nrows();
  Vector retval(dim);
  for (i=0 ; i<dim ; ++i)
    retval[i] = gam_random(nums[i], 1.0);
  retval /= retval.sum();
  return retval;
}

double logdiri_density(Vector& x, Vector& parms)
{
  int i,n=x.size();
  double tx=0.0,s=parms.sum();
  for (i=0 ; i<n ; ++i)
    {
      tx += (parms[i]-1.0)*log(x[i]);
      tx -= gsl_sf_lngamma(parms[i]);
    }
  tx += gsl_sf_lngamma(s);
  return tx;
}

static double logtable[500];
static int logsfilled=0;

double logpois_density(double x, double lambda)
{
  int i,ix = (int)(x+0.5);

  if (!logsfilled)
    {
      logtable[0] = 0.0;
      for (i=1 ; i<500 ; ++i)
	logtable[i] = log((double)i)+logtable[i-1];
      logsfilled = 1;
    }

  // returns log poisson pmf (x should be a non-neg. integer)
  double ld = -lambda + x*log(lambda);
  ld -= logtable[ix<500?ix:499];
  for (i=500 ; i<ix ; ++i)
    ld -= log((double)i);
  return ld;
}

double logpois_density_logl(double x, double loglambda)
{
  int i,ix = (int)(x+0.5);

  if (!logsfilled)
    {
      logtable[0] = 0.0;
      for (i=1 ; i<500 ; ++i)
	logtable[i] = log((double)i)+logtable[i-1];
      logsfilled = 1;
    }
  
  // returns log poisson pmf (x should be a non-neg. integer)
  double ld = -exp(loglambda) + (x*loglambda);
  ld -= logtable[ix<500?ix:499];
  for (i=500 ; i<ix ; ++i)
    ld -= log((double)i);
  return ld;
}


Vector int_to_mvb(long int li, int dim)
{
  int i;
  Vector retval(dim);
  for (i=0 ; i<dim ; ++i)
    {
      retval[dim-1-i] = (li%2);
      li = li >> 1;
    }
  return retval;
}

long int mvb_to_int(Vector& x)
{
  int i,dim = x.size();
  long int place=1,total=0;
  for (i=0 ; i<dim ; ++i)
    {
      total += x[dim-1-i]*place;
      place = place << 1;
    }
  return total;
}

double quad_form(Vector& tv, Matrix& zetas)
  // returns 0.5*tv.dot(zetas*tv), assuming all elements of tv=1 or 0,
  // and assuming zetas is symmetric
{
  int i,j,dim = tv.nrows(),nones;
  double retval=0.0;
  int *ones = new int[dim];
  for (i=0,j=0 ; i<dim ; ++i)
    if (tv[i]==1)
      ones[j++] = i;
  nones = j;

  for (i=0 ; i<nones ; ++i)
    for (j=i+1 ; j<nones ; ++j)
      retval += zetas[ones[i]][ones[j]];

  delete[] ones;

  return retval;
}

double build_mvb_table(double *table, double *cumtable,
		     Vector& betas, Matrix& zetas)
{
  // puts all the exp(...) pmfs into table,
  // along with cumulative values into cumtable,
  // and the total is returned

  int i,dim=betas.size();
  long int li, len;
  Vector logbetas(dim), tv;
  double tx,nconst;

  len = 1 << dim;

  for (i=0 ; i<dim ; ++i)
    logbetas[i] = log(betas[i]);

  for (li = 0 ; li<len ; ++li)
    {
      tv = int_to_mvb(li, dim);
      tx = tv.dot(logbetas) + quad_form(tv,zetas);
      table[li] = exp(tx);
      if (li==0)
	cumtable[li] = table[li];
      else
	cumtable[li] = cumtable[li-1] + table[li];
    }
  nconst = cumtable[len-1];
  return nconst;
}

Vector mvbern_random(Vector& betas, Matrix& zetas, Vector& groups)
  // (1) betas = log(p/(1-p)) when all other x's=0
  // (2) zetas should be a symmetric matrix
  //     containing interaction terms (0=independence)
  //     it should have zeros down the diagonal
{
  int i,j,k,dim=betas.size(),ng = groups.size();
  Vector retval(dim), marginal;

  if (ng>1)
    {
      int i1,i2;
      Matrix subzeta;
      Vector subbeta,subgroups(1),tv;
      for (i=0 ; i<ng ; ++i)
	{
	  i1 = (int)(0.5+groups[i]);
	  i2 = i<ng-1 ? groups[i+1] : dim;
	  subzeta = zetas.submatrix(i1,i1,i2,i2);
	  subbeta = betas.subvector(i1,i2);
	  subgroups[0] = 0;
	  tv = mvbern_random(subbeta,subzeta,subgroups);
	  if (i==0)
	    retval = tv;
	  else
	    retval.append_bottom(tv);
	}
    }
  else
    {
      // do this one by tabulation
      Vector tv(dim);
      long int li, len = (1 << dim);
      double tx, *table = new double[len], *cumtable = new double[len],
	nconst;

      nconst = build_mvb_table(table, cumtable, betas, zetas);

      double u = unif_random()*nconst;
      long int lower=0,upper=len-1,mid;
      while (lower!=upper)
	{
	  mid = (lower+upper)/2;
	  if (u>cumtable[mid])
	    lower = mid+1;
	  else
	    upper = mid;
	}

      retval = int_to_mvb(lower, dim);

      delete[] table;   delete[] cumtable;
    }

  return retval;
}

double logmvbern_pmf(Vector& x, Vector& betas, Matrix& zetas, Vector& groups)
{
  int i,i1,i2;
  double retval=0;
  int dim=x.size(),ngrps=groups.size();
  Matrix subzeta;
  Vector subbeta, subgroups(1), subx;
  subgroups[0] = 0.0;

  if (ngrps==1)
    {
      long int li, len = (1 << dim);
      double tx, *table = new double[len], *cumtable = new double[len],
	nconst;

      nconst = build_mvb_table(table, cumtable, betas, zetas);
      li = mvb_to_int(x);
      retval = log(table[li]) - log(nconst);

      delete[] table;   delete[] cumtable;
    }
  else
    for (i=0 ; i<ngrps ; ++i)
      {
	i1 = (int)(0.5+groups[i]);
	i2 = i<ngrps-1 ? (int)(0.5+groups[i+1]) : dim;
	subzeta = zetas.submatrix(i1,i1,i2,i2);
	subbeta = betas.subvector(i1,i2);
	subx = x.subvector(i1,i2);
	retval += logmvbern_pmf(subx, subbeta, subzeta, subgroups);
      }

  return retval;
}

Vector logmvbern_pmf(Matrix& x, Vector& betas, Matrix& zetas, Vector& grps)
  // note: now is horribly inefficient!  should build table only once!
{
  int i,dim = x.nrows();
  Vector retval(dim),tv;

  for (i=0 ; i<dim ; ++i)
    {
      tv = x.extract_row(i);
      retval[i] = logmvbern_pmf(tv, betas, zetas, grps);
    }

  return retval;
}

Matrix mvbern_cov(Vector& betas, Matrix& zetas, Vector& grps)
{
  int i,i1,i2,j,k;
  double tx;
  int dim=betas.size(),ngrps=grps.size();
  Matrix subzeta;
  Matrix retcov(dim,dim),sr;
  retcov.zeroes();
  Vector subbeta, subgroups(1), subx, tv;
  subgroups[0] = 0.0;

  if (ngrps==1)
    {
      long int li, len = (1 << dim);
      double tx, *table = new double[len], *cumtable = new double[len],
	nconst;

      nconst = build_mvb_table(table, cumtable, betas, zetas);

       for (li=0 ; li<len ; ++li)
	{
	  tx = table[li];
	  tv = int_to_mvb(li, dim);
	  
	  retcov = retcov + tv.outer(tv) * (table[li]/nconst);
	}

      delete[] table;   delete[] cumtable;
    }
  else
    for (i=0 ; i<ngrps ; ++i)
      {
	i1 = (int)(0.5+grps[i]);
	i2 = i<ngrps-1 ? (int)(0.5+grps[i+1]) : dim;
	subzeta = zetas.submatrix(i1,i1,i2,i2);
	subbeta = betas.subvector(i1,i2);
	sr = mvbern_cov(subbeta,subzeta,subgroups);
	for (j=i1 ; j<i2 ; ++j)
	  for (k=i1 ; k<i2 ; ++k)
	    retcov[j][k] = sr[j-i1][k-i1];
      }

  return retcov;
}

Vector mvbern_marginal(Vector& betas, Matrix& zetas, Vector& grps)
{
  int i,i1,i2,j;
  double tx;
  int dim=betas.size(),ngrps=grps.size();
  Matrix subzeta;
  Vector subbeta, subgroups(1), subx, margin(dim), tv;
  subgroups[0] = 0.0;

  if (ngrps==1)
    {
      long int li, len = (1 << dim);
      double tx, *table = new double[len], *cumtable = new double[len],
	nconst;

      nconst = build_mvb_table(table, cumtable, betas, zetas);

      margin.zeroes();
      for (li=0 ; li<len ; ++li)
	{
	  tx = table[li];
	  tv = int_to_mvb(li, dim);
	  for (j=0 ; j<dim ; ++j)
	    if (tv[j]==1)
	      margin[j]+=tx;
	}
      margin /= nconst;

      delete[] table;   delete[] cumtable;
    }
  else
    for (i=0 ; i<ngrps ; ++i)
      {
	i1 = (int)(0.5+grps[i]);
	i2 = i<ngrps-1 ? (int)(0.5+grps[i+1]) : dim;
	subzeta = zetas.submatrix(i1,i1,i2,i2);
	subbeta = betas.subvector(i1,i2);
	tv = mvbern_marginal(subbeta,subzeta,subgroups);
	for (j=i1 ; j<i2 ; ++j)
	  margin[j] = tv[j-i1];
      }

  return margin;
}

}
