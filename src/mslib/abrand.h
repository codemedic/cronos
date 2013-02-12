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
 */


#ifndef ABRAND_HPP
#define ABRAND_HPP

#include "matrix.h"

namespace mslib {

void init_generator(long int sd);

// sampling routines
double unif_random();                           // Uniform on [0,1)
double norm_random();                           // N(0,1)
double norm_random(double mu, double sig2);     // N(mu,sig2)
int bino_random(int, double);                   // Binonial (n,p)
double exp_random(double lambda=1.0);          // Exponential mean 1/lambda
double t_random(int degreesoffreedom);       // t-distribution
double gam_random(double alpha, double beta);  // Gamma, mean=alpha/beta, var=alpha/beta^2
double chi2_random(double degreesoffreedom); // Chi-squared r.v.
int pois_random(double lambda);                 // Poisson r.v. with mean lambda
double beta_random(double a, double b);      // Beta, mean=a/(a+b)
double cauchy_random(double location=0);     // Cauchy + location
int geo_random(double p);                    // Geometric, mean = 1/p
Vector diri_random(Vector& nums);            // Dirichlet
Vector mvbern_random(Vector& betas, Matrix& zetas, Vector& groups); // Multivariate Bernoulli

// along with some mvn sampling routines
Vector mvn_random(Vector& mean, Matrix& sig2, int *error=NULL);
Vector mvn_with_sqrtvar(Vector& mean, Matrix& sqrtsig2);
Vector cond_mvn_random(Vector& mean, Matrix& var, 
		Matrix& H, Vector& rel, Vector& k,
		Matrix *suppliedP = NULL, Matrix *suppliedPinv = NULL,
		int fillinP = 0);

// miscellaneous routines
Matrix cov_msqrt(Matrix& covm, int reduce, int& error);


// density, cdf routines
double norm_density(double x, double mu=0, double sig2=1);
double lognorm_density(double x, double mu=0, double sig2=1);
double norm_cdf(double);        // (mu=0, sig2=1)
double norm_cdf(double x, double mu, double sig2);
double norm_invcdf(double);     // (mu=0, sig2=1)
double logpois_density(double x, double lambda);
double logpois_density_logl(double x, double loglambda);
double logmvn_density(Vector& x, Vector *mn, Matrix *Sig,
		     Matrix *SigInv = NULL, double Sigdet=0);
double exp_density(double x, double lambda=1.0);
double logexp_density(double x, double lambda=1.0);
double logdiri_density(Vector& x, Vector& parms);
double exp_cdf(double x, double lambda=1.0);
double exp_invcdf(double x, double lambda=1.0);
double logbeta_density(double x, double a, double b);
double logmvbern_pmf(Vector& x, Vector& betas, Matrix& zetas, Vector& groups);
Vector logmvbern_pmf(Matrix& x, Vector& betas, Matrix& zetas, Vector& groups);
Vector mvbern_marginal(Vector& betas, Matrix& zetas, Vector& groups);
Matrix mvbern_cov(Vector& betas, Matrix& zetas, Vector& groups); 

} // end of namespace

#endif
