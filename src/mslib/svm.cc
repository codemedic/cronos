/*
 * -------------------------------------------------------------------
 *
 * Copyright 2005 Anthony Brockwell
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
 * This module contains the source code for
 * stochastic volatility model analysis.
 *
 */

#include <iomanip.h>
#include <strstream>
#include <complex>
#include <string>
#include <functional>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include "matrix.h"
#include "abrand.h"
#include "poly.h"
#include "timeseries.h"
#include "svm.h"

const double  zero=1e-10;

SSVModel::SSVModel()
  // constructor
{
  estimation_mask.resize(4);  // it has 4 parameters
  estimation_mask[0] =        // default is to estimate all parms
  estimation_mask[1] = estimation_mask[2] = estimation_mask[3] = 1.0;
  mean = 0;
  phi = 0;
  nu = 1;
  mux = 0;
  loglikelihood_is_valid = false;
}

int SSVModel::FitModel(TimeSeries *ts, const int method, const int ilimit,
			void (*callback)(int,void *), void *cb_parms,
			ostringstream& msg, ostringstream& supplemental,
			bool get_parameter_cov)
{
  int i,j;
  Vector parms = ParameterBundle(), mask;
  Matrix parmcovs;
  double corr, tx;

  // create mask
  mask = GetMask();
  iterationlimit = ilimit;

  if (method==METHOD_MLE)
    {      
      parmcovs = MLEFit(ts, mask, callback, cb_parms, get_parameter_cov);
      StreamMLESummary(supplemental, parmcovs, get_parameter_cov);

      msg << "Fitted model computed using " << iterationcounter
	  << " likelihood evaluations." << ends;
      return SUCCESS;
    }
  if (method==METHOD_BAYESIAN)
    {

      supplemental << "Bayesian analysis of SSVMs not yet implemented." << ends;
      return UNABLE;
    }
}

void SSVModel::SimulateInto(TimeSeries *ts, int n, bool ovw)
{
  int i,nn=n;
  double tx;

  // If overwrite is true, we don't use n, and we just
  // use the existing time series of whatever length it is.
  // Otherwise, we clear it out and create a new ts of length n.
  if (!ovw)
    {
      ts->ClearOut();
      for (i=0 ; i<nn ; ++i)
        ts->Append(0);
    }
  else
    nn = ts->GetN();

  // 1. simulate from volatility process
  Vector logvols(nn);
  tx = nu*nu/(1-phi*phi);
  logvols[0] = norm_random(0,tx);
  for (i=1 ; i<nn ; ++i)
    logvols[i] = phi*logvols[i-1] + norm_random(0,nu*nu);

  // 2. compute observations (log-returns)
  for (i=0 ; i<nn ; ++i)
    (*ts)[i] = norm_random(mean, exp(logvols[i]+mux)); 
}

Matrix SSVModel::SimulatePredictive(TimeSeries& ts, int nsteps, int nsims)
{
  Matrix retval(nsims, nsteps);
  retval.zeroes();

  // step 1: use approximation to compute one-step predictive d-ns of volatilities
  int i,j,n=ts.GetN(),sim_num;
  double tx, tw, q;
  Vector volmeans(n),volstds(n);  // one-step predictive means and variances

  TimeSeries transformed;
  for (i=0 ; i<n ; ++i)
    {
      tx = log((ts[i]-mean)*(ts[i]-mean)+zero) - mux + 1.27;
      transformed.Append(tx, ts.IsMissing(i));
    }

  // this is y'_t = log[(y_t - mu)^2 + small] - mux + 1.27;     1.27 is mean of log(chi^2)
  // y_t' \simeq  w_t + N(0,4.94)
  // w_t = phi w_{t-1} + \epsilon_t,    \epsilon_t \sim iid(0,\nu^2)

  double curmean=0,curvar=nu*nu/(1-phi*phi),p,k;

  // now use Kalman recursions to get one-step predictive d-ns
  for (i=0 ; i<n ; ++i)
    {
      // curmean and curvar are filtering mean and variance from previous step (-1 to begin with)
      p = phi*phi*curvar + nu*nu;        // one step predictive state variance
      volstds[i] = sqrt(p);
      volmeans[i] = phi*curmean;
      k = p/(p+4.93);                    // 4.93 is var of log(chi^2)
      if (ts.IsMissing(i))
	{
	  curmean = phi*curmean;
	  curvar = p;
	}
      else
	{
	  curmean = phi*curmean + k*(transformed[i] - phi*curmean);
	  curvar = (1.0-k)*p;
	}
    }

  // now simulate beyond the end of the original time series, nsims times
  for (sim_num=0 ; sim_num<nsims ; ++sim_num)
    {
      tx = norm_random(curmean, curvar);  // log-vol

      for (i=0 ; i<nsteps ; ++i)
	{ 
	  // simulate the log-return itself
	  retval[sim_num][i] = norm_random(mean, exp(tx+mux)); 
	  // then update log-vol
	  tx = phi*tx + norm_random(0, nu*nu);
	}
    }

  return retval;
}

void SSVModel::RenderInto(ostringstream& os)
{
  os << "\\bStochastic Volatility Model\\n" << endl;
  os << "X(t) = \u03c3(t) Z(t) ";
  if (mean>=0)
    os << " + " << mean;
  else
    os << " - " << (-mean);
  os << endl << "\u03c3(t)\\s2\\n = exp[Y(t)";
  if (mux>=0)
    os << " + " << mux << "]" << endl;
  else
    os << " - " << (-mux) << "]" << endl;
  os << "Y(t) = " << phi << " Y(t-1) + \u03b5(t)" << endl;
  os << "{\u03b5(t)} ~ iidn(0," << (nu*nu) << ")" << endl;

  if (loglikelihood_is_valid)
    {
      os << endl << endl << "Log Likelihood = " << GetLogLikelihood() << endl;
      os << "AICC = " << GetAICC() << endl;
    }
  else
    {
      os << endl << endl << "Log Likelihood = <Not Computed>" << endl;
      os << "AICC = <Not Computed>" << endl;
    }
}

bool SSVModel::RenderCtsVersionInto(ostringstream &output, ostringstream &error, double delta)
{
  double alpha, beta, muprime, gamma;

  if (phi>0)
    {
      // translate parms (see cronos-scrawls1.pdf for details)
      alpha = log(phi)/delta;
      beta = mux/2 - log(delta)/2;
      gamma = sqrt(nu*nu/4 / ((exp(2*alpha*delta)-1)/(2*alpha)) );
      muprime = (mean + delta/2 * exp(2*beta + nu*nu/(2*(1-phi*phi)))) / delta;

      output << "\\bAssumptions:\\n" << endl;
      output << "X(t) is the log-return series defined by" << endl;
      output << "X(t) = log( S(t\u03b4\\n)/S((t-1)\u03b4) )" << endl;
      output << "\u03b4 = " << delta << endl << endl;

      output << "\\b(Approximate) Equivalent Continuous-Time Model:\\n" << endl;
      output << "dS(t) = \u03bc S(t) dt + \u03c4(t) S(t) dW\\d1\\n(t)" << endl;
      output << "\u03bc = " << muprime << endl;
      output << "\u03c4(t) = exp(V(t))" << endl;
      output << "dV(t) = " << alpha << "( V(t) ";
      if (beta<0)
	output << "+ " << (-beta) << " ) dt + ";
      else
	output << "- " << beta << " ) dt + ";
      output << gamma << " dW\\d2\\n(t)" << endl;
      output << ends;

      return true;
    }
  else
    {
      error << "The standard continuous-time representation is not valid when phi<=0." << ends;
      return false;
    }
}

void SSVModel::ComputeACFPACF(Vector &acf, Vector *pacf, bool normalizeacf)
{
  acf.zeroes();
  if (normalizeacf)
    acf[0] = 1.0;     
  else
    {
      double gamma = nu*nu/(1-phi*phi);
      acf[0] = exp(mux + gamma*gamma/2);
    }
  if (pacf!=NULL)
    {
      pacf->zeroes();
      (*pacf)[0] = 1.0;
    }
}


double SSVModel::Cdf(double y, double mu, double sig2)
{
  // returns cdf of RV with mean mean, and mean and var of log-vol = mu+mux and sig2
  int i,j,n;
  double u,tx,ty;

  for (n=0,u=0.01,ty=0.0 ; u<=0.99 ; u+=0.01)
    {
      tx = norm_invcdf(u) * sqrt(sig2) + mu;
      ++n;
      ty += norm_cdf(y, mean, exp(mux + tx));
    }
  ty /= n;
  return ty;
}


double SSVModel::InvCdf(double u, double mu, double sig2)
{
  double retval, tx = exp(mu + sig2/2);  // tx is expected volatility for mu, sig2
  bool done=false;

  // do a direct numerical optimization, holding mu and sig2 fixed, find y such that Cdf(y,mu,sig2)=u
  // use simple bisection algorithm
  tx = exp(mu + sig2/2);
  double upper = mean + 7*sqrt(tx), lower = mean - 7*sqrt(tx), mid;
  double cdfupper = Cdf(upper, mu, sig2), cdflower = Cdf(lower, mu, sig2), cdfmid;
  
  while (!done)
    {
      mid = (lower + upper)/2;
      cdfmid = Cdf(mid, mu, sig2);
      
      if (u > cdfmid)
	{
	  lower = mid;
	  cdflower = cdfmid;
	}
      else
	{
	  upper = mid;
	  cdfupper = cdfmid;
	}

      done = (fabs(upper-lower)<0.01*sqrt(tx));
    }

  retval = (upper + lower)/2;
  return retval;
}

void SSVModel::Forecast(TimeSeries *tsp, int nsteps, double *fret, double *fmse)
{
  // step 1: use approximation to compute one-step predictive d-ns of volatilities
  int i,j,n=tsp->GetN();
  double retval, tx, tw, q;
  Vector volmeans(n),volstds(n);  // one-step predictive means and variances

  TimeSeries transformed;
  for (i=0 ; i<n ; ++i)
    {
      tx = log(((*tsp)[i]-mean)*((*tsp)[i]-mean)+zero) - mux + 1.27;
      transformed.Append(tx, tsp->IsMissing(i));
    }

  // this is y'_t = log[(y_t - mu)^2 + small] - mux + 1.27;     1.27 is mean of log(chi^2)
  // y_t' \simeq  w_t + N(0,4.94)
  // w_t = phi w_{t-1} + \epsilon_t,    \epsilon_t \sim iid(0,\nu^2)

  double curmean=0,curvar=nu*nu/(1-phi*phi),p,k;

  // now use Kalman recursions to get one-step predictive d-ns
  for (i=0 ; i<n ; ++i)
    {
      // curmean and curvar are filtering mean and variance from previous step (-1 to begin with)
      p = phi*phi*curvar + nu*nu;        // one step predictive state variance
      volstds[i] = sqrt(p);
      volmeans[i] = phi*curmean;
      k = p/(p+4.93);                    // 4.93 is var of log(chi^2)
      if (tsp->IsMissing(i))
	{
	  curmean = phi*curmean;
	  curvar = p;
	}
      else
	{
	  curmean = phi*curmean + k*(transformed[i] - phi*curmean);
	  curvar = (1.0-k)*p;
	}
    }

  // now generate forecasts for the original time series
  for (i=0 ; i<nsteps ; ++i)
    {
      // generate 0.025 and 0.975 quantiles of predictor
      tx = InvCdf(0.025, curmean, curvar);  // pass current mean and variance of log-vol
      fret[i] = mean;
      fmse[i] = fabs(mean-tx)/1.96;
      fmse[i] = fmse[i]*fmse[i];
      
      // update predictors in AR process for log-volatility
      curmean = phi*curmean;
      curvar = phi*phi*curvar + nu*nu;
    }
}

void SSVModel::ComputeStdResiduals(TimeSeries *tsp, TimeSeries *resids)
{
  // step 1: use approximation to compute one-step predictive d-ns of volatilities
  int i,j,n=tsp->GetN();
  double retval, tx, tw, q;
  Vector volmeans(n),volstds(n);  // one-step predictive means and variances

  TimeSeries transformed;
  for (i=0 ; i<n ; ++i)
    {
      tx = log(((*tsp)[i]-mean)*((*tsp)[i]-mean)+zero) - mux + 1.27;
      transformed.Append(tx, tsp->IsMissing(i));
    }

  // this is y'_t = log[(y_t - mu)^2 + small] - mux + 1.27;     -1.27 is mean of log(chi^2)
  // y_t' \simeq  w_t + N(0,4.94)
  // w_t = phi w_{t-1} + \epsilon_t,    \epsilon_t \sim iid(0,\nu^2)

  double curmean=0,curvar=nu*nu/(1-phi*phi),p,k;

  // now use Kalman recursions to get one-step predictive d-ns
  for (i=0 ; i<n ; ++i)
    {
      // curmean and curvar are filtering mean and variance from previous step (-1 to begin with)
      p = phi*phi*curvar + nu*nu;        // one step predictive state variance
      volstds[i] = sqrt(p);
      volmeans[i] = phi*curmean;
      k = p/(p+4.93);                    // 4.93 is var of log(chi^2)
      if (tsp->IsMissing(i))
	{
	  curmean = phi*curmean;
	  curvar = p;
	}
      else
	{
	  curmean = phi*curmean + k*(transformed[i] - phi*curmean);
	  curvar = (1.0-k)*p;
	}
    }

  // now find residuals
  resids->ClearOut();
  for (i=0 ; i<n ; ++i)
    {
      if (tsp->IsMissing(i))
	tx = 0.0;
      else
	{
	  // generate 0.025 and 0.975 quantiles of predictor
	  tx = Cdf((*tsp)[i], volmeans[i], volstds[i]*volstds[i]);
	  tx = norm_invcdf(tx);  // transform to N(0,1)
	}
      resids->Append(tx, tsp->IsMissing(i));
    }
}

double SSVModel::ComputeLogLikelihood(TimeSeries *tsp)
{
  int i,j,n=tsp->GetN();
  double retval, tx, tw, q;
  Vector volmeans(n),volstds(n);  // one-step predictive means and variances

  if (n==0)
    return 0.0;

  // step 1: use approximation to compute one-step predictive d-ns of volatilities
  TimeSeries transformed;
  for (i=0 ; i<n ; ++i)
    {
      tx = log(((*tsp)[i]-mean)*((*tsp)[i]-mean)+zero) - mux + 1.27;
      transformed.Append(tx, tsp->IsMissing(i));
    }

  // this is y'_t = log[(y_t - mu)^2 + small] - mux + 1.27;     1.27 is mean of log(chi^2)
  // y_t' \simeq  w_t + N(0,4.94)
  // w_t = phi w_{t-1} + \epsilon_t,    \epsilon_t \sim iid(0,\nu^2)

  double curmean=0,curvar=nu*nu/(1-phi*phi),p,k;

  // now use Kalman recursions to get one-step predictive d-ns
  for (i=0 ; i<n ; ++i)
    {
      // curmean and curvar are filtering mean and variance from previous step (-1 to begin with)
      p = phi*phi*curvar + nu*nu;        // one step predictive state variance
      volstds[i] = sqrt(p);
      volmeans[i] = phi*curmean;
      k = p/(p+4.93);                    // 4.93 is var of log(chi^2)
      if (tsp->IsMissing(i))
	{
	  curmean = phi*curmean;
	  curvar = p;
	}
      else
	{
	  curmean = phi*curmean + k*(transformed[i] - phi*curmean);
	  curvar = (1.0-k)*p;
	}
    }

  // now work out one-step predictive densities
  for (i=0,retval = 0.0 ; i<n ; ++i)
    {
      if (!tsp->IsMissing(i))
	{
	  // approximate the integral of p(y_t | w_t) p(w_t | y_{1:t-1})
	  for (tx=0, q=0.01, j=0 ; q<=0.99 ; q+=0.02)
	    {
	      tw = norm_invcdf(q) * volstds[i] + volmeans[i];
	      tw = exp(tw + mux);
	      tx += norm_density((*tsp)[i], mean, tw);
	      ++j;
	    }
	  retval += log(tx/j);
	}
    }

  loglikelihood = retval;
  AICC = -2*loglikelihood + 8;  // 4 parameters
  loglikelihood_is_valid = true;

  //cout << "Log-like(" << mean << "," << mux << "," << phi << "," << nu << ") -> " << loglikelihood << endl;

  return loglikelihood;
}

bool SSVModel::CheckModel(bool flipit)
{
  bool fixed;
  Vector tv  = ValidityConstraint();
  double tx = tv.min();
  fixed = tx >= 0;
  if (!fixed) // is invalid model
    if (flipit)
      {
	// fix it temporarily!
	if (tv[0]<zero) // nu negative is easily fixed
	  {
	    nu = fabs(nu);
	    if (nu<zero) nu=zero;
	  }
	if (tv[1]<zero) // non-causal AR(1)
	  {
	    phi = (1.0-zero)*(1.0-zero) / phi;
	  }
	fixed = true;
      }
  return fixed;
}

Vector SSVModel::ParameterBundle()
{
  Vector parms(4);
  parms[0] = mean;
  parms[1] = mux;
  parms[2] = phi;
  parms[3] = nu;
  return parms;
}

void SSVModel::UnbundleParameters(Vector& v)
{
  mean = v[0];
  mux = v[1];
  phi = v[2];
  nu = v[3];
}

Vector SSVModel::ValidityConstraint()
  // this function returns all positive elements of a vector
  // if the parameters are OK
  // nu has to be positive and |phi|<1
{
  Vector retval(2);
  retval[0] = nu;                   // nu negative is a problem
  retval[1] = 1.0-fabs(phi)-zero;   // |phi|>1 is a problem -> returns negative value then
  return retval;
}

Vector SSVModel::GetMask()
  // returns estimation_mask, padded with defaults
  // the only role this function plays is the padding with defaults
{
  Vector retval;
  retval = estimation_mask;
  return retval;                 // never any padding since you can't change number of params with this model
}

Vector SSVModel::GetDefaultParameters(Vector& mask, TimeSeries *tsp)
  // returns initial guesses for non-held model parameters; other parameters are fixed at current values
{
  int i;
  Vector parms = ParameterBundle(), acf(1);
  double tx, tm = tsp->SampleMean();
  TimeSeries tempts;

  tsp->ComputeSampleACF(&acf, NULL, false);  // ACF is un-normalized cov of time series

  // mean is handled differently
  if (mask[0]) // mean is to be estimated
    {
      mask[0] = 0;
      parms[0] = tm;
    }

  // initialize all non-masked parameters now

  // start by estimating variance of volatility
  for (i=0 ; i<tsp->GetN() ; ++i)
    {
      tx = (*tsp)[i] - tm;
      tx = tx*tx;
      tempts.Append(tx);
    }
  Vector vacf(1);
  double vm;
  tempts.ComputeSampleACF(&vacf,NULL,false);  // now vacf[0] = variance of volatility
  vm = tempts.SampleMean();

  if (mask[1])
    parms[1] = log(vm);

  if (mask[2])
    parms[2] = 0.0;  // phi

  if (mask[3])
    parms[3] = 0.0;  // nu

  return parms;
}

void SSVModel::StreamParameterName(ostringstream& strm, int parmnum)
{
  switch (parmnum) {
  case 0:
    strm << "Mean";   break;
  case 1:
    strm << "Mean of log-volatility";   break;
  case 2:
    strm << "\u03c6";   break;   // greek mu
  case 3:
    strm << "\u03bd";   break;   // \nu of log-volatility process
  }
}

Vector SSVModel::ComputeSpectralDensity(Vector& omegas)
// SSVM is special case of WN->hence has flat spectrum
{
  int i,j;
  Vector retval(omegas.size()),acf(1);
  complex tc1,tc2,eil,tsum,thetapart,phipart;

  ComputeACFPACF(acf,NULL,false);  // get the variance
  for (i=0 ; i<omegas.size() ; ++i)
    retval[i] = acf[0]/(2*M_PI);

  return retval;
}
