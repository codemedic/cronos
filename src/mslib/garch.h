/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004, 2005, 2006 Anthony Brockwell
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

#ifndef GARCH_H
#define GARCH_H

#include "timeseries.h"

using namespace std;
using namespace mslib;


const double root_margin = 0.001;


class GARCHModel : public TimeSeriesModel {
  protected:
    int p,q;
    double mean;
    Vector as, bs;

    virtual double OneStepVariance(double *xt, double *st, int t, double gamma0); 
    // returns one-step predictive volatility, t starts at 0, goes up to T-1.
 
  public:
    GARCHModel();
    ~GARCHModel();

    int FitModel(TimeSeries *ts, const int method, int numits,
		 void (*itercallback)(int,void *), void *cb_parms, 
		 ostringstream& msg, ostringstream& supplementalinfo, 
		 bool get_parameter_cov);
    void BayesianFit(int p, int q, TimeSeries *ts, 
		     void (*callback)(int,void *) = NULL, void *cb_parms = NULL);
    void SimulateInto(TimeSeries *, int n, bool ovw);
    Matrix SimulatePredictive(TimeSeries &ts, int nsteps, int nsims);
    void RenderInto(ostringstream&);
    bool RenderCtsVersionInto(ostringstream &output, ostringstream &error, double delta);
    double GetVar();
    void ComputeACFPACF(Vector &acf, Vector *pacf, bool normalizeacf=true);
    void Forecast(TimeSeries *, int nsteps, double *fret, double *fmse);
    void ComputeStdResiduals(TimeSeries *data, TimeSeries *resids);
    double LogPrior(bool& isok);
    double ComputeLogLikelihood(TimeSeries *);

    Vector ParameterBundle();
    void UnbundleParameters(Vector& v);

    void SetParameters(int p, int q, Vector& as, Vector& bs, double mu);
    void GetParameters(int& p, int& q, Vector& as, Vector& bs, double& mu);
    Vector GetMask();          // returns estimation_mask, padded with defaults
    Vector GetDefaultParameters(Vector& mask, TimeSeries *tsp);     // returns initial guesses
    void StreamParameterName(ostringstream& os, int pnum);
    Vector ComputeSpectralDensity(Vector& omegas);


    bool IsCausal(bool causalify=false);
    bool CheckModel(bool flipit);
    Vector ValidityConstraint();
};


class EGARCHModel : public GARCHModel
{
 protected:
  Vector gs;             // has as and bs inherited from GARCH, gs are new
  
  double OneStepVariance(double *xt, double *st, int t, double gamma0); 
  // returns one-step predictive volatility, t starts at 0, goes up to T-1.

  void StabilizeBs();    // makes sure b polynomial is well-behaved

  
 public:
  EGARCHModel();
  ~EGARCHModel();

  int FitModel(TimeSeries *ts, const int method, int numits,
	       void (*itercallback)(int,void *), void *cb_parms, 
	       ostringstream& msg, ostringstream& supplementalinfo, 
	       bool get_parameter_cov);

  void ComputeACFPACF(Vector &acf, Vector *pacf, bool normalizeacf=true);
  void RenderInto(ostringstream&);
  bool CheckModel(bool flipit);
  Vector ValidityConstraint();


  Vector ParameterBundle();
  void UnbundleParameters(Vector& v);
  
  void SetParameters(int p, int q, Vector& as, Vector& gs, Vector& bs, 
		     double mu);
  void GetParameters(int& p, int& q, Vector& as, Vector& gs, Vector& bs, 
		     double& mu);
  Vector GetMask();          // returns estimation_mask, padded with defaults
  Vector GetDefaultParameters(Vector& mask, TimeSeries *tsp);     // returns initial guesses
  void StreamParameterName(ostringstream& os, int pnum);
};


#endif
