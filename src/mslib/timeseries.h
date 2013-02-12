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

/*
** July 2006 update: adapted TimeSeries class to handle
** multivariate time series (but not yet TimeSeriesModel).
** is backwards-compatible, i.e. ts[i] returns i'th element of
** first component TS if multivariate
*/

#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <string>
#include <strstream>
#include "matrix.h"
#include <gsl/gsl_vector.h>

using namespace std;
using namespace mslib;

const int maxorder = 100;
const int maxtransforms = 50;
const int METHOD_MLE = 1;
const int METHOD_BAYESIAN = 2;


class TimeSeriesModel;
class ARMAModel;

// The first class is simply a data structure containing
// observed values of a time series.

class TimeSeries {
protected:
  void ParseDescriptor(char *);
  string title, description;

public:
  Matrix data, missing;  // data contains data, each component TS in a row, missing contains 1/0s
  bool anymissing;
  int n, nallocated, n_components, current_component;
  
public:
  TimeSeries();
  ~TimeSeries();
  
  void ClearOut(int dimension=0);              // clears it out and sets dimension (0 means leave alone)
  void Append(double x, bool msng=false);
  void Append(Vector& x, Vector& msng);
  int GetN() {return n;}
  int GetDim() {return n_components;}
  bool IsMissing(int t) {return missing[current_component][t]; }
  bool IsMissing(int component, int t) {return missing[component][t];}
  bool HasAnyMissing();
  double *GetData() {return &data[current_component][0];}
  Vector GetDataVector();
  Matrix GetDataMatrix();
  Vector GetMissingVector();
  Matrix GetMissingMatrix();
  void SetMissingVector(Vector&);
  void SetMissing(int t, bool msng);
  void Clip(int n0, int n1);
  void RenderInto(ostringstream &);
  string& GetDescription() {return description;}
  string& GetTitle() {return title;}
  void SetTitle(string& t) {title=t;}
  double& operator[](int index) {return data[current_component][index];}
  
  // useful stuff
  double SampleMean();
  double SampleMin();
  void ComputeSampleACF(Vector *acf, Vector *pacf,
			bool normalizeacf=true);
  Matrix *GetInnovations(double *vhat, int m);
  Vector ComputePeriodogram(Vector *omegas = NULL);

  TimeSeries& operator=(const TimeSeries &other);
  bool operator==(const TimeSeries &other);

  friend ostream& operator<<(ostream&, TimeSeries&);
  friend istream& operator>>(istream&, TimeSeries&);
};


// The following class defines a model.

class TimeSeriesModel {
  protected:
    double loglikelihood, AICC;
    bool loglikelihood_is_valid;
    Vector estimation_mask;            // mask used in estimation (1=estimate, 0=hold fixed)
    int iterationcounter, iterationlimit, error_code;

  public:
    TimeSeriesModel();
    virtual ~TimeSeriesModel();

    Matrix MLEFit(TimeSeries *ts, Vector& mask, 
		  void (*callback)(int,void *) = NULL, void *cb_parms = NULL, 
		  bool get_parameter_cov=true);    // returns covariance 
                                                   // matrix of parameter estimates

    //---------------------------------------------------------------------+
    // Here are the functions which MUST be overridden in derived classes. |
    //---------------------------------------------------------------------+

    virtual int FitModel(TimeSeries *ts, const int method, const int numits,  // typically this f-n calls MLEFit
			 void (*itercallback)(int,void *), void *cb_parms,
			 ostringstream& msg, ostringstream& supplementalinfo,
			 bool get_parameter_cov=true)=0;

    virtual void SimulateInto(TimeSeries *, int n, bool ovw)=0;               // generates a simulated realization
    virtual void RenderInto(ostringstream &output)=0;                         // streams a description of the model
                                                                              // ends i NOT appended
    virtual bool RenderCtsVersionInto(ostringstream &output, ostringstream &error, double delta)=0;
                                                                              // constructs cts-time representation
                                                                              // and streams it
    virtual void ComputeACFPACF(Vector &acf, Vector *pacf,                    // autocorrelation/covariance+PACF
      bool normalizeacf=true)=0;
    virtual void Forecast(TimeSeries *, int nsteps, double *fret, double *fmse)=0;  // fill in predictive means and MSEs
 
    virtual Matrix SimulatePredictive(TimeSeries&, int nsteps, int nsims)=0;             // simulate nsteps into the future,
                                                                              // returning future values in a vector 
    virtual void ComputeStdResiduals(TimeSeries *data, TimeSeries *resids)=0; // if model is correct, these should be
                                                                              // realizations of iid std. normals
    virtual double ComputeLogLikelihood(TimeSeries *)=0;                      // compute log-likelihood with curr. parms
    virtual Vector ParameterBundle() = 0;                                     // returns all parms bundled into vector
    virtual void UnbundleParameters(Vector& v) = 0;                           // takes parms from vector
                                                                              // and copies into model
    virtual Vector ValidityConstraint()=0;                                    // checks validity of model parms:
                                                                              // all vec. components must be >= 0
    virtual Vector GetMask()=0;                                               // returns estimation_mask, 
                                                                              // padded with defaults (1=estimate)
    virtual Vector GetDefaultParameters(Vector& mask, TimeSeries *tsp)=0;     // returns initial guesses for 
                                                                              // non-held model parameters; 
                                                                              // other parameters are fixed at current values
    virtual void StreamParameterName(ostringstream& strm, int parmnum)=0;     // text descrip. of parm (no ends!)
    virtual Vector ComputeSpectralDensity(Vector& omegas) = 0;                // model spec. dens. eval'd at omegas


    // These remaining functions are defined for TimeSeriesModel, and do not necessarily need to be overridden
    void SetMask(Vector& m) {estimation_mask=m;}
    virtual bool CheckModel(bool flipit=false)=0;
    static Vector param_mask(Vector& parms, Vector& mask);
    static Vector param_expand(const gsl_vector *x, Vector& curparms, Vector& mask);  // useful internal functions
    void StreamMLESummary(ostringstream&, Matrix& parmcovs, bool get_parameter_cov);

    double GetLogLikelihood() {return loglikelihood;}
    double GetAICC() {return AICC;}

    static const int SUCCESS=1, NONCONVERGENT=0, UNABLE=-1;
};

// Here we define ARMAModel as a special case of the more general TimeSeriesModel.

class ARMAModel : public TimeSeriesModel {
  protected:
    int p,q;
    double sigma;
    double mean;
    double kappa(int i, int j, int m, Vector &acfs);	// used for prediction
    double fracdiff;
    double phi[maxorder], theta[maxorder];
    bool using_whittle;

    // cached info for whittle likelihood
    bool whittle_cached;
    TimeSeries whittle_ts;
    Vector whittle_pdgm;

    // some private functions
     double LogPrior(bool& isvalid);
    void BayesianFit(TimeSeries *ts, Vector& mask, 
		     void (*callback)(int,void *) = NULL, void *cb_parms = NULL);
    bool is_short_mem();

    void ComputeSpecialResiduals(TimeSeries *tp, TimeSeries *rret, double *rstore);
    void ComputeSpecialResidualsWithMissing(TimeSeries *tp, TimeSeries *rret, double *rstore);
    void ForecastWithMissing(TimeSeries *, int nsteps, double *fret, double *fmse);
    double ComputeReducedLogLikelihood(TimeSeries *tp, double *bestsig=NULL);
    Matrix cfuncsfor(double d, double rho_real, double rho_complex,
		     int p, int h, int extent);
 
    Vector ValidityConstraint();

  public:
    ARMAModel();
    ~ARMAModel();
    ARMAModel& operator=(const ARMAModel &other);


    int FitModel(TimeSeries *, const int method, const int numits,
		 void (*itercallback)(int,void *), void *cb_parms, 
		 ostringstream& msg, ostringstream& supplementalinfo, 
		 bool get_parameter_cov);
    void SimulateInto(TimeSeries *, int, bool);
    Matrix SimulatePredictive(TimeSeries&, int nsteps, int nsims);
    void Forecast(TimeSeries *, int nsteps, double *fret, double *fmse);
    void ComputeStdResiduals(TimeSeries *data, TimeSeries *resids);
    void RenderInto(ostringstream&);
    bool RenderCtsVersionInto(ostringstream &output, ostringstream &error, double delta);
    void ComputeACFPACF(Vector& acf, Vector* pacf, bool normalizeacf);
    double ComputeLogLikelihood(TimeSeries *tp);

    Vector ParameterBundle();
    void UnbundleParameters(Vector& v);
    Vector GetMask();
    Vector GetDefaultParameters(Vector& mask, TimeSeries *tsp);
    void StreamParameterName(ostringstream& os, int pnum);
    Vector ComputeSpectralDensity(Vector& omegas);

    // functions for getting/setting model
    void SetOrder(int arorder, int maorder)
      { p=arorder;   q=maorder; }
    void SetSigma2(double sigma2)
      { sigma = sqrt(sigma2); }
    void SetCoefficients(double *phis, double *thetas);
    void SetMean(double m)
      { mean=m; }
    void SetFracDiff(double d)
      { fracdiff=d; }
    int GetP() {return p;}
    int GetQ() {return q;}
    double GetFracDiff() {return fracdiff;}
    double GetMean() {return mean;}
    double GetSigma() {return sigma;}
    void GetCoefficients(Vector& phi, Vector& theta);
    bool IsCausal(bool causalify=false);
    bool IsInvertible(bool invertibilify=false);
    void SetWhittleMode(bool wm) {using_whittle=wm;}
    bool GetWhittleMode() {return using_whittle;}

    bool CheckModel(bool flipit);
    double MinAbsRoot();
};


class JulianTime {
  protected:
    int julianday;
    double dayfraction;

  public:
    JulianTime() {julianday=0;}
    ~JulianTime() {};
    void GetYearMonthDay(int *y, int *m, int *d);
    void SetYearMonthDay(int y, int m, int d);
    void IncreaseBy(double amount);
};

// here's a generic function used in MLE
double minus_log_like_for(const gsl_vector *x, void *min_info_container);

class MinimizerInfoContainer {
public:
  TimeSeriesModel *tsmp;
  TimeSeries *tsp;
  Vector themask;
};

#endif
