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
 */

#ifndef SVM_H
#define SVM_H

// simple stochastic volatility model first
//   Y_t \sim N(\mu_Y,\sigma_t^2)
//   \sigma_t = \exp(X_t + \mu_X)
//   X_t = \phi X_{t-1} + \epsilon_t, {\epsilon_t} \sim IIDN(0,\nu^2)

class SSVModel : public TimeSeriesModel {
  protected:
    double mean, mux, phi, nu;


  public:
    SSVModel();    // constructor

    int FitModel(TimeSeries *ts, const int method, const int ilimit,
		 void (*itercallback)(int,void *), void *cb_parms,
		 ostringstream& msg, ostringstream& supplementalinfo,
		 bool get_parameter_cov=true);
    void SimulateInto(TimeSeries *, int n, bool ovw);
    Matrix SimulatePredictive(TimeSeries&, int nsteps, int nsims);
    void RenderInto(ostringstream&);
    bool RenderCtsVersionInto(ostringstream &output, ostringstream &error, double delta);
    void ComputeACFPACF(Vector &acf, Vector *pacf, bool normalizeacf=true);
    void Forecast(TimeSeries *, int nsteps, double *fret, double *fmse);
    void ComputeStdResiduals(TimeSeries *data, TimeSeries *resids);
    double ComputeLogLikelihood(TimeSeries *);
    bool CheckModel(bool flipit);
    Vector ParameterBundle();
    void UnbundleParameters(Vector& v);
    Vector ValidityConstraint();
    Vector GetMask();                         // returns estimation_mask, padded with defaults
    Vector GetDefaultParameters(Vector& mask, TimeSeries *tsp);     // returns initial guesses for non-held model parameters; other parameters are fixed at current values
    void StreamParameterName(ostringstream& strm, int parmnum);        // nms[0],... must already be valid allocated strings with length >= 20
    Vector ComputeSpectralDensity(Vector& omegas);

    double Cdf(double y, double mu, double sig2);
    double InvCdf(double u, double mu, double sig2);
};

#endif
