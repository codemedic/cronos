/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004,2005,2006,2007 Anthony Brockwell
 *
 * This program is free software; you can redistribute it and/or
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
 * Non-dynamically linked windows for additional info,
 * including residual analysis windows, forecast windows,
 * spectral analysis windows, etc.
 * Nov, 2003, Anthony Brockwell.
 *
 */

#include "includes.h"
#include "textrender.h"
#include "tsmrender.h"
#include "otherpopups.h"
#include <gsl/gsl_cdf.h>

using namespace std;
using namespace Gtk;
using namespace SigC;

ResidualWin::ResidualWin(TimeSeries *resids, TimeSeriesModel *model, TimeSeries *tseries)
  : Table(2,2,true), residuals(*resids)
{
  set_name("MITSMGraphWidget");
  string s1 = "Residuals", s2 = "Residual ACF", s3 = "Residual PACF";

  tsp = manage(new TSPlot(*resids));
  tsp->SetTitle(s1);
  attach(*tsp, 0,1,0,1);

  //PatchPlot *ppp = manage(new PatchPlot(*resids));
  //attach(*ppp, 2,3,0,1);

  acfpp = manage(new ACFPlot(*resids));
  acfpp->SetTitle(s2);
  //  pacfpp = manage(new PACFPlot(*resids));
  //  pacfpp->SetTitle(s3);
  qqpp = manage(new QQPlot(*resids));
  //qqpp->SetTitle("QQ Plot (Normal)");
  attach(*acfpp, 0,1,1,2);
  attach(*qqpp, 1,2,1,2);

  tr = manage(new TextRenderer());
  RunTests();
  ostringstream os;
  os << "\\bIID test p-values:\\n" << endl;

  if (pvals[0]<0.05)
    os << "\\0";
  else
    os << "\\n";
  os << "  Ljung-Box (40 lags) = " << pvals[0] << endl;

  if (pvals[1]<0.05)
    os << "\\0";
  else
    os << "\\n";
  os << "  McLeod-Li (40 lags) = " << pvals[1] << endl;

  if (pvals[2]<0.05)
    os << "\\0";
  else
    os << "\\n";
  os << "  Turning-point       = " << pvals[2] << endl;

  os << endl;
  model->RenderInto(os);
  os << endl;
  tseries->RenderInto(os);
  os << endl << ends;
  tr->Update(os.str().c_str());
  ScrolledWindow *sw = manage(new ScrolledWindow());
  sw->set_policy(POLICY_AUTOMATIC, POLICY_AUTOMATIC);  
  sw->add(*tr);
  Frame *frame = manage(new Frame());
  frame->add(*sw);
  attach(*frame, 1,2,0,1, FILL, FILL, 10, 10);
}

ResidualWin::~ResidualWin()
{
  delete tsp;
  delete acfpp;
  delete pacfpp;
}

void ResidualWin::RunTests()
{
  int i;
  double T;

  // start by computing acf and pacf
  int n = residuals.GetN();
  int h = min(n-1,40);
  acf.resize(h+1);   pacf.resize(h+1);

  // now tests
  Ts.resize(3);
  pvals.resize(3);

  // Ljung-Box test
  residuals.ComputeSampleACF(&acf, &pacf);
  for (T=0.0,i=1 ; i<=h ; ++i)
    T += acf[i]*acf[i] / (n-i);
  T *= n*(n+2); // now T should be
  Ts[0] = T;
  pvals[0] = 1.0-gsl_cdf_chisq_P(T, h);

  // McLeod-Li test
  TimeSeries res2;
  res2 = residuals;
  for (i=0 ; i<n ; ++i)
    res2[i] = res2[i]*res2[i];   // residuals squared
  res2.ComputeSampleACF(&acf, &pacf);
  for (T=0.0,i=1 ; i<=h ; ++i)
    T += acf[i]*acf[i] / (n-i);
  T *= n*(n+2); // now T should be
  Ts[1] = T;
  pvals[1] = 1.0-gsl_cdf_chisq_P(T, h);

  // Turning-Point test
  for (i=1,T=0.0 ; i<n-1 ; ++i)
    if (((residuals[i] > residuals[i-1]) && (residuals[i] > residuals[i+1]))
	|| ((residuals[i] < residuals[i-1]) && (residuals[i] < residuals[i+1])))
	++T;
  Ts[2] = T;
  pvals[2] = gsl_cdf_gaussian_P(fabs(Ts[2]-(2.0*n/3.0)),sqrt((16.0*n-29.0)/90.0));
  pvals[2] = 2*(1.0 - pvals[2]);  // want right tail, and two-sided test
}

ModelACFWin::ModelACFWin(TimeSeriesModel *m)
  : Table(2,2,true)
{
  string s1 = "Model ACF", s2 = "Model PACF";
  acfpp = manage(new ACFPlot(*m));
  acfpp->SetTitle(s1);
  pacfpp = manage(new PACFPlot(*m));
  pacfpp->SetTitle(s2);
  TSMRenderer *mr = manage(new TSMRenderer(m,NULL,NULL));
  //mr->SetTitle("Blah");

  attach(*acfpp, 0, 1, 1, 2);
  attach(*pacfpp, 1, 2, 1, 2);
  attach(*mr, 0, 2, 0, 1);
}

void ForecastWin::RebuildData()
{
  int i;

  // create matrices for xyplot
  int len = local_ts_copy.GetN(), nr = 4, tlen = len+nf;
  x.resize(tlen,nr);
  y.resize(tlen,nr);
  missing.resize(tlen,nr);

  // The idea here is two superimpose two plots: the data and the forecasts,
  // where data are missing at times > n and forecasts are missing at times <= n.
  // first col of x,y,missing is the data
  // second col of x,y,missing is the forecasts
  
  Vector pmeans = all_sims.mean(), range_lower(nf), range_upper(nf);

  for (i=0 ; i<len+nf ; ++i)
    x[i][0] = x[i][1] = i+1;

  for (i=0 ; i<len ; ++i)
    {
      y[i][0] = local_ts_copy[i];
      y[i][1] = 0;
      missing[i][0] = local_ts_copy.IsMissing(i); 
      missing[i][1] = 1;
    }
  for (i=0 ; i<nf ; ++i)
    {
      y[i+len][0] = 0;
      y[i+len][1] = pmeans[i];
      missing[i+len][0] = 1;
      missing[i+len][1] = 0;
    }

  // add two more series to plot: forecs + 1.96 mse and forecs - 1.96 mse
  for (i=0 ; i<len ; ++i)
    {
      x[i][2] = i+1;   x[i][3] = i+1;
      y[i][2] = 0;     y[i][3] = 0;
      missing[i][2] = 1;
      missing[i][3] = 1;
    }
  range_lower = all_sims.quantiles(0.025, true);
  range_upper = all_sims.quantiles(0.975, true);
  for (i=0 ; i<nf ; ++i)
    {
      x[i+len][2] = x[i+len][3] = i+len+1;
      y[i+len][2] = range_lower[i];
      y[i+len][3] = range_upper[i];
      missing[i+len][2] = 0;
      missing[i+len][3] = 0;
    }
}

ForecastWin::ForecastWin(TimeSeries& ts, TimeSeriesModel& model, int horizon,
			 vector<Transformation*> *tlistcopy, int num_sims)
  : Table(1,1,true), trans(tlistcopy), nf(horizon)
{
  int i, j;
  set_name("MITSMGraphWidget");
  bool show_mse = true;
  Vector tv;


  // perform simulation to get multiple paths into the future
  all_sims = model.SimulatePredictive(ts, horizon, num_sims);

  ntrans = trans->size();
  transpos = ntrans;     // can be from 0 ... ntrans

  local_ts_copy = ts;

  RebuildData();  // create the x,y and missing matrices and fill them in

  VBox *vbox = manage(new VBox());

  xyplot = manage(new XYPlot(&x,&y,&missing));
  xyplot->SetIntegerxaxis(true);
  string t = ts.GetTitle() + " Forecasts";
  xyplot->SetTitle(t);
  // fix color_map too
  tv.resize(8);
  tv[0]=0;
  tv[1]=1;
  tv[2]=tv[3]=2;
  xyplot->SetColorMap(tv);

  HBox *hbox = manage(new HBox(false, 2));
  Button *b1 = manage(new Button(Gtk::Stock::GO_BACK)),
    *b2 = manage(new Button(Gtk::Stock::GO_FORWARD));
  b1->signal_clicked().connect(mem_fun(*this, &ForecastWin::UndoTransform));
  b2->signal_clicked().connect(mem_fun(*this, &ForecastWin::RedoTransform));
  hbox->pack_start(*b1, false, false);
  hbox->pack_start(*b2, false, false);
  Frame *frame = manage(new Frame());
  frame->add(*hbox);
  frame->set_label("Transformations");

  vbox->pack_start(*frame, false, false);
  vbox->pack_start(*xyplot, true, true);

  attach(*vbox,0,1,0,1);
}

ForecastWin::~ForecastWin()
{
  // clean up transformlist
  Transformation *tp;
  int i,n = trans->size();
  // iterate through and delete
  for (i=0 ; i<n ; ++i)
    {
      tp = (*trans)[i];
      delete tp;
    }
  // clear it out!
  trans->clear();
}

void ForecastWin::UndoTransform()
{
  int i,j,k, nobs = x.nrows()-nf, num_to_do = all_sims.nrows(), new_len;
  Matrix nx(x), ny(y), nmissing(missing), new_sims;
  ostrstream errst;
  bool success;

  if (transpos==0)
    return;  // can't go back before first transform!

  TimeSeries tempts;

  // we either need to do reverse transformation on each predictive simulation
  new_sims = all_sims;
  success = true;
  for (i=0 ; (i<num_to_do) && success ; ++i)
    {
      tempts = local_ts_copy;

      for (j=0 ; j<nf ; ++j)
	tempts.Append(all_sims[i][j]);

      success &= (*trans)[transpos-1]->Reverse(&tempts, errst, false);
      new_len = tempts.GetN();
      
      for (j=0 ; j<nf ; ++j)
        new_sims[i][j] = tempts[new_len-nf+j];
    }
  if (success)
    {
      all_sims = new_sims;

      // we also need to update local_ts_copy
      local_ts_copy.ClearOut();
      for (i=0 ; i<new_len-nf ; ++i)
	local_ts_copy.Append(tempts[i], tempts.IsMissing(i));
    }
  else
    {
      errst << endl << endl << "Undo transform aborted." << ends;
      char *temps = errst.str();
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window*) wp;
      MessageDialog ed(*winp, temps);
      ed.run();
      return;
    }

  // unbundle into new x,y,missing matrices
  RebuildData();

  xyplot->Update(&x,&y,&missing);
  --transpos;
}

void ForecastWin::RedoTransform()
{
  int i,j,k,num_to_do = all_sims.nrows(), nobs = x.nrows()-nf, new_len;
  Matrix nx(x),ny(y),nmissing(missing);
  ostrstream errst;

  if (transpos==ntrans)
    return;  // can't go forward

  TimeSeries tempts;

  // we need to do (local_ts_copy + all_sims[i][:]) for each i=1,...,all_sims.nrows()
  for (i=0 ; i<num_to_do ; ++i)
    {
      tempts = local_ts_copy;
      for (j=0 ; j<nf ; ++j)
	tempts.Append(all_sims[i][j]);

      (*trans)[transpos]->Apply(&tempts, errst, false);

      new_len = tempts.GetN();
      for (j=0 ; j<nf ; ++j)
	all_sims[i][j] = tempts[new_len-nf+j];
    }

  // we also need to update local_ts_copy
  local_ts_copy.ClearOut();
  for (i=0 ; i<new_len-nf ; ++i)
    local_ts_copy.Append(tempts[i], tempts.IsMissing(i));

  RebuildData();

  xyplot->Update(&x,&y,&missing);
  ++transpos;
}

SpectrumWin::SpectrumWin(Matrix& f, Matrix& p, string& ttl)
  : Table(1,1,true), freqs(f), pdgm(p)
{
  int i;
  set_name("MITSMGraphWidget");
  logged = false;

  VBox *vbox = manage(new VBox());

  xyplot = manage(new XYPlot(&f, &p));
  string t = ttl;
  xyplot->SetTitle(t);
  xyplot->SetBasicUnits(M_PI,1.0);
  xyplot->SetMargins(0,0.05);

  HBox *hbox = manage(new HBox(false, 2));
  Button *b1 = manage(new Button("Log")),
    *b2 = manage(new Button("No Log"));
  b1->signal_clicked().connect(mem_fun(*this, &SpectrumWin::LogTransform));
  b2->signal_clicked().connect(mem_fun(*this, &SpectrumWin::ExpTransform));
  hbox->pack_start(*b1, false, false);
  hbox->pack_start(*b2, false, false);

  vbox->pack_start(*hbox, false, false);
  vbox->pack_start(*xyplot, true, true);

  attach(*vbox,0,1,0,1);
}

SpectrumWin::~SpectrumWin()
{
}

void SpectrumWin::LogTransform()
{
  int i, j, n=freqs.nrows(),nc=freqs.ncols();
  if (logged)
    return;
  for (i=0 ; i<n ; ++i)
    for (j=0 ; j<nc ; ++j)
      pdgm[i][j] = log(pdgm[i][j]);
  logged = true;
  xyplot->Update(&freqs, &pdgm, NULL, true);
}

void SpectrumWin::ExpTransform()
{
  int i, j, n=freqs.nrows(),nc=freqs.ncols();
  if (!logged)
    return;
  for (i=0 ; i<n ; ++i)
    for (j=0 ; j<nc ; ++j)
      pdgm[i][j] = exp(pdgm[i][j]);
  logged = false;
  xyplot->Update(&freqs, &pdgm, NULL, true);
}

