/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004 Anthony Brockwell
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
 */

#ifndef RESIDWIN_H
#define RESIDWIN_H

class TextRenderer;

class ModelACFWin : public Gtk::Table {
 protected:
  ACFPlot *acfpp;
  PACFPlot *pacfpp;

 public:
  ModelACFWin(TimeSeriesModel *m);
};

class ResidualWin : public Gtk::Table {
 protected:
  TSPlot *tsp;
  ACFPlot *acfpp;
  PACFPlot *pacfpp;
  TextRenderer *tr;
  TimeSeries residuals;
  Vector acf, pacf;
  QQPlot *qqpp;

  void RunTests();
  Vector Ts, pvals;

public:
  ResidualWin(TimeSeries *resids, TimeSeriesModel *model, TimeSeries *ts);
  ~ResidualWin();
};

class SpectrumWin : public Gtk::Table {
 protected:
  Matrix freqs;    // Nx1 vector or Nx2 vector
  Matrix pdgm;     // often just a vector, but if smoothed as well, is Nx2
  XYPlot *xyplot;
  bool logged;

 public:
  SpectrumWin(Matrix& f, Matrix& p, string& title);
  ~SpectrumWin();
  void LogTransform();
  void ExpTransform();
};

class ForecastWin : public Gtk::Table {
 protected:
  vector<Transformation*> *trans;
  int ntrans,transpos;
  Matrix x,y,missing;   // data fed to xyplot class
  int nf;
  XYPlot *xyplot;
  Matrix all_sims;
  void RebuildData();
  TimeSeries local_ts_copy;

 public:
  ForecastWin(TimeSeries&, TimeSeriesModel&, int horizon, 
	      vector<Transformation*> *tlistcopy, int num_sims);
  ~ForecastWin();
  void UndoTransform();
  void RedoTransform();
};

#endif
