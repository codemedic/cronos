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

#ifndef OTHERPLOTS_HPP
#define OTHERPLOTS_HPP

#include "xyplot.h"
#include "gtkmm.h"

using namespace Gtk;

class TSPlot : public XYPlot
{
protected:

public:
  TSPlot(TimeSeries& ts);
  void Update(TimeSeries&);
  void Update(TimeSeries& ts1, TimeSeries& ts2);

protected:
};

class QQPlot : public XYPlot
{
 protected:
 public:
  QQPlot(TimeSeries& ts);
  void Update(TimeSeries&);
};

class ACFPlot : public XYPlot
{
protected:
  Vector localacf, locallag;
  double gamma0;
  Adjustment *lagadjust;
  SpinButton *lagspin;
  CheckButton *cb1,*cb2;
  int nlags;
  bool normalizing, includinglag0;

public:
  ACFPlot(TimeSeries& ts);
  ACFPlot(TimeSeriesModel& tm);
  ~ACFPlot();
  void Update(TimeSeries&);
  void Update(TimeSeriesModel&);
  void RenderStuff();

  virtual void InitPropDialog(CenteringDialog& cd);
  virtual void ClosePropDialog(CenteringDialog& cd);
};

class PACFPlot : public XYPlot
{
protected:
  Vector localacf, localpacf, locallag;
public:
  PACFPlot(TimeSeries& ts);
  PACFPlot(TimeSeriesModel& tm);
  ~PACFPlot();
  void Update(TimeSeries&);
  void Update(TimeSeriesModel&);
};

class PatchPlot : public Gtk::DrawingArea
{
 protected:
  string title;
  Glib::RefPtr<Gdk::GC> gc;
  Gdk::Color blue, lightblue, red, green, black, white, grey, darkgrey;
  Gdk::Color *shades;

  int n, nshades;
  Vector *u;
  Matrix pseudomap, baselinemap;

 public:
  PatchPlot(TimeSeries& resids);
  ~PatchPlot();
  void Update();
  void SetTitle(string&);  

  virtual void on_realize();
  virtual bool on_expose_event(GdkEventExpose* e);
  void ComputePseudomap(int nr, int nc);
  double  Intensity(double x, double y, double centerx, double centery, double sig);
};

#endif
