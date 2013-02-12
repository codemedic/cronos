/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004,2005 Anthony Brockwell
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


#ifndef XYPLOT_HPP
#define XYPLOT_HPP

#include <gdkmm/colormap.h>
#include <gdkmm/color.h>
#include <gdkmm/gc.h>
#include "dialogs.h"

using namespace std;
using namespace mslib;

const int GMODE_XY = 1, GMODE_ACF = 2;

class XYPlot : public Gtk::DrawingArea
{
protected:
  // local copies of data
  Matrix x,y,missing;
  Vector minx, maxx, miny, maxy, color_map;
  int nobs, actualnobs;  // actualnobs used for ACF/PACF bounds

  string title;
  Glib::RefPtr<Gdk::GC> gc;
  Gdk::Color blue, lightblue, red, green, black, white, grey, darkgrey;
  Gdk::Color *ColorNumber(int i);

  // some properties
  int integerxaxis, drawinglines;
  double basic_xunit, basic_yunit;
  double margin_fraction_x, margin_fraction_y;
  bool common_x_axis, common_y_axis;

public:
  XYPlot(Matrix *x, Matrix *y, Matrix *missing = NULL);
  ~XYPlot();
  void Update(Matrix *x, Matrix *y, Matrix *missing = NULL, bool redraw=true);
  void SetTitle(string&);
  void SetBasicUnits(double x, double y)
  { basic_xunit = x;   basic_yunit = y; }
  void SetMargins(double mfx, double mfy)
  { margin_fraction_x = mfx;  margin_fraction_y = mfy; }
  void SetMultiMode(bool common_x, bool common_y, bool redraw=true);
  void SetIntegerxaxis(bool val) {integerxaxis=val;}
  void SetColorMap(Vector& c) {color_map=c;}

  int mode;

  virtual void Oop1_Pressed(), Oop2_Pressed(), Oop3_Pressed(),
    Oop4_Pressed(), Oop5_Pressed();

  virtual void InitPropDialog(CenteringDialog& cd);
  virtual void ClosePropDialog(CenteringDialog& cd);

protected:
  double GetNiceUnitFor(double x1, double x2, double basic_unit=1.0);
  void LittleBox(Glib::RefPtr<Gdk::Window> win, Glib::RefPtr<Gdk::GC> gc, 
		 int x1, int y1, int r, Gdk::Color &fillcolor);

  virtual void on_realize();
  virtual bool on_expose_event(GdkEventExpose* e);
  virtual bool on_button_press_event(GdkEventButton* event);
  virtual void on_button_copy();

  void RedrawIt();

  void eps_draw(string &fnm);
};

#endif
