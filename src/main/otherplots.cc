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
 *
 *
 * This module contains classes derived from XYPlot,
 * including TSPlot (time series plot),
 * ACFPlot, PACFPlot (obvious)
 *
 */


#include "includes.h"
#include "rect.h"

using namespace Gtk;
using namespace std;

const int default_maxlag = 50, max_possible_lag = 2000;

//------- TSPlot ---------

TSPlot::TSPlot(TimeSeries& ts)
  : XYPlot(NULL, &ts.GetDataVector(), &ts.GetMissingVector())
{
  integerxaxis = true;
}

void TSPlot::Update(TimeSeries& ts)
{
  Matrix y,m;
  y = ts.GetDataMatrix();
  m = ts.GetMissingMatrix();
  string s1 = ts.GetTitle();
  SetTitle(s1);
  SetMultiMode(true, false, false);  // common x-axis, no common y-axis, no redraw 
  XYPlot::Update(NULL, &y, &m);
}

void TSPlot::Update(TimeSeries& ts1, TimeSeries& ts2)
{
  Matrix y,m,y1,m1,y2,m2;
  y1 = ts1.GetDataVector();
  m1 = ts1.GetMissingVector();
  y2 = ts2.GetDataVector();
  m2 = ts2.GetMissingVector();
  y = hblockmatrix(y1,y2);
  m = hblockmatrix(m1,m2);
  string s1 = ts1.GetTitle();
  SetTitle(s1);
  XYPlot::Update(NULL, &y, &m);
}

//------ QQPlot -------------

QQPlot::QQPlot(TimeSeries& ts)
  : XYPlot(NULL, NULL, NULL)
{
  int i,j,nnonmissing;
  Vector tx,ty;
  double alpha;

  drawinglines = false;

  tx = ts.GetMissingVector();
  nnonmissing = ts.GetN() - tx.sum();
  tx.resize(nnonmissing);
  ty.resize(nnonmissing);

  for (i=0,j=0 ; i<nnonmissing ; ++i)
    {
      while (ts.IsMissing(j))
	++j;
      tx[i] = ts[j];
      ++j;
    }

  tx.sort();
  for (i=0 ; i<nnonmissing ; ++i)
    {
      alpha = (i+0.5)/nnonmissing;
      ty[i] = norm_invcdf(alpha);
    }

  string s1 = "Q-Q Plot (Normal)";
  SetTitle(s1);
  XYPlot::Update(&ty, &tx, NULL, is_realized());
}

void QQPlot::Update(TimeSeries& ts)
{
}

//------ ACFPlot ------------

void ACFPlot::RenderStuff()
{
  int i,uselags=nlags,startlag=includinglag0 ? 0 : 1;
  Vector templag, tempacf;

  if (localacf.nrows()<nlags)
    uselags = localacf.nrows();
  if (uselags - startlag <= 0)
    return;
  templag.resize(uselags - startlag);
  tempacf = localacf.subvector(startlag,uselags);
  for (i=startlag ; i<uselags ; ++i)
    templag[i-startlag] = i;
  XYPlot::Update(&templag, &tempacf, NULL, is_realized());
}

ACFPlot::ACFPlot(TimeSeries& ts)
  : XYPlot(NULL, NULL, 0)
{
  int i,actlag = min(max_possible_lag, ts.GetN());

  mode = GMODE_ACF;
  integerxaxis = true;

  nlags = default_maxlag;
  normalizing = true;
  includinglag0 = true;
  localacf.resize(actlag);
  locallag.resize(actlag);
  for (i=0; i<actlag ; ++i)
    locallag[i] = i;
  if (actlag>0)
    {
      ts.ComputeSampleACF(&localacf, NULL, 0);
      gamma0 = localacf[0];
      localacf = localacf/gamma0;
      actualnobs = ts.GetN();
      RenderStuff(); 
    }
}

ACFPlot::ACFPlot(TimeSeriesModel& tm)
  : XYPlot(NULL, NULL, 0)
{
  int i,actlag = max_possible_lag;
  mode = GMODE_ACF;
  nlags = default_maxlag;
  normalizing = true;
  includinglag0 = true;
  localacf.resize(actlag);
  locallag.resize(actlag);
  for (i=0; i<actlag ; ++i)
    locallag[i] = i;
  tm.ComputeACFPACF(localacf, NULL, 0);  // not normalized
  gamma0 = localacf[0];
  localacf = localacf/gamma0;
  actualnobs = 0;
  RenderStuff();
}

ACFPlot::~ACFPlot()
{
}

void ACFPlot::Update(TimeSeries& ts)
{
  int i, actlag = min(max_possible_lag, ts.GetN());
  localacf.resize(actlag);
  locallag.resize(actlag);
  for (i=0 ; i<actlag ; ++i)
    locallag[i] = i;
  if (actlag>0) {
    ts.ComputeSampleACF(&localacf, NULL, 0);
    //for (i=0 ; i<actlag ; ++i)
    //localacf[i] = actlag-i;
    gamma0 = localacf[0];
    if (normalizing)
      localacf = localacf/gamma0;
  }
  actualnobs = ts.GetN();
  RenderStuff();
}

void ACFPlot::Update(TimeSeriesModel& tm)
{
  int i;
  localacf.resize(max_possible_lag);
  locallag.resize(max_possible_lag);
  for (i=0 ; i<max_possible_lag ; ++i)
    locallag[i] = i;
  tm.ComputeACFPACF(localacf, NULL, 0);
  gamma0 = localacf[0];
  if (normalizing)
    localacf = localacf/gamma0;
  RenderStuff();
}

void ACFPlot::InitPropDialog(CenteringDialog& cd)
{
  Table *table = manage(new Table(2,2,true));
  table->attach(*(manage(new Label("# Lags"))), 
		0, 1, 0, 1, EXPAND, EXPAND, 10, 2);
  lagadjust = manage(new Adjustment(nlags,10,localacf.nrows(),5,25));
  lagspin = manage(new SpinButton(*lagadjust,0,0));
  table->attach(*lagspin, 0,1,1,2, EXPAND, EXPAND, 10, 2);
  cb1 = manage(new CheckButton("Normalize"));
  cb1->set_active(normalizing);
  cb2 = manage(new CheckButton("Include Lag 0"));
  cb2->set_active(includinglag0);
  table->attach(*cb1, 1,2,0,1, EXPAND, EXPAND);
  table->attach(*cb2, 1,2,1,2, EXPAND, EXPAND);
  cd.get_vbox()->pack_start(*table);
}

void ACFPlot::ClosePropDialog(CenteringDialog& cd)
{
  bool nn,i0;
  Vector templag, tempacf;

  // get num lags
  nlags = lagspin->get_value_as_int();
  // get new normalization
  nn = cb1->get_active();
  if (nn!=normalizing)
    {
      normalizing = nn;
      if (normalizing)
	localacf /= gamma0;
      else
	localacf *= gamma0;
    }
  includinglag0 = cb2->get_active();
  RenderStuff();
}

//------ PACFPlot ------------

PACFPlot::PACFPlot(TimeSeries& ts)
  : XYPlot(NULL, NULL, 0)
{
  int i;

  integerxaxis = true;
  mode = GMODE_ACF;

  localpacf.resize(default_maxlag);
  localacf.resize(default_maxlag);
  locallag.resize(default_maxlag);
  for (i=0; i<default_maxlag ; ++i)
    locallag[i] = i;
  ts.ComputeSampleACF(&localacf, &localpacf);
  actualnobs = ts.GetN();
  if (actualnobs>0)
    XYPlot::Update(&locallag, &localpacf, NULL, false);
}

PACFPlot::PACFPlot(TimeSeriesModel& tm)
  : XYPlot(NULL, NULL, 0)
{
  int i;
  mode = GMODE_ACF;
  localpacf.resize(default_maxlag);
  localacf.resize(default_maxlag);
  locallag.resize(default_maxlag);
  for (i=0; i<default_maxlag ; ++i)
    locallag[i] = i;
  tm.ComputeACFPACF(localacf, &localpacf);
  actualnobs = 0;
  XYPlot::Update(&locallag, &localpacf, NULL, false);
}

PACFPlot::~PACFPlot()
{
}

void PACFPlot::Update(TimeSeries& ts)
{
  ts.ComputeSampleACF(&localacf, &localpacf);
  actualnobs = ts.GetN();
  XYPlot::Update(&locallag, &localpacf, NULL);
}

void PACFPlot::Update(TimeSeriesModel& tm)
{
  tm.ComputeACFPACF(localacf, &localpacf);
  XYPlot::Update(&locallag, &localpacf, NULL);
}

//------- PatchPlot --------------------

PatchPlot::PatchPlot(TimeSeries& resids)
 : Gtk::DrawingArea()
{
  int i,j,nseries;

  // some defaults
  set_name("GraphWidget");

  // make local copy of data

  // initialize data
  n = resids.GetN()/2;
  u = new Vector[n];
  for (i=0 ; i<n ; ++i)
    {
      u[i].resize(2);
      u[i][0] =  norm_cdf(resids[2*i]);
      u[i][1] =  norm_cdf(resids[2*i+1]);
      cout << u[i][0] << " " << u[i][1] << endl;
    }

  title = "Patch Plot";

  nshades = 50;
  shades = new Gdk::Color[nshades];
}

PatchPlot::~PatchPlot()
{
  delete[] shades;
}

void PatchPlot::SetTitle(string& s)
{
  title = s;
}

void PatchPlot::Update()
{
}

void PatchPlot::on_realize()
{
  int i;

  // we need to do the default realize
  Gtk::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::Window> window = get_window();
  gc = Gdk::GC::create(window);

  // and allocate colors
  Glib::RefPtr<Gdk::Colormap> colormap = get_default_colormap();

  blue = Gdk::Color("blue");
  lightblue = Gdk::Color("lightblue");
  red = Gdk::Color("red");
  green = Gdk::Color("green");
  black = Gdk::Color("black");
  white = Gdk::Color("white");
  grey = Gdk::Color("grey");
  darkgrey = Gdk::Color("darkgrey");
  colormap->alloc_color(blue);
  colormap->alloc_color(lightblue);
  colormap->alloc_color(red);
  colormap->alloc_color(green);
  colormap->alloc_color(black);
  colormap->alloc_color(white);
  colormap->alloc_color(grey);
  colormap->alloc_color(darkgrey);

  for (i=0 ; i<nshades ; ++i)
    {
      shades[i].set_grey(65536*i/nshades);
      colormap->alloc_color(shades[i]);
    }

  window->clear();

  ComputePseudomap(80,80);
}

double PatchPlot::Intensity(double x, double y, double centerx, double centery, double sig)
// intensity always integrates to one, sig is a range of effect (like std. dev.)
{
  return (1/(sqrt(2*M_PI)*sig)) * exp(-((centerx-x)*(centerx-x)+(centery-y)*(centery-y))/(2*sig*sig));
}

void PatchPlot::ComputePseudomap(int nr, int nc)
{
  int i,j,k;
  double total_intensity = 1.0, range = pow(n,-0.3);
  double x0,x1,y0,y1,deltax=1.0/nc,deltay=1.0/nr,tx,ty;
  pseudomap.resize(nr,nc);
  baselinemap.resize(nr,nc);
  pseudomap.zeroes();
  baselinemap.zeroes();
  for (i=0 ; i<n ; ++i)
    {
      // take the residual and put it into pseudomap:
      // we do this by storing integral of intensity centered on residual over block
      
      for (j=0 ; j<nr ; ++j)
	for (k=0 ; k<nc ; ++k)
	  {
	    // get 4 corners of block
	    x0 = k*deltax;  x1 = x0+deltax;
	    y0 = j*deltay;  y1 = y0+deltay;
	    pseudomap[j][k] += Intensity((x0+x1)/2, (y0+y1)/2, u[i][0], u[i][1], range);
	  }
    }

  // compute baseline map (of biases due to area constraints)
  for (i=0 ; i<20*20 ; ++i)
    {
      tx = ((i%20)+0.5)/20;
      ty = ((i/20)+0.5)/20;
      for (j=0 ; j<nr ; ++j)
	for (k=0 ; k<nc ; ++k)
	  {
	    // get 4 corners of block
	    x0 = k*deltax;  x1 = x0+deltax;
	    y0 = j*deltay;  y1 = y0+deltay;
	    baselinemap[j][k] += Intensity((x0+x1)/2, (y0+y1)/2, tx, ty, range);
	  }
    }

  for (j=0 ; j<nr ; ++j)
    for (k=0 ; k<nc ; ++k)
      pseudomap[j][k] /= baselinemap[j][k];
}

bool PatchPlot::on_expose_event(GdkEventExpose *e)
{
  // get window size
  int ix0,iy0,ix1,iy1,w,h,d,i1,i2,nr=pseudomap.nrows(),nc=pseudomap.ncols();
  double deltax=1.0/nc,deltay=1.0/nr;
  int i, j, k, shade_index, ix, iy,  efsz = 8;
  double lx0,linelim,x0,x1,y0,y1;
  const double aspect = 0.5;
  int szx,szy;

  Glib::RefPtr<Gdk::Window> window = get_window();
  Glib::RefPtr<Gdk::GC> whitegc = get_style()->get_white_gc();

  window->clear();

  // here is where we draw on the window      
  window->get_geometry(ix0,iy0,w,h,d);
  ix0 = iy0 = 0;


  TRect outer(ix0+10,iy0+10,w-20,h-20);
  TRect plotrect(ix0+20,iy0+20,w-40,h-40);
  gc->set_foreground(shades[0]);
  window->draw_rectangle(gc, true, outer.left, outer.top,
			 outer.Width(), outer.Height());

  // go through and plot the pseudomap
  double m0 = pseudomap.min(), m1 = pseudomap.max();
  m0 = 0.0;

  Gdk::Color tcol;
  for (j=0 ; j<nr ; ++j)
    for (k=0 ; k<nr ; ++k)
      {
	// get 4 corners of block
	x0 = k*deltax;  x1 = x0+deltax;
	y0 = j*deltay;  y1 = y0+deltay;
	// scale to pixel coordinates
	ix0 = x0*plotrect.Width() + plotrect.left;
	ix1 = x1*plotrect.Width() + plotrect.left;
	iy0 = y0*plotrect.Height() + plotrect.top;
	iy1 = y1*plotrect.Height() + plotrect.top;

	// fill in rectangle
	shade_index = (pseudomap[j][k] - m0)*nshades/(m1-m0);
	if (shade_index==nshades)
	  --shade_index;
	gc->set_foreground(shades[shade_index]);
	window->draw_rectangle(gc, true, ix0, iy0, ix1-ix0, iy1-iy0);
      }
}
