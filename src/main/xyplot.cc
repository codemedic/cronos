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
 *
 * Updated Nov, 2003, Anthony Brockwell.
 * Now structure is in place for multivariate time series plots.
 * Popup menus for exporting data/graphics are implemented.
 *
 */

#include "includes.h"
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include "dialogs.h"
#include "rect.h"
#include "epsdc.h"

using namespace mslib;
using namespace Gtk;
using namespace SigC;


// The following variable is a pointer to an EpsPlot variable.
// If non-null, then the expose event should draw into it.

static TEpsPlot *eplot = NULL;


double min(double a, double b)
{
  return (a < b ? a : b);
}


XYPlot::XYPlot(Matrix *tx, Matrix *ty, Matrix *m)
  : Gtk::DrawingArea()
{
  int i,j,n,nseries;

  // some defaults
  common_x_axis = common_y_axis = true;
  drawinglines = true;
  integerxaxis = false;
  set_name("GraphWidget");
  basic_xunit = basic_yunit = 1.0;
  margin_fraction_x = margin_fraction_y = 0.05;
  color_map.resize(8);
  for (i=0 ; i<8 ; ++i)
    color_map[i] = i;

  if (ty!=NULL)
    {
      n=ty->nrows();
      nseries=ty->ncols();
    }
  else
    { n=0; nseries=1; }

  // make local copy of data
  nobs = n;
  if (n==0)
    n=1;
  x.resize(n,nseries);
  y.resize(n,nseries);
  missing.resize(n,nseries);
  if (nobs>0) {
    y = *ty;
    if (m==NULL)
      missing.zeroes(); // all false = no missing data
    else
      missing = *m;
    if (tx==NULL)
      for (i=0 ; i<n ; ++i)
	for (j=0 ; j<nseries ; ++j)
	  x[i][j] = i+1;
    else
      x = *tx;
  }

  title = "Time Series Plot";
  mode = GMODE_XY;

  // also enable relevant events
  Gdk::EventMask flags = get_events();
  flags |= (Gdk::BUTTON_PRESS_MASK);
  set_events(flags);
}

void XYPlot::SetTitle(string& s)
{
  title = s;
}

void XYPlot::SetMultiMode(bool cx, bool cy, bool redraw)
{
  common_x_axis = cx;
  common_y_axis = cy;
  if (redraw)
    RedrawIt();
}

void XYPlot::Update(Matrix *tx, Matrix *ty, Matrix *m, bool redraw)
{
  int i,j,n=ty->nrows(),nseries=ty->ncols();
  nobs = n;
  x.resize(n,nseries);   y.resize(n,nseries);   
  missing.resize(n,nseries);
  if (nobs>0) {
    y = *ty;
    if (m!=NULL)
      missing = *m;
    else
      missing.zeroes();
    if (tx!=NULL)
      x = *tx;
    else
      for (i=0 ; i<n ; ++i)
	for (j=0 ; j<nseries ; ++j)
	  x[i][j] = i+1;
  }

  if (redraw)
      RedrawIt();
}

void XYPlot::RedrawIt()
{
  Glib::RefPtr<Gdk::Window> window = get_window();
  int x1,y1,w,h,d;
  window->get_geometry(x1,y1,w,h,d);
  x1 = y1 = 0;
  Gdk::Rectangle* rect = new Gdk::Rectangle(x1,y1,w,h);
  window->invalidate_rect(*rect, false);
  delete rect;
}

XYPlot::~XYPlot()
{
}

void XYPlot::Oop1_Pressed()
{
  // get file name
  FileSelection fs;
  if (fs.run()==RESPONSE_OK)
    {
      string s1 = fs.get_filename();
      eps_draw(s1);
    }
}

void XYPlot::Oop5_Pressed()
  // allow user to change title
{
  string s = title;
  Widget *wp = get_toplevel();
  Gtk::Window *winp = (Gtk::Window*) wp;

  GeneralTextDialog gtd(*winp,"Plot Title: ", s);
  if (gtd.run()==RESPONSE_OK)
    // update the title
    {
      s = gtd.get_text();
      SetTitle(s);
      RedrawIt();
    }
}

void XYPlot::Oop2_Pressed()
// export to file option
{
  // get file name
  int i,j;
  FileSelection fs;
  if (fs.run()==RESPONSE_OK)
    {
      string s1 = fs.get_filename();
      Matrix z;
      z = hblockmatrix(x,y);

      // write manually, so missing values can be stored!
      ofstream outfile(s1.c_str());
      for (i=0 ; i<y.nrows() ; ++i)
	{
	  for (j=0 ; j<y.ncols() ; ++j)
	    if (missing[i][j])
	      outfile << "NA ";
	    else
	      outfile << y[i][j] << " ";
	  outfile << endl;
	}
    }
}

void XYPlot::Oop3_Pressed()
{
  // create a new plot, zoomed in on some portion
  double x0,x1; // new range
}

void XYPlot::Oop4_Pressed()
{
  Glib::ustring nm="Plot Properties";
  // allow user to modify properties
  Widget *wp = get_toplevel();
  Gtk::Window *winp = (Gtk::Window*) wp;

  CenteringDialog pd(nm, *winp);
  // action area
  pd.add_button("Cancel", RESPONSE_CANCEL);
  pd.add_button("O.K.", RESPONSE_OK);
  InitPropDialog(pd);
  pd.show_all();
  if (pd.run()==RESPONSE_OK)
    ClosePropDialog(pd);
}

void XYPlot::InitPropDialog(CenteringDialog& cd)
{
}

void XYPlot::ClosePropDialog(CenteringDialog& cd)
{
}

bool XYPlot::on_button_press_event(GdkEventButton* event)
{
  // do something!
  if( (event->type == GDK_BUTTON_PRESS) && (event->button == 3) )
  {
    Menu *local_menu = manage(new Menu());
    Menu_Helpers::MenuList& list1 = local_menu->items();
    
    Menu_Helpers::MenuElem m1("Export to EPS File",
			      mem_fun(*this, &XYPlot::Oop1_Pressed)), 
      m2("Export Data to Text File", mem_fun(*this, &XYPlot::Oop2_Pressed)),
      m6("Export Data to Clipboard", mem_fun(*this, &XYPlot::on_button_copy)),
      m3("Create Zoomed View", mem_fun(*this, &XYPlot::Oop3_Pressed)),
      m5("Edit Title", mem_fun(*this, &XYPlot::Oop5_Pressed)),
      m4("Properties",mem_fun(*this,&XYPlot::Oop4_Pressed));
    list1.push_back(m1);
    list1.push_back(m2);
    list1.push_back(m6);
    list1.push_back(m3);
    list1.push_back(m5);
    list1.push_back(m4);
    
    local_menu->popup(event->button, event->time);
    return true; // signify that it has been handled
  }
  else
    return false;
}

void XYPlot::LittleBox(Glib::RefPtr<Gdk::Window> win, Glib::RefPtr<Gdk::GC> gc, int x1, int y1, int r, Gdk::Color &fillcolor)
{
  if (eplot==NULL)
    win->draw_rectangle(gc, true, x1-r, y1-r, 2*r+1, 2*r+1);
  else
    {
      if (r>1)
	eplot->FillRect(x1-r,y1-r,x1+r,y1+r,255,255,255);
      if (r>=1)
	eplot->DrawRect(x1-r,y1-r,x1+r,y1+r);
    }
}

double XYPlot::GetNiceUnitFor(double x1, double x2, double prelim_scale)
//
// This function returns a "nice" subunit of the range (x1,x2) to use
// in between graph divisions.
// By choosing prelim_scale, for instance, to be M_PI,
// it gives nice divisions regarding the unit as M_PI (useful for plotting spectral densities)
//
{
  x1 /= prelim_scale;
  x2 /= prelim_scale;
  double w=(x2-x1), pw=floor (log10(w)), tens=pow(10,pw), units;
  w /= tens;		// now we have a number from 1.0 to 9.999...
  if (w<1.4)
    units = 0.2;
  else if (w<2.0)
    units = 0.25;
  else if (w<4)
    units = 0.5;
  else
    units = 1.0;
  return (units*tens*prelim_scale);
}

Gdk::Color *XYPlot::ColorNumber(int i)
{
  Gdk::Color *tref;
  switch ((int)(color_map[i])) {
    case 0 : tref=&black;  break;
    case 1 : tref=&red;    break;
    case 2 : tref=&green;  break;
    case 3 : tref=&blue;   break;
    case 4 : tref=&lightblue;  break;
    case 5 : tref=&grey;   break;
    case 6 : tref=&darkgrey;   break;
    case 7 : tref=&white;  break;
  }
  return tref;
}

void XYPlot::on_realize()
{
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

  window->clear();
}

void XYPlot::on_button_copy()
{
  //Build a string representation of the stuff to be copied:
  Glib::ustring strData;
  ostringstream os;
  int i,j;

  Matrix z;
  z = hblockmatrix(x,y);
  
  // write manually, so missing values can be stored!
  for (i=0 ; i<y.nrows() ; ++i)
    {
      for (j=0 ; j<y.ncols() ; ++j)
	if (missing[i][j])
	  os << "NA ";
	else
	  os << y[i][j] << " ";
      os << endl;
    }
  os << ends;

  strData += os.str();
   
  Glib::RefPtr<Gtk::Clipboard> refClipboard = Gtk::Clipboard::get();
  refClipboard->clear();
  refClipboard->set_text(strData);
}

void XYPlot::eps_draw(string& fnm)
{
  TRect r1(0,0,1500,1000);
  TEpsPlot e1(r1, fnm.c_str());

  eplot = &e1;
  on_expose_event(NULL);  // render to eps
  eplot = NULL;
  RedrawIt();             // and then redraw to screen
}

bool XYPlot::on_expose_event(GdkEventExpose* e)
{
  // get window size
  int ix0,iy0,w,h,d,i1,i2,szx,szy;
  int i, j, ix, iy, ix2, iy2, efsz = 8, series, nseries = y.ncols();
  double linelim, margin, tx, ty;
  Vector lx0, x0,x1,y0,y1,xunit,yunit;
  const double aspect = 0.5;
  int epsmode = (eplot!=NULL);
  Glib::RefPtr<Gdk::Window> window = get_window();
  Glib::RefPtr<Gdk::GC> whitegc = get_style()->get_white_gc();

  int active_series = 0;

  if (mode==GMODE_ACF)
    if (nseries>1)
      cout << "Internal error: cannot use ACF plots when num series > 1." << endl;

  // blank it out
  window->clear();  

  if (nobs==0)
    return true;  // and return if no data to plot
  
  if (!epsmode) // drawing on the window
    {
      window->get_geometry(ix0,iy0,w,h,d);
      ix0 = iy0 = 0;
    }
  else // for encapsulated postscript files, use a default pixel size of 1500x1000
    {
      ix0 = 0;   iy0 = 0;
      w = 1500;   h = 1000;
    }

  // 1. Subdivide the rectangle
  TRect clientrect(ix0,iy0,w,h),innerrect(ix0,iy0,w,h), ylabelrect, xlabelrect, graphrect, titlerect;
  i = clientrect.Height()/30;
  innerrect.top += i;   innerrect.bottom -= i;
  i = clientrect.Width()/30;
  innerrect.left += i;  innerrect.right -= i;

  i = innerrect.Height();   j= innerrect.Width();
  xlabelrect = ylabelrect = graphrect = titlerect = innerrect;
  if (title.length()>0)
    titlerect.bottom = titlerect.top + (int)(i*0.15);
  else
    titlerect.bottom = titlerect.top + (int)(i*0.03); // small padding if no title
  graphrect.top = titlerect.bottom;
  ylabelrect.right = ylabelrect.left + (int)(j*0.1);
  graphrect.left = ylabelrect.right;
  titlerect.left = graphrect.left;
  xlabelrect.top = xlabelrect.bottom - (int)(i*0.08);
  graphrect.bottom = xlabelrect.top;

  // 2. Determine limits for each time series
  minx.resize(nseries);
  maxx.resize(nseries);
  miny.resize(nseries);
  maxy.resize(nseries);
  x0.resize(nseries);
  x1.resize(nseries);
  y0.resize(nseries);
  y1.resize(nseries);
  lx0.resize(nseries);
  xunit.resize(nseries);
  yunit.resize(nseries);

  // find minima and maxima
  for (series=0 ; series<nseries ; ++series)
    {
      minx[series]=miny[series]=1e100;
      maxx[series]=maxy[series]=-1e100;
      for (i=0 ; i<nobs ; ++i)
	  if (!missing[i][series])
	    {
	      if (x[i][series]<minx[series])
		minx[series] = x[i][series];
	      if (x[i][series]>maxx[series])
		maxx[series] = x[i][series];
	      if (y[i][series]<miny[series])
		miny[series] = y[i][series];
	      if (y[i][series]>maxy[series])
		maxy[series] = y[i][series];
	    }
      if (fabs(minx[series]-maxx[series])<1e-10) {
	minx[series] -= 1.0;
	maxx[series] += 1.0;
      }
      if (fabs(miny[series]-maxy[series])<1e-10) {
	miny[series] -= 1.0;
	maxy[series] += 1.0;
      }
    }

  // adjust if there are common axes
  if (common_x_axis)
    {
      tx = minx.min();
      for (i=0 ; i<nseries ; ++i)
	minx[i] = tx;
      tx = maxx.max();
      for (i=0 ; i<nseries ; ++i)
	maxx[i] = tx;
    }
  if (common_y_axis)
    {
      tx = miny.min();
      for (i=0 ; i<nseries ; ++i)
	miny[i] = tx;
      tx = maxy.max();
      for (i=0 ; i<nseries ; ++i)
	maxy[i] = tx;
    }

  // determine "good" units
  for (series=0 ; series<nseries ; ++series)
    {
      if (mode==GMODE_ACF)  // note: GMODE_ACF should only be used when nseries==1
	{
	  if (nobs>0)
	    {
	      linelim = 1.96*y[0][0]/sqrt((double)actualnobs);
	      maxx[series] = x[0][nobs-1]+1;
	    }
	  else
	    linelim = 0;
	  
	  // and symmetrize miny, maxy
	  maxy[series] = max(-miny[series],maxy[series]);
	  miny[series] = -maxy[series];
	}

      xunit[series] = GetNiceUnitFor(minx[series],maxx[series],basic_xunit);
      yunit[series] = GetNiceUnitFor(miny[series],maxy[series],basic_yunit);

      margin=(maxx[series]-minx[series])*margin_fraction_x;
      lx0[series]=x0[series]=floor((minx[series]-margin)/xunit[series])*xunit[series];
      if ((minx[series]>=0) && (x0[series]<0))
	lx0[series]=x0[series]=0;
      x1[series]=ceil((maxx[series]+margin)/xunit[series])*xunit[series];
      if (integerxaxis)
	{
	  lx0[series] = x0[series] = 0;
	  x1[series] = maxx[series];
	}
      
      margin=(maxy[series]-miny[series])*margin_fraction_y;
      y0[series]=floor((miny[series]-margin)/yunit[series])*yunit[series];
      if ((miny[series]>=0) && (y0[series]<0))
	y0[series] = 0;
      y1[series]=ceil((maxy[series]+margin)/yunit[series])*yunit[series];
      
      if (mode==GMODE_ACF)
	{
	  x0[series] = -1;                   // acfs start at lag 0 (making left xlim = -1 separates lag 0 from the y axis)
	  // also fix scaling on normalized ACF plots
	  if (maxy[series]==1.0)
	    { yunit[series] = 0.25;   y0[series] = -1.0;   y1[series] = 1.0;   }
	}
    } // end of series loop


  // 3. DRAWING BEGINS HERE

  // 3a. Do the frame amd title first
  if (epsmode)
    eplot->DrawRect(graphrect);
  else
    {
      gc->set_foreground(white);
      window->draw_rectangle(gc, true, graphrect.left, graphrect.top,
			     graphrect.Width(), graphrect.Height());
      gc->set_foreground(darkgrey);
      window->draw_rectangle(gc, false, graphrect.left, graphrect.top,
			     graphrect.Width()-1, graphrect.Height()-1);
      gc->set_foreground(black);
    }

  // Setup pango context for displaying text
  int ht = (int)(titlerect.Height()*0.3),
    wd = (int)(titlerect.Width()/48);
  if (ht*aspect<wd)
    wd = (int)(ht*aspect);
  else
    ht = (int)(wd/aspect); // desired font size
  char title_fontnm[30], little_fontnm[30];
  sprintf(title_fontnm, "Sans Bold %d", ht);
  sprintf(little_fontnm, "Sans %d", (int)(ht*0.6));

  Pango::FontDescription title_fdescrip(title_fontnm),
    little_fdescrip(little_fontnm);
  Glib::RefPtr<Pango::Context> pc = Gtk::Widget::create_pango_context();
  Glib::RefPtr<Pango::Layout> play = Pango::Layout::create(pc);

  // Put title in titlebox
  if (title.length()>0)
    {
      play->set_font_description(title_fdescrip);
      play->set_text(title);
      play->get_pixel_size(w,h);
      i1=(titlerect.left+titlerect.right)/2;
      i2=(titlerect.bottom+titlerect.top)/2;
      if (!epsmode)
	window->draw_layout(gc, i1-w/2, i2-h/2, play);
      else
	{
	  eplot->SetFontSize(efsz*5,efsz*7.5);
	  eplot->TextOut(i1, i2, title.c_str(), true);
	  eplot->SetFontSize(efsz*2.5,efsz*3.75);
	}
      play->set_font_description(little_fdescrip);
    }

  // ACF bars if in ACF mode (this is nseries==1 only! so linelim doesn't need to be vector)
  if (mode==GMODE_ACF)
    {
      // also plot \pm 1.96/sqrt(n)
      int l1,l2;
      l1 = -(linelim-y0[0])/(y1[0]-y0[0])*graphrect.Height() + graphrect.bottom;
      l2 = -(-linelim-y0[0])/(y1[0]-y0[0])*graphrect.Height() + graphrect.bottom;
       gc->set_foreground(lightblue);
      if (!epsmode)
	window->draw_rectangle(gc, true, graphrect.left+1, l1, graphrect.Width()-2, l2-l1+1);
      else
	eplot->FillRect(graphrect.left+1, l1, graphrect.right-1, l2, 190, 190, 255);
    }

  // 3b. Draw x-axis labels and lines
  char temps[40];
  ostrstream os(temps, 40);
  bool inepspath;
  int xlabelprecision = -floor(log10(xunit[active_series]))+1;
  if (xlabelprecision<0)
    xlabelprecision = 0;

  j = graphrect.Height()/60;
  for (tx=lx0[active_series] ; tx<=x1[active_series] ; tx+=xunit[active_series])
    {
      ix = (int)((tx-x0[active_series])/(x1[active_series]-x0[active_series])*graphrect.Width() + graphrect.left);
      if ((ix>graphrect.left+1) && (ix<graphrect.right))
	if (!epsmode)
	  {
	    gc->set_foreground(black);
	    window->draw_line(gc, ix, graphrect.top, ix, graphrect.top-j); 
	    window->draw_line(gc, ix, graphrect.bottom, ix, graphrect.bottom+j);
	    gc->set_foreground(grey);
	    window->draw_line(gc, ix, graphrect.top+1, ix, graphrect.bottom); 
	  }
        else
	  {
	    eplot->SetLineMode(0);
	    eplot->MoveTo(ix,graphrect.top); eplot->LineTo(ix, graphrect.top-j);
	    eplot->FinishPath();
	    eplot->MoveTo(ix,graphrect.bottom); eplot->LineTo(ix, graphrect.bottom+j);
	    eplot->FinishPath();
	    eplot->SetLineMode(1);
	    eplot->MoveTo(ix,graphrect.top+1);
	    eplot->LineTo(ix,graphrect.bottom);
	    eplot->FinishPath();
	  }

      gc->set_foreground(black);
      iy = (xlabelrect.top + xlabelrect.bottom)/2;
      os.seekp(0);
      os << setiosflags(ios::fixed) << setprecision(xlabelprecision) << tx << ends;

      play->set_text(temps);
      play->get_pixel_size(szx,szy);
      if (epsmode)
	eplot->TextOut(ix,iy,temps,true);
      else
	window->draw_layout(gc, ix-szx/2, iy-szy/2, play);
    }
  if (epsmode)
    eplot->SetLineMode(0);

  // 3c. Draw y-axis labels and lines
  int localp = (int)(3-floor(log10(yunit[active_series])));
  if (localp<0)
    localp=0;
  if (epsmode)
    eplot->SetLineMode(1);
  for (ty=y0[active_series] ; ty<=y1[active_series] ; ty+=yunit[active_series])
    {
      gc->set_foreground(black);
      iy = (int)(-(ty-y0[active_series])/(y1[active_series]-y0[active_series])*graphrect.Height() + graphrect.bottom);
      ix = (int)(ylabelrect.right * 0.75 + ylabelrect.left * 0.25);
      os.seekp(0);
      os << setiosflags(ios::fixed) << setprecision(localp) << ty << ends;

      play->set_text(temps);
      play->get_pixel_size(szx,szy);
      if (epsmode)
	eplot->TextOut(ix,iy+efsz*1.5,temps,false);
      else
	window->draw_layout(gc, ix-szx, iy-szy/2, play);

      gc->set_foreground(grey);
      if ((iy>graphrect.top+1) && (iy<graphrect.bottom)) 
	if (!epsmode)
	  window->draw_line(gc, graphrect.left, iy, graphrect.right, iy);
	else {
	  eplot->MoveTo(graphrect.left, iy);
	  eplot->LineTo(graphrect.right, iy);
	}
    }

  // draw the y=0 line in red 
  gc->set_foreground(red);
  iy = (int)(y0[active_series]/(y1[active_series]-y0[active_series])*graphrect.Height() + graphrect.bottom);
  if ((iy>graphrect.top+1) && (iy<graphrect.bottom)) 
    if (!epsmode)
      window->draw_line(gc, graphrect.left, iy, graphrect.right, iy);
    else {
      eplot->MoveTo(graphrect.left, iy);
      eplot->LineTo(graphrect.right, iy);
    }  
  gc->set_foreground(grey);

  if (epsmode)
    eplot->SetLineMode(0);

  // 3d. Draw Graph Contents
  int boxsize = (int)(min(graphrect.Width(), graphrect.Height())/100),
    otherlimit = graphrect.Width()/nobs/2, snum;
  if (boxsize>otherlimit)
    boxsize = otherlimit;
  if (boxsize<1)
    boxsize = 1;

  for (snum=0 ; snum<nseries ; ++snum) // go through multiple series if necessary
    {
      if (snum==active_series)
	gc->set_foreground(black);
      else
	gc->set_foreground(*ColorNumber(snum));

      ix = (int)((x[0][snum]-x0[snum])/(x1[snum]-x0[snum])*graphrect.Width() + graphrect.left);
      iy = (int)(-(y[0][snum]-y0[snum])/(y1[snum]-y0[snum])*graphrect.Height() + graphrect.bottom);
      if (mode==GMODE_XY)
	{
	  inepspath = false;
	  if (drawinglines)
	    for (i=1 ; i<nobs ; ++i)
	      {
		ix2 = (int)((x[i][snum]-x0[snum])/(x1[snum]-x0[snum])*graphrect.Width() + graphrect.left);
		iy2 = (int)(-(y[i][snum]-y0[snum])/(y1[snum]-y0[snum])*graphrect.Height() + graphrect.bottom);
		if (!missing[i-1][snum] && !missing[i][snum])
		  {
		    if (!epsmode)
		      window->draw_line(gc,ix,iy,ix2,iy2);
		    else
		      {
			if (!inepspath)
			  eplot->MoveTo(ix,iy);
			eplot->LineTo(ix2,iy2);
		      }
		  }
		else if (epsmode)
		  inepspath = false;
		ix = ix2;   iy = iy2;
	      }
	  if (epsmode)
	    eplot->FinishPath();
	  for (i=0 ; i<nobs ; ++i)
	    {
	      ix = (int)((x[i][snum]-x0[snum])/(x1[snum]-x0[snum])*graphrect.Width() + graphrect.left);
	      iy = (int)(-(y[i][snum]-y0[snum])/(y1[snum]-y0[snum])*graphrect.Height() + graphrect.bottom);
	      if (!missing[i][snum])
		LittleBox(window, gc, ix, iy, boxsize, white);
	    }
	}
      
      if (mode==GMODE_ACF) // nseries==1 case only!	
	{
	  // draw vertical bars instead of lines
	  int barwidth = 1.0/(x1[0]-x0[0])*graphrect.Width() * 0.3, baroffset,iy2;
	  if (barwidth<1)
	    barwidth = 1;

	  baroffset = -barwidth/2;
	  for (i=0 ; i<nobs ; ++i)
	    {
	      ix = (x[i][snum]-x0[0])/(x1[0]-x0[0])*graphrect.Width() + graphrect.left;
	      iy = -(y[i][snum]-y0[0])/(y1[0]-y0[0])*graphrect.Height() + graphrect.bottom;
	      iy0 = (y0[0])/(y1[0]-y0[0])*graphrect.Height() + graphrect.bottom;
	      h = (iy-iy0+1);
	      if (h>0)
		if (!epsmode)
		  window->draw_rectangle(gc, 1, ix+baroffset+1, iy0, barwidth, h);
		else
		  eplot->FillRect(ix+baroffset+1, iy0, ix+baroffset+1+barwidth,
				  iy0+h, 0, 0, 0);
	      else
		if (!epsmode)
		  window->draw_rectangle(gc, 1, ix+baroffset+1, iy0+h, barwidth, -h);
		else
		  eplot->FillRect(ix+baroffset+1, iy0+h, ix+baroffset+1+barwidth,
				  iy0, 0, 0, 0);
	    }
	}
    } // end of for snum loop

  gc->set_foreground(black);  // restore color
  return true;
}

