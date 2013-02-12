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
 * Code for CoreWin class in Cronos.
 * This is the top-level window, which keeps track of 
 * a time series and a time series model, and keeps the
 * corresponding displays updated.
 */


#include "includes.h"
#include <gtkmm/notebook.h>
#include <gtkmm/scrolledwindow.h>
#include <fstream>
#include "textrender.h"
#include "tsmrender.h"
#include "dialogs.h"
#include <string>
#include <vector>
#include "otherpopups.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_multimin.h>
#include "tsinterface.h"
#include <sigc++/bind.h>
#include <sigc++/slot.h>
#include "mainwin.h"

#include "newabout.xpm"
#include "flogs.xpm"
#include "ffore.xpm"
#include "fdiff.xpm"
#include "ftreverse.xpm"
#include "farma.xpm"
#include "fgarch.xpm"
#include "fegarch.xpm"
#include "fsvm.xpm"
#include "fresids.xpm"
#include "fmle.xpm"
#include "fintegrate.xpm"

using namespace SigC;
using namespace Gtk;
using namespace std;


CoreWin::CoreWin()
  // main constructor
  : Window()
{
  set_name("MITSMGraphWidget");
  set_default_size(960,700);
  set_title("Cronos");

  timeseries = new TimeSeries();

  ARMAModel *tap = new ARMAModel();
  model = tap;
  interface = new ARMAInterface(tap);
  modeltype = 1;  // for ARMA models

  //
  // Construct the menu bar.
  //
  MenuBar *menubar = manage(new MenuBar());
  Menu_Helpers::MenuList& menubarlist = menubar->items();
  
  // File menu
  Menu *file_menu = manage(new Menu());
  Menu_Helpers::MenuList& list = file_menu->items();  
  Menu_Helpers::MenuElem m0("_New", mem_fun(*this, &CoreWin::file_new_cb)),
    m1("_Open", mem_fun(*this, &CoreWin::file_open_cb)),
    m1b("_Save", mem_fun(*this, &CoreWin::file_save_cb)),
    m2("E_xit", mem_fun(*this, &CoreWin::file_exit_cb));
  list.push_back(m0);
  list.push_back(m1);
  list.push_back(m1b);
  list.push_back(m2);
  menubarlist.push_front(Menu_Helpers::MenuElem("_File", *file_menu));

  // Transform menu comes under a data menu now
  Menu *data_menu = manage(new Menu());
  list = data_menu->items();

  Menu_Helpers::MenuElem m21("_Import from Clipboard", mem_fun(*this,&CoreWin::data_import_clipboard));
  list.push_back(m21);

  Menu *transform_menu = manage(new Menu()), *ls_submenu = manage(new Menu());
  Menu_Helpers::MenuList& sublist = transform_menu->items(),
    ls_sublist = ls_submenu->items();
  Menu_Helpers::MenuElem m4("_Box-Cox/Log"), m5("_Difference"), m5b("_Integrate"),
    m7("_Aggregate"), m8("A_bsolute  Value"), mls("_Least Squares", *ls_submenu),
    m9("_Remove Spline Trend"), m9b("De_seasonalize"), m9c("_Clip"),
    m9d("Custom"), m9e("Reverse _Time");

  sublist.push_back(m4);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),1,true));
  sublist.push_back(m5);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),2,true));
  sublist.push_back(m5b);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),10,true));
  sublist.push_back(m7);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),3,true));
  sublist.push_back(m8);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),4,true));
  ls_sublist.push_back(m9);
  ls_sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),5,true));
  ls_sublist.push_back(m9b);
  ls_sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),6,true));
  sublist.push_back(mls);
  sublist.push_back(m9c);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),7,true));
  sublist.push_back(m9d);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),8,true));
  sublist.push_back(m9e);
  sublist.back().signal_activate().connect(bind<int>(mem_fun(*this,&CoreWin::transform_cb),9,true));

  list.push_back(Menu_Helpers::MenuElem("_Transform", *transform_menu));
  Menu_Helpers::MenuElem m20("_Periodogram",mem_fun(*this,&CoreWin::data_periodogram));
  list.push_back(m20);
  menubarlist.push_back(Menu_Helpers::MenuElem("_Data", *data_menu));


  // Model menu
  Menu *model_menu = manage(new Menu()), *submodel_menu = manage(new Menu()),
    *subfit_menu = manage(new Menu()), *subview_menu = manage(new Menu());
  Menu_Helpers::MenuList& viewsublist = subview_menu->items();
  list = model_menu->items();
  sublist = submodel_menu->items();
  Menu_Helpers::MenuList& subfitlist = subfit_menu->items();
  Menu_Helpers::MenuElem m3("_Specify", *submodel_menu),
    m3a("ARMA",bind<int>(mem_fun(*this,&CoreWin::model_specify_cb),1)),
    m3b("GARCH",bind<int>(mem_fun(*this,&CoreWin::model_specify_cb),2)),
    m3d("EGARCH",bind<int>(mem_fun(*this,&CoreWin::model_specify_cb),4)),
    m3c("SVM",bind<int>(mem_fun(*this,&CoreWin::model_specify_cb),3)),
    m60("_Fit", *subfit_menu),
    m60a("MLE", bind<int>(mem_fun(*this,&CoreWin::model_fit_cb),METHOD_MLE)),
    m60b("Bayesian", bind<int>(mem_fun(*this,&CoreWin::model_fit_cb),METHOD_BAYESIAN)),
    m6a("Model _ACF/PACF", mem_fun(*this,&CoreWin::model_acf_pacf_cb)),
    m6b("_Residual Analysis", mem_fun(*this,&CoreWin::model_residuals_cb)),
    m6c("Simulate",mem_fun(*this,&CoreWin::model_simulate)),
    m6d("_Forecast",mem_fun(*this,&CoreWin::model_forecast_cb)),
    m61("_Spectral Density",mem_fun(*this,&CoreWin::model_spectral_density)),
    m62("_View", *subview_menu),
    m62a("Continuous-time representation", mem_fun(*this,&CoreWin::model_view_cts_version_cb));
  sublist.push_back(m3a);
  sublist.push_back(m3b);
  sublist.push_back(m3d);
  sublist.push_back(m3c);
  viewsublist.push_back(m62a);
  subfitlist.push_back(m60a);
  subfitlist.push_back(m60b);
  list.push_back(m3);
  list.push_back(m60);
  list.push_back(m62);
  list.push_back(m6a);
  list.push_back(m61);
  list.push_back(m6b);
  list.push_back(m6c);
  list.push_back(m6d);
  menubarlist.push_back(Menu_Helpers::MenuElem("_Model", *model_menu));

  // Help menu
  Menu *help_menu = manage(new Menu());
  list = help_menu->items();
  Menu_Helpers::MenuElem m10("About");
  list.push_back(m10);
  list.back().signal_activate().connect(mem_fun(*this,&CoreWin::help_about_cb));

  menubarlist.push_back(Menu_Helpers::MenuElem("_Help", *help_menu));

  // Set up the main window.
  plotwin = manage(new TSPlot(*timeseries));
  acfplot = manage(new ACFPlot(*timeseries));
  pacfplot = manage(new PACFPlot(*timeseries));
  string mytitle = "Sample ACF";
  acfplot->SetTitle(mytitle);
  mytitle = "Sample PACF";
  pacfplot->SetTitle(mytitle);
  modelrenderer = manage(new TSMRenderer(model, timeseries, &transformlist));

  Table *table = manage(new Table(3,3,true));
  table->set_name("MITSMGraphWidget");
  table->attach(*plotwin, 1,3,0,2);
  table->attach(*acfplot, 1,2,2,3);
  table->attach(*pacfplot, 2,3,2,3);

  ScrolledWindow *sw = manage(new ScrolledWindow());
  sw->set_policy(POLICY_AUTOMATIC, POLICY_AUTOMATIC);  
  sw->add(*modelrenderer);
  Frame *frame = manage(new Frame());
  frame->add(*sw);
  table->attach(*frame, 0,1,0,2);

  // Setup the shortcut buttons
  HBox *hbox = manage(new HBox(false, 2));

  Button *b1 = manage(new Button(Gtk::Stock::OPEN)), 
    *b2 = manage(new Button()), *b3 = manage(new Button(Gtk::Stock::GO_BACK)), 
    *b4 = manage(new Button()), *b5 = manage(new Button()), 
    *b6 = manage(new Button()), *b7 = manage(new Button()),
    *b8 = manage(new Button()), *b9 = manage(new Button()),
    *b10 = manage(new Button()), *b11 = manage(new Button()),
    *b12 = manage(new Button()), *b13 = manage(new Button());

  Gtk::Image *p1 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(ffore_xpm))), 
    *p2 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(flogs_xpm))),
    *p3 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(fdiff_xpm))),
    *p4 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(ftreverse_xpm))),
    *p7 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(farma_xpm))),
    *p8 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(fresids_xpm))),
    *p9 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(fmle_xpm))),
    *p10 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(fgarch_xpm))),
    *p11 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(fsvm_xpm))),
    *p12 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(fegarch_xpm))),
    *p13 = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(fintegrate_xpm)));

  Frame *data_shortcuts = manage(new Frame()), *model_shortcuts = manage(new Frame()),
    *model_and_data_shortcuts = manage(new Frame());
  HBox *data_hbox = manage(new HBox(false, 2)), *model_hbox = manage(new HBox(false, 2)),
    *model_and_data_hbox = manage(new HBox(false, 2));

  b1->signal_clicked().connect(mem_fun(*this,&CoreWin::file_open_cb));
  b2->add(*p2);
  b2->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::transform_cb), 1, false));
  b5->add(*p3);
  b5->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::transform_cb),2, false));  
  b3->signal_clicked().connect(mem_fun(*this,&CoreWin::undo_transform_cb));
  b4->add(*p1);
  b4->signal_clicked().connect(mem_fun(*this,&CoreWin::model_forecast_cb));
  b6->add(*p4);
  b6->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::transform_cb), 9, false));
  b7->add(*p7);
  b7->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::model_specify_cb), 1));
  b8->add(*p8);
  b8->signal_clicked().connect(mem_fun(*this,&CoreWin::model_residuals_cb));
  b9->add(*p9);
  b9->signal_clicked().connect(bind<int>(mem_fun(*this,&CoreWin::model_fit_cb),METHOD_MLE));
  b10->add(*p10);
  b10->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::model_specify_cb), 2));
  b11->add(*p11);
  b11->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::model_specify_cb), 3));
  b12->add(*p12);
  b12->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::model_specify_cb), 4));
  b13->add(*p13);
  b13->signal_clicked().connect(bind<int>(mem_fun(*this, &CoreWin::transform_cb), 10, false));

  data_hbox->pack_start(*b1, false, false);  // open
  data_hbox->pack_start(*b2, false, false);  // log
  data_hbox->pack_start(*b5, false, false);  // diff
  data_hbox->pack_start(*b13, false, false); // integrate
  data_hbox->pack_start(*b6, false, false);  // reverse time
  data_hbox->pack_start(*b3, false, false);  // reverse previous transform

  data_shortcuts->add(*data_hbox);
  data_shortcuts->set_label("Data");

  model_hbox->pack_start(*b7, false, false);  // ARMA select
  model_hbox->pack_start(*b10, false, false);  // GARCH select
  model_hbox->pack_start(*b12, false, false);  // EGARCH select
  model_hbox->pack_start(*b11, false, false);  // SVM select

  model_shortcuts->add(*model_hbox);
  model_shortcuts->set_label("Model Select");

  model_and_data_hbox->pack_start(*b9, false, false);  // MLE fit
  model_and_data_hbox->pack_start(*b8, false, false);  // residual analysis
  model_and_data_hbox->pack_start(*b4, false, false);  // forecast

  model_and_data_shortcuts->add(*model_and_data_hbox);
  model_and_data_shortcuts->set_label("Data+Model");

  hbox->pack_start(*data_shortcuts,false,false);
  hbox->pack_start(*model_shortcuts,false,false,20);
  hbox->pack_start(*model_and_data_shortcuts,false,false,20);


  VBox *vbox = manage(new VBox());
  vbox->pack_start(*menubar, false, false);
  vbox->pack_start(*hbox, false, false);
  vbox->pack_start(*table, true, true);
  add(*vbox);
}

void corewin_callback(int stage, void *parms)
{
  // stage: 0=optimization, 1=partial deriv evaluation, 2=MCMC
  CoreWin *cp = (CoreWin *) parms;
  if (cp!=NULL)
    cp->model_fit_iteration(stage);
}

void CoreWin::file_new_cb()
{
  // 2. build window to put them in
  Window *extra_win = new CoreWin();
  extra_win->set_title("Cronos Subsidiary Window");
  extra_win->show_all();
}

void CoreWin::file_open_cb()
{
  fs = new FileSelection();
  fs->set_filename("*.dat");

  if (fs->run()==RESPONSE_OK)
    {
      string s1 = fs->get_filename();
      const char *c = s1.c_str();

      ifstream infile(c);
      if (!infile.good())
	{
	  char temps[40]="Problem opening file.";
	  // give error message and return
	  MessageDialog ed(*this, temps);
	  ed.run();
	  delete fs;
          return;
	}

      clear_transforms();
      infile >> (*timeseries);
      infile.close();

      if (timeseries->GetN()==0)
	{
	  char temps[80]="No data read.  This could be a problem with the file format.";
	  // give error message and return
	  MessageDialog ed(*this, temps);
	  ed.run();
	  delete fs;
	  return;
	}

      model->ComputeLogLikelihood(timeseries);
      update_all(true);
    }
  delete fs;
}

void CoreWin::data_import_clipboard()
{
  timeseries->ClearOut();
  Glib::RefPtr<Gtk::Clipboard> refClipboard = Gtk::Clipboard::get();
  local_data = refClipboard->wait_for_text();
  string str1 = local_data;
  istringstream ins(str1);
  ins >> (*timeseries);
  model_specified();  // recompute likelihood, etc.
  update_all(true);
}


void CoreWin::on_clipboard_received(const Gtk::SelectionData& selection_data)
{
  const std::string target = selection_data.get_target();

  if(target == local_data) //It should always be this, because that' what we asked for when calling request_contents().
  {
    cout << "Got it damn it!" << endl;
  }
  else
    {
      cout << "Notunderstood." << endl;
    }
  update_all(true);
}


void CoreWin::file_save_cb()
{
  fs = new FileSelection();
  if (fs->run()==RESPONSE_OK)
    {
      string s1 = fs->get_filename();
      ofstream outfile(s1.c_str());
      outfile << *timeseries;
      outfile.close();
    }
  delete fs;
}

void CoreWin::file_exit_cb()
{
  hide();
}

void CoreWin::data_periodogram()
{
  if (timeseries->GetN()==0)
    {
      char temps[80]="No data available for periodogram.  Load or simulate data and retry.";
      // give error message and return
      MessageDialog ed(*this, temps);
      ed.run();
      delete fs;
      return;
    }

  // 1. compute the periodogram
  int i,j,nn,na,bw,tsn = timeseries->GetN();
  Vector pdgm = timeseries->ComputePeriodogram();
  Matrix freqs;
  nn = (tsn+1)/2;
  bw = (int)(floor(sqrt(tsn))+0.5);
  freqs.resize(nn,2);
  for (i=0 ; i<nn ; ++i)
    freqs[i][0] = freqs[i][1] = (i+1)*2*M_PI/tsn;

  // and the smoothed periodogram (half-bandwidth = 0.5*sqrt(n))
  Matrix withsmoothed(nn,2);
  for (i=0 ; i<nn ; ++i)
    {
      withsmoothed[i][0] = pdgm[i+1];     // skip constant mean term
      withsmoothed[i][1] = 0.0;
      for (na=0, j=i-bw+1 ; j<=i+bw+1 ; ++j)
	if (j>=1 && j<=nn)
	  { ++na;  withsmoothed[i][1] += pdgm[j]; }
      withsmoothed[i][1] /= na; 
    }


  // 2. build window to put them in
  Window *pd_win = new Window;
  pd_win->set_default_size(600,400);
  string t=timeseries->GetTitle() + " Periodogram";
  SpectrumWin* xyp =  manage(new SpectrumWin(freqs, withsmoothed, t));
  pd_win->add(*xyp);
  pd_win->show_all();
}

void CoreWin::clear_transforms()
{
  Transformation *tp;
  int i,n = transformlist.size();
  // iterate through and delete
  for (i=0 ; i<n ; ++i)
    {
      tp = transformlist[i];
      delete tp;
    }
  // clear it out!
  transformlist.clear();
}

void CoreWin::update_all(bool acf)
{
  if (acf)
    {
      acfplot->Update(*timeseries);
      pacfplot->Update(*timeseries);
    }
  // now update plotwin
  plotwin->Update(*timeseries);
  modelrenderer->Update(model,timeseries,&transformlist);
}

void CoreWin::undo_transform_cb()
{
  int i,nt=transformlist.size();
  ostringstream err;
  if (nt==0)
    return;
  trans_loc = transformlist[nt-1];
  bool wasok = trans_loc->Reverse(timeseries, err);
  if (!wasok)
    {
      MessageDialog ed(*this, err.str());
      ed.run();
    }
  else
    {
      // clean up memory
      delete trans_loc;
      transformlist.pop_back();
    }
  model->ComputeLogLikelihood(timeseries);
  update_all(true);
}

void CoreWin::transform_specified(bool getparms)
{
  ostringstream err;
  bool wasok = trans_loc->Apply(timeseries, err, getparms);
  if (!wasok)
    {
      delete trans_loc;
      MessageDialog ed(*this, err.str());
      ed.run();
    }
  else
    transformlist.push_back(trans_loc);
}

void CoreWin::transform_cb(int tt, bool getparms)
{
  TransformDialog *tdp;
  RemoveSplineTransform *rttp;
  RemoveSeasonal *rsp;
  ClipTransform *ctp;
  CustomTransform *cstp;
  TimeReverseTransform *trtp;

  if (timeseries->GetN()==0)
    {
      char temps[80] = "You cannot apply a data transformation when there is no data.";
      MessageDialog ed(*this, temps);
      ed.run();
      return;
    }

  switch (tt)
    {
    case 1: // take logs
      trans_loc = new LogTransform();
      break;
    case 2: // difference
      trans_loc = new DiffTransform();
      break;
    case 10: // integrate
      trans_loc = new IntegrateTransform();
      break;
    case 3: // aggregate
      trans_loc = new AggTransform();
      break;
    case 4: // abs. value
      trans_loc = new AbsTransform();
      break;
    case 5: // trend removal
      rttp = new RemoveSplineTransform();
      rttp->SetPlotWin(plotwin);
      trans_loc = rttp;
      break;
    case 6: // seasonal comp. removal
      rsp = new RemoveSeasonal();
      rsp->SetPlotWin(plotwin);
      trans_loc = rsp;
      break;
    case 7: // clipping transformation
      ctp = new ClipTransform();
      ctp->SetPlotWin(plotwin);
      trans_loc = ctp;
      break;
    case 8: // custom transformation (using mu_parser)
      cstp = new CustomTransform();
      trans_loc = cstp;
      break;
    case 9: // time reversal
      trtp = new TimeReverseTransform();
      trans_loc = trtp;
      break;
    }


  trans_loc->GetParList();
  tdp = NULL;
  if ((trans_loc->npars>0) && (getparms)) {
    tdp = new TransformDialog(trans_loc, *this);
    int result = tdp->run();
    delete tdp;
    if (result==RESPONSE_OK)
      transform_specified(getparms);
   }
  else
    transform_specified(getparms);

  model->ComputeLogLikelihood(timeseries);
  update_all(true);
}

void CoreWin::help_about_cb()
{
  Dialog *about_win = new Dialog("About", *this);
  about_win->set_position(Gtk::WIN_POS_CENTER_ON_PARENT);
  Image *im = manage(new Image(Gdk::Pixbuf::create_from_xpm_data(newabout_xpm)));
  Frame *frm = manage(new Frame());
  frm->add(*im);
  about_win->get_vbox()->set_border_width(0);
  about_win->get_vbox()->pack_start(*frm);
  about_win->set_border_width(0);
  about_win->show_all();
}

void CoreWin::model_specified()
{
  // re-compute things
  model->ComputeLogLikelihood(timeseries);
  modelrenderer->Update(model,timeseries,&transformlist);
}

void CoreWin::model_specify_cb(int mtype)
{
  ARMAModel *tap;
  GARCHModel *gar;
  EGARCHModel *egar;
  SSVModel* ssvm;

  TimeSeriesModel *new_model;
  TSMInterface *new_interface;

  if (mtype!=modeltype)
    {
      switch (mtype)
	{
	case 1:
	  tap = new ARMAModel();
	  new_model = tap;
	  new_interface = new ARMAInterface(tap);
	  break;
	case 2:
	  gar = new GARCHModel();
	  new_model = gar;
	  new_interface = new GARCHInterface(gar);
	  break;
	case 3:
          ssvm = new SSVModel();
	  new_model = ssvm;
	  new_interface = new SSVMInterface(ssvm);
	  break;
	case 4:
	  egar = new EGARCHModel();
	  new_model = egar;
	  new_interface = new EGARCHInterface(egar);
	  break;
	}

      if (new_interface->Specify(this))  // is successful (i.e. not cancelled)
	{
	  delete model;
	  model = new_model;
	  interface = new_interface;
	  model_specified();
	  modeltype = mtype;
	}
    }
  else // same family, just change order, parameters, etc.
    if (interface->Specify(this))
      model_specified();
}

void CoreWin::model_fit_iteration(int stage)
{
  Glib::TimeVal nt, diff;

  nt.assign_current_time();
  diff = nt - ittime;

  itcount->increment(stage);

  if (diff.as_double()>0.1) // seconds
    {
      while (gtk_events_pending())
	gtk_main_iteration();
      ittime = nt;
    }
}

void CoreWin::model_fit_cb(int methodtype)
{
  ostringstream res, supplemental;

  if (timeseries->GetN()==0)
    {
      res << "Cannot fit model when there is no data." << ends;
      MessageDialog md(*this, res.str());
      md.run();
      return;
    }

  ittime.assign_current_time();
  itcount = new IterationCounter(*this);
  itcount->show_all();
  int i, result = model->FitModel(timeseries, methodtype, 1000,
				  &corewin_callback, this, res, supplemental);
  delete itcount;

  // now display result
  MessageDialog md(*this, res.str());
  md.add_button("Show Details", 17);
  i = md.run();

  if (result==TimeSeriesModel::SUCCESS)
    {
      string localstr = supplemental.str();
      // create new textview to display supplemental results
      if (i==17)
	if (localstr.size()>0)
	  {
	    Window *supp_win = new Window;
	    supp_win->set_default_size(600,400);
	    TextRenderer  *rw = manage(new TextRenderer());
	    rw->Update(localstr.c_str());
	    supp_win->add(*rw);
	    supp_win->show_all();
	  }

      update_all();
    }
}

void CoreWin::model_view_cts_version_cb()
{
  bool success;
  int i;
  string s1;
  double delta;
  // get delta
  GeneralTextDialog gtd(*this, 
			"Enter the value of Delta.\nDelta is the continuous-time model time difference corresponding to one step in discrete-time.",
			"0.00396");  // delta = cts time corresponding to one discrete time unit
  i=gtd.run();
  if (i==RESPONSE_OK)
    {
      s1 = gtd.get_text();
      delta = atof(s1.c_str());

      // now render
      ostringstream modelstring, errorstring;
      success = model->RenderCtsVersionInto(modelstring, errorstring, delta);
      if (success)
	{
	  // 1. build window to display conversion in
	  Window *supp_win = new Window;
	  supp_win->set_default_size(600,400);
	  TextRenderer  *rw = manage(new TextRenderer());
	  rw->Update(modelstring.str().c_str());
	  supp_win->add(*rw);
	  supp_win->show_all();
	}
      else
	{
	  MessageDialog ed(*this, errorstring.str());
	  ed.run();
	}
    }
}

void CoreWin::model_residuals_cb()
{
  // 1. compute the residuals
  TimeSeries resids;
  model->ComputeStdResiduals(timeseries, &resids);

  // 2. build window to put them in
  Window *res_win = new Window;
  res_win->set_title("Cronos: Residuals");
  res_win->set_default_size(800,600);
  ResidualWin *rw = manage(new ResidualWin(&resids, model, timeseries));
  res_win->add(*rw);
  res_win->show_all();
}

void CoreWin::model_acf_pacf_cb()
{
  // build window for model acf/pacf
  Window *mod_acf_win = new Window;
  mod_acf_win->set_default_size(640,400);
  ModelACFWin *mw = manage(new ModelACFWin(model));
  mod_acf_win->add(*mw);
  mod_acf_win->show_all();
}

void CoreWin::model_spectral_density()
{
  int i;
  const int res = 300;

  // 1. compute the spec. density
  Vector lambdas(res), results;
  for (i=0 ; i<res ; ++i)
    lambdas[i] = (i+1)*M_PI/res;

  ARMAModel *ap = (ARMAModel *)((void *) model);
  results = ap->ComputeSpectralDensity(lambdas);


  // 2. build window to put them in
  Window *pd_win = new Window;
  pd_win->set_default_size(640,400);
  string t = "Model Spectral Density";
  SpectrumWin *xyp = manage(new SpectrumWin(lambdas, results, t));
  pd_win->add(*xyp);
  pd_win->show_all();
}


void CoreWin::model_simulate()
{
  SimulateDialog sd(*this, timeseries->GetN());
  sd.show_all();
  if (sd.run()==RESPONSE_OK)
    {
      int n = sd.get_n();
      model->SimulateInto(timeseries, n, false);
      model->ComputeLogLikelihood(timeseries);
      update_all(true);
    }
  return;
}

void CoreWin::model_forecast_cb()
{
  // 0. come up with a reasonable forecast horizon
  int default_horizon = timeseries->GetN()/10;

  // 1. create dialog to get forecast parameters
  ForecastDialog fd(*this, default_horizon);
  fd.show_all();
  if (fd.run()==RESPONSE_OK)
    {
      int i,n,h = fd.get_n(), ns = fd.get_sims();

      // 1.5. create copy of transform list
      vector<Transformation *> *transform_list_copy
	= new vector<Transformation *>;
      Transformation *tp,*tcopy;
      n = transformlist.size();
      // iterate through and copy
      for (i=0 ; i<n ; ++i)
	{
	  tp = transformlist[i];
	  tcopy = tp->CreateCopy();
	  transform_list_copy->push_back(tcopy);
	}

      // 2. build window to put them in
      Window *f_win = new Window;
      f_win->set_default_size(800,600);
      f_win->set_title("Cronos: Forecasts");
      ForecastWin *fw = manage(new ForecastWin(*timeseries, *model, h,
					       transform_list_copy, ns));
      f_win->add(*fw);
      f_win->show_all();
    }
}
