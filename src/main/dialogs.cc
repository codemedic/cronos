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
 * This module contains definitions of
 * the dialog windows used in ostsap.
 *
 */

#include "includes.h"
#include "dialogs.h"

using namespace Gtk;
using namespace std;
using namespace SigC;

//----------------------------------------------------------------

CenteringDialog::CenteringDialog(Glib::ustring title, Window &wi)
  : Dialog(title, wi)
{
  cd_parent_win = &wi;
  set_position(Gtk::WIN_POS_CENTER_ON_PARENT);
}


//----------------------------------------------------------------


TransformDialog::TransformDialog(Transformation *t, Window& wi)
  : CenteringDialog("Transform", wi), trans(t)
{
  int i,np = trans->npars;

  add_button("Cancel", RESPONSE_CANCEL);
  add_button("O.K.", RESPONSE_OK);
  signal_response().connect(mem_fun(*this, &TransformDialog::terminate));

  // Then vbox area
  table = manage(new Table(np,2,true));
  es = new Gtk::Entry*[np];
  ls = new Gtk::Label*[np];
  for (i=0 ; i<np ; ++i)
    {
      ls[i] = manage(new Label(trans->parnames[i]));
      es[i] = manage(new Entry());
      es[i]->set_text(t->pardefaults[i]);
      es[i]->signal_activate().connect(mem_fun(*this, &TransformDialog::enter_pressed));
      table->attach(*ls[i],0,1,i,i+1,Gtk::EXPAND,Gtk::EXPAND,0,5);
      table->attach(*es[i],1,2,i,i+1,Gtk::EXPAND,Gtk::EXPAND,0,5);
    }

  get_vbox()->pack_start(*table);

  show_all();
}

void TransformDialog::enter_pressed()
{
  // we really should validate text entered here ...
  response(RESPONSE_OK);
}

void TransformDialog::terminate(int code)
{
  // extract parameters
  int i, np = trans->npars;
  string s1;
  for (i=0 ; i<np ; ++i)
    {
      // extract parameter i
      s1 = es[i]->get_text();
      istrstream inst(s1.c_str());
      switch (trans->partypes[i])
	{
	case 0: // integer
	  inst >> trans->ipars[i];
	  break;
	case 1: // real
	  inst >> trans->dpars[i];
	  break;
	case 2: // string
	  inst >> trans->spars[i];
	  break;
	}

      delete es[i];
      delete ls[i];
    }
  delete[] es;
  delete[] ls;
}

IterationCounter::IterationCounter(Window& wp)
  : CenteringDialog("Progress", wp)
{
  localcount = 0;
  Table *table = manage(new Table(3,2,true));
  pbar = manage(new ProgressBar());
  label1 = manage(new Label("# Iterations: "));
  label2 = manage(new Label("0"));
  stage_label = manage(new Label());
  table->attach(*stage_label, 0, 2, 0, 1, EXPAND, EXPAND, 10, 2);
  table->attach(*pbar, 0, 2, 2, 3, EXPAND, EXPAND, 10, 2);
  table->attach(*label1, 0, 1, 1, 2, EXPAND, EXPAND, 10, 2);
  table->attach(*label2, 1, 2, 1, 2, EXPAND, EXPAND, 10, 2);
  curstage = -1; // no stage to begin with
  get_vbox()->pack_start(*table);
}

int IterationCounter::increment(int stage)
{
  char temps[80];
  ostrstream os(temps, 80);

  ++localcount;
  os << localcount << ends;
  label2->set_label(temps);

  if (stage!=curstage)
    {
      curstage = stage;
      os.seekp(0);
      switch (curstage) {
      case 0: os << "<b>Likelihood Maximization</b>" << ends;   break;
      case 1: os << "<b>Hessian Matrix Calculation</b>" << ends;  break;
      case 2: os << "<b>Markov Chain Monte Carlo</b>" << ends;  break;
      default: os << "<b>Bad Stage</b>" << ends; break;
      };
      stage_label->set_markup(temps);
    }

  pbar->pulse();
  return localcount;
}

//-------------------------------------------------------------

CoefficientDialog::CoefficientDialog(Window& wi, int n,
				     string *labs, double *vals, double *mask)
  : CenteringDialog("Coefficient Entry", wi), nentries(n)
{
  int i;
  char temps[30];

  ScrolledWindow *sw = manage(new ScrolledWindow());
  sw->set_size_request(300,min(30*n+10,160));
  sw->set_policy(POLICY_AUTOMATIC, POLICY_AUTOMATIC);

  Table *table = manage(new Table(nentries,1,false));

  Label *lp;
  HBox *vb;
  eps = new Entry*[nentries];
  checks = new CheckButton*[nentries];
  for (i=0 ; i<nentries ; ++i)
    {
      vb = manage(new HBox());

      lp = manage(new Label(labs[i]));
      vb->pack_start(*lp);

      eps[i] = manage(new Entry());
      ostrstream os(temps, 30);
      os << vals[i] << ends;
      eps[i]->set_text(temps);
      vb->pack_start(*eps[i]);

      lp = manage(new Label("  Hold:"));
      vb->pack_start(*lp);

      checks[i] = manage(new CheckButton());
      checks[i]->set_active(mask != NULL ? !mask[i] : 0);
      vb->pack_start(*checks[i]);

      table->attach(*vb, 0, 1, i, i+1, EXPAND, EXPAND);
    }

  sw->add(*table);
  get_vbox()->pack_start(*sw);
  add_button("OK",RESPONSE_OK);
  add_button("Cancel",RESPONSE_CANCEL);

  show_all_children();
}

CoefficientDialog::~CoefficientDialog()
{
  delete[] eps;
  delete[] checks;
}

void CoefficientDialog::on_response(int rid)
{
  string s1;
  if (rid==RESPONSE_OK)
    {
      retvals.resize(nentries);
      maskvals.resize(nentries);
      for (int i=0 ; i<nentries ; ++i)
	{
	  s1 = eps[i]->get_text();
	  istrstream inst(s1.c_str());
	  inst >> retvals[i];
	  maskvals[i] = !checks[i]->get_active();
	}
    }
}

//-------------------------------------------------------------------

SimulateDialog::SimulateDialog(Window& wi, int defaultn)
  : CenteringDialog("Simulation", wi)
{
  int i,nn = defaultn==0 ? 250 : defaultn;
  char temps[30];

  // action area
  add_button("Cancel", RESPONSE_CANCEL);
  add_button("O.K.", RESPONSE_OK);

  // Then vbox area
  table = manage(new Table(1,2,false));
  table->attach(*(manage(new Label("# Obs"))), 0, 1, 0, 1, EXPAND, EXPAND, 10, 2);
  nadjust = manage(new Adjustment(nn,0,100000,10,10,10));
  nspin = manage(new SpinButton(*nadjust,0,0));
  table->attach(*nspin, 0,1,1,2, EXPAND, EXPAND, 10, 2);
  get_vbox()->pack_start(*table);

  show_all_children();
}

int SimulateDialog::get_n()
{
  return nspin->get_value_as_int();
}

//-------------------------------------------------------------------

ForecastDialog::ForecastDialog(Window& wi, int defaultn)
  : CenteringDialog("Forecast", wi)
{
  int i,nn = defaultn==0 ? 10 : defaultn;
  char temps[30];

  // action area
  add_button("Cancel", RESPONSE_CANCEL);
  add_button("O.K.", RESPONSE_OK);

  // Then vbox area
  table = manage(new Table(2,3,false));
  table->attach(*(manage(new Label("Forecast Horizon"))), 0,1,1,2, EXPAND, EXPAND, 10, 2);
  table->attach(*(manage(new Label("Number of Trajectories"))), 1,2,1,2, EXPAND, EXPAND, 10, 2);

  nadjust = manage(new Adjustment(nn,0,100000,1,10,10));
  nspin = manage(new SpinButton(*nadjust,0,0));

  simsadjust = manage(new Adjustment(2500, 100, 100000, 100, 10, 10));
  simsspin = manage(new SpinButton(*simsadjust,0,0));

  table->attach(*nspin, 0,1,2,3, EXPAND, EXPAND, 10, 2);
  table->attach(*simsspin, 1,2,2,3, EXPAND, EXPAND, 10, 2);
  get_vbox()->pack_start(*table);

  Label *label = manage(new Label("Forecasts are based on simulation of (future) trajectories.\nSome Monte-Carlo simulation error can be expected in predictive distributions."));
  table->attach(*label, 0,2,0,1, EXPAND, EXPAND);
  
  show_all_children();
}

int ForecastDialog::get_n()
{
  return nspin->get_value_as_int();
}

int ForecastDialog::get_sims()
{
  return simsspin->get_value_as_int();
}

//-------------------------------------------------------------------

GeneralTextDialog::GeneralTextDialog(Window& wi, string label, string thetext)
  : CenteringDialog("", wi)
{
  int i;

  // action area
  add_button("Cancel", RESPONSE_CANCEL);
  add_button("O.K.", RESPONSE_OK);

  // Then vbox area
  table = manage(new Table(1,2,false));
  table->attach(*(manage(new Label(label.c_str()))), 
		0, 1, 0, 1, EXPAND, EXPAND, 10, 2);
  entry = manage(new Gtk::Entry());
  entry->set_text(thetext.c_str());
  table->attach(*entry, 0,1,1,2, EXPAND, EXPAND, 10, 2);
  get_vbox()->pack_start(*table);

  show_all_children();
}

string GeneralTextDialog::get_text()
{
  string ts = entry->get_text();
  return ts;
}
