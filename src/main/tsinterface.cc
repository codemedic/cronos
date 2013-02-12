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
 * Time Series Model Interfaces:
 *
 * This module contains dialogs which allow the specification
 * of various time series models.  For the sake of making
 * TimeSeriesModel and TimeSeries classes self-contained and not
 * dependent on a particular GUI, the code for the interface is
 * contained here rather than being included in timeseries.cc/h.
 *
 */

#include "includes.h"
#include "tsinterface.h"
#include <sigc++/bind.h>

using namespace SigC;

ARMAModelDialog::ARMAModelDialog(ARMAModel *a, Window& wi)
  : CenteringDialog("ARMA Model Parameters",  wi), amp(a)
{
  int i;
  localphi.resize(100);   localtheta.resize(100);
  localphi.zeroes();
  localtheta.zeroes();
  localmask1.resize(3);  // for mu,d,sigma
  localmask2.resize(50); // for phis
  localmask3.resize(50); // for thetas

  // create 4 default masks which use current mask if available
  Vector tmask = amp->GetMask();
  for (i=0 ; i<3 ; ++i)
    localmask1[i] = tmask[i];
  for (i=3 ; i<3+amp->GetP() ; ++i)
    localmask2[i-3] = tmask[i];
  for ( ; i<53 ; ++i)
    localmask2[i-3] = 1;
  for (i=3+amp->GetP() ; i<3+amp->GetP()+amp->GetQ() ; ++i)
    localmask3[i-3-amp->GetP()] = tmask[i];
  for ( ; i<53+amp->GetP() ; ++i)
    localmask3[i-3-amp->GetP()] = 1;

  Vector tphi, ttheta;
  amp->GetCoefficients(tphi, ttheta);
  for (i=0 ; i<amp->GetP() ; ++i)
    localphi[i] = tphi[i];
  for (i=0 ; i<amp->GetQ() ; ++i)
    localtheta[i] = ttheta[i];

  BaseConstruction();
  show_all();
}

void ARMAModelDialog::EnterCoefficients(int id)
{
  int i;
  string *sp;
  if (id==1) // phi coefficients
    {
      int localn = arspin->get_value_as_int();
      if (localn==0)
	return;
      sp = new string[localn];
      for (i=0 ; i<localn ; ++i)
	{
	  ostrstream os;
	  os << "Phi(" << (i+1) << "): " << ends;
	  sp[i] = os.str();
	  os.rdbuf()->freeze(0);
	}
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window *) wp;
      CoefficientDialog cd(*winp, localn, sp, &localphi[0], &localmask2[0]);
      cd.show_all();
      if (cd.run()==RESPONSE_OK)
	for (i=0 ; i<localn ; ++i)
	  {
	    localphi[i] = cd.retvals[i];
	    localmask2[i] = cd.maskvals[i];
	  }
	
      delete[] sp;
    }
  if (id==2) // thetas
    {
      int localn = maspin->get_value_as_int();
      if (localn==0)
	return;
      sp = new string[localn];
      for (i=0 ; i<localn ; ++i)
	{
	  ostrstream os;
	  os << "Theta(" << (i+1) << "): " << ends;
	  sp[i] = os.str();
	  os.rdbuf()->freeze(0);
	}
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window *) wp;
      CoefficientDialog cd(*winp, localn, sp, &localtheta[0],&localmask3[0]);
      cd.show_all();
      if (cd.run()==RESPONSE_OK)
	for (i=0 ; i<localn ; ++i)
	  {
	    localtheta[i] = cd.retvals[i];
	    localmask3[i] = cd.maskvals[i];
	  }

      delete[] sp;
    }
}

void ARMAModelDialog::BaseConstruction()
{
  // Action area first
  add_button("Cancel", RESPONSE_CANCEL);
  add_button("OK", RESPONSE_OK);
  signal_response().connect(mem_fun(*this, &ARMAModelDialog::terminate));

  // Then vbox area
  table = manage(new Table(2,3,true));
  // buttons for coefficient entry
  Button *b1 = manage(new Button("Phi(j)s")),
    *b2 = manage(new Button("Theta(j)s"));

  // AR stuff
  VBox *arvbox = manage(new VBox());
  arvbox->pack_start(*(manage(new Label("AR Order"))));
  aradjust = manage(new Adjustment(amp->GetP(),0,100,1,10,10));
  arspin = manage(new SpinButton(*aradjust,0,0));
  arvbox->pack_start(*arspin);
  b1->signal_clicked().connect(bind<int>(mem_fun(*this,&ARMAModelDialog::EnterCoefficients),1));
  arvbox->pack_start(*b1);
  table->attach(*arvbox, 0, 1, 0, 1, EXPAND, EXPAND, 10, 2);

  // MA stuff
  VBox *mavbox = manage(new VBox());
  mavbox->pack_start(*(manage(new Label("MA Order"))));
  b2->signal_clicked().connect(bind<int>(mem_fun(*this,&ARMAModelDialog::EnterCoefficients),2));
  maadjust = manage(new Adjustment(amp->GetQ(),0,100,1,10,10));
  maspin = manage(new SpinButton(*maadjust,0,0));
  mavbox->pack_start(*maspin);
  mavbox->pack_start(*b2);
  table->attach(*mavbox, 1, 2, 0, 1, EXPAND, EXPAND, 10, 2);

  // Whittle or not
  VBox *whittlebox = manage(new VBox());
  //whittlebox->pack_start(*(manage(new Label("Whittle Approx."))));
  whittlecheck = manage(new CheckButton("Use Whittle Approx.", false));
  whittlebox->pack_start(*whittlecheck);
  table->attach(*whittlebox, 2,3,0,1, EXPAND, EXPAND, 10, 2);

  // frac diff stuff
  VBox *fdvbox = manage(new VBox());
  fdvbox->pack_start(*(manage(new Label("Frac. Diff"))));
  fracadjust = manage(new Adjustment(amp->GetFracDiff(),-0.5,0.5,0.01,10,10));
  fracspin = manage(new SpinButton(*fracadjust,0,0));
  fracspin->set_digits(3);
  fraccheck = manage(new CheckButton("Hold fixed", false));
  fdvbox->pack_start(*fracspin);
  fdvbox->pack_start(*fraccheck);
  table->attach(*fdvbox, 2,3,1,2, EXPAND, EXPAND, 10, 2);

  // mean
  VBox *mubox = manage(new VBox());
  muadjust = manage(new Adjustment(amp->GetMean(),-100000,100000,0.01,10,10));
  muspin = manage(new SpinButton(*muadjust,0,0));
  muspin->set_digits(4);
  mucheck = manage(new CheckButton("Hold fixed", false));
  mubox->pack_start(*(manage(new Label("Mean"))));
  mubox->pack_start(*muspin);
  mubox->pack_start(*mucheck);
  table->attach(*mubox, 0, 1, 1, 2, EXPAND, EXPAND, 10, 2);

  // sigma
  VBox *sigbox = manage(new VBox());
  sigmaadjust = manage(new Adjustment(amp->GetSigma(),0,100000,0.01,10,10));
  sigmaspin = manage(new SpinButton(*sigmaadjust,0,0));
  sigmaspin->set_digits(4);  
  sigmacheck = manage(new CheckButton("Hold fixed", false));
  sigbox->pack_start(*(manage(new Label("Sigma"))));
  sigbox->pack_start(*sigmaspin);
  sigbox->pack_start(*sigmacheck);
  table->attach(*sigbox, 1, 2, 1, 2, EXPAND, EXPAND, 10, 2);
  
  // set checkbutton toggle states
  Vector mask = amp->GetMask();
  mucheck->set_active(!localmask1[0]);
  fraccheck->set_active(!localmask1[1]);
  sigmacheck->set_active(!localmask1[2]);
  whittlecheck->set_active(amp->GetWhittleMode());

  get_vbox()->pack_start(*table);
}



void ARMAModelDialog::terminate(int code)
{
  int i;
  if (code==RESPONSE_OK)
    { // was OK
      int newp=arspin->get_value_as_int(), newq=maspin->get_value_as_int(), np=newp+newq+3;
      amp->SetOrder(newp, newq);
      amp->SetFracDiff(fracspin->get_value());
      amp->SetCoefficients(&localphi[0], &localtheta[0]);
      amp->SetSigma2(sigmaspin->get_value()*sigmaspin->get_value());
      amp->SetMean(muspin->get_value());

      // now set constraints
      Vector mask(np);
      for (i=0 ; i<np ; ++i)
	mask[i] = 1;
      mask[0] = !mucheck->get_active();
      mask[1] = !fraccheck->get_active();
      mask[2] = !sigmacheck->get_active();
      for (i=0 ; i<newp ; ++i)
	mask[i+3] = localmask2[i];
      for (i=0 ; i<newq ; ++i)
	mask[i+3+newp] = localmask3[i];
      amp->SetMask(mask);

      amp->SetWhittleMode(whittlecheck->get_active());
    }
  else
    { // was cancelled
    }
}

int ARMAModelDialog::run()
  // we override the run() function to add checks for model validity
{
  int retval;
  retval = CenteringDialog::run();
  if (retval==RESPONSE_OK)
    {
      // do a causality/invertibility check
      if (!amp->IsCausal() || !amp->IsInvertible()) // fix if necessary
	{
	  amp->IsCausal(true);
	  amp->IsInvertible(true);
	  // notify that model is changed
	  char temps[160]="The specified model was non-causal and/or non-invertible.  It has been altered by reflecting polynomial roots outside the unit circle.";
	  MessageDialog md(*this, temps);
	  md.run();
	}
    }
  return retval;
}


bool ARMAInterface::Specify(Gtk::Window *mainwin)
{
  bool retval;
  ARMAModelDialog *ad;
  ad = new ARMAModelDialog(mp, *mainwin);
  retval =  (ad->run()==RESPONSE_OK);
  delete ad;
  return retval;
}

// +---------------------------------+
// |          GARCH STUFF            |
// +---------------------------------+

GARCHModelDialog::GARCHModelDialog(GARCHModel *a, Window& wi)
  : CenteringDialog("GARCH Model Parameters",  wi), gmp(a)
{
  int i,p,q;
  double mu;
  Vector as,bs;
  gmp->GetParameters(p,q,as,bs,mu);
  localas.resize(100);   localbs.resize(100);
  localas.zeroes();
  localbs.zeroes();

  // get default mask and copy to local stuff
  Vector tmask = a->GetMask();
  localmask1.resize(1);  // mean
  localmask2.resize(100);  // as
  localmask3.resize(100);  // bs
  for (i=0 ; i<100 ; ++i)
    { localmask2[i]=1;   localmask3[i]=1; }
  localmask1[0] = tmask[0];
  for (i=0 ; i<=p ; ++i)
    localmask2[i] = tmask[i+1];
  for (i=0 ; i<q ; ++i)
    localmask3[i] = tmask[i+p+2];

  // copy current param values to local
  for (i=0 ; i<=p ; ++i)
    localas[i] = as[i];
  for (i=0 ; i<q ; ++i)
    localbs[i] = bs[i];
  BaseConstruction();
  show_all();
}

void GARCHModelDialog::EnterCoefficients(int id)
{
  int i;
  string *sp;
  if (id==1) // a coefficients
    {
      int localn = aspin->get_value_as_int();
      if (localn<0)
	return;
      sp = new string[localn+1];
      for (i=0 ; i<=localn ; ++i)
	{
	  ostrstream os;
	  os << "a(" << i << ")" << ends;
	  sp[i] = os.str();
	  os.rdbuf()->freeze(0);
	}
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window *) wp;
      CoefficientDialog cd(*winp, localn+1, sp, &localas[0], &localmask2[0]);
      cd.show_all();
      if (cd.run()==RESPONSE_OK)
	{
          for (i=0 ; i<=localn ; ++i)
	    {
	      localas[i] = cd.retvals[i];
	      localmask2[i] = cd.maskvals[i];
	    }
	}
      delete[] sp;
    }
  if (id==2) // thetas
    {
      int localn = bspin->get_value_as_int();
      if (localn==0)
	return;
      sp = new string[localn];
      for (i=0 ; i<localn ; ++i)
	{
	  ostrstream os;
	  os << "b(" << (i+1) << ")" << ends;
	  sp[i] = os.str();
	  os.rdbuf()->freeze(0);
	}
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window *) wp;
      CoefficientDialog cd(*winp, localn, sp, &localbs[0], &localmask3[0]);
      cd.show_all();
      if (cd.run()==RESPONSE_OK)
	{
          for (i=0 ; i<localn ; ++i)
	    {
	      localbs[i] = cd.retvals[i];
	      localmask3[i] = cd.maskvals[i];
	    }
	}
      delete[] sp;
    }
}

void GARCHModelDialog::BaseConstruction()
{
  int i,p,q;
  double mu;
  Vector as,bs;
  gmp->GetParameters(p,q,as,bs,mu);

  // Action area first
  add_button("Cancel", RESPONSE_CANCEL);
  add_button("OK", RESPONSE_OK);
  signal_response().connect(mem_fun(*this, &GARCHModelDialog::terminate));

  // Then vbox area
  table = manage(new Table(1,3,true));

  VBox *vbox1 = manage(new VBox());
  vbox1->pack_start(*manage(new Label("P")));
  aadjust = manage(new Adjustment(p,0,100,1,10,10));
  aspin = manage(new SpinButton(*aadjust,0,0));
  vbox1->pack_start(*aspin);
  Button *b1 = manage(new Button("Alpha coeffs."));
  b1->signal_clicked().connect(bind<int>(mem_fun(*this,&GARCHModelDialog::EnterCoefficients),1));
  vbox1->pack_start(*b1);
  table->attach(*vbox1, 1, 2, 0, 1, Gtk::EXPAND,Gtk::EXPAND,10,2);

  VBox *vbox2 = manage(new VBox());
  vbox2->pack_start(*manage(new Label("Q")));
  badjust = manage(new Adjustment(q,0,100,1,10,10));
  bspin = manage(new SpinButton(*badjust,0,0));
  vbox2->pack_start(*bspin);
  Button *b2 = manage(new Button("Beta coeffs."));
  b2->signal_clicked().connect(bind<int>(mem_fun(*this,&GARCHModelDialog::EnterCoefficients),2));
  vbox2->pack_start(*b2);
  table->attach(*vbox2, 2, 3, 0, 1, Gtk::EXPAND,Gtk::EXPAND,10,2);

  VBox *vbox3 = manage(new VBox());
  vbox3->pack_start(*(manage(new Label("Mean"))));
  muadjust = manage(new Adjustment(mu,-100000,100000,0.01,10,10));
  muspin = manage(new SpinButton(*muadjust,0,0));
  muspin->set_digits(6);
  vbox3->pack_start(*muspin);
  table->attach(*vbox3, 0, 1, 0, 1, EXPAND, EXPAND, 10, 2);

  get_vbox()->pack_start(*table);
}


void GARCHModelDialog::terminate(int code)
{
  int i;
  if (code==RESPONSE_OK)
    { 
      // was OK
      int newp=aspin->get_value_as_int(), newq=bspin->get_value_as_int();
      double newmu=muspin->get_value();
      gmp->SetParameters(newp,newq,localas,localbs,newmu);

      // now set constraints
      int np = newp+newq+2;
      Vector mask(np);
      for (i=0 ; i<np ; ++i)
	mask[i] = 1;
      mask[0] = localmask1[0];
      for (i=1 ; i<=newp+1 ; ++i)
	mask[i] = localmask2[i-1];
      for (i=0 ; i<newq ; ++i)
	mask[i+newp+2] = localmask3[i];
      gmp->SetMask(mask);
    }
  else
    { // was cancelled
    }
}


int GARCHModelDialog::run()
  // we override the run() function to add checks for model validity
{
  int retval;
  retval = CenteringDialog::run();
  if (retval==RESPONSE_OK)
    {
      // do a validity check
      if (!gmp->CheckModel(false))
	{
	  gmp->CheckModel(true);
	  // notify that model is changed
	  char temps[160]="The specified GARCH model had invalid parameters.  It has been altered to a model with valid parameters.";
	  MessageDialog md(*this, temps);
	  md.run();
	}
    }
  return retval;
}


bool GARCHInterface::Specify(Gtk::Window *mainwin)
{
  bool retval;
  GARCHModelDialog *ad;
  ad = new GARCHModelDialog(mp, *mainwin);
  retval =  (ad->run()==RESPONSE_OK);
  delete ad;
  return retval;
}

// +---------------------------------+
// |          EGARCH STUFF           |
// +---------------------------------+

EGARCHModelDialog::EGARCHModelDialog(EGARCHModel *a, Window& wi)
  : CenteringDialog("EGARCH Model Parameters",  wi), egmp(a)
{
  int i,p,q;
  double mu;
  Vector as,bs,gs;

  egmp->GetParameters(p,q,as,gs,bs,mu);
  localas.resize(100);   localbs.resize(100);   localgs.resize(100);
  localas.zeroes();   localbs.zeroes();   localgs.zeroes();

  // get default mask and copy to local stuff
  Vector tmask = a->GetMask();
  localmask1.resize(1);  // mean
  localmask2.resize(100);  // as
  localmask3.resize(100);  // bs
  localmask4.resize(100);  // gs
  for (i=0 ; i<100 ; ++i)
    { localmask2[i]=localmask3[i]=localmask4[i]=1; }
  localmask1[0] = tmask[0];
  for (i=0 ; i<=p ; ++i)
    localmask2[i] = tmask[i+1];
  for (i=p+2 ; i<2*p+2 ; ++i)
    localmask4[i-p-2] = tmask[i];
  for (i=0 ; i<q ; ++i)
    localmask3[i] = tmask[i+2*p+2];

  // copy current param values to local
  localas[0] = as[0];
  for (i=1 ; i<=p ; ++i)
    {
      localas[i] = as[i];
      localgs[i] = gs[i];
    }
  for (i=0 ; i<q ; ++i)
    localbs[i] = bs[i];

  BaseConstruction();
  show_all();
}

void EGARCHModelDialog::EnterCoefficients(int id)
{
  int i;
  string *sp;
  if (id==1) // alpha coefficients
    {
      int localn = aspin->get_value_as_int();
      if (localn<0)
	return;
      sp = new string[localn+1];
      for (i=0 ; i<=localn ; ++i)
	{
	  ostrstream os;
	  os << "a(" << i << ")" << ends;
	  sp[i] = os.str();
	  os.rdbuf()->freeze(0);
	}
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window *) wp;
      CoefficientDialog cd(*winp, localn+1, sp, &localas[0], &localmask2[0]);
      cd.show_all();
      if (cd.run()==RESPONSE_OK)
	{
          for (i=0 ; i<=localn ; ++i)
	    {	    
	      localas[i] = cd.retvals[i];
	      localmask2[i] = cd.maskvals[i];
	    }
	    
	}
      delete[] sp;
    }
  if (id==3) // gamma coefficients
    {
      int localn = aspin->get_value_as_int();
      if (localn==0)
	return;
      sp = new string[localn+1];
      for (i=1 ; i<=localn ; ++i)
	{
	  ostrstream os;
	  os << "g(" << i << ")" << ends;
	  sp[i] = os.str();
	  os.rdbuf()->freeze(0);
	}
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window *) wp;
      CoefficientDialog cd(*winp, localn, &sp[1], &localgs[1], &localmask4[0]);
      cd.show_all();
      if (cd.run()==RESPONSE_OK)
	{
          for (i=1 ; i<=localn ; ++i)	    
	    {
	      localgs[i] = cd.retvals[i-1];
	      localmask4[i-1] = cd.maskvals[i-1];
	    }
	}
      delete[] sp;
    }
  if (id==2) // beta coefficients
    {
      int localn = bspin->get_value_as_int();
      if (localn==0)
	return;
      sp = new string[localn];
      for (i=0 ; i<localn ; ++i)
	{
	  ostrstream os;
	  os << "b(" << (i+1) << ")" << ends;
	  sp[i] = os.str();
	  os.rdbuf()->freeze(0);
	}
      Widget *wp = get_toplevel();
      Gtk::Window *winp = (Gtk::Window *) wp;
      CoefficientDialog cd(*winp, localn, sp, &localbs[0], &localmask3[0]);
      cd.show_all();
      if (cd.run()==RESPONSE_OK)
	{
          for (i=0 ; i<localn ; ++i)
	    {
	      localbs[i] = cd.retvals[i];
	      localmask3[i] = cd.maskvals[i];
	    }
	}
      delete[] sp;
    }
}

void EGARCHModelDialog::BaseConstruction()
{
  int i,p,q;
  double mu;
  Vector as,bs,gs;
  egmp->GetParameters(p,q,as,gs,bs,mu);

  // Action area first
  add_button("Cancel", RESPONSE_CANCEL);
  add_button("OK", RESPONSE_OK);
  signal_response().connect(mem_fun(*this, &EGARCHModelDialog::terminate));

  // Then vbox areas
  table = manage(new Table(1,3,true));

  // vbox for P and associated alpha and gamma coefficients
  VBox *vbox1 = manage(new VBox());
  vbox1->pack_start(*manage(new Label("P")));
  aadjust = manage(new Adjustment(p,0,100,1,10,10));
  aspin = manage(new SpinButton(*aadjust,0,0));
  vbox1->pack_start(*aspin);

  Button *b1 = manage(new Button("Alpha coeffs."));
  b1->signal_clicked().connect(bind<int>(mem_fun(*this,&EGARCHModelDialog::EnterCoefficients),1));
  vbox1->pack_start(*b1);

  Button *b4 = manage(new Button("Gamma coeffs."));
  b4->signal_clicked().connect(bind<int>(mem_fun(*this,&EGARCHModelDialog::EnterCoefficients),3));
  vbox1->pack_start(*b4);

  table->attach(*vbox1, 1, 2, 0, 1, Gtk::EXPAND,Gtk::EXPAND,10,2);

  
  // vbox for Q and associated beta coefficients
  VBox *vbox2 = manage(new VBox());
  vbox2->pack_start(*manage(new Label("Q")));
  badjust = manage(new Adjustment(q,0,100,1,10,10));
  bspin = manage(new SpinButton(*badjust,0,0));
  vbox2->pack_start(*bspin);
  Button *b2 = manage(new Button("Beta coeffs."));
  b2->signal_clicked().connect(bind<int>(mem_fun(*this,&EGARCHModelDialog::EnterCoefficients),2));
  vbox2->pack_start(*b2);
  table->attach(*vbox2, 2, 3, 0, 1, Gtk::EXPAND,Gtk::EXPAND,10,2);

  VBox *vbox3 = manage(new VBox());
  vbox3->pack_start(*(manage(new Label("Mean"))));
  muadjust = manage(new Adjustment(mu,-100000,100000,0.01,10,10));
  muspin = manage(new SpinButton(*muadjust,0,0));
  muspin->set_digits(6);
  vbox3->pack_start(*muspin);
  table->attach(*vbox3, 0, 1, 0, 1, EXPAND, EXPAND, 10, 2);

  get_vbox()->pack_start(*table);
}


void EGARCHModelDialog::terminate(int code)
{
  int i;
  if (code==RESPONSE_OK)
    { // was OK
      int newp=aspin->get_value_as_int(), newq=bspin->get_value_as_int();
      double newmu=muspin->get_value();
      egmp->SetParameters(newp,newq,localas,localgs,localbs,newmu);

      // now set constraints
      int np = 2*newp+newq+2;
      Vector mask(np);
      for (i=0 ; i<np ; ++i)
	mask[i] = 1;
      mask[0] = localmask1[0];
      for (i=1 ; i<=newp+1 ; ++i)
	mask[i] = localmask2[i-1];
      for (i=newp+2 ; i<2*newp+2 ; ++i)
	mask[i] = localmask4[i-newp-2];
      for (i=0 ; i<newq ; ++i)
	mask[i+2*newp+2] = localmask3[i];
      egmp->SetMask(mask);
    }
  else
    { // was cancelled
    }
}

int EGARCHModelDialog::run()
  // we override the run() function to add checks for model validity
{
  int retval;
  retval = CenteringDialog::run();
  if (retval==RESPONSE_OK)
    {
      // do a validity check
      if (!egmp->CheckModel(false))
	{
	  egmp->CheckModel(true);
	  // notify that model is changed
	  char temps[160]="The specified EGARCH model had invalid parameters.  It has been altered to a model with valid parameters.";
	  MessageDialog md(*this, temps);
	  md.run();
	}
    }
  return retval;
}

bool EGARCHInterface::Specify(Gtk::Window *mainwin)
{
  bool retval;
  EGARCHModelDialog *ad;
  ad = new EGARCHModelDialog(mp, *mainwin);
  retval =  (ad->run()==RESPONSE_OK);
  delete ad;
  return retval;
}

// +---------------------------------+
// |          SVM STUFF              |
// +---------------------------------+


SSVMDialog::SSVMDialog(SSVModel *mp, Window& wi)
  : CenteringDialog("SSVM Model",  wi), svmp(mp)
{
  localmask = svmp->GetMask();
  BaseConstruction();
  show_all();
}

void SSVMDialog::BaseConstruction()
{
  int i,p,q;
  double mu;
  Vector parms = svmp->ParameterBundle();


  // Action area first
  add_button("Cancel", RESPONSE_CANCEL);
  add_button("OK", RESPONSE_OK);
  signal_response().connect(mem_fun(*this, &SSVMDialog::terminate));

  // Then vbox area
  Table *table = manage(new Table(1,2,true));

  VBox *vbox1 = manage(new VBox());
  vbox1->pack_start(*(manage(new Label("Mean"))));
  meanadjust = manage(new Adjustment(parms[0], -1000,1000,0.01,10,10));
  meanspin = manage(new SpinButton(*meanadjust,0,0));
  meanspin->set_digits(4);
  meancheck = manage(new CheckButton("Hold fixed", false));
  vbox1->pack_start(*meanspin);
  vbox1->pack_start(*meancheck);


  vbox1->pack_start(*(manage(new Label("Mean of Log-Vol."))));
  muxadjust = manage(new Adjustment(parms[1], -1000,1000,0.01,10,10));
  muxspin = manage(new SpinButton(*muxadjust,0,0));
  muxspin->set_digits(4);
  muxcheck = manage(new CheckButton("Hold fixed", false));
  vbox1->pack_start(*muxspin);
  vbox1->pack_start(*muxcheck);

  VBox *vbox2 = manage(new VBox());
  vbox2->pack_start(*(manage(new Label("Phi"))));
  phiadjust = manage(new Adjustment(parms[2], -1, 1, 0.01, 10, 10));
  phispin = manage(new SpinButton(*phiadjust,0,0));
  phispin->set_digits(3);
  phicheck = manage(new CheckButton("Hold fixed", false));
  vbox2->pack_start(*phispin);
  vbox2->pack_start(*phicheck);

  vbox2->pack_start(*(manage(new Label("Nu"))));
  nuadjust = manage(new Adjustment(parms[3], -1000,1000,0.01, 10, 10));
  nuspin = manage(new SpinButton(*nuadjust,0,0));
  nuspin->set_digits(3);
  nucheck = manage(new CheckButton("Hold fixed", false));
  vbox2->pack_start(*nuspin);
  vbox2->pack_start(*nucheck);
		  

  table->attach(*vbox1, 0, 1, 0, 1, EXPAND, EXPAND, 10, 2);
  table->attach(*vbox2, 1, 2, 0, 1, EXPAND, EXPAND, 10, 2);

  meancheck->set_active(!localmask[0]);
  muxcheck->set_active(!localmask[1]);
  phicheck->set_active(!localmask[2]);
  nucheck->set_active(!localmask[3]);

  get_vbox()->pack_start(*table);
}

void SSVMDialog::terminate(int code)
{
  // if code was OK, copy values
  if (code==RESPONSE_OK)
    {
      Vector newparms(4);
      newparms[0] = meanspin->get_value();
      newparms[1] = muxspin->get_value();
      newparms[2] = phispin->get_value();
      newparms[3] = nuspin->get_value();

      svmp->UnbundleParameters(newparms);

      // now set constraints
      localmask[0] = !meancheck->get_active();
      localmask[1] = !muxcheck->get_active();
      localmask[2] = !phicheck->get_active();
      localmask[3] = !nucheck->get_active();
      svmp->SetMask(localmask);
    }
}

bool SSVMInterface::Specify(Gtk::Window *mainwin)
{
  bool retval = true;
  SSVMDialog *svmdp;
  svmdp = new SSVMDialog(svmmp, *mainwin);
  retval =  (svmdp->run()==RESPONSE_OK);
  delete svmdp;
  return retval;
}

int SSVMDialog::run()
  // we override the run() function to add checks for model validity
{
  int retval;
  retval = CenteringDialog::run();
  if (retval==RESPONSE_OK)
    {
      // do a causality/invertibility check
      if (!svmp->CheckModel(false)) // fix if necessary
	{
	  char temps[160]="The specified model had invalid parameter values.  It has been altered to make parameters valid.";
	  MessageDialog md(*this, temps);
	  md.run();
	  svmp->CheckModel(true);
	}
    }
  return retval;
}
