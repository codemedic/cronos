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

#ifndef tsinterface_h
#define tsinterface_h

#include "timeseries.h"
#include "garch.h"
#include "svm.h"

class TSMInterface
{
 protected:
  TimeSeriesModel *tsmp;

 public:
  TSMInterface() {}
  virtual bool Specify(Gtk::Window *mainwin) = 0;
};

class ARMAModelDialog : public CenteringDialog {
protected:
  ARMAModel *amp;

  Gtk::Button *ok_button, *cancel_button;
  Gtk::Table *table;
  Gtk::SpinButton *arspin, *maspin, *fracspin, *sigmaspin, *muspin;
  Gtk::Adjustment *aradjust, *maadjust, *fracadjust, *sigmaadjust,
    *muadjust;
  Gtk::CheckButton *mucheck, *fraccheck, *sigmacheck, *whittlecheck;

  void terminate(int);
  void BaseConstruction();
  void EnterCoefficients(int id);

  Vector localphi, localtheta, localmask1,localmask2,localmask3;

public:
  ARMAModelDialog(ARMAModel *am, Window& w); // centered over w
  int run();  // local version
};

class ARMAInterface : public TSMInterface
{
 protected:
  ARMAModel *mp;
 public:
  ARMAInterface(ARMAModel *amp) :
    TSMInterface(), mp(amp) {}

  bool Specify(Gtk::Window *mainwin);
};



class GARCHModelDialog : public CenteringDialog {
protected:
  GARCHModel *gmp;

  Gtk::Button *ok_button, *cancel_button;
  Gtk::Table *table;
  Gtk::SpinButton *aspin, *bspin, *muspin;
  Gtk::Adjustment *aadjust, *badjust, *muadjust;

  void terminate(int);
  void BaseConstruction();
  void EnterCoefficients(int id);

  Vector localas, localbs, localmask1,localmask2,localmask3;

public:
  GARCHModelDialog(GARCHModel *gm, Window& w); // centered over w
  int run();  // local version
};

class GARCHInterface : public TSMInterface
{
 protected:
  GARCHModel *mp;
 public:
  GARCHInterface(GARCHModel *gmp) :
    TSMInterface(), mp(gmp) {}

  bool Specify(Gtk::Window *mainwin);
};



class EGARCHModelDialog : public CenteringDialog {
protected:
  EGARCHModel *egmp;

  Gtk::Button *ok_button, *cancel_button;
  Gtk::Table *table;
  Gtk::SpinButton *aspin, *bspin, *muspin;
  Gtk::Adjustment *aadjust, *badjust, *muadjust;

  void terminate(int);
  void BaseConstruction();
  void EnterCoefficients(int id);

  Vector localas, localbs, localgs, 
    localmask1, localmask2, localmask3, localmask4;

public:
  EGARCHModelDialog(EGARCHModel *egm, Window& w); // centered over w
  int run();   // local version checks validity of params
};


class EGARCHInterface : public TSMInterface
{
 protected:
  EGARCHModel *mp;
 public:
  EGARCHInterface(EGARCHModel *gmp) :
    TSMInterface(), mp(gmp) {}

  bool Specify(Gtk::Window *mainwin);
};



class SSVMDialog : public CenteringDialog {
 protected:
  SSVModel *svmp;
  Vector localmask;

  Gtk::Adjustment *meanadjust, *muxadjust, *phiadjust, *nuadjust;
  Gtk::SpinButton *meanspin, *muxspin, *phispin, *nuspin;
  Gtk::CheckButton *meancheck, *muxcheck, *phicheck, *nucheck;

  void BaseConstruction();
  void terminate(int);

 public:
  SSVMDialog(SSVModel *mp, Window& w);   
  int run();         // local version checks for validity of params
};

class SSVMInterface : public TSMInterface
{
 protected:
  SSVModel *svmmp;
 public:
  SSVMInterface(SSVModel *mp)
    : TSMInterface(), svmmp(mp) {}
  
  bool Specify(Gtk::Window *mainwin);             // called in response to model->specify menu option
};

#endif
