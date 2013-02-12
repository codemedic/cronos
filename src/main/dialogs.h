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

#ifndef DIALOGS_HPP
#define DIALOGS_HPP

#include "timeseries.h"

class Transformation;
class TransformDialog;

class CenteringDialog : public Gtk::Dialog {
public:
  Window *cd_parent_win;
  CenteringDialog(Glib::ustring title, Window& w);
};


class TransformDialog : public CenteringDialog {
protected:
  Transformation *trans;
  Gtk::Table *table;
  Gtk::Entry **es;
  Gtk::Label **ls;

public:
  TransformDialog(Transformation *t, Window& wi);

  void terminate(int);
  void enter_pressed();
};

class IterationCounter : public CenteringDialog {
protected:
  int localcount;
  Gtk::ProgressBar *pbar;
  Gtk::Label *label1, *label2, *stage_label;
  int curstage;

public:
  IterationCounter(Window& wp);

  int increment(int stage);
};

class CoefficientDialog : public CenteringDialog {
protected:
  int nentries;
  Gtk::Entry **eps;
  Gtk::CheckButton **checks;

public:
  CoefficientDialog(Window& wi, int n, string *, double *val, double *mask);
  ~CoefficientDialog();
  void on_response(int response_id);

  Vector retvals, maskvals;
};

class SimulateDialog : public CenteringDialog {
protected:
  Gtk::Table *table;
  Gtk::SpinButton *nspin;
  Gtk::Adjustment *nadjust;

public:
  SimulateDialog(Window& wi, int defaultn);
  int get_n();
};

class ForecastDialog : public CenteringDialog {
protected:
  Gtk::Table *table;
  Gtk::SpinButton *nspin, *simsspin;
  Gtk::Adjustment *nadjust, *simsadjust;

public:
  ForecastDialog(Window& wi, int defaultn);
  int get_n();
  int get_sims();
};

class GeneralTextDialog : public CenteringDialog {
protected:
  Gtk::Table *table;
  Gtk::Entry *entry;

public:
  GeneralTextDialog(Window& wi, string label, string thetext);
  string get_text();
};

#endif
