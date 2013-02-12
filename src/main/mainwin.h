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

#ifndef MAINWIN_H
#define MAINWIN_H

class CoreWin : public Window {
 protected:
  FileSelection *fs;
  TimeSeriesModel *model;
  TSMInterface *interface;
  int modeltype;
  TimeSeries *timeseries;
  vector<Transformation*> transformlist;
  TSPlot *plotwin;
  ACFPlot *acfplot;
  PACFPlot *pacfplot;
  TSMRenderer *modelrenderer;
  Transformation *trans_loc;
  IterationCounter *itcount;
  Glib::TimeVal ittime;
  Glib::ustring local_data;

  //void on_clipboard_get(Gtk::SelectionData& selection_data, guint);
 public:
  CoreWin();  // constructor (empty TS, etc.)

  void clear_transforms();
  void update_all(bool doacf=false);

  void file_new_cb();
  void file_open_cb();
  void file_save_cb();
  void file_exit_cb();
  void undo_transform_cb();
  void data_periodogram();
  void data_import_clipboard();
  void transform_specified(bool getparms);
  void transform_cb(int tt, bool getparms);
  void help_about_cb();
  void model_specify_cb(int mtype);
  void model_specified();
  void model_fit_iteration(int stage);
  void model_fit_cb(int methodtype);
  void model_residuals_cb();
  void model_acf_pacf_cb();
  void model_spectral_density();
  void model_simulate();
  void model_forecast_cb();
  void model_view_cts_version_cb();

  void on_clipboard_received(const Gtk::SelectionData& selection_data);


};

#endif
