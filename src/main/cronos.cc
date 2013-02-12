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
 * This is the main module of the open-source time series analysis
 * package (tsap).
 *
 */

#include "includes.h"
#include <gtkmm/notebook.h>
#include <gtkmm/scrolledwindow.h>
#include <fstream.h>
#include "tsmrender.h"
#include "dialogs.h"
#include "tsinterface.h"
#include "mainwin.h"
#include <string>
#include <vector>
#include "otherpopups.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_multimin.h>


using namespace SigC;
using namespace Gtk;
using namespace std;


Window *mainwin;
Window *get_main_win() // used by other modules
{ return mainwin; }


int main(int argc, char** argv)
{
  Gtk::Main kit(argc,argv);

  init_generator(0);  // seed the random number generator

  mainwin = new CoreWin(); // a corewin keeps track of a time series and a time series model
  //mainwin = new Window();
  mainwin->show_all();

  kit.run(*mainwin);    // the argument ensures that closing the window terminates the program
  
  return 0;
}
