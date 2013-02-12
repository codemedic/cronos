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
 */

#ifndef CUBESPLINE_H
#define CUBESPLINE_H

using namespace mslib;

class cubespline {
 protected:
  double clip(double);
  double cube(double);
  static double cost(Vector& knotgaps, void *vp); // used for free knot fit
  static double cost2(Vector& gapsandparms, void *vp); // used for Pn version
  static double cost3(double *xpos, void *vp);     // used for finding argmax
  static double cost4(Vector& coeffs, void *vp);  // used for Pn MLE fit

 public:
  int numknots;
  Vector beta, knotpos;
  double sig2;
  Vector fitted;
  Matrix X, betacov;

  double minx, maxx;
  Matrix *tempdata;

  cubespline(int numknots);
  ~cubespline();

  void lsfit(Matrix& data);
  void fitPn(Matrix& data, int fileit=0);
  double loglikelihood(Matrix& data, int mode=0, double dispersion=1.0); // mode=1 for Poisson
  void display(FILE *gnupipe=NULL, double minx=0.0, double maxx=1.0,
	       double miny=0.0, double maxy=1.0,
	       Matrix *data=NULL);
  void setknots(Vector& knots);
  double getfitfor(double x);          // get fitted y(x)

  Vector getbeta(Matrix& cov);
  Vector knotsforgaps(double min, double max, Vector dirch);
  Vector gapsforknots(double min, double max, Vector knots);
  Matrix getdesignmatrix(Matrix& data);
};

#endif
