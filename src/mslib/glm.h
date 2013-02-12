/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004 Anthony Brockwell
 *
 * This library is free software; you can redistribute it and/or
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

#ifndef GLM_H
#define GLM_H

#include "matrix.h"

namespace mslib {

class GenLinearModel {
 protected:
  virtual double linkfn(double)=0;
  virtual double invlinkfn(double)=0;
  virtual double linkderiv(double)=0;
  virtual double variancefn(double etai)=0;
  virtual double dispersion(int i)=0;

  Vector betahat;
  Matrix betacov;
  Vector wts;

 public:
  GenLinearModel();

  Vector Estimate(Vector& y, Matrix& X, Vector *initbetahat=NULL,
		   int maxits=40);
  Vector Fitted(Matrix& X);
  Matrix *GetBetahatCovInv() {return &betacov;}
  void SetBetaHat(Vector& bh) {betahat=bh;}
};

class BinomialGLM : public GenLinearModel {
 protected:
  int binm;

 public:
  BinomialGLM(int m);

  double linkfn(double);
  double invlinkfn(double);
  double linkderiv(double);
  double variancefn(double etai);
  double dispersion(int i);
};

class PoissonGLM : public GenLinearModel {
 public:
  PoissonGLM();

  double linkfn(double);
  double invlinkfn(double);
  double linkderiv(double);
  double variancefn(double etai);
  double dispersion(int i);
  double scaleddeviance(Vector& y, Matrix& X);
};

}

#endif
