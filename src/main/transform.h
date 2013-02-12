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

#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

#include "timeseries.h"
#include "cubespline.h"
#include "otherplots.h"

class Transformation {
protected:
public:
  Transformation();
  virtual ~Transformation() {};
  virtual bool Apply(TimeSeries *, ostream&, bool getparms=true)=0;
  virtual bool Reverse(TimeSeries *, ostream&, bool discard=true)=0;
  virtual void DescribeIn(ostream& os);
  virtual void GetParList() {};
  virtual bool GetParameters() {return true;};
  virtual Transformation *CreateCopy()=0;

  double dpars[10];
  int ipars[10];
  string spars[10];
  string parnames[10], pardefaults[10], transname;
  int partypes[10];  // 0=int, 1=real, 2=string
  int npars;
};

class RemoveSplineTransform : public Transformation
{
protected:
  int numknots, nforfit;  
  cubespline *cs;
  TSPlot *pp;

public:
  RemoveSplineTransform();
  virtual ~RemoveSplineTransform(); 
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList()
  {
    npars = 1;
    partypes[0] = 0;  // non-negative integer
    parnames[0] = "Number of Knots";
    pardefaults[0] = "3";
    transname = "Cubic Spline Trend Removal";
  };
  void SetPlotWin(TSPlot *tpp) { pp=tpp; };
  Transformation *CreateCopy();
};

class RemoveSeasonal : public Transformation
{
protected:
  int period;
  TSPlot *pp;
  Vector scomp;

public:
  RemoveSeasonal();
  virtual ~RemoveSeasonal(); 
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList()
  {
    npars = 1;
    partypes[0] = 0;  // non-negative integer
    parnames[0] = "Period";
    pardefaults[0] = "12";
    transname = "Seasonal Component Removal";
  };
  void SetPlotWin(TSPlot *tpp) { pp=tpp; };
  Transformation *CreateCopy();
};

class DiffTransform : public Transformation
{
protected:
  int lag;
  double *beginning;
  bool isreversible;

public:
  DiffTransform();
  DiffTransform(DiffTransform& other);
  virtual ~DiffTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList()
  {
    npars = 1;
    partypes[0] = 0;  // non-negative integer
    parnames[0] = "Lag";
    pardefaults[0] = "1";
    transname = "Difference";
  };
  Transformation *CreateCopy();
};

class IntegrateTransform : public Transformation
{
protected:
  int lag;
  double *beginning;
  bool isreversible;

public:
  IntegrateTransform();
  IntegrateTransform(IntegrateTransform& other);
  virtual ~IntegrateTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList()
  {
    npars = 0;
    transname = "Integrate";
  };
  Transformation *CreateCopy();
};

class AbsTransform : public Transformation
{
public:
  AbsTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList() { npars = 0; };
  Transformation *CreateCopy();
};

class TimeReverseTransform : public Transformation
{
public:
  TimeReverseTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList() { npars = 0; };
  Transformation *CreateCopy();
};

class SMeanTransform : public Transformation
{
protected:
  double subtracted;
  
public:
  SMeanTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void GetParList()
  {
    npars = 1;
    partypes[0] = 1;  // double
    parnames[0] = "Constant";
    pardefaults[0] = "0";
    transname = "Subtract";
  };
};

class LogTransform : public Transformation
{
protected:
  bool isreversible;
  double lambda;

public:
  LogTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList()
  {
    npars = 1;
    partypes[0] = 1;  // Box-Cox param is double-valued
    pardefaults[0] = "0.0";  // default=0 for log
    parnames[0] = "Lambda";
    transname = "Box-Cox Parameter";
  };
  Transformation *CreateCopy();
};

class AggTransform : public Transformation
{
protected:
  int length;
public:
  AggTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList()
  {
    npars = 1;
    partypes[0] = 0;  // non-negative integer
    parnames[0] = "Interval Length";
    pardefaults[0] = "1";
    transname = "Aggregate";
  };
  Transformation *CreateCopy();
};

class PolyTransform : public Transformation
{
protected:
  int power;
public:
  PolyTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
};

class ClipTransform : public Transformation
{
protected:
  int n1, n2;
  TSPlot *pp;

public:
  ClipTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters(TimeSeries *);
  void GetParList()
  {
    npars = 2;
    partypes[0] = 0;  // non-negative integer
    parnames[0] = "First obs.";
    pardefaults[0] = "1";
    partypes[1] = 0;
    parnames[1] = "Last obs." ;
    pardefaults[1] = "100";
    transname = "Clip";
  };
  void DescribeIn(ostream&);
  void SetPlotWin(TSPlot *tpp) { pp=tpp; };
  Transformation *CreateCopy();
};

class CustomTransform : public Transformation
{
protected:
  int lag;
  double *beginning;
  bool isreversible;

public:
  CustomTransform();
  CustomTransform(CustomTransform& other);
  virtual ~CustomTransform();
  bool Apply(TimeSeries *, ostream&, bool getparms=true);
  bool Reverse(TimeSeries *, ostream&, bool discard=true);
  bool GetParameters();
  void DescribeIn(ostream&);
  void GetParList()
  {
    npars = 2;
    partypes[0] = 2;  // string
    parnames[0] = "Function of x (and t=0,1,...)";
    pardefaults[0] = "x";
    partypes[1] = 2;
    parnames[1] = "Inverse function (may be blank)";
    pardefaults[1] = "";
    transname = "Custom";
  };
  Transformation *CreateCopy();
};


#endif
