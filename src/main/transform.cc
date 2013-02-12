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
 * This module contains code which applies
 * various transformations to timeseries.
 */

#include "includes.h"
#include "timeseries.h"
#include "transform.h"
#include "cubespline.h"
#include "muParser.h"

using namespace Gtk;
using namespace mu;


Window *get_main_win();  // defined in main mitsm2.cc, used by other modules


Transformation::Transformation()
{ }

void Transformation::DescribeIn(ostream &os)
{
  os << "(unknown)";
}

double boxcox(double x, double lambda)
{
  if (lambda>0)
    return (pow(x,lambda)-1.0)/lambda;
  else
    return log(x);
}

double iboxcox(double x, double lambda)
{
  if (lambda>0)
    return pow(x*lambda+1.0,1/lambda);
  else
    return exp(x);
}

// --- Differencing Transformation ---

DiffTransform::DiffTransform()
  : Transformation()
{
  beginning = NULL;
  isreversible = true;
  lag = 1;   // default lag
}

DiffTransform::DiffTransform(DiffTransform& other)
  : Transformation()
{
  lag = other.lag;
  isreversible = other.isreversible;
  if (lag>0)
    {
      beginning = new double[lag];
      for (int i=0 ; i<lag ; ++i)
	beginning[i] = other.beginning[i];
    }
}

DiffTransform::~DiffTransform()
{
  delete[] beginning;
}

Transformation *DiffTransform::CreateCopy()
{
  DiffTransform *retval = new DiffTransform(*this);  
  return (Transformation *) retval;
}

void DiffTransform::DescribeIn(ostream &os)
{
  os << "Differenced at Lag " << lag;
}

bool DiffTransform::GetParameters()
{
  lag = ipars[0];
  return true;
}

bool DiffTransform::Apply(TimeSeries *tp, ostream &err, bool getparms)
{
  int i,n=tp->GetN();
  Vector mv;

  if (getparms)
    if (!GetParameters())
      return false;

  if (lag>=n)
    {
      err << "There are only " << n << " data points.  You cannot difference ";
      err << "at lag " << lag << ends;
      return false;
    }

  delete[] beginning;
  beginning = new double[lag];
  double *data = tp->GetData();
  mv = tp->GetMissingVector();
  for (i=0 ; i<lag ; ++i)
    beginning[i] = data[i];
  for (i=n-1 ; i>=lag ; --i)
    {
      data[i] -= data[i-lag];
      mv[i] += mv[i-lag];
      if (mv[i]>1)
	mv[i] = 1;
      if (mv[i])
	isreversible = false;
    }
  tp->SetMissingVector(mv);
  tp->Clip(lag,n-1);

  return true;
}

bool DiffTransform::Reverse(TimeSeries *tp, ostream &err, bool discard)
{
  int i,n;
  double *data,*newdata;

  if (!isreversible)
    {
      err << "You can't undo differencing when there are missing values.";
      return false;
    }

  for (i=0 ; i<lag ; ++i)
    tp->Append(0);

  data = tp->GetData();
  n = tp->GetN();
  newdata = new double[n];
  for (i=0 ; i<lag ; ++i)
    newdata[i] = beginning[i];
  for (i=lag ; i<n ; ++i)
    newdata[i] = newdata[i-lag] + data[i-lag];
  for (i=0 ; i<n ; ++i)
    data[i] = newdata[i];
  delete[] newdata;

  if (discard) {
    delete[] beginning;
    beginning = NULL;
  }

  return true;
}


// --- Integration Transformation ---

IntegrateTransform::IntegrateTransform()
  : Transformation()
{
  isreversible = true;
}

IntegrateTransform::IntegrateTransform(IntegrateTransform& other)
  : Transformation()
{
  isreversible = other.isreversible;
}

IntegrateTransform::~IntegrateTransform()
{
}

Transformation *IntegrateTransform::CreateCopy()
{
  IntegrateTransform *retval = new IntegrateTransform(*this);  
  return (Transformation *) retval;
}

void IntegrateTransform::DescribeIn(ostream &os)
{
  os << "Integrated";
}

bool IntegrateTransform::GetParameters()
{
  return true;
}

bool IntegrateTransform::Apply(TimeSeries *tp, ostream &err, bool getparms)
{
  int i,j,n=tp->GetN();
  Vector mv;

  if (getparms)
    if (!GetParameters())
      return false;

  double *data = tp->GetData();
  mv = tp->GetMissingVector();

  for (i=1 ; i<n ; ++i)
    {
      j=i-1;
      while (j>=0 && mv[j])
	--j;
      if (j!=-1)
	data[i] += data[j];  // last non-missing
    }

  if (tp->HasAnyMissing())
    isreversible = false;

  return true;
}

bool IntegrateTransform::Reverse(TimeSeries *tp, ostream &err, bool discard)
{
  int i,n;
  double *data,*newdata;

  if (!isreversible)
    {
      err << "You can't undo integration when there were missing values" << endl;
      err << "in the original data.";    // this is the only possible cause for irreversibility
      return false;
    }


  data = tp->GetData();
  n = tp->GetN();
  for (i=n-1 ; i>=1 ; --i)
    data[i] -= data[i-1];

  return true;
}


// --- Absolute Value Transformation ---

AbsTransform::AbsTransform()
{
}

bool AbsTransform::GetParameters()
{
  return true;
}

Transformation *AbsTransform::CreateCopy()
{
  AbsTransform *retval = new AbsTransform(*this);
  return (Transformation *) retval;
}

bool AbsTransform::Apply(TimeSeries *tp, ostream& err, bool getparms)
{
  int i,j,n=tp->GetN(),newn;
  double *data = tp->GetData(), sum;

  for (i=0 ; i < n ; ++i)
    data[i] = fabs(data[i]);

  return true;
}

bool AbsTransform::Reverse(TimeSeries *tp, ostream &os, bool discard)
{
  os << "Absolute value transformations are not reversible." << ends;
  return false;
}

void AbsTransform::DescribeIn(ostream &os)
{
  os << "Absolute value taken";
}

// --- TimeReverse Transformation ---

TimeReverseTransform::TimeReverseTransform()
{
}

bool TimeReverseTransform::GetParameters()
{
  return true;
}

Transformation *TimeReverseTransform::CreateCopy()
{
  TimeReverseTransform *retval = new TimeReverseTransform(*this);
  return (Transformation *) retval;
}

bool TimeReverseTransform::Apply(TimeSeries *tp, ostream& err, bool getparms)
{
  // getparms is irrelevant since nparms=0

  int i,n=tp->GetN();
  double *data = tp->GetData();
  Vector tcopy(n);

  for (i=0 ; i < n ; ++i)
    tcopy[i] = data[i];
  for (i=0 ; i < n ; ++i)
    data[i] = tcopy[n-1-i];

  return true;
}

bool TimeReverseTransform::Reverse(TimeSeries *tp, ostream &os, bool discard)
{
  Apply(tp, os, false);
  return true;
}

void TimeReverseTransform::DescribeIn(ostream &os)
{
  os << "Time reversed";
}


// --- Subtract mean transformation ---

SMeanTransform::SMeanTransform()
{
}

bool SMeanTransform::GetParameters()
{
  return true;
}

bool SMeanTransform::Apply(TimeSeries *tp, ostream& err, bool getparms)
{
  int n=tp->GetN();
  if (getparms)
    if (!GetParameters())
      return false;
  double *data = tp->GetData(), mn = tp->SampleMean();
  for (int i=0 ; i<n ; ++i)
    data[i] -= mn;
  subtracted = mn;
  return true;
}

bool SMeanTransform::Reverse(TimeSeries *tp, ostream &err, bool discard)
{
  int n=tp->GetN();
  double *data = tp->GetData();
  for (int i=0 ; i<n ; ++i)
    data[i] += subtracted;
  return true;
}

// --- Log Transform ---

LogTransform::LogTransform()
{
  lambda = 0.0;  // default value
}

bool LogTransform::GetParameters()
{
  return true;
}

Transformation *LogTransform::CreateCopy()
{
  LogTransform *retval = new LogTransform(*this);
  return (Transformation *) retval;
}

void LogTransform::DescribeIn(ostream &os)
{
  if (fabs(lambda)<1e-8)
    os << "Log";
  else
    os << "Box-Cox(" << lambda << ")";
}

bool LogTransform::Apply(TimeSeries *tp, ostream& err, bool getparms)
{
  if (getparms)
    {
      if (!GetParameters())
	return false;
      lambda = dpars[0];  // first double-valued param. is lambda
    }

  if (tp->SampleMin()<=0.0)
    {
      Window *mw = get_main_win();
      MessageDialog dialog(*mw, "There are non-positive elements in the time series.  Do you wish to replace them by missing values?", false,
				MESSAGE_QUESTION,
				BUTTONS_YES_NO, true);
      int result = dialog.run();
      if (result==RESPONSE_NO)
	{
	  err << "Log transformation cancelled." << ends;
	  return false;
	}
    }

  isreversible = true;
  int n=tp->GetN();
  double *data = tp->GetData(), tx;
  for (int i=0 ; i<n ; ++i)
    {
      tx = data[i];
      if (tx<=0.0)
	{ tp->SetMissing(i,true);  isreversible = false; }
      else
	(*tp)[i] = boxcox(tx, lambda);
    }
  return true;
}

bool LogTransform::Reverse(TimeSeries *tp, ostream& err, bool discard)
{
  if (!isreversible)
    {
      err << "Transformation is not reversible since it introduced missing values." << ends;
      return false;
    }

  int n=tp->GetN();
  double *data = tp->GetData();
  for (int i=0 ; i<n ; ++i)
    data[i] = iboxcox(data[i], lambda);
  return true;
}

// --- Aggregation ---

AggTransform::AggTransform()
{
}

Transformation *AggTransform::CreateCopy()
{
  AggTransform *retval = new AggTransform(*this);
  return (Transformation *)retval;
}

bool AggTransform::GetParameters()
{
  length = ipars[0];
  return true;
}

bool AggTransform::Apply(TimeSeries *tp, ostream& err, bool getparms)
{
  int i,j,n=tp->GetN(),newn;
  double *data = tp->GetData(), sum;
  bool tbool;

  if (getparms)
    if (!GetParameters())
      return false;

  newn = n/length;

  for (i=0 ; i*length < n ; ++i)
    {
      sum = 0.0;
      for (j=0,tbool=false ; j<min(n-i*length, length) ; ++j)
	{
	  sum += data[i*length + j];
	  tbool |= tp->IsMissing(i*length+j);
	}
      data[i] = sum;
      tp->SetMissing(i,tbool);
    }

  tp->Clip(0, newn-1);
  return true;
}

bool AggTransform::Reverse(TimeSeries *tp, ostream &os, bool discard)
{
  os << "Aggregation transformations are not reversible." << ends;
  return false;
}

void AggTransform::DescribeIn(ostream &os)
{
  os << "Aggregated in lots of " << length;
}

//--------------------------------------------------------

PolyTransform::PolyTransform()
{
}

bool PolyTransform::GetParameters()
{
  power = 2;
  return true;
}

bool PolyTransform::Apply(TimeSeries *tp, ostream& err, bool getparms)
{
  int i,j,n=tp->GetN(),newn;
  double *data = tp->GetData(), tx;

  if (getparms)
    if (!GetParameters())
      return false;

  for (i=0 ; i<n ; ++i)
    {
      tx = data[i];
      for (j=1 ; j<power ; ++j)
        data[i] *= tx;
    }

  return true;
}

bool PolyTransform::Reverse(TimeSeries *tp, ostream& err, bool discard)
{
  err << "Polynomial transformations are not reversible." << ends;
  return false;
}

void PolyTransform::DescribeIn(ostream &os)
{
  if (power==2)
    os << "Squared";
  else if (power==3)
    os << "Cubed";
  else
    os << "Raised to power " << power;
}

//--------------------------------------------------------

RemoveSeasonal::RemoveSeasonal()
{
}

RemoveSeasonal::~RemoveSeasonal()
{
}

Transformation *RemoveSeasonal::CreateCopy()
{
  RemoveSeasonal *retval = new RemoveSeasonal(*this);
  return (Transformation *) retval;
}

void RemoveSeasonal::DescribeIn(ostream &os)
{
  os << "Period " << period << " seasonal component removed";
}

bool RemoveSeasonal::GetParameters()
{
  period = ipars[0];
  return true;
}

bool RemoveSeasonal::Apply(TimeSeries *tp, ostream &err, bool getparms)
{
  double tx;
  int i,j,p,n=tp->GetN();
  if (getparms)
    {
      if (!GetParameters())
	return false;

      scomp.resize(period);

      for (p=0 ; p<period ; ++p)
	{
	  scomp[p]=0;
	  for (i=p,j=0 ; i<n ; i+=period)
	    if (!tp->IsMissing(i)) {
	      scomp[p] += (*tp)[i];
	      ++j;
	    }
	  scomp[p] /= j;
	}

      // display in plotwin
      TimeSeries theseason;
      theseason = *tp;
      for (i=0 ; i<n ; ++i)
	theseason[i] = scomp[i % period];
      
      pp->SetMultiMode(true,true,false);
      pp->Update(*tp, theseason);
      
      // get OK/Cancel
      Window *mw = get_main_win();
      Gtk::MessageDialog dialog(*mw, "The estimated seasonal component is dispayed on the time series plot in red.  Do you wish to remove it?", false,
				Gtk::MESSAGE_QUESTION,
				Gtk::BUTTONS_YES_NO, false);
      int result = dialog.run();

      if (result==Gtk::RESPONSE_NO)
	{
	  err << "Deseasonalization cancelled." << ends;
	  return false;
	}
    }

  // go ahead and remove it
  for (i=0 ; i<n ; ++i)
    {
      tx=((double)i)/n;
      (*tp)[i] -= scomp[i % period];
    }

  return true;
}

bool RemoveSeasonal::Reverse(TimeSeries *tp, ostream &err, bool discard)
{
  double tx;
  int i,n=tp->GetN();
  for (i=0 ; i<n ; ++i)
    (*tp)[i] += scomp[i%period];
  return true;
}

//--------------------------------------------------------

RemoveSplineTransform::RemoveSplineTransform()
  : Transformation()
{
  cs = NULL;
}

RemoveSplineTransform::~RemoveSplineTransform()
{
  delete cs;
}

Transformation *RemoveSplineTransform::CreateCopy()
{
  RemoveSplineTransform *retval = new RemoveSplineTransform(*this);
  if (cs!=NULL)
    {
      cubespline *localcs = new cubespline(*cs);
      cs = localcs;   // we keep the copy of cubespline, the copy of transform gets the old one
    }
  return (Transformation *) retval;
}

void RemoveSplineTransform::DescribeIn(ostream &os)
{
  os << numknots << "-knot cubic spline removed";
}

bool RemoveSplineTransform::GetParameters()
{
  numknots = ipars[0];
  return true;
}

bool RemoveSplineTransform::Apply(TimeSeries *tp, ostream &err, bool getparms)
{
  double tx;
  int i,j,n=tp->GetN(),nummissing;

  if (getparms)
    {
      if (!GetParameters())
	return false;

      // now fit it
      for (i=0,nummissing=0 ; i<n ; ++i)
	nummissing += tp->IsMissing(i);
      
      Matrix data(n-nummissing, 2);
      for (i=0,j=0 ; i<n ; ++i)
	{
	  if (!tp->IsMissing(i))
	    {
	      tx = ((double)i)/n; 
	      data[j][0] = tx;
	      data[j++][1] = (*tp)[i];
	    }
	}
      
      cs = new cubespline(numknots);
      Vector knots(numknots);
      for (i=0 ; i<numknots ; ++i)
	knots[i] = ((double)(i+1))/(numknots+1);
      cs->setknots(knots);
      cs->lsfit(data);
      
      // display in plotwin
      TimeSeries thetrend;
      thetrend = *tp;
      for (i=0 ; i<n ; ++i)
	{
	  tx = ((double)i)/n;
	  thetrend[i] = cs->getfitfor(tx);
	}
      nforfit = n;
      
      pp->SetMultiMode(true,true,false);
      pp->Update(*tp, thetrend);
      
      // get OK/Cancel
      Window *mw = get_main_win();
      Gtk::MessageDialog dialog(*mw, "The estimated trend is displayed on the time series plot in grey.  Do you wish to remove it?", false,
				Gtk::MESSAGE_QUESTION,
				Gtk::BUTTONS_YES_NO, true);
      int result = dialog.run();
      if (result==Gtk::RESPONSE_NO)
	{
	  err << "Trend removal cancelled." << ends;
	  return false;
	} 
    }

  // now subtract the fit
  for (i=0 ; i<n ; ++i)
    {
      tx=((double)i)/nforfit;
      (*tp)[i] -= cs->getfitfor(tx);
    }

  return true;
}

bool RemoveSplineTransform::Reverse(TimeSeries *tp, ostream &err, bool discard)
{
  double tx;
  int i,n=tp->GetN();
  for (i=0 ; i<n ; ++i)
    {
      tx=((double)i)/nforfit;
      (*tp)[i] += cs->getfitfor(tx);
    }
  return true;
}


//--------------------------------------------------------

ClipTransform::ClipTransform()
{
}

bool ClipTransform::GetParameters(TimeSeries *tp)
{
  n1 = ipars[0];
  n2 = ipars[1];
  return true;
}

bool ClipTransform::Apply(TimeSeries *tp, ostream& err, bool getparms)
{
  int i,j,n=tp->GetN(),newn;
  double *data = tp->GetData(), tx;

  if (getparms)
    if (!GetParameters(tp))
      return false;

  if ((n1<1) || (n2>n) || (n1>=n2))
    {
      err << "Clipping bounds must be in the range from 1 to " << n << ", with n(1)<n(2)." << ends;
      return false;
    }

  tp->Clip(n1-1,n2-1);

  return true;
}

bool ClipTransform::Reverse(TimeSeries *tp, ostream& err, bool discard)
{
  err << "Clipping transformations are not reversible." << ends;
  return false;
}

void ClipTransform::DescribeIn(ostream &os)
{
  os << "Clipped to Observations " << n1 << " ... " << n2;
}

Transformation *ClipTransform::CreateCopy()
{
  ClipTransform *retval = new ClipTransform(*this);
  return (Transformation *) retval;
}


// --- Custom Transformation ---

CustomTransform::CustomTransform()
  : Transformation()
{
}

CustomTransform::CustomTransform(CustomTransform& other)
  : Transformation()
{
}

CustomTransform::~CustomTransform()
{
}

Transformation *CustomTransform::CreateCopy()
{
  CustomTransform *retval = new CustomTransform(*this);  
  return (Transformation *) retval;
}

void CustomTransform::DescribeIn(ostream &os)
{
  os << "Custom transformation: " << spars[0];
}

bool CustomTransform::GetParameters()
{
  return true;
}

bool CustomTransform::Apply(TimeSeries *tp, ostream &err, bool getparms)
{
  int i,n=tp->GetN();
  Vector mv;
  double *data = tp->GetData(), tx, tx0;
  mu::Parser parser, invparser;
  value_type vVarVal[2];  // array containing variables defined for parser (x and t)

  parser.DefineVar("x", &vVarVal[0]);
  invparser.DefineVar("x", &vVarVal[0]);
  parser.DefineVar("t", &vVarVal[1]);
  invparser.DefineVar("t", &vVarVal[1]);
 
  parser.SetExpr(spars[0]);

 
  // check for validity of expression
  vVarVal[0]=data[0];
  vVarVal[1]=0;
  try { parser.Eval(); }
  catch (mu::Parser::exception_type &e)
    {
      err << "Transformation expression non-interpretable.";
      return false;
    }

  // check for validity of inverse (if there)
  if (spars[1].length()>0)
    {
      // check validity of inverse (holding t=vVarVal[1]=1)
      vVarVal[1]=1;
      invparser.SetExpr(spars[1]);
      vVarVal[0] = parser.Eval();
      try  { invparser.Eval(); }
      catch (mu::Parser::exception_type &e)
	{
	  err << "Inverse transformation expression non-interpretable.";
	  return false;
	}
      // check invertibility
      for (i=0 ; i<n ; ++i)
	if (!tp->IsMissing(i))
	  {
	    vVarVal[0] = data[i];
	    vVarVal[1] = i;
	    tx0 = parser.Eval();
	    vVarVal[0] = tx0;
	    tx = invparser.Eval();
	    if ((fabs(tx-data[i])>1e-10) || (isnan(tx) && !isnan(data[i])))
	      {
		err << "Inverse transformation does not yield original values." << endl;
		err << "time " << i << ": x=" << data[i] << ", f(x)=" << tx0 << ", inv(f(x))=" << tx << endl; 
		return false;
	      }
	  }
    }

  for (i=0 ;i<n ; ++i)
    {
      vVarVal[0] = data[i];
      vVarVal[1] = i;

      data[i] = parser.Eval();
      if (isnan(data[i]))
	tp->SetMissing(i, true);
    }

  return true;
}

bool CustomTransform::Reverse(TimeSeries *tp, ostream &err, bool discard)
{
  int i,n = tp->GetN();
  double *data = tp->GetData();
  value_type tx, tt;

  if (spars[1].length()==0)
    {
      err << "Custom transforms are only reversible if you specify the inverse transformation." << endl;
      return false;
    }

  // otherwise reverse it all
  mu::Parser invparser;
  invparser.DefineVar("x", &tx);
  invparser.DefineVar("t", &tt);
  invparser.SetExpr(spars[1]);

  for (i=0 ; i<n ; ++i)
    if (!tp->IsMissing(i))
      {
	tx = data[i];
	tt = i;
	data[i] = invparser.Eval();
	if (isnan(data[i]))
	  tp->SetMissing(i, true);
      }

  return true;
}
