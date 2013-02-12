/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004,2005,2006 Anthony Brockwell
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
 *
 * This module contains the source code defining a
 * TimeSeries class.  Operations for reading from a file,
 * writing to a file, etc. are included.
 * The module  also contains the
 * base TimeSeriesModel and derived ARMAModel class.
 */

#include <iomanip.h>
#include <strstream>
#include <complex>
#include <string>
#include <functional>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "matrix.h"
#include "abrand.h"
#include "poly.h"
#include "timeseries.h"

using namespace std;

const int alloc_unit = 512, KFburnin=12;
double accuracy = 0.00005;  // used in FindMinimum
const double zero=1e-10;


//--------------------------------------------------+
// Time stuff                                       :
//--------------------------------------------------+
void JulianTime::GetYearMonthDay(int *y, int *m, int *d)
{
  int l,n,i,j;
  l = julianday + 68569;
  n = ( 4 * l ) / 146097;
  l = l - ( 146097 * n + 3 ) / 4;
  i = ( 4000 * ( l + 1 ) ) / 1461001;
  l = l - ( 1461 * i ) / 4 + 31;
  j = ( 80 * l ) / 2447;
  *d = l - ( 2447 * j ) / 80;
  l = j / 11;

  *m = j + 2 - ( 12 * l );
  *y = 100 * ( n - 49 ) + i + l;
}

void JulianTime::SetYearMonthDay(int y, int m, int d)
{
  julianday = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +
          ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -
          ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 +
          d - 32075;
}

void JulianTime::IncreaseBy(double amount)
{
  int f = (int)(floor(amount));
  julianday += f;
  dayfraction += amount-f;
}


//---------------------------+
// Miscellaneous             :
//---------------------------+

int min(int a, int b)
{
  return a < b ? a : b;
}

int max(int a, int b)
{
  return a < b ? b : a;
}

int abs(int x)
{
  return x < 0 ? -x : x;
}

//--------------------------------------------------+
// TimeSeries stuff                                 :
//--------------------------------------------------+

TimeSeries::TimeSeries()
{
  anymissing = false;
  n = 0;
  nallocated = 0;
  description = "";   title = "";
  current_component = 0;
  n_components = 1;
  data.resize(n_components, 1);
  missing.resize(n_components, 1);
}

TimeSeries::~TimeSeries()
{
  ClearOut();
}

bool TimeSeries::HasAnyMissing()
{
  return anymissing;
}

Vector TimeSeries::GetDataVector()
{
  int i;
  Vector temp;
  temp.resize(n);
  for (i=0 ; i<n ; ++i)
    temp[i] = data[current_component][i];
  return temp;
}

Matrix TimeSeries::GetDataMatrix()
{
  int i,j;
  Matrix tm(n,n_components);
  for (i=0 ; i<n ; ++i)
    for (j=0 ; j<n_components ; ++j)
      tm[i][j] = data[j][i];
  return tm;
}

Vector TimeSeries::GetMissingVector()
{
  int i;
  Vector temp;
  temp.resize(n);
  for (i=0 ; i<n ; ++i)
    temp[i] = missing[current_component][i];
  return temp;
}

Matrix TimeSeries::GetMissingMatrix()
{
  int i,j;
  Matrix tm(n,n_components);
  for (i=0 ; i<n ; ++i)
    for (j=0 ; j<n_components ; ++j)
      tm[i][j] = missing[j][i];
  return tm;
}

void TimeSeries::SetMissingVector(Vector& mv)
{
  int i;
  anymissing = false;
  for (i=0 ; i<n ; ++i)
    {
      missing[current_component][i] = (mv[i]!=0);
      anymissing |= (mv[i]!=0);
    }
}

void TimeSeries::SetMissing(int t, bool msng)
{
  missing[current_component][t] = msng;
  anymissing |= msng;
}

void TimeSeries::ClearOut(int dim)     // wipes it out and sets dimension
{
  n=0;
  anymissing = false;
  nallocated = 0;
  description = "";   title = "";
  if (dim!=0)
    n_components = dim;
  data.resize(n_components, 1);
  missing.resize(n_components, 1);
}

double TimeSeries::SampleMean()
{
  int i,n=GetN(),m;
  double total=0.0;
  for (i=0,m=0 ; i<n ; ++i)
    if (!missing[current_component][i])
      {
	total += data[current_component][i];
	++m;
      }
  if (m>0)
    total /= m;
  return total;
}

double TimeSeries::SampleMin()
{
  int i,n=GetN();
  double min=1e100;
  for (i=0 ; i<n ; ++i)
    if (!missing[current_component][i])
      if (data[current_component][i]<min)
	min = data[current_component][i];
  return min;
}

void ACFtoPACF(Vector& acvf, Vector& pacfstorage)
{
  int i,j,maxlag=acvf.nrows()-1;
  if (maxlag<=0)
    {
      pacfstorage[0] = 0;
      return;
    }

  // compute the sample PACF using the Durbin-Levinson algorithm (p.169 B&D)
  Vector phis(maxlag), phis2(maxlag);
  double vi, phinn;
  phis[0] = acvf[1]/acvf[0];
  pacfstorage[0] = 1.0;
  pacfstorage[1] = phis[0];
  vi = acvf[0];
  vi = vi*(1-phis[0]*phis[0]);
  for (i=2 ; i<=maxlag ; ++i) // i=iteration number
    {
      for (j=0 ; j<i-1 ; ++j)
        phis2[j] = phis[i-j-2];
      phinn = acvf[i];
      for (j=1 ; j<i ; ++j)
        phinn -= phis[j-1]*acvf[i-j];
      phinn /= vi;
      for (j=0 ; j<i-1 ; ++j)
        phis[j] -= phinn*phis2[j];
      vi = vi*(1-phinn*phinn);
      pacfstorage[i] = phis[i-1] = phinn;
    }
}

Matrix *Innovations(double *retvhat, int maxm, double *acfs)
{
  int j,k,m;
  Matrix *thetas = new Matrix(maxm,maxm);

  double *vhats = new double[maxm+1];
  vhats[0] = acfs[0];
  for (m=1 ; m<=maxm ; ++m)
    {
      // compute thetahat(i,1,...,i)
      for (k=0 ; k<m ; ++k)
        {
          (*thetas)[m-1][m-k-1] = acfs[m-k];
          for (j=0 ; j<k ; ++j)
            (*thetas)[m-1][m-k-1] -= (*thetas)[m-1][m-j-1]*(*thetas)[k-1][k-j-1]
              *vhats[j];
          (*thetas)[m-1][m-k-1] /= vhats[k];
        }
      vhats[m] = acfs[0];
      for (j=0 ; j<m ; ++j)
        vhats[m] -= (*thetas)[m-1][m-j-1]*(*thetas)[m-1][m-j-1]*vhats[j];
    }

  *retvhat = vhats[maxm];
  delete[] vhats;
  return thetas;
}

void TimeSeries::ComputeSampleACF(Vector *acvf, Vector *pacfstorage,
  bool normalizeacf)
{
  int i,j,n=GetN(),maxlag=acvf->nrows()-1,nonmissing=0;
  if (n==0)
    return;
  double mean=SampleMean(), total;
  double *tdata = &data[current_component][0];
  double *tmissing = &missing[current_component][0];

  // First compute the sample autocovariance function.
  if (!HasAnyMissing())
    for (i=0 ; i<=maxlag ; ++i)
      {
	total = 0;
	for (j=i ; j<n ; ++j)
	  total += (tdata[j]-mean)*(tdata[j-i]-mean);
	(*acvf)[i] = total/n;
      }
  else
    {
      for (i=0,nonmissing=0 ; i<n ; ++i)
	nonmissing += floor(1.5-tmissing[i]);
      for (i=0 ; i<=maxlag ; ++i)
	{
	  total = 0;
	  for (j=i ; j<n ; ++j)
	    if ((tmissing[j]==0) && (tmissing[j-i]==0))
	      total += (tdata[j]-mean)*(tdata[j-i]-mean);
	  (*acvf)[i] = total/nonmissing;
	}
    }

  // then the sample PACF using the Durbin-Levinson algorithm (p.169 B&D)
  if (pacfstorage!=NULL)
    ACFtoPACF(*acvf, *pacfstorage);

  // now normalize the sample ACFs
  if (normalizeacf)
    for (i=maxlag ; i>=0 ; --i)
      (*acvf)[i] /= (*acvf)[0];
}

void TimeSeries::Append(double x, bool miss)
{
  int i;
  double *tx;
  bool *tb;

  if (n_components>1)
    {
      cout << "Error: appending scalar to multivariate TS." << endl;
      return;
    }

  if (nallocated<n+1)
    {
      nallocated += alloc_unit;
      data.resize(n_components, nallocated);
      missing.resize(n_components, nallocated);
    }

  data[0][n]=x;
  missing[0][n++]=miss;
  anymissing |= miss;
}


void TimeSeries::Append(Vector& x, Vector& miss)
{
  int i,k;
  double *tx;
  bool *tb;

  if (n_components!=x.nrows())
    {
      cout << "Error: Append() dimension mismatch." << endl;
      return;
    }

  if (nallocated<n+1)
    {
      nallocated += alloc_unit;
      data.resize(n_components, nallocated);
      missing.resize(n_components, nallocated);
    }

  for (k=0 ; k<x.nrows() ; ++k)
    {
      data[k][n]=x[k];
      missing[k][n]=miss[k];
      anymissing != (miss[k]==1);
    }
  ++n;
}


void TimeSeries::Clip(int n0, int n1)
{
  int i;
  data = data.submatrix(0,n0,n_components,n1+1);
  missing = missing.submatrix(0,n0,n_components,n1+1);
  n = (n1-n0+1);
  nallocated = n;
}


TimeSeries& TimeSeries::operator=(const TimeSeries& other)
{
  int i;

  ClearOut();
  n = other.n;
  nallocated = other.nallocated;
  if (n>0) {
    data = other.data;
    missing = other.missing;
    nallocated = n;
  }
  anymissing = other.anymissing;
  current_component = other.current_component;

  return (*this);
}

bool TimeSeries::operator==(const TimeSeries& other)
{
  int i,j;
  bool retval=true;
  double tx;
  if (n != other.n)
    retval = false;
  for (j=0 ; j<n_components ; ++j)
    for (i=0 ; i<n ; ++i)
      {
	tx = other.missing[j][i];
	if (missing[j][i]!=tx)
	  retval=false;
	if (!missing[j][i])
	  if (data[j][i]!=other.data[j][i])
	    retval=false;
      }
  return retval;
}

void TimeSeries::ParseDescriptor(char *cp)
{
  int i, b,e;
  // first determine the command
  istrstream is(cp);
  char command[40], temps[1024], cr=13;
  string white;

  white = "\t\n ";
  white = white + cr;

  command[0]=0;
  is >> command;

  if (strcmp(command, "Title:")==0)
    {
      is.getline(temps, 1024, '\n');
      title = temps;
      b = title.find_first_not_of(white);
      e = title.find_last_not_of(white);
      title = string(title, b, e-b+1);
    }

  if (strcmp(command, "Description:")==0)
    {
      is.getline(temps, 1024, '\n');
      description = temps;
      b = description.find_first_not_of(white);
      e = description.find_last_not_of(white);
      description = string(description, b, e-b+1);
      // now replace %s with newlines
      i = description.find('%');
      while (i!=-1)
	{
	  description[i] = '\n';
	  i = description.find('%');
	}
    }
}

ostream& operator<<(ostream& os, TimeSeries& ts)
{
  int i,j;
  for (i=0 ; i<ts.n ; ++i)
    for (j=0 ; j<ts.n_components ; ++j)
      {
	if (!ts.IsMissing(j,i))
	  os << ts.data[j][i] << endl;
	else
	  os << "NA" << endl;
      }
  return(os);
}

istream& operator>>(istream& is, TimeSeries& ts)
{
  int i=0,j,k,nfields;
  const int maxlinelength = 65536;
  char temps[maxlinelength], temps2[maxlinelength],c;
  double tx;
  bool isok;
  int streampos;
  string s1;

  // First check to see if file contains lines which are
  // numeric, with endl at end, or NA (missing).  
  // If not, return with timeseries length=0.
  streampos = is.tellg();
  isok = true;
  for (i=0 ; i==0 ; )
    {
      is.getline(temps, maxlinelength);
      istrstream inst(temps);
      double tx;
      inst >> ws;             // skip over white space
      c = inst.peek();        // and see what's next

      // check it
      if ((c=='N') || (c>='0' && c<='9') 
	  || (c=='n') || (c=='+') || (c=='-')
	  || (c=='.') )
	{
	  i=1;  // found a line of data, so break and assume the remaining lines are OK
	  // but first count the number of fields
	  nfields = 0;
	  bool done=false;
	  while (!done)
	    {
	      inst >> tx;
	      if (inst)
		++nfields;
	      else
		{
		  if (inst.eof())
		    done=true;
		  else if (inst.fail()) // it's probably a NaN! still counts as a field
		    {
		      inst.clear();
		      inst >> s1;
		      ++nfields;
		    }
		}
	      if (inst.eof())
		done = true;
	    }
	}
      else if (c=='%')
	{ } // do nothing if it's just a comment line
      else
	{ // if it's not numeric and it's not a comment, then break and return
	  cout << "Error, character " << ((int)c) << " found." << endl;
	  isok = false;
	  i = 1;
	}
    }

  is.seekg(streampos);  // put it back to where it was
  if (!isok)
    return is;          // read failed for some reason

  // read data
  ts.ClearOut(nfields);
  Vector tv(nfields),tm(nfields);

  temps[0]=1;
  while ((!is.eof()) && (temps[0]!=0))
    {
      is.getline(temps, maxlinelength);
      istrstream inst(temps);

      // remove white space at beginning of line
      inst >> ws;
      c=inst.peek();

      // now read contents of line
      if (c!=-1) // -1 appears to happen at eof
        if (c=='%')
	  {
	    inst >> s1;
	    char *cp = (char *) (void *) s1.c_str();
            ts.ParseDescriptor(cp);
	  }
        else // read in nfields fields
          {
	    for (k=0 ; k<nfields ; ++k)
	      {
		inst >> tv[k];
		if (inst)
		  tm[k] = 0;          // not missing
		else
		  {
		    inst.clear();
		    inst >> s1;       // skip over  it
		    tm[k] = 1;
		    tv[k] = 0;
		  }
	      }
	    ts.Append(tv, tm);
	  }     
    }
  return is;
}

Vector TimeSeries::ComputePeriodogram(Vector *omegas)
{
  complex sum,consti(0,1),tc;
  double omega_j;
  int i,j,t;
  Vector periodogram(n);
  double *tdata = &data[current_component][0];

  if (n==0)
    return periodogram;  // can't compute periodogram of nothing.


  // new efficient way (same but enormously faster)
  for (i=0 ; i<n ; ++i)
    periodogram[i] = tdata[i];

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  work = gsl_fft_real_workspace_alloc (n);
  real = gsl_fft_real_wavetable_alloc (n);
  gsl_fft_real_transform (&periodogram[0], 1, n,
			  real, work);

  // now we have it stored in the efficient halfcomplex format
  // check gsl documentation (info gsl in linux) for details
  // so we unpack it into standard format
  double complex_coeffs[2*n+2];
  gsl_fft_halfcomplex_unpack (&periodogram[0], complex_coeffs, 1, n);

  // and free up memory
  gsl_fft_real_wavetable_free (real);
  gsl_fft_real_workspace_free (work);

  for (i=0 ; i<n ; ++i)
    periodogram[i] = ((complex_coeffs[2*i]*complex_coeffs[2*i]) 
		      + (complex_coeffs[2*i+1]*complex_coeffs[2*i+1])) / n;
    
  // fill in frequencies if necessary
  if (omegas!=NULL)
    {
      omegas->resize(n);
      for (i=0 ; i<n ; ++i)
	(*omegas)[i] = 2*M_PI*i/n;
    }

  return periodogram;
}

void TimeSeries::RenderInto(ostringstream &os)
{
  int i;

  os << "\\b";  // make it bold
  if (title.length()==0)
    os << "Time Series";
  else
    os << title;
  os << ":" << endl;
  os << "\\n";  // make it normal

  if (description.length()>0)
    os << description << endl << endl;

  os << "n = " << GetN() << endl;
  os << endl;
}

Matrix *TimeSeries::GetInnovations(double *retvhat, int maxm)
// this procedure returns a Matrix containing thetahat_{m,k}
// \hat{v} is returned in vhat
{
  int nacfs=maxm+1;

  Vector acfs(nacfs+1);
  ComputeSampleACF(&acfs, NULL, 0);
  Matrix *tm = Innovations(retvhat, maxm, &acfs[0]);

  return tm;
}

//-----------------------------------------------+
// TimeSeriesModel stuff                         :
//-----------------------------------------------+

TimeSeriesModel::TimeSeriesModel()
{
  loglikelihood = 0;
  AICC = 0;
  loglikelihood_is_valid = false;
  estimation_mask.resize(0);  // default is an empty mask
}

TimeSeriesModel::~TimeSeriesModel()
{
}

// the following two functions apply a mask to the
// parameter vector: the mask contains 1s if the respective
// component is to be included, 0s otherwise

Vector TimeSeriesModel::param_mask(Vector& parms, Vector& mask)
{
  int i, j, sz = (int) (mask.sum()+0.5); // round it
  Vector sub(sz);
  for (i=0,j=0 ; i<parms.nrows() ; ++i)
    if (mask[i])
      sub[j++] = parms[i];
  return sub;
}

Vector TimeSeriesModel::param_expand(const gsl_vector *x, Vector& curparms, Vector& mask)
{
  Vector retval = curparms;
  int i,j;
  for (i=0,j=0 ; i<curparms.size() ; ++i)
    if (mask[i])
      retval[i] = gsl_vector_get(x, j++);
  return retval;
}

// --- ARMAModel stuff ---

ARMAModel& ARMAModel::operator=(const ARMAModel& other)
{
  int i;
  p=other.p;   q=other.q;
  mean=other.mean;
  sigma=other.sigma;
  fracdiff=other.fracdiff;
  estimation_mask = other.estimation_mask;
  for (i=0 ; i<maxorder ; ++i)
    {
      phi[i] = other.phi[i];
      theta[i] = other.theta[i];
    }
  using_whittle = other.using_whittle;
  whittle_cached = other.whittle_cached;
  whittle_ts = other.whittle_ts;
  whittle_pdgm = other.whittle_pdgm;
  return *this;
}

ARMAModel::ARMAModel()
  : TimeSeriesModel()
{
  p=q=0;
  mean=0;
  sigma=1;		// default is white noise with mean 0, variance 1
  fracdiff = 0;
  int i;
  for (i=0 ; i<maxorder ; ++i)
    phi[i] = theta[i] = 0;
  estimation_mask.resize(3);
  estimation_mask[0]=1;
  estimation_mask[1]=0;  // default: don't estimate fracdiff parm
  estimation_mask[2]=1;  
  using_whittle = false;
  whittle_cached = false;
}

ARMAModel::~ARMAModel()
{
}

void ARMAModel::SetCoefficients(double *phis, double *thetas)
{
  int i;
  for (i=0 ; i<p ; ++i)
    phi[i] = phis[i];
  for (i=0 ; i<q ; ++i)
    theta[i] = thetas[i];
}

void ARMAModel::GetCoefficients(Vector& phis, Vector& thetas)
{
  int i;
  phis.resize(p);
  for (i=0 ; i<p ; ++i)
    phis[i] = phi[i];
  thetas.resize(q);
  for (i=0 ; i<q ; ++i)
    thetas[i] = theta[i];
}

bool ARMAModel::is_short_mem()
{
  return !(fabs(fracdiff) > zero);
}

void ARMAModel::SimulateInto(TimeSeries *ts, int nn, bool ovw)
// This function is perfectly fine.  I no longer
// have complaints about it.
{
  int i,j,k,q=GetQ(),p=GetP();
  double sum;

  if (!ovw)
    {
      ts->ClearOut();
      for (i=0 ; i<nn ; ++i)
        ts->Append(0);
    }
  else
    nn = ts->GetN();

  // uses Durbin-Levinson recursions to simulate from
  // successive one-step predictive d-ns
  // (works for long-memory T-S, unlike the obvious
  //  constructive approach)

  Vector acf(nn), nu(nn), olda(nn), a(nn);
  ComputeACFPACF(acf, NULL, false);  // fill in ACF
 
  nu[0] = acf[0];  // nu(0) = 1-step pred. variance of X(0)
  ts->data[0][0] = norm_random(0,nu[0]);
  a.zeroes();
  for (i=1 ; i<nn ; ++i)
    {
      for (j=0 ; j<nn ; ++j)
	olda[j] = a[j];

      // compute the new a vector
      for (j=1,sum=0.0 ; j<i ; ++j)
	sum += olda[j-1]*acf[i-j];
      a[i-1] = 1/nu[i-1]*(acf[i]-sum);
      for (k=0 ; k<i-1 ; ++k)
	a[k] = olda[k] - a[i-1]*olda[i-2-k];

      // update nu
      nu[i] = nu[i-1]*(1-a[i-1]*a[i-1]);

      // compute xhat
      for (j=0,sum=0.0 ; j<i ; ++j)
	sum += a[j]*ts->data[0][i-1-j];
      ts->data[0][i] = norm_random(sum, nu[i]);
    }
  for (i=0 ; i<nn ; ++i)
    ts->data[0][i] += mean;
}

Matrix ARMAModel::SimulatePredictive(TimeSeries &ts, int nsteps, int nsims)
{
  // note: handling of missing values is not optimal.
  //       they are just filled in with draws from one-step predictive d-ns
  //       (instead of the full-conditional distributions)

  int i,j,k,sim_num,q=GetQ(),p=GetP(), n0=ts.GetN(), nn=n0 + nsteps;
  double sum;

  // uses Durbin-Levinson recursions to simulate from
  // successive one-step predictive d-ns
  // (works for long-memory T-S, unlike the obvious
  //  constructive approach)

  Vector acf(nn), nu(nn), olda(nn), a(nn);
  ComputeACFPACF(acf, NULL, false);  // fill in ACF
  Vector new_ts(nn);
  Matrix new_sims(nsims, nsteps);

  for (i=0 ; i<n0 ; ++i)
    new_ts[i] = ts.data[0][i] - mean;

  nu[0] = acf[0];  // nu(0) = 1-step pred. variance of X(0)
  if (n0>0)
    if (ts.IsMissing(0))
      new_ts[0] = norm_random(0,nu[0]);
  if (n0==0)
    if (nsteps>0)
      for (sim_num=0 ; sim_num<nsims ; ++sim_num)
	new_sims[sim_num][0] = norm_random(0,nu[0]);
  a.zeroes();
  for (i=1 ; i<nn ; ++i)
    {
      for (j=0 ; j<nn ; ++j)
	olda[j] = a[j];

      // compute the new a vector
      for (j=1,sum=0.0 ; j<i ; ++j)
	sum += olda[j-1]*acf[i-j];
      a[i-1] = 1/nu[i-1]*(acf[i]-sum);
      for (k=0 ; k<i-1 ; ++k)
	a[k] = olda[k] - a[i-1]*olda[i-2-k];

      // update nu
      nu[i] = nu[i-1]*(1-a[i-1]*a[i-1]);

      // compute xhat
      if (i<n0) // nothing but filling in missing values
	{
	  for (j=0,sum=0.0 ; j<i ; ++j)
	    sum += a[j]*new_ts[i-1-j];

	  if (ts.IsMissing(i))
	    new_ts[i] = norm_random(sum, nu[i]);
	}
      else // we need multiple sims
	for (sim_num=0 ; sim_num<nsims ; ++sim_num)
	  {
	    for (j=0,sum=0.0 ; j<i ; ++j)
	      sum += a[j]*(i-j<=n0 ? new_ts[i-1-j] : new_sims[sim_num][i-1-j-n0]);
	    new_sims[sim_num][i-n0] = norm_random(sum, nu[i]);
	  }
    }

  for (sim_num=0 ; sim_num<nsims ; ++sim_num)
    for (i=0 ; i<nsteps ; ++i)
      new_sims[sim_num][i] += mean;

  return new_sims;
}

void ARMAModel::RenderInto(ostringstream &os)
// Create text description of the current model.
// The special characters (bold, greek, etc.)
// are for use in the "TextRenderer" class.
{
  int i;

  os << "\\b";  // make it bold
  os << "AR";
  if (!is_short_mem())
    os << "FI";
  os << "MA(" << p << ",";
  if (!is_short_mem())
    os << fracdiff << ",";
  os << q << ") Model" << endl;
  os << "\\n";  // make it normal
  os << "X(t)";
  for (i=1 ; i<=p ; ++i)
    {
      if (phi[i-1]>0)
        os << setiosflags(ios::fixed) << setprecision(3) << " - " << phi[i-1] << "X";
      else
        os << setiosflags(ios::fixed) << setprecision(3) << " + " << -phi[i-1] << "X";
      os << "(t-" << i << ")";
    }
  os << " = Z(t)";
  for (i=1 ; i<=q ; ++i)
    {
      if (theta[i-1]>0)
        os << setiosflags(ios::fixed) << setprecision(3) << " + " << theta[i-1] << "Z";
      else
        os << setiosflags(ios::fixed) << setprecision(3) << " - " << -theta[i-1] << "Z";
      os << "(t-" << i << ")";
    }

  os << endl << "Mean = " << setprecision(5) << mean << ",  " << endl;
  os << "\u03c3";  // UTF-8 unicode for greek lower case sigma
  os << "\\s2\\n = " << setprecision(5) << sigma*sigma;
  if (loglikelihood_is_valid)
    {
      os << endl << endl << "Log Likelihood = " << GetLogLikelihood() << endl;
      os << "AICC = " << GetAICC() << endl;
    }
  else
    {
      os << endl << endl << "Log Likelihood = <Not Computed>" << endl;
      os << "AICC = <Not Computed>" << endl;
    }
}

bool ARMAModel::RenderCtsVersionInto(ostringstream &output, ostringstream &error, double delta)
{
  error << "Feature not yet supported for ARMA models." << ends;
  return false;
}

double ARMAModel::MinAbsRoot()
{
  int i, p2=0, q2=0;
  double min, tx;
  double *cp;
  complex *roots;

  // first construct autoregressive polynomial
  for (i=0 ; i<p ; ++i)
    if (fabs(phi[i])>zero)
      p2 = i+1;					// p2 is the effective AR order, is < p if phi_p=0
  if (p2==0)
    min = 1e6;                                  // a bad approximation to \infty
  else
    {
      cp = new double[p2+1];
      cp[0] = 1.0;
      for (i=1 ; i<=p2 ; ++i)
	cp[i] = -phi[i-1];
      
      roots = rootsof(p2, cp);
      
      // now find minimum magnitude
      min = abs(roots[0]);
      for (i=1 ; i<=p2 ; ++i)
	{
	  tx = abs(roots[i-1]);
	  if (tx<min)
	    min = tx;
	}
      delete[] cp;
      delete[] roots;
    }

  // then construct moving average polynomial
  for (i=0 ; i<q ; ++i)
    if (fabs(theta[i])>zero)
      q2 = i+1;					// q2 is the effective MA order
  if (q2>0)
    {
      cp = new double[q2+1];
      cp[0] = 1.0;
      for (i=1 ; i<=q2 ; ++i)
	cp[i] = theta[i-1];

      roots = rootsof(q2, cp);

      // now find minimum magnitude
      for (i=1 ; i<=q2 ; ++i)
	{
	  tx = abs(roots[i-1]);
	  if (tx<min)
	    min = tx;
	}
      delete[] cp;
      delete[] roots;
    }

  return min; // this is min of roots both AR and MA polynomials (or 1e6)
}

bool ARMAModel::IsCausal(bool causalify)
{
  int i, p2=0;
  double min, tx, ty, tz;
  bool nonstationary = false;
  complex *roots;
  double *cp;

  // first construct autoregressive polynomial
  for (i=0 ; i<p ; ++i)
    if (fabs(phi[i])>zero)
      p2 = i+1;					// p2 is the effective AR order, is < p if phi_p=0
  if (p2==0)
    return true;

  cp = new double[p+1];
  cp[0] = 1.0;
  for (i=1 ; i<=p ; ++i)
    cp[i] = -phi[i-1];

  roots = rootsof(p2, cp);

  // now find minimum magnitude
  min = abs(roots[0]);
  for (i=1 ; i<=p2 ; ++i)
    {
      tx = abs(roots[i-1]);
      if (tx<min)
        min = tx;
      if (fabs(tx-1.0)<zero)
        nonstationary = true;
    }
  delete[] cp;
  if ((min>1.0001) && (!nonstationary))
    {
      delete[] roots;
      return true;
    }

  // otherwise it is non-causal,
  // we may now force it to be causal if causalify is true

  if (!causalify)
    {
      delete[] roots;
      return false;
    }

  complex *recon;

  for (i=1 ; i<=p2 ; ++i)
    {
      tx = abs(roots[i-1]);
      ty = arg(roots[i-1]);
      if (tx<1.0002)
	// map it into the interval (1.0001,1.0002)
	{
	  tz = (1.0002 - tx)/1.0002;
	  tx = 1.0001 + 0.0001*pow((1-tz),2.0);
	}
      roots[i-1] = polar(tx,ty);
    }

  recon = coeffsof(p2, roots);
  for (i=0 ; i<p2 ; ++i)
    phi[i] = -real(recon[i+1]/recon[0]);

  delete[] recon;
  delete[] roots;

  return true;
}

bool ARMAModel::IsInvertible(bool invertibilify)
{
  int i,q2=0;
  double min, tx, ty, tz;
  complex *roots;
  double *cp;

  // first construct moving average polynomial
  for (i=0 ; i<q ; ++i)
    if (fabs(theta[i])>zero)
      q2 = i+1;					// p2 is the effective MA order
  if (q2==0)
    return true;

  cp = new double[q+1];
  cp[0] = 1.0;
  for (i=1 ; i<=q ; ++i)
    cp[i] = theta[i-1];

  roots = rootsof(q2, cp);

  // now find minimum magnitude
  min = abs(roots[0]);
  for (i=1 ; i<=q2 ; ++i)
    if (abs(roots[i-1])<min)
      min = abs(roots[i-1]);
  delete[] cp;
  if (min>1.0001)
    {
      delete[] roots;
      return true;
    }

  // otherwise we may force it to be causal!
  if (!invertibilify)
    {
      delete[] roots;
      return false;
    }

  complex *recon;

  for (i=1 ; i<=q2 ; ++i)
    {
      tx = abs(roots[i-1]);
      ty = arg(roots[i-1]);

      if (tx<1.0002)
	// map it into the interval (1.0001,1.0002)
	{
	  tz = (1.0002 - tx)/1.0002;
	  tx = 1.0001 + 0.0001*pow((1-tz),2.0);
	}

      roots[i-1] = polar(tx,ty);
    }
  recon = coeffsof(q2, roots);
  for (i=0 ; i<q2 ; ++i)
    theta[i] = real(recon[i+1]/recon[0]);

  delete[] recon;
  delete[] roots;

  return true;
}

bool ARMAModel::CheckModel(bool flipit)
{
  int i;
  // AR and MA polynomials first
  bool causal = IsCausal(flipit), invertible = IsInvertible(flipit);
  if ((!causal) || (!invertible))
    return false;

  // then frac diff parm d
  double d = GetFracDiff(), intpart;
  if (flipit)
    {
      if (d>0.499)
	d=0.499;
      if (d<-0.499)
	d=-0.499;
      SetFracDiff(d);
    }
  else if ((d>=0.5) || (d<=-0.5))
    return false;

  // then sigma >= 0
  double s = GetSigma();
  if (flipit)
    {
      if (s<zero)
	s=zero;
      SetSigma2(s*s);
    }
  else
    if (s<zero)
      return false;

  return true;
}


double penaltymultiplier;
void (*iter_call)(int,void *);


Vector ARMAModel::ValidityConstraint()
{
  double tx;
  Vector retval(3);

  // 3 constraints
  // (a) d in (-1/2,1/2)
  tx = GetFracDiff();
  tx = 0.2499 - tx*tx;
  retval[0] = tx;

  // (b) sigma > 0
  retval[1] = GetSigma();

  // (c) min |root| > 1
  retval[2] = MinAbsRoot()-1.0001;

  return retval;
}


double minus_log_like_for(const gsl_vector *x, void *minparms)
{
  int i;
  double retval, ppenalty;
  MinimizerInfoContainer *ic = (MinimizerInfoContainer *) minparms;

  // take vector x and expand it into a proper parameter vector
  TimeSeriesModel *ap = (TimeSeriesModel *) ic->tsmp;
  Vector curp = ap->ParameterBundle();
  Vector p = TimeSeriesModel::param_expand(x, curp, ic->themask);

  // set parameters in current model
  ap->UnbundleParameters(p);

  // compute penalty for invalid parms
  Vector penalty = ap->ValidityConstraint();
  for (i=0,ppenalty=0.0 ; i<penalty.size() ; ++i)
    ppenalty += penaltymultiplier*(penalty[i] < 0 ? -penalty[i] : 0);

  // now map parameters to valid parms to likelihood evaluation works
  ap->CheckModel(true);  // make sure parameters are OK

  // and compute -log like + penalty for constraints
  retval = -ic->tsmp->ComputeLogLikelihood(ic->tsp) + ppenalty;

  p.set_iomode(DISP_MODE);
  //cout << p << " " << retval << endl;

  // now restore parameters (just in case CheckModel changed them)
  ap->UnbundleParameters(p);

  return retval;
}


//-------------------------------------
// Here's a lot of messy stuff for
// calculating 2nd derivatives of
// the log-likelihood. 
//
// Basic idea: use gsl_deriv_central
//             function, recursively
//-------------------------------------

class mll_parm {
public:
  int i,j;
  const gsl_vector *x;
  MinimizerInfoContainer *minparms;
};

double other_partial(double x, void *parms)
{
  double retval;
  mll_parm *mp = (mll_parm *) parms;
  int np = mp->x->size;
  gsl_vector *z = gsl_vector_alloc(np);
  gsl_vector_memcpy(z, mp->x);
  gsl_vector_set(z, mp->i, x);
  retval = minus_log_like_for(z, mp->minparms);

  gsl_vector_free(z);
  return retval;
}

double ith_partial(gsl_vector *y, void *parms)
{
  double pd, abserr, h, tx;
  const gsl_vector *backup_copy;
  gsl_function F;
  F.function = &other_partial;
  F.params = parms;
  mll_parm *mp = (mll_parm *) parms;
  int i = mp->i, np = y->size;
  tx = gsl_vector_get(y, i);
  gsl_vector *z = gsl_vector_alloc(np);
  gsl_vector_memcpy(z, y);
  backup_copy = mp->x;
  mp->x = z;
  h = fabs(tx/100);
  if (h<0.01)
    h=0.01;
  gsl_deriv_central(&F, tx, h, &pd, &abserr);
  mp->x = backup_copy;
  gsl_vector_free(z);
  return pd;
}

double partial_minus_log_like_for(double x, void *parms)
  // parms is a pointer to a block of 2 pointers: one to an int = index, one to a gsl_vector
  // this should return the i-th partial derivative!
{
  double retval;
  mll_parm *mp = (mll_parm *) parms;
  int np = mp->x->size;
  gsl_vector *y = gsl_vector_alloc(np);
  gsl_vector_memcpy(y, mp->x);
  gsl_vector_set(y, mp->j, x);
  retval = ith_partial(y, parms);
  gsl_vector_free(y);
  return retval;
}

//--- end of messy 2nd derivative calculation stuff ---
//---  (apart from section below at end of MLEFit)  ---


Vector ARMAModel::GetDefaultParameters(Vector& mask, TimeSeries *tsp)
  // can modify mask if necessary
{
  int i;
  Vector parms = ParameterBundle();
  // mean is handled differently
  if (mask[0]) // mean is to be estimated
    {
      mask[0] = 0;
      parms[0] = tsp->SampleMean();
    }

  // initialize all non-masked parameters now
  Vector tacf(1);
  tsp->ComputeSampleACF(&tacf,NULL,0);
  for (i=0 ; i<mask.size() ; ++i)
    if (mask[i])
      if (i==2)
	parms[i] = sqrt(tacf[0]);
      else
	parms[i] = 0.0;
  return parms;
}


Matrix TimeSeriesModel::MLEFit(TimeSeries *ts, Vector& mask,
			       void (*callback)(int, void *), void *cb_parms, 
			       bool get_parameter_cov)
{
  int i, j, nump, total_parameters = mask.size();
  double epsabs=1e-4, size;
  int status;
  bool isvalid;
  gsl_vector *res, *res_copy;
  Vector parms, initialparms;
  MinimizerInfoContainer ic;
  Matrix bigcov(total_parameters, total_parameters);  // used to return cov. matrix of estimator
  bigcov.zeroes();

  iter_call = callback;
  iterationcounter = 0;

  //
  // Basic idea: Use a standard miminization scheme on -(log likelihood), but
  // add a penalty function for degree of constraint violation
  // Successively: minimize, then check to see if constraints are violated,
  //               if they are, increase the penalty function and try again, otherwise stop.
  //
  //

  parms = GetDefaultParameters(mask,ts);
  nump = mask.sum();
  initialparms = parms;

  ic.tsp = ts;
  ic.tsmp = this;
  ic.themask = mask;

  res_copy = gsl_vector_alloc(nump);

  if (nump>0) // i.e. at least one param to estimate
    {
      penaltymultiplier = ts->GetN()/10.0;
      isvalid = false;
      iterationcounter = 0;
      
      while (!isvalid)
	{
	  penaltymultiplier *= 3;
	  UnbundleParameters(initialparms);
	  parms = initialparms;
	  
	  Vector subparms=param_mask(parms,mask);
	  gsl_vector *x = gsl_vector_alloc(nump), *step_size =gsl_vector_alloc(nump);
	  
	  gsl_multimin_function ff;
	  gsl_multimin_fminimizer *solver = 
	    gsl_multimin_fminimizer_alloc
	    (gsl_multimin_fminimizer_nmsimplex, nump);
	  
	  ff.n = nump;
	  ff.params = &ic;
	  ff.f = &minus_log_like_for;
	  
	  for (i=0 ; i<nump ; ++i)
	    gsl_vector_set(x, i, subparms[i]);  // copy in subparms
	  gsl_vector_set_all(step_size, 0.2);
	  gsl_multimin_fminimizer_set(solver,&ff,x,step_size);
	  
	  for (i=0 ; i<iterationlimit ; ++i) { // maximum iterations is constant defined somewhere else
	    status = gsl_multimin_fminimizer_iterate(solver);
	    
	    size = gsl_multimin_fminimizer_size (solver);
	    status = gsl_multimin_test_size (size, epsabs);
	    //cout << "Status = " << status << endl;
	    if (status==GSL_SUCCESS)
	      i=iterationlimit;
	    
	    ++iterationcounter;
	    if (callback!=NULL)
	      (*callback)(0, cb_parms);
	  }

	  error_code = status;
	  
	  res = gsl_multimin_fminimizer_x(solver);
	  gsl_vector_memcpy(res_copy, res);

	  // copy results back into the current model
	  parms = param_expand(res, parms, mask);
	  UnbundleParameters(parms);

	  // free up memory
	  gsl_multimin_fminimizer_free(solver);
	  gsl_vector_free(x);

	  isvalid = CheckModel(false);
	}
      
      // evaluate log-likelihood !!!
      ComputeLogLikelihood(ts);


      // also compute partial derivatives of log-likelihood, 
      // for information matrix to get parameter uncertainty
      Matrix parmcov(nump,nump), sub;
      Vector badones(nump);

      if (get_parameter_cov)
	{
	  badones.zeroes();            // start with no problems in Hessian
	  mll_parm temp_holder;
	  temp_holder.x = res_copy;
	  temp_holder.minparms = &ic;

	  double pd, abserr, h, tx;
	  gsl_function F;
	  F.function = &partial_minus_log_like_for;
	  F.params = &temp_holder;

	  for (i=0 ; i<nump ; ++i)
	    {
	      temp_holder.i = i;
	      for (j=0 ; j<=i ; ++j)
		{
		  if (callback!=NULL)
		    (*callback)(1,cb_parms);
		  temp_holder.j = j;
		  tx = gsl_vector_get(res_copy, temp_holder.j);
		  h = fabs(tx/100);
		  if (h<0.01)
		    h = 0.01;
		  gsl_deriv_central(&F, tx, h, &pd, &abserr);
		  parmcov[temp_holder.i][temp_holder.j] = 
		    parmcov[temp_holder.j][temp_holder.i] = pd;
		}
	      // now check to see if the i'th parameter is "ok", i.e.
	      // its row doesn't make the matrix non-positive-definite
	      sub = parmcov.submatrix(0,0,i+1,i+1);
	      sub.singular(&tx);  // compute determinant
	      if (tx<=0) // we have a problem!  mark this parameter as bad!
		{
		  for (j=0 ; j<i ; ++j)
		    parmcov[i][j] = parmcov[j][i] = 0.0;
		  parmcov[i][i] = 1.0;
		  badones[i] = 1.0;   // mark this as bad!
		}
	    }

	  // translate back into full matrix of parameter estimate covariances
	  parmcov = parmcov.inverse();

	  // and finally, map into full size matrix by padding
	  // with zeroes for non-estimated parameters
	  Matrix tm(total_parameters, nump);
	  tm.zeroes();
	  for (i=0,j=0 ; i<total_parameters ; ++i)
	    if (mask[i])
	      {
		tm[i][j] = 1.0 - badones[j];  // if it's a bad parameter, zero out the row/col of cov. matrix
		++j;
	      }
	  bigcov = tm*parmcov*tm.transpose();
	  for (i=0,j=0 ; i<total_parameters ; ++i)
	    if (mask[i])
	      {
		if (badones[j])
		  bigcov[i][i] = -1.0;       // in fact, mark bad ones with Var=-1
		// this differentiates from Var=0 case
		// which simply means parm. unestimated
		++j;
	      }
	} // end of if get_parameter_cov
      else
	bigcov = parmcov;                // just dummy return matrix with meaningless values in it
    }
  else // WN model
    cout << "Error: no parameters to estimate!" << endl;

  gsl_vector_free(res_copy);

  UnbundleParameters(parms); // restore parms

  return bigcov;
}


// The following functions use the
// general parameter vector format:
//
// parm[0,...,p-1] = phi_1,...phi_p
// parm[p,...,p+q-1] = theta_1,...,theta_q 
// parm[p+q] = sigma
// parm[p+q+1] = d
// parm[p+q+2] = mean
//

const double epsilon=1e-10;  // used to keep roots away from unit circle

double logit(double x)
{
  double tx = exp(x);
  return (1.0 / (1.0+tx));
}

double invlogit(double x)
{
  return log(1.0/x - 1.0);
}

void RCTranslate(double x, double y, complex& c1, complex& c2)
{
  double r,r1,r2,theta;
  if (y<0) // make a complex-conjugate pair
    {
      r = logit(x); // magnitude
      theta = y;
      c1 = polar(r,theta);
      c2 = polar(r,-theta);
    }
  else // make a real pair
    {
      c1 = logit(x+y)*2-1;
      c2 = logit(x-y)*2-1;
    }
}
  
Vector ARMAModel::ParameterBundle()
{
  int i,p=GetP(),q=GetQ();
  Vector pvec(p+q+3);
  for (i=0 ; i<p ; ++i)
    pvec[i+3] = phi[i];
  for (i=0 ; i<q ; ++i)
    pvec[i+p+3] = theta[i];
  pvec[0] = GetMean();
  pvec[1] = GetFracDiff();
  pvec[2] = GetSigma();
  return pvec;
}

void ARMAModel::UnbundleParameters(Vector& pvec)
{
  int p=GetP(),q=GetQ();
  SetMean(pvec[0]);
  SetFracDiff(pvec[1]);
  SetSigma2(pvec[2]*pvec[2]);
  SetCoefficients(&pvec[3], &pvec[p+3]);
}

void ARMAModel::StreamParameterName(ostringstream& os, int pnum)
{
  int p = GetP(), q=GetQ();
  switch (pnum) {
  case 0: 
    os << "\u03bc = mean"; break;
  case 1:
    os << "d = frac. diff. parm"; break;
  case 2: 
    os << "\u03c3"; break;   // sigma
  default:
    if (pnum<3+p)
      os << "\u03c6" << "(" << (pnum-2) << ")";
    else
      os << "\u03b8(" << (pnum-2-p) << ")";
  }
}

Vector ARMAModel::GetMask()
{
  int i,p=GetP(),q=GetQ(),np=p+q+3;
  Vector retval(np);

  for (i=0 ; i<min(estimation_mask.size(),np) ; ++i)
    retval[i] = estimation_mask[i];
  for ( ; i<np ; ++i)
    {
      switch (i)
	{
	case 0:  case 2: retval[i] = 1;   break;
	case 1: retval[i] = 0;   break;
	default: retval[i] = 1;   break;
	}
    }
  return retval;
}

double ARMAModel::LogPrior(bool &isvalid)
{
  isvalid = true;
  double tx = 0.0;

  // check for causality, invertibility
  if ((!IsCausal()) || (!IsInvertible()))
    isvalid = false;

  tx = lognorm_density(log(GetSigma()),0,10000);

  return tx;
}

void ARMAModel::BayesianFit(TimeSeries *ts, Vector& mask, void (*callback)(int,void *), void *cb_parms)
{
  int i,j;
  bool isok, accept;
  Vector pvec, acf(1), localmask;
  Matrix allparms;
  double delta;
  double tx,ty;

  // initialize parameter vector pvec=(\mu,d,\sigma,\phi,\theta)
  // - only initializing unmasked elements
  localmask = GetMask();
  pvec = GetDefaultParameters(mask, ts);
  UnbundleParameters(pvec);

  tx = ComputeLogLikelihood(ts);
  tx += LogPrior(isok);
  if (!isok)
    cout << "ERROR: Initial model not valid." << endl;
  allparms = pvec.transpose();

  // go through the MCMC simulation process
  for (iterationcounter = 0 ; iterationcounter < iterationlimit ;
       ++ iterationcounter)
    {
      // register the iteration
      if (callback!=NULL)
	(*callback)(2,cb_parms);

      // carry out phi and theta updates
      for (i=0 ; i<p+q ; ++i)
	if (localmask[i+3])
	  {
	    delta = norm_random(0,0.004);
	    pvec[i+3] += delta;
	    
	    // unbundle params
	    UnbundleParameters(pvec);
	    isok = CheckModel(false);
	    accept = false;
	    if (isok)
	      {
		ty = LogPrior(isok);
		accept = isok;
		if (isok)
		  {
		    ty += ComputeLogLikelihood(ts);
		    accept = log(unif_random())<(ty-tx);
		  }
	      }
	    if (accept)
	      tx = ty;
	    else
	      pvec[i+3] -= delta;
	  }
      

      // carry out sigma update
      double scale = exp(norm_random(0,0.0001));
      if (localmask[2])
	{
	  pvec[2] *= scale;
	  UnbundleParameters(pvec);
	  ty = LogPrior(isok);
	  accept = isok;
	  if (isok)
	    {
	      ty += ComputeLogLikelihood(ts);
	      accept = (log(unif_random())<(ty-tx));
	    }
	  if (accept)
	    tx = ty;
	  else
	    pvec[2] /= scale;
	}
      
      // carry out mu update
      
      // carry out d update
      
      // make sure parameters are correct
      UnbundleParameters(pvec);
      
      // update matrix
      allparms.append_bottom(pvec, 1);
    } // end of iteration loop

  error_code = GSL_SUCCESS;
  
  // now fitted model is posterior mean
  pvec.zeroes();
  for (i=100 ; i<iterationlimit ; ++i)
    for (j=0 ; j<pvec.size() ; ++j)
      pvec[j] += allparms[i][j];
  pvec = pvec*(1.0/(iterationlimit-100));
  
  UnbundleParameters(pvec);
}

void TimeSeriesModel::StreamMLESummary(ostringstream& supplemental, Matrix& parmcovs, bool get_parameter_cov)
  // fill in info about fit
{
  int i,j;
  Vector ps = ParameterBundle();
  int np = ps.size();
  double corr,tx;

  supplemental << "\\bMaximum Likelihood Estimation Results\\n" << endl << endl;
  supplemental << "\\uParameter estimates\\n" << endl << endl;
  for (i=0 ; i<np ; ++i)
    {
      supplemental << "Param. #" << (i+1) << ": ";
      StreamParameterName(supplemental, i);
      if (get_parameter_cov)
	{
	  supplemental << "\\n = " << ps[i] << " (";
	  if (parmcovs[i][i]>0)
	    supplemental << sqrt(parmcovs[i][i]);
	  if (parmcovs[i][i]==0)
	    supplemental << "NA";
	  if (parmcovs[i][i]==-1.0)
	    supplemental << "US";
	  supplemental << ")";
	}
      supplemental << endl;
    }

  if (get_parameter_cov)
    {
      supplemental << endl << "Standard errors are given in parentheses." << endl;
      supplemental << "NA indicates parameter was held fixed or not estimated by maximizing the likelihood." << endl;
      supplemental << "US indicates that numerical differentation of the log-likelihood failed for the parameter." << endl;

      supplemental << endl << "\\uCorrelations of Parameter Estimates\\n" << endl << endl;

      supplemental << "\\m            ";  // switch to monospacing mode so table entries line up
      for (j=0 ; j<np ; ++j)
	supplemental << setw(2) << (j+1) << "|";
      supplemental << endl;
      for (i=0 ; i<np ; ++i)
	{
	  supplemental << "Param. #" << setw(2) << (i+1) <<" |";
	  for (j=0 ; j<np ; ++j)
	    {
	      corr = parmcovs[i][j];
	      if (corr==-1.0)
		corr = 0;
	      tx = (parmcovs[i][i]*parmcovs[j][j]);
	      if ((parmcovs[i][i]==-1.0) || (parmcovs[j][j]==-1.0))
		tx = 0.0;
	      if (fabs(tx)<zero)
		supplemental << "--|";
	      else
		{
		  tx = sqrt(tx);
		  corr /= tx;
		  if (corr<-0.7)
		    supplemental << "\\0  \\m|";
		  else if (corr<-0.4)
		    supplemental << "\\1  \\m|";
		  else if (corr<-0.1)
		    supplemental << "\\2  \\m|";
		  else if (corr<0.1)
		    supplemental << "  |";
		  else if (corr<0.4)
		    supplemental << "\\3  \\m|";
		  else if (corr<0.7)
		    supplemental << "\\4  \\m|";
		  else
		    supplemental << "\\5  \\m|";
		}
	    }
	  supplemental << endl;
	}

      supplemental << endl << "\\nCorrelation Legend:\\m" << endl;
      supplemental << "|\\0  \\m| = [-1.0,-0.7]" << endl << "|\\1  \\m| = [-0.7,-0.4]" << endl << "|\\2  \\m| = [-0.4,-0.1]" << endl;
      supplemental << "|\\m  \\m| = [-0.1, 0.1]" << endl;
      supplemental << "|\\3  \\m| = [ 0.1, 0.4]" << endl << "|\\4  \\m| = [ 0.4, 0.7]" << endl << "|\\5  \\m| = [ 0.7, 1.0]" << endl << endl;
    }

  supplemental << ends;
}

int ARMAModel::FitModel(TimeSeries *ts, const int method, const int ilimit,
			void (*callback)(int,void *), void *cb_parms, 
			ostringstream& msg, ostringstream& supplemental, bool get_parameter_cov)
{
  int i,j,retval;
  Vector parms = ParameterBundle(), mask;
  Matrix parmcovs;
  double corr, tx;

  // create mask
  mask = GetMask();
  iterationlimit = ilimit;

  if (method==METHOD_MLE)
    {
      if (ts->HasAnyMissing())
	if (!is_short_mem())
	  {
	    msg << "Unable to fit fractionally differenced model by maximum likelihood estimation when time series has missing observations."
		<< ends;
	    return UNABLE;
	  }
      
      parmcovs = MLEFit(ts, mask, callback, cb_parms, get_parameter_cov);
      StreamMLESummary(supplemental, parmcovs, get_parameter_cov);

      msg << "Fitted model computed using " << iterationcounter
	  << " likelihood evaluations." << endl;
      if (error_code!=GSL_SUCCESS)
	msg << "Warning: Optimizer did not converge within specified tolerance." 
	    << endl;
      msg << ends;
    }
  if (method==METHOD_BAYESIAN)
    {
      BayesianFit(ts, mask, callback, cb_parms);
      msg << "MCMC estimation carried out with " << ilimit
	  << " iterations." << endl << endl;
      supplemental << "Analysis of MCMC output not yet implemented." << ends;
    }

  retval = SUCCESS;
  if (error_code!=GSL_SUCCESS)
    retval = NONCONVERGENT;
  return retval;
}


void ARMAModel::ComputeACFPACF(Vector &acf, Vector *pacf, bool normalizeacf)
// Method 3: B\&D, p. 95
{
  int i,h,j,k,maxlag = acf.nrows()-1;
  double tx, ty;
  double *psi = new double[q+1];

  if (is_short_mem()) // i.e. if it's not fractionally integrated
    {
      // first step: compute $\psi_0,\ldots,\psi_q$.
      psi[0] = 1.0;
      for (i=1 ; i<=q ; ++i)
	{
	  if (i<=q)
	    psi[i] = theta[i-1];
	  else
	    psi[i] = 0.0;
	  for (k=1 ; k<=i ; ++k)
	    if (k<=p)
	      psi[i] += phi[k-1]*psi[i-k];
	}

      // second step: solve Y-W for $\gamma(0),\ldots,\gamma(p)$
      Matrix phistuff(p+1,p+1);
      Vector gammas(p+1), psistuff(p+1);
      phistuff.zeroes();
      for (i=0 ; i<=p ; ++i)
	for (j=0 ; j<=p ; ++j)
	  {
	    tx = 1.0;
	    if (j>0)
	      if (j<=p)
		tx = -phi[j-1];
	      else
		tx = 0;
	    phistuff[i][abs(i-j)] += tx;
	  }
      for (i=0 ; i<=p ; ++i)
	{
	  psistuff[i] = 0;
	  for (j=i ; j<=q ; ++j)
	    {
	      tx = 1.0;
	      if (j>0)
		tx = theta[j-1];
	      if (j>q)
		tx = 0.0;
	      psistuff[i] += tx*psi[j-i];
	    }
	  psistuff[i] *= sigma*sigma;
	}

      gammas = psistuff;
      phistuff.solve_system(&gammas);

      // copy into return array
      for (i=0 ; i<=min(p,maxlag) ; ++i)
	acf[i] = gammas[i];

      // and then compute the rest recursively
      for (i=p+1 ; i<=maxlag ; ++i)
	{
	  for (j=1, tx=0.0 ; j<=p ; ++j)
	    tx += phi[j-1]*acf[i-j];
	  if (i<max(p,q+1))
	    for (j=i ; j<=q ; ++j)
	      {
		ty = 1.0;
		if (j>0)
		  ty = theta[j-1];
		if (j>q)
		  ty = 0.0;
		tx += sigma*sigma*ty*psi[j-i];
	      }
	  acf[i] = tx;
	}
    } // end of normal fracdiff==0.0 case
  else
    {
      // fractionally differenced case
      Vector psis(q+1);
      Matrix zetas(p+1,2);
      complex *roots, tc, tc2, tc3;
      double *cp;
      double tx, tx2;
      int p2,l;
      p2 = 0;
      for (i=0 ; i<p ; ++i)
	if (fabs(phi[i])>zero)
	  p2 = i+1;				// p2 is the effective AR order, is < p if phi_p

      // compute psis
      for (i=0 ; i<=q ; ++i)
	for (psis[i]=0, j=i ; j<=q ; ++j)
	  {
	    tx = j>0 ? theta[j-1] : 1.0;
	    tx2 = j>i ? theta[j-1-i] : 1.0;
	    psis[i] += tx*tx2;
	  }

      // get AR polynomial roots
      cp = new double[p2+1];
      cp[0] = 1.0;
      for (i=1 ; i<=p2 ; ++i)
	cp[i] = -phi[i-1];
      if (p2>0)
	roots = rootsof(p2, cp);
      else
	roots = NULL;
      delete[] cp;

      // transform to inverse roots
      for (i=0 ; i<p2 ; ++i)
	roots[i] = 1.0/roots[i];

      // get zetas
      for (i=0 ; i<p2 ; ++i)
	{
	  for (j=0,tc=1.0 ; j<p2 ; ++j)
	    tc *= (1.0-roots[i]*roots[j]);
	  for (j=0 ; j<p2 ; ++j)
	    if (j!=i)
	      tc *= (roots[i]-roots[j]);
	  tc = 1.0/tc;
	  zetas[i][0] = tc.real();
	  zetas[i][1] = tc.imag();
	}

      // compute gamma
      Matrix C;
      if (p2==0)
	{
	  acf.zeroes();
	  for (l=-q ; l<=q ; ++l)
	    {
	      int al = l<0 ? -l : l;
	      tx = gsl_sf_gamma(1-2*fracdiff)*gsl_sf_gamma(fracdiff+l) /
		( gsl_sf_gamma(fracdiff)*gsl_sf_gamma(1.0-fracdiff)*gsl_sf_gamma(1.0-fracdiff+l) );
	      acf[0] += psis[al]*tx;
	      for (h=1 ; h<=maxlag ; ++h)
		{
		  tx *= (1.0-fracdiff-h+l)/(fracdiff-h+l);
		  acf[h] += psis[al]*tx;
		}
	      acf = acf*sigma*sigma;
	    }
	}
      else
	{
	  acf.zeroes();
	  for (j=0 ; j<p2 ; ++j)
	    {
	      C = cfuncsfor(fracdiff, roots[j].real(), roots[j].imag(), 
			    p2, p2+q, 2*q+maxlag+1);
	      for (l=-q ; l<=q ; ++l)
		for (h=0 ; h<=maxlag ; ++h)
		  {
		    tc = complex(zetas[j][0],zetas[j][1]);
		    tc2 = complex(C[h+q-l][0], C[h+q-l][1]);
		    tc3 = sigma*sigma * psis[l<0 ? -l : l] * tc * tc2;
		    acf[h] += tc3.real();  // we know imag. parts will cancel!
		  }
	    }
	}

      delete[] roots;  // free mem when finished with them
    }

  if (normalizeacf)
    for (i=maxlag ; i>=0 ; --i)
      acf[i] /= acf[0];

  // Second part: compute the model partial autocorrelation function
  if (pacf!=NULL)
    ACFtoPACF(acf, *pacf);

  delete[] psi;
}

double ARMAModel::kappa(int i, int j, int m, Vector& acf)
// this is autocovariance of the modified X_t process (Ansley 1979)
{
  int m1 = min(i,j), m2 = max(i,j), ti, r;
  double sum, tx, t1, t2;

  if (m2 <= m)
    return acf[abs(i-j)]/(sigma*sigma);
  if (m1 > m)
    {
      sum = 0.0;
      for (r=0 ; r <= q ; ++r)
        {
          if (r==0)
            t1 = 1.0;
          else
            t1 = theta[r-1];
          ti = r+abs(i-j);
          if (ti==0)
            t2 = 1.0;
          else if (ti<=q)
            t2 = theta[ti-1];
          else
            t2 = 0.0;
          sum += t1*t2;
        }
      return sum;
    }
  if ((m1 <= m) && (m2 > m) && (m2 <= 2*m))
    {
      sum = acf[abs(i-j)];
      for (r=1 ; r<=p ; ++r)
        {
          ti = r-abs(i-j);
          sum -= phi[r-1]*acf[abs(ti)];
        }
      return sum/(sigma*sigma);
    }
  return 0;
}

void ARMAModel::ForecastWithMissing(TimeSeries *tp, int nsteps, 
				    double *forecret, double *fmse)
{
  // This function uses the state-space model formulation of the ARMA model
  // to get one-step predictors, residuals, and additional forecasts.
  // See comments under Forecast(...) below for details on parameters
  // If fracdiff!=0 then this routine fails!
  int i,j,t,nobs=tp->GetN();

  if (is_short_mem())
    {
      // set up standard SS model representation for ARMA model
      // (note: ARMA model is zero-mean version)
      Matrix F,G,Q;
      Vector resids(nobs),rs(nobs+nsteps); 
      double yhat,ysig2,ri,ss,slri;
      int nonmissing = 0;

      int r=max(p,q+1);
      F.resize(r,r);
      F.zeroes();
      for (i=0 ; i<r-1 ; ++i)
	F[i][i+1] = 1.0;
      for (i=0 ; i<p ; ++i)
	F[r-1][r-1-i] = phi[i];
      G.resize(1,r);
      G.zeroes();
      for (i=1 ; i<=q ; ++i)
	G[0][r-1-i] = theta[i-1];
      G[0][r-1] = 1.0;
      Q.resize(r,r);
      Q.zeroes();
      Q[r-1][r-1] = sigma*sigma;

      // now run the Kalman Filter
      Matrix Sigmat(r,r);  // filtering d-n variance
      Vector mut(r);       // filtering d-n mean 
      Matrix Ktp1,Ptp1;
      
      // initialize them
      mut.zeroes();
      Sigmat = Q;
      for (i=0 ; i<KFburnin ; ++i)
	{
	  mut =  F*mut;
	  Sigmat = F*Sigmat*F.transpose() + Q;
	}
      // now we have mut and Sigmat for t==-1
 
      loglikelihood = ss = slri = 0.0;  // really this will be reduced log-likelihood
      for (t=0 ; t<nobs ; ++t)
	// update for t
	{
	  // 1. add in contribution to reduced log-likelihood 
	  // (see B&D, eq. 8.7.4, from
	  //    one-step predictive based on current mut, Sigmat
	  yhat = (G*F*mut)[0];
	  ysig2 = (G*(F*Sigmat*F.transpose()+Q)*G.transpose())[0][0];
	  ri = ysig2/(sigma*sigma);
	  rs[t] = ri;
	  if (!tp->IsMissing(t))
	    {
	      ++nonmissing;
 	      slri += log(ri);
	      resids[t] = (yhat-(tp->data[0][t]-mean))/sqrt(ri);
	      ss += resids[t]*resids[t];
	    }
	  else
	    // else just add the constant 0.0 to loglikelihood, and set resid=0
	    resids[t] = 0.0;
	  
	  // 2. iterate to obtain next mut, Sigmat
	  if (!tp->IsMissing(t))
	    {
	      Ptp1 = F*Sigmat*F.transpose() + Q;
	      double g = (G*Ptp1*G.transpose())[0][0];
	      Ktp1 = Ptp1*G.transpose()*(1.0/g);
	      Sigmat = Ptp1 - Ktp1*G*Ptp1;  // Sigmat<-Sigma_{t+1}
	      mut = F*mut + Ktp1*((tp->data[0][t]-mean) - (G*F*mut)[0]);
	    }
	  else
	    {
	      mut =  F*mut;
	      Sigmat = F*Sigmat*F.transpose() + Q;
	    }
	}

      for (i=0 ; i<nsteps ; ++i) // fill in forecasts if necessary
	{
	  yhat = (G*F*mut)[0];
	  ysig2 = (G*(F*Sigmat*F.transpose()+Q)*G.transpose())[0][0];
	  forecret[i] = yhat+mean;
	  fmse[i] = ysig2;
	  mut =  F*mut;
	  Sigmat = F*Sigmat*F.transpose() + Q;
	}
    }
  else
    {
      cout << "Missing data not supported when d not equal to 0." << endl;
    }
}

void ARMAModel::ComputeSpecialResidualsWithMissing(TimeSeries *tp, TimeSeries *residreturn, double *rstore)
{
  // This function uses the state-space model formulation of the ARMA model
  // to get one-step predictors, residuals, and additional forecasts.
  // See comments under Forecast(...) below for details on parameters
  // If fracdiff!=0 then this routine fails!
  int i,j,t,nobs=tp->GetN();

  if (is_short_mem())
    {
      // set up standard SS model representation for ARMA model
      // (note: ARMA model is zero-mean version)
      Matrix F,G,Q;
      Vector resids(nobs),rs(nobs); 
      double yhat,ysig2,ri,ss,slri;
      int nonmissing = 0;

      int r=max(p,q+1);
      F.resize(r,r);
      F.zeroes();
      for (i=0 ; i<r-1 ; ++i)
	F[i][i+1] = 1.0;
      for (i=0 ; i<p ; ++i)
	F[r-1][r-1-i] = phi[i];
      G.resize(1,r);
      G.zeroes();
      for (i=1 ; i<=q ; ++i)
	G[0][r-1-i] = theta[i-1];
      G[0][r-1] = 1.0;
      Q.resize(r,r);
      Q.zeroes();
      Q[r-1][r-1] = sigma*sigma;

      // now run the Kalman Filter
      Matrix Sigmat(r,r);  // filtering d-n variance
      Vector mut(r);       // filtering d-n mean 
      Matrix Ktp1,Ptp1;
      
      // initialize them
      mut.zeroes();
      Sigmat = Q;
      for (i=0 ; i<KFburnin ; ++i)
	{
	  mut =  F*mut;
	  Sigmat = F*Sigmat*F.transpose() + Q;
	}
      // now we have mut and Sigmat for t==-1
 
      loglikelihood = ss = slri = 0.0;  // really this will be reduced log-likelihood
      for (t=0 ; t<nobs ; ++t)
	// update for t
	{
	  // 1. add in contribution to reduced log-likelihood 
	  // (see B&D, eq. 8.7.4, from
	  //    one-step predictive based on current mut, Sigmat
	  yhat = (G*F*mut)[0];
	  ysig2 = (G*(F*Sigmat*F.transpose()+Q)*G.transpose())[0][0];
	  ri = ysig2/(sigma*sigma);
	  rs[t] = ri;
	  if (!tp->IsMissing(t))
	    {
	      ++nonmissing;
 	      slri += log(ri);
	      resids[t] = (yhat-(tp->data[0][t]-mean))/sqrt(ri);
	      ss += resids[t]*resids[t];
	    }
	  else
	    // else just add the constant 0.0 to loglikelihood, and set resid=0
	    resids[t] = 0.0;
	  
	  // 2. iterate to obtain next mut, Sigmat
	  if (!tp->IsMissing(t))
	    {
	      Ptp1 = F*Sigmat*F.transpose() + Q;
	      double g = (G*Ptp1*G.transpose())[0][0];
	      Ktp1 = Ptp1*G.transpose()*(1.0/g);
	      Sigmat = Ptp1 - Ktp1*G*Ptp1;  // Sigmat<-Sigma_{t+1}
	      mut = F*mut + Ktp1*((tp->data[0][t]-mean) - (G*F*mut)[0]);
	    }
	  else
	    {
	      mut =  F*mut;
	      Sigmat = F*Sigmat*F.transpose() + Q;
	    }
	}
      
      if (residreturn!=NULL)
	{
	  residreturn->ClearOut();
	  for (i=0 ; i<nobs ; ++i)
	    residreturn->Append(resids[i], tp->IsMissing(i));
	}
      if (rstore!=NULL)
	for (i=0 ; i<nobs ; ++i)
	  rstore[i] = rs[i];
      loglikelihood_is_valid = true;
    }
  else
    {
      // fill in garbage
      if (residreturn!=NULL)
	{
	  residreturn->ClearOut();
	  for (i=0 ; i<nobs ; ++i)
	    residreturn->Append(0, false);
	}
      if (rstore!=NULL)
	for (i=0 ; i<nobs ; ++i)
	  rstore[i] = 1.0;
      loglikelihood_is_valid = false;
    }
}

void ARMAModel::ComputeStdResiduals(TimeSeries *tp, TimeSeries *resids)
{
  int i,n = tp->GetN();
  double *rret = new double[n];

  ComputeSpecialResiduals(tp, resids, rret);
  for (i=0 ; i<n ; ++i)
    (*resids)[i] /= sigma;

  delete[] rret;
}

void ARMAModel::ComputeSpecialResiduals(TimeSeries *tp, TimeSeries *rret, double *rstore)
{
  // This function uses the approach described
  // in Section 8.7 of B&D if fracdiff parm==0
  // to get one-step predictors, residuals, and additional forecasts

  // If fracdiff!=0 then it uses the D-L algorithm
  // with ACF computed by Sowell's formula.

  // 1.  resids are stored in rret if it's non-NULL
  // 2.  rs[...] are one-step predictive MSEs / sigma^2,
  //     based on sigma when this routine was called
  //     rs[...] are put in rstore if it's non-NULL
  // 3.  this function returns (*rret) =  [ (x_i - \hat{x}_i)/\sqrt{rs_{i-1}} ]

  if (tp->HasAnyMissing())
    {
      ComputeSpecialResidualsWithMissing(tp,rret,rstore);
      return;
    }

  int h,i,j,k,n,nobs = tp->GetN();
  Vector acf, xhat(max(nobs+1,1)), residuals(nobs);
  double t1,t2,sum,tx;
  Vector rs(nobs+1);
  double *tsdata = tp->GetData();

  if (is_short_mem())
    {
      int m=max(p,q);
      acf.resize(m+1);
      // first we need to compute Gamma(0..m)
      ComputeACFPACF(acf, NULL, 0);

      // then apply the innovations algorithm to compute $theta_{n,j}$s
      int minwidth = max(max(q,m-1),1);
      Matrix thetas(nobs, minwidth);
      
      rs[0] = kappa(1,1,m,acf);
      for (n=1 ; n<=nobs ; ++n)
	{
	  // first find $\theta_{n,q},\ldots,\theta_{n,1}$
	  for (k=max(n-minwidth,0) ; k<n ; ++k)
	    {
	      sum = kappa(n+1,k+1,m,acf);
	      for (j=max(n-minwidth,0) ; j<k ; ++j)
		{
		  if ((k-j-1)<minwidth)
		    t1 = thetas[k-1][k-j-1];
		  else
		    t1 = 0.0;
		  t2 = thetas[n-1][n-j-1];
		  sum -= t1*t2*rs[j];
		}
	      thetas[n-1][n-k-1] = sum/rs[k];
	    }
	  
	  // then find $r_n$
	  rs[n] = kappa(n+1,n+1,m,acf);
	  for (j=max(n-minwidth,0) ; j<n ; ++j)
	    {
	      t1 = thetas[n-1][n-j-1];
	      rs[n] -= t1*t1*rs[j];
	    }
	}
      
      // now compute all the one-step predictors and the residuals    
      xhat[0] = mean;
      for (n=1 ; n<nobs; ++n)
	{
	  // work out $\hat{X}_{n+1}$
	  sum = 0.0;
	  if (n<m)
	    {
	      for (j=1 ; j<=min(n,minwidth) ; ++j)
		sum += thetas[n-1][j-1] * (tsdata[n-j] - xhat[n-j]);
	    }
	  else
	    {
	      for (j=1 ; j<=p ; ++j)
		sum += phi[j-1] * (tsdata[n-j] - mean);
	      for (j=1 ; j<=minwidth ; ++j)
		sum += thetas[n-1][j-1] * (tsdata[n-j] - xhat[n-j]);
	    }
	  xhat[n] = sum + mean;
	}
    } // end of fracdiff==0.0 case
  else
    {
      // USE the standard Durbin-Levison Algorithm
      acf.resize(nobs);
      // first we need to compute Gamma(0..m)
      ComputeACFPACF(acf, NULL, 0);
      
      // Now apply the Durbin-Levinson algorithm to compute 
      // all the one-step predictors and their variances
      int ncols=nobs;
      Vector a(ncols), olda(ncols);

      rs[0] = acf[0];  // nu(0)
      xhat[0] = mean;
      for (i=1 ; i<=ncols ; ++i)
	{
	  for (j=0 ; j<ncols ; ++j)
	    olda[j] = a[j];

	  // compute the new a vector
	  for (j=1,sum=0.0 ; j<i ; ++j)
	    sum += olda[j-1]*acf[i-j];
	  a[i-1] = 1/rs[i-1]*(acf[i]-sum);
	  for (k=0 ; k<i-1 ; ++k)
	    a[k] = olda[k] - a[i-1]*olda[i-2-k];

	  // update nu
	  rs[i] = rs[i-1]*(1-a[i-1]*a[i-1]);

	  for (j=0,tx=0.0 ; j<i ; ++j)
	    tx += a[j]*
	      (i-1-j<nobs ? (tsdata[i-1-j]-mean) : (xhat[i-1-j]-mean));
	  xhat[i] = tx + mean;
	}

      for (i=0 ; i<ncols ; ++i)
	rs[i] = rs[i]/(sigma*sigma);
   }

  // store the results
  for (n=0 ; n<nobs ; ++n)
    residuals[n] = (tsdata[n] - xhat[n])/sqrt(rs[n]); 

  // construct new time series for the residuals
  if (rret != NULL)
    {
      rret->ClearOut();
      for (i=0 ; i<nobs ; ++i)
        rret->Append(residuals[i]);
    }

  if (rstore!=NULL)
    for (i=0 ; i<nobs ; ++i)
      rstore[i] = rs[i];

  loglikelihood_is_valid = true;
}

void ARMAModel::Forecast(TimeSeries *tp, int nsteps, double *fret, double *fmse)
{
  // This function uses the approach described
  // in Section 8.7 of B&D if fracdiff parm==0
  // to get one-step predictors, residuals, and additional forecasts

  // If fracdiff!=0 then it uses the D-L algorithm
  // with ACF computed by Sowell's formula.

  // 1.  forecasts (1..nsteps predictive means) are stored in fret if it's non-NULL
  // 2.  predictive variances are returned in fmse if it's non-NULL


  if (tp->HasAnyMissing())
    {
      ForecastWithMissing(tp,nsteps,fret,fmse);
      return;
    }

  int h,i,j,k,n,nobs = tp->GetN();
  Vector acf, xhat(max(nobs+nsteps+1,1)), residuals(nobs);
  double t1,t2,sum,tx;
  Vector rs(nobs+1+nsteps), local_fmse(nsteps+1);
  double *tsdata = tp->GetData();

  if (is_short_mem())
    {
      int m=max(p,q);

      // first we need to compute Gamma(0..m)
      acf.resize(m+1);
      ComputeACFPACF(acf, NULL, 0);

      // then apply the innovations algorithm to compute $theta_{n,j}$s
      int minwidth = max(max(q,m-1),1);
      Matrix thetas(nobs+nsteps, minwidth);
      
      rs[0] = kappa(1,1,m,acf);
      for (n=1 ; n<=nobs+nsteps ; ++n)
	{
	  // first find $\theta_{n,q},\ldots,\theta_{n,1}$
	  for (k=max(n-minwidth,0) ; k<n ; ++k)
	    {
	      sum = kappa(n+1,k+1,m,acf);
	      for (j=max(n-minwidth,0) ; j<k ; ++j)
		{
		  if ((k-j-1)<minwidth)
		    t1 = thetas[k-1][k-j-1];
		  else
		    t1 = 0.0;
		  t2 = thetas[n-1][n-j-1];
		  sum -= t1*t2*rs[j];
		}
	      thetas[n-1][n-k-1] = sum/rs[k];
	    }
	  
	  // then find $r_n$
	  rs[n] = kappa(n+1,n+1,m,acf);
	  for (j=max(n-minwidth,0) ; j<n ; ++j)
	    {
	      t1 = thetas[n-1][n-j-1];
	      rs[n] -= t1*t1*rs[j];
	    }
	}
      
      // now compute all the one-step predictors and the residuals    
      xhat[0] = mean;
      for (n=1 ; n<nobs; ++n)
	{
	  // work out $\hat{X}_{n+1}$
	  sum = 0.0;
	  if (n<m)
	    {
	      for (j=1 ; j<=min(n,minwidth) ; ++j)
		sum += thetas[n-1][j-1] * (tsdata[n-j] - xhat[n-j]);
	    }
	  else
	    {
	      for (j=1 ; j<=p ; ++j)
		sum += phi[j-1] * (tsdata[n-j] - mean);
	      for (j=1 ; j<=minwidth ; ++j)
		sum += thetas[n-1][j-1] * (tsdata[n-j] - xhat[n-j]);
	    }
	  xhat[n] = sum + mean;
	}
      
      // now compute the h-step predictors and their MSEs: 
      //   store the MSEs(/sig^2) in local_fmse[0,...,h-1]
      Vector psis(nsteps);
      psis[0] = 1.0;            // other components are updated in loop below
      for (h=1 ; h<=nsteps ; ++h)
	{
	  // first get P_n X_{n+h}
	  sum = 0.0;
	  if (h<=m-nobs)
	    {
	      for (j=h ; j<=nobs+h-1 ; ++j)
		if (j-1<minwidth)
		  sum += thetas[nobs+h-2][j-1] *
		    (tsdata[nobs+h-j-1] - xhat[nobs+h-j-1]);
	    }
	  else
	    {
	      for (i=1 ; i<=p ; ++i)
		{
		  if (h>i)
		    sum += phi[i-1]*(xhat[nobs+h-i-1]-mean);
		  else
		    sum += phi[i-1]*(tsdata[nobs+h-i-1]-mean);
		}
	      for (j=h ; j<=minwidth ; ++j)
		sum += thetas[nobs+h-2][j-1] *
		  (tsdata[nobs+h-j-1] - xhat[nobs+h-j-1]);
	    }
	  xhat[nobs+h-1] = sum + mean;

	  // then get its MSE: using approximation (5.3.24) from B&D
	  if (h==1)
	    local_fmse[h-1] = 1.0;
	  else
	    {
	      // compute psis[h-1]
	      j=h-1;
	      psis[j] = (j-1 < q ? theta[j-1] : 0);
	      for (k=1 ; k<=p ; ++k)
		psis[j] += phi[k-1]*(j>=k ? psis[j-k] : 0);
	      // update mse
	      local_fmse[h-1] = local_fmse[h-2] + psis[h-1]*psis[h-1];
	    }
	}
      // now fix local_fmse
      local_fmse = local_fmse * sigma*sigma;

    } // end of fracdiff==0.0 case
  else
    {
      // USE the standard Durbin-Levison Algorithm
      acf.resize(nobs+nsteps);
      // first we need to compute Gamma(0..m)
      ComputeACFPACF(acf, NULL, 0);
      
      // Now apply the Durbin-Levinson algorithm to compute 
      // all the one-step predictors and their variances
      int ncols=nobs+nsteps;
      Vector a(ncols), olda(ncols);

      rs[0] = acf[0];  // nu(0)
      xhat[0] = mean;
      for (i=1 ; i<=ncols ; ++i)
	{
	  for (j=0 ; j<ncols ; ++j)
	    olda[j] = a[j];

	  // compute the new a vector
	  for (j=1,sum=0.0 ; j<i ; ++j)
	    sum += olda[j-1]*acf[i-j];
	  a[i-1] = 1/rs[i-1]*(acf[i]-sum);
	  for (k=0 ; k<i-1 ; ++k)
	    a[k] = olda[k] - a[i-1]*olda[i-2-k];

	  // update nu
	  rs[i] = rs[i-1]*(1-a[i-1]*a[i-1]);

	  for (j=0,tx=0.0 ; j<i ; ++j)
	    tx += a[j]*
	      (i-1-j<nobs ? (tsdata[i-1-j]-mean) : (xhat[i-1-j]-mean));
	  xhat[i] = tx + mean;
	}

      for (i=0 ; i<ncols ; ++i)
	rs[i] = rs[i]/(sigma*sigma);

      // now compute MSEs: use same approx. as above, but
      // compute psis by phi(B) (1-B)^d psi(B) = theta(B) and matching coefficients!
      Vector psis(nsteps),xis(nsteps);
      psis[0] = 1.0;
      xis[0] = 1.0;

      // first compute xis = expansion (1-B)^d = \sum B^j xis[j]
      for (i=1 ; i<nsteps ; ++i)
	xis[i] = -xis[i-1]*(fracdiff-(i-1))/((double)i);

      // then do psis: the recursion I derived is
      //
      //  \psi_k = \theta_k + \sum_{m=0}^{k-1} \psi_m \sum_{i=0}^{\max(p,k-m)} \phi_i \xi_{k-m-i},
      //  with the convention that \phi_0 = -1.0
      for (k=1 ; k<nsteps ; ++k)
	{
	  psis[k] = k-1 < q ? theta[k-1] : 0;
	  for (j=0 ; j<k ; ++j)
	    {
	      sum = 0.0;
	      for (i=0 ; i<=max(p,k-j) ; ++i)
		sum += (i==0 ? -1.0 : phi[i-1]) * xis[k-j-i];
	      psis[k] += psis[j] * sum;
	    }
	} 

      // finally, use approx. (B&D eqn (5.3.24)) as before
      local_fmse[0] = 1.0;
      for (i=1 ; i<nsteps ; ++i)
	local_fmse[i] = local_fmse[i-1] + psis[i]*psis[i];
      local_fmse = local_fmse * sigma*sigma;
   }

  // actual forecasts
  if (fret != NULL)
    for (i=1 ; i<=nsteps ; ++i)
      fret[i-1] = xhat[nobs-1+i];

  // fill in the predictive variances
  if (fmse!=NULL)
    for (i=0 ; i<nsteps ; ++i)
      fmse[i] = local_fmse[i];
    
}

double ARMAModel::ComputeLogLikelihood(TimeSeries *tp)
  // computes loglikelihood and AICC, possibly using data with
  // missing values (fracdiff not yet supported for missing vals),
  // returns loglikelihood, also sets loglikelihood and AICC
  // members of ARMAModel class
{
  TimeSeries resids;
  int i,nobs = tp->GetN(),nonmissing;
  double *rs = new double[nobs+1], *res;
  loglikelihood = 0.0;
  if (nobs==0)
    return 0.0;  // no data!


  if (!using_whittle)
    {
      ComputeSpecialResiduals(tp, &resids, rs);
      res = resids.GetData();

      // now it is easy to compute the likelihood
      // this is a straightforward implementation of eqn (8.7.4) in Brockwell & Davis,
      // Time Series: Theory and Methods (2nd edition)
      loglikelihood = 0;
      for (i=0,nonmissing=0 ; i<nobs ; ++i)
	if (!tp->IsMissing(i))
	  {
	    ++nonmissing;
	    loglikelihood -= log(rs[i])/2;
	    loglikelihood -= res[i]*res[i]/(2*sigma*sigma);
	  }
      loglikelihood -= log(2*M_PI*sigma*sigma)*nonmissing/2.0;
    }
  else
    {
      // use Whittle's approximation to likelihood!
      loglikelihood = nobs*log(2*M_PI) + 2*nobs*log(sigma);
      Vector omegas, pdgm, actualspecdens;
      bool using_cache = false;
      if (whittle_cached)
	if ((*tp)==whittle_ts)
	  using_cache = true;
	else
	  whittle_cached = false;  // reset cache when time series doesnt match
      if (using_cache)
	{
	  pdgm = whittle_pdgm;
	  omegas.resize(nobs);
	  for (i=0 ; i<nobs ; ++i)
	    omegas[i] = 2*M_PI*i/nobs;
	}
      else  
	{
	  pdgm = tp->ComputePeriodogram(&omegas);
	  whittle_pdgm = pdgm;
	  whittle_cached = true;
	  whittle_ts = (*tp);
	}
	  
      omegas[0] = omegas[1];  // no zero frequency!
      actualspecdens = ComputeSpectralDensity(omegas);

      actualspecdens *= (2*M_PI/(sigma*sigma));
      for (i=1 ; i<pdgm.size() ; ++i)
	{
	  loglikelihood += pdgm[i]/actualspecdens[i] / (sigma*sigma);
	  loglikelihood += log(actualspecdens[i]);
	}
      loglikelihood /= -2.0;
      loglikelihood_is_valid = true;
    }

  if (nonmissing-p-q>2)
    {
      AICC = -2*loglikelihood + 2*(p+q+1.0)*nonmissing/(nonmissing-p-q-2.0);
      if (fabs(fracdiff)>zero)
	AICC += 2.0;
    }
  else
    AICC = 0.0;   // probably should flag an error instead!

  delete[] rs;
  return loglikelihood;
}

double my_hyperg1(double t, void *parms)
{
  Vector *p = (Vector *) parms;
  double a = (*p)[0], c = (*p)[1], rho_re = (*p)[2], rho_im = (*p)[3];

  complex tc(1-t*rho_re, -rho_im);
  complex rval = (c-1)*pow(tc,-a)*pow(1-t,c-2);
  return rval.real();
}

double my_hyperg2(double t, void *parms)
{
  Vector *p = (Vector *) parms;
  double a = (*p)[0], c = (*p)[1], rho_re = (*p)[2], rho_im = (*p)[3];

  complex tc(1-t*rho_re, -rho_im);
  complex rval = (c-1)*pow(tc,-a)*pow(1-t,c-2);
  return rval.imag();
}

void my_err_handler(const char *reason, const char *file,
		     int line, int gsl_err_num)
{
  cout << "GSL Err #" << gsl_err_num << ": ";
  cout << reason << endl;
}

complex sf_hyperg_2F1(double a, double b, double c, double rho_real, double rho_imag)
{
  // this one sums the series directly, assuming b=1
  int i,k,m;
  double paonpc,magrho = sqrt(rho_real*rho_real+rho_imag*rho_imag);
  double rhostar = 1 - (1-magrho)/2.0;
  complex total,chro(rho_real,rho_imag),chroprod;
  double tolerance = 0.0001;

  // determine how many terms to use
  k = (int) (ceil((fabs(a)*magrho + fabs(c)*rhostar) / (magrho-rhostar))+0.5);
  if (k<2)
    k = 2;
  m = (int) (ceil(log(tolerance) / log(rhostar)) + 0.5);
  if (m<3) 
    m = 3;
 
  total=1.0;
  paonpc=1.0;
  chroprod=1.0;

  for (i=0 ; i<k+m ; ++i)
    {
      paonpc *= (a+i)/(c+i);
      chroprod *= chro;
      total += paonpc*chroprod;
    }
  return total;
}


complex old_sf_hyperg_2F1(double a, double b, double c, double rho_real, double rho_imag)
{
  const int workspacelen=65536;
  complex ret_val;

  gsl_error_handler_t *oldhandler = gsl_set_error_handler(&my_err_handler);

  if (rho_imag==0.0)
    {
      ret_val = complex(gsl_sf_hyperg_2F1(a,b,c,rho_real),0);
    }
  else
    {
      int err;
      double result1, result2, abserr;
      // b is assumed to be equal to one!
      Vector parms(4);
      parms[0] = a;
      parms[1] = c;
      parms[2] = rho_real;  parms[3] = rho_imag;
      gsl_function F1,F2;
      F1.function = &my_hyperg1;
      F2.function = &my_hyperg2;
      F1.params = F2.params = &parms;
      gsl_integration_workspace *ws = gsl_integration_workspace_alloc(workspacelen);
      err = gsl_integration_qags((const gsl_function *) &F1, 0.0, 1.0,
				 5e-6, 2e-3, workspacelen, ws, &result1, &abserr);
      err = gsl_integration_qags((const gsl_function *) &F2, 0.0, 1.0,
				 5e-6, 2e-3, workspacelen, ws, &result2, &abserr);
      gsl_integration_workspace_free(ws);
      ret_val = complex(result1, result2);
    }
  return ret_val;
}

Matrix ARMAModel::cfuncsfor(double d, double rho_real, double rho_imag, 
		 int p, int h, int extent)
  // returns vector C^*(d,h,rho)...C^*(d,h-extent,rho)
{
  int i,j;
  double c0,c1,c0onc1;
  complex tc, rho(rho_real, rho_imag), tc2, tc3;
  Matrix result(extent,2);
  result.zeroes();

  // Deal with numerical problems as in Doornik & Ooms
  int glen = 2*(extent-2)+1, gmid=extent-2;
  Matrix gval(glen,2);
  double a,c;

  // Step 1: compute gval vector
  //    gval(j) = G(d+i,1-d+i,rho), i=j-gmid
  c = 1-d+extent-2;
  a = d+extent-2;
  tc = sf_hyperg_2F1(a,1,c,rho_real,rho_imag);
  tc -= complex(1,0);
  tc /= rho; // fix as rho->0 !
  gval[glen-1][0] = tc.real();
  gval[glen-1][1] = tc.imag();

  // then the rest are computed recursively
  for (i=extent-3,j=2 ; i>=-extent+2 ; --i,++j)
    {
      c = 1-d+i;
      a = d+i;
      tc = complex(gval[glen-j+1][0], gval[glen-j+1][1]);
      tc = (a/c)*(1.0+rho*tc);
      gval[glen-j][0] = tc.real();
      gval[glen-j][1] = tc.imag();
    } 

  // Step 2: compute C function
  c0 = gsl_sf_gamma(1-2*d)*gsl_sf_gamma(d+h);
  c1 = gsl_sf_gamma(1-d+h)*gsl_sf_gamma(1-d)*gsl_sf_gamma(d);
  c0onc1 = c0/c1;
  tc2 = complex(gval[gmid+h][0], gval[gmid+h][1]);
  tc3 = complex(gval[gmid-h][0], gval[gmid-h][1]);
  tc = c0onc1 * (pow(rho,2*p)*tc2 + pow(rho,2*p-1) + tc3);
  result[0][0] = tc.real();
  result[0][1] = tc.imag();
  for (i=1 ; i<extent ; ++i)
    {
      int lh=h-i;
      c0onc1 *= (1-d+lh)/(d+lh);
      tc2 = complex(gval[gmid+h-i][0], gval[gmid+h-i][1]);
      tc3 = complex(gval[gmid-h+i][0], gval[gmid-h+i][1]);
      tc = c0onc1 * (pow(rho,2*p)*tc2 + pow(rho,2*p-1) + tc3);
      result[i][0] = tc.real();
      result[i][1] = tc.imag();
    }

  return result;
}

Vector ARMAModel::ComputeSpectralDensity(Vector& omegas)
{
  int i,j;
  Vector retval(omegas.size());
  complex tc1,tc2,eil,tsum,thetapart,phipart;

  for (i=0 ; i<omegas.size() ; ++i)
    {
      complex ilambda(0,-omegas[i]);
      eil = exp(ilambda);

      tc1 = tsum = 1.0;
      for (j=0 ; j<q ; ++j)
	{
	  tc1 *= eil;
	  tsum += theta[j]*tc1;
	}
      // now tsum is theta(exp(-i lambda)))
      thetapart = tsum;

      tc1 = tsum = 1.0;
      for (j=0 ; j<p ; ++j)
	{
	  tc1 *= eil;
	  tsum -= phi[j]*tc1;
	}
      // now tsum is phi(exp(-i lambda)))
      phipart = tsum;

      retval[i] = sigma*sigma/(2*M_PI)
	*(abs(thetapart)*abs(thetapart))/(abs(phipart)*abs(phipart))
	*pow(abs(1.0-eil), -2*fracdiff);
    }

  return retval;
}

