#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <string>
#include "matrix.h"
using namespace std;
using namespace mslib;

const int maxorder = 100;
const int maxtransforms = 50;

class TimeSeriesModel;
class ARMAModel;
class Transformation;
class matrix;

int min(int,int);
int max(int,int);
int abs(int);

class JulianTime {
  protected:
    int julianday;
    double dayfraction;

  public:
    JulianTime() {julianday=0;}
    ~JulianTime() {};
    void GetYearMonthDay(int *y, int *m, int *d);
    void SetYearMonthDay(int y, int m, int d);
    void IncreaseBy(double amount);
};

class TimeSeries {
protected:
  void ParseDescriptor(char *);
  string title, description;
  
public:
  double *data;
  bool *missing, anymissing;
  int n, nallocated;
  int desiredforecasts;  // This supersedes the realitypoint variable above.
  int ntransforms;
  Transformation *transforms[maxtransforms];
  
public:
  TimeSeries();
  ~TimeSeries();
  
  void ClearOut();
  void Append(double x, bool msng=false);
  int GetN() {return n;}
  bool IsMissing(int t) {return missing[t];}
  bool HasAnyMissing();
  double *GetData() {return data;}
  Vector GetDataVector();
  Vector GetMissingVector();
  void Clip(int n0, int n1);
  void RenderInto(ostrstream &);
  string& GetDescription() {return description;}
  string& GetTitle() {return title;}
  
  // useful stuff
  double SampleMean();
  double SampleMin();
  void ComputeSampleACF(Vector *acf, Vector *pacf,
			int normalizeacf=1);
  Matrix *GetInnovations(double *vhat, int m);
  void ComputePeriodogram();
  double *periodogram;
  TimeSeries& operator=(const TimeSeries &other);
  bool ApplyTransformation(Transformation *, ostrstream&);
  bool ReverseATransformation(ostrstream&);
  
  friend ostream& operator<<(ostream&, TimeSeries&);
  friend istream& operator>>(istream&, TimeSeries&);
};

class TimeSeriesModel {
  protected:
    double loglikelihood, AICC;

  public:
    TimeSeriesModel();
    virtual ~TimeSeriesModel();

    virtual void Specify()=0;
    virtual int FitModel(TimeSeries *ts, void (*itercallback)(),
			 ostrstream& msg)=0;
    virtual void SimulateInto(TimeSeries *, int n, bool ovw)=0;
    virtual void RenderInto(ostrstream&)=0;
    virtual void ComputeACFPACF(Vector &acf, Vector *pacf,
      int normalizeacf=1)=0;
    virtual void Forecast(TimeSeries *, int nsteps,
      double *fret, TimeSeries *resids, double *rstore=NULL)=0;
    virtual void ForecastWithMissing(TimeSeries *, int nsteps,
      double *fret, TimeSeries *resids, double *rstore=NULL)=0;
    virtual double ComputeLogLikelihood(TimeSeries *)=0;
    virtual void AutoFitModel(TimeSeries *ts)=0;
    double GetLogLikelihood() {return loglikelihood;}
    double GetAICC() {return AICC;}

    static const int SUCCESS=1, NONCONVERGENT=0, UNABLE=-1;
};

class ARMAModel : public TimeSeriesModel {
  protected:
    int p,q;
    double sigma;
    double mean;
    double kappa(int i, int j, int m, Vector &acfs);	// used for prediction
    double fracdiff;
 
  public:
    ARMAModel();
    ~ARMAModel();
    ARMAModel& operator=(const ARMAModel &other);

    double phi[maxorder], theta[maxorder];

    void InteriorFit(int p, int q, TimeSeries *ts, bool usemle=false,
		     void (*callback)() = NULL);
    void SetOrder(int arorder, int maorder)
      { p=arorder;   q=maorder; }
    void SetSigma2(double sigma2)
      { sigma = sqrt(sigma2); }
    void SetCoefficients(double *phis, double *thetas);
    void SetMean(double m)
      { mean=m; }
    void SetFracDiff(double d)
      { fracdiff=d; }

    int GetP() {return p;}
    int GetQ() {return q;}
    double GetFracDiff() {return fracdiff;}
    double GetMean() {return mean;}
    double GetSigma() {return sigma;}
    int CreateStateSpaceRepresentation(Matrix& F, Matrix& G, Matrix& Q);

    bool IsCausal(bool causalify=false);
    bool IsInvertible(bool invertibilify=false);
    bool CheckModel(bool flipit=false);

    void Specify();
    int FitModel(TimeSeries *, void (*itercallback)(), ostrstream& msg);
    void AutoFitModel(TimeSeries *ts);
    void SimulateInto(TimeSeries *, int, bool);
    void Forecast(TimeSeries *, int nsteps,
      double *fret, TimeSeries *resids, double *rstore=NULL);
    void ForecastWithMissing(TimeSeries *, int nsteps,
      double *fret, TimeSeries *resids, double *rstore=NULL);
    void RenderInto(ostrstream&);
    void ComputeACFPACF(Vector& acf, Vector* pacf, int normalizeacf);
    double ComputeReducedLogLikelihood(TimeSeries *tp, double *bestsig=NULL);
    double ComputeLogLikelihood(TimeSeries *tp);

    bool hasmissingflag;
};

Matrix cfuncsfor(double d, double rho_real, double rho_complex,
		 int p, int h, int extent);

#endif
