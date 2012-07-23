/*********************************************************
* family.c	Last modified 02/04/11                   *
* Imported from BAS 0.93                                 *
*                                                        *
* Description: Support functions for bayesglm_fit().     *
*                                                        *
* Author: Merlise Clyde                                  *
*********************************************************/

#include <Rinternals.h>
#include <Rconfig.h>
#include <R_ext/Constants.h>
#include "sampling.h"
#include "family.h"
#include <float.h>
#include <R_ext/BLAS.h>

static const double THRESH = 30.;
static const double MTHRESH = -30.;
static const double INVEPS = 1/DOUBLE_EPS;

/**
 * Evaluate x/(1 - x). An inline function is used so that x is
 * evaluated once only.
 *
 * @param x input in the range (0, 1)
 *
 * @return x/(1 - x)
 */
static R_INLINE double x_d_omx(double x) {
    if (x < 0 || x > 1)
	error(_("Value %d out of range (0, 1)"), x);
    return x/(1 - x);
}

/**
 * Evaluate x/(1 + x). An inline function is used so that x is
 * evaluated once only. [but inlining is optional!]
 *
 * @param x input
 *
 * @return x/(1 + x)
 */
static R_INLINE double x_d_opx(double x) {return x/(1 + x);}

void logit_variance(double *mu, double *var, int n) {

  int i;

  for (i = 0; i<n; i++) {
    var[i] = mu[i]*(1.0 - mu[i]);
  }
}

void logit_link(double *rmu, double *rans, int n)
{
    int i;

    for (i = 0; i < n; i++)
	rans[i] = log(x_d_omx(rmu[i]));
}

void logit_linkinv(double *reta, double *rans, int n)
{
    int i;

    for (i = 0; i < n; i++) {
	double etai = reta[i], tmp;
	tmp = (etai < MTHRESH) ? DOUBLE_EPS :
	    ((etai > THRESH) ? INVEPS : exp(etai));
	rans[i] = x_d_opx(tmp);
    }
}

void logit_mu_eta(double *reta, double *rans, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	double etai = reta[i];
	double opexp = 1 + exp(etai);

	rans[i] = (etai > THRESH || etai < MTHRESH) ? DOUBLE_EPS :
	    exp(etai)/(opexp * opexp);
    }
}

static R_INLINE
double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

void binomial_dev_resids(double *ry, double *rmu, double *rwt, double *rans, int n)
{
  int i;
  double mui, yi;
  	for (i = 0; i < n; i++) {
	    mui = rmu[i];
	    yi = ry[i];
	    rans[i] = 2 * rwt[i] *
		(y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
	}
}

double binomial_dispersion(double *resid,  double *weights, int n, int rank) {
  return(1.0);
}

void binomial_initialize(double *Y, double *mu,  double *weights, int n) {
  int i;
  for (i = 0; i < n; i++) {
    if (weights[1] == 0) Y[i] = 0.0;
    mu[i] = (weights[i] * Y[i] + 0.5)/(weights[i] + 1.0) ;
  }
}

/* Poisson */

double poisson_dispersion(double *resid,  double *weights, int n, int rank) {
  return(1.0);
}

/* Gaussian */
double Gaussian_dispersion(double *resid,  double *weights, int n, int rank) {
  double dispersion = 0.0;
  int i, dgrf = 0; // Changed from df; df is defined in Rmath.h. gl37

  for (i = 0; i<n; i++) {
    if (weights[i] > 0) dgrf += 1;
    dispersion += weights[i]*resid[i]*resid[i];
  }
   return(dispersion/(double) (dgrf - rank));
}


/* generic functions */

double deviance(double *res, int n) {
  int i;
  double   dev = 0;

  for (i=0; i<n; i++) {
    dev += res[i];
  }
  return dev;
}


double quadform (double *bwork, double *R,  int p) {

  double Q = 0.0;
  // double  alpha=1.0, beta=0.0;
  int inc = 1;
  char uplo[] = "U", trans[]="T", diag[]="N";
  //  F77_NAME(dcopy)(&p, &b[0], &inc,  &bwork[0], &inc); 
  F77_NAME(dtrsv)(uplo, trans, diag, &p, &R[0], &p, &bwork[0], &inc);  	
  Q = F77_NAME(dnrm2)(&p, &bwork[0], &inc);
  Q *=Q;
  return(Q);
}

void chol2se(double *qr, double *se, double *R, double *covwork, int p, int n) {

  int i, j, l;
  // int info;

  for (j=0, l=0; j < p; j++) {

    for (i = 0; i <p; i++, l++) { 
      R[l] = 0;
      if (i < (j+1))   R[j*p+i] = qr[j*n + i];
    }	
  }	
  //  F77_NAME(ch2inv)(&R[0], &p, &p, &covwork[0], &info);
  Lapack_chol2inv(R, p, covwork);

for (j=0; j < p; j++) {
  se[j] = sqrt(covwork[j*p + j]);
}

 return;
}

void QR2cov(double *qr,  double *R, double *covwork, int p,  int n) {

  int i, j, l;
  // int info;

  for (j=0, l=0; j < p; j++) {

    for (i = 0; i <p; i++, l++) { 
      R[l] = 0;
      if (i < (j+1))   R[j*p+i] = qr[j*n + i];
    }	
  }	
  //  F77_NAME(ch2inv)(&R[0], &p, &p, &covwork[0], &info);
  Lapack_chol2inv(R, p, covwork);
  return;
}


void  Lapack_chol2inv(double *A, int sz, double *ans)
{
  // int inc = 1;
  int i, j;
  //	F77_NAME(dcopy)(&sz, &A[0], &inc,  &ans[0], &inc); 
	for (j = 0; j < sz; j++) {
	    for (i = 0; i <= j; i++)
		ans[i + j * sz] = A[i + j * sz];
	}

	F77_CALL(dpotri)("Upper", &sz, &ans[0], &sz, &i);
	if (i != 0) {
	    if (i > 0)
		error(_("element (%d, %d) is zero, so the inverse cannot be computed"),
		      i, i);
	    error(_("argument %d of Lapack routine %s had invalid value"),
		  -i, "dpotri");
	}

	for (j = 0; j < sz; j++) {
	    for (i = j+1; i < sz; i++)
		ans[i + j * sz] = ans[j + i * sz];
	}
}


