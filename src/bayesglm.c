#include "sampling.h"
#include "family.h"

/*********************************************************
* bayesglm.c	Last modified 02/04/11                   *
* Imported from BAS 0.93                                 *
*                                                        *
* Description: A plug-in replacement for stats::glm.fit. *
* Also returns marginal likelihoods for Bayesian model   *
* comparison.                                            *
*                                                        *
* Authors: Merlise Clyde, Gary Lipton                    *
*********************************************************/
typedef struct glmfamilystruc {
  const char *family;
  const char *link;
  void (*mu_eta_fn)(double *eta, double *mu, int n);
  void (*linkfun)(double *mu, double *eta, int n);
  void (*variance)(double * mu, double *var, int n);
  void (*dev_resids)(double *y, double *mu, double *weights, double *resids, 
		     int n);
  void (*linkinv)(double *eta, double *mu, int n);
  void (*initialize)(double *Y, double *mu, double *weights, int n);
  double (*dispersion)(double *resid,  double *weights, int n, int rank);
} glmstptr;


typedef struct coefpriorstruc {
  const char *family;
  const char *class;
  double *hyper;
  double (*log_marginal_likelihood)(double dev, double regSS, int n, int p, 
				    int pgamma, double g, double *hyper);
  double (*shrinkage)(double dev,  double regSS, int n, int p, int pgamma, 
		      double g, double *hyper);
  double (*g)(double dev,  double regSS, int n, int p, int pgamma, 
	      double *hyper);
} coefdistptr;

typedef struct glmfitstruc {
   double *coef;
   double *se;
   double *mu;
   double *deviance;
   int    *rank;
   double *g;
   double *shrinkage;
   double *regSS;
   double *log_marg_lik;
} glmfitptr;

SEXP bayesglm_fit(SEXP RX, SEXP RY, SEXP Rfamily, SEXP Roffset, SEXP Rweights, 
		  SEXP Rpriorcoef, SEXP Rcontrol);


double no_shrinkage(double dev,  double regSS, int n, int p, int pgamma, 
		    double g,  double *hyper) 
{
  return 1.0;
}

double shrinkage_gprior(double dev, double regSS, int n, int p, int pgamma,  
			double g, double *hyper) 
{
  return g/(1.0 + g) ;
}

double g_EB_local(double dev, double regSS, int n, int p, int pgamma, 
		  double *hyper) 
{ 
  double g = regSS/pgamma - 1.0;
  return g > 0.0 ? g : 0.0;
}

double g_gprior(double dev, double regSS, int n, int p, int pgamma, 
		double *hyper) 
{ 
  return hyper[0];
}

double no_g(double dev,  double regSS, int n, int p, int pgamma, double *hyper) 
{ 
  return 1.0;
}

double log_marginal_likelihood_IC(double dev, double regSS, int n, int p, 
				  int pgamma, double g, double *hyper) 
{
  return -.5*(dev +  pgamma*hyper[0]);
}

double log_marginal_likelihood_gprior(double dev, double regSS, int n, int p, 
				      int pgamma, double g, double *hyper) 
{
  return -.5 * (dev  + pgamma * log(g + 1.0) + regSS / (g + 1.0));
}

/* Static structures to avoid reallocating memory with each call. */
static double *DEV     = NULL;
static double *COEF    = NULL;
static double *SE      = NULL;
static double *MU      = NULL;
static int    *RANKP   = NULL;
static double *REGSS   = NULL;
static double *G       = NULL;
static double *SHRINKP = NULL;
static double *LOGMLIK = NULL;
static int    NMODELS  = 0;
static int    NOBS     = 0;
static int    P        = 0;


/* Allocate structures for glm_fit_C(). As long as the arguments don't 
   change, allocations are done for only the first call to
   bayesglm_fit(). If arguments change, collect garbage & reallocate. */   
void allocStaticStructs(int nmodels, int nobs, int p) {
  if (nmodels > NMODELS) {
    if (DEV)     Free(DEV);
    if (RANKP)   Free(RANKP);
    if (REGSS)   Free(REGSS);
    if (G)       Free(G);
    if (SHRINKP) Free(SHRINKP);
    if (LOGMLIK) Free(LOGMLIK);
    DEV      = (double *)Calloc((size_t)(nmodels), double);
    RANKP    = (int    *)Calloc((size_t)(nmodels), int);
    REGSS    = (double *)Calloc((size_t)(nmodels), double);
    G        = (double *)Calloc((size_t)(nmodels), double);
    SHRINKP  = (double *)Calloc((size_t)(nmodels), double);
    LOGMLIK  = (double *)Calloc((size_t)(nmodels), double);
    NMODELS  = nmodels;
  }
  if (p > P) {
    if (p < 100) p = 100;
    if (COEF) Free(COEF);
    if (SE)   Free(SE);
    COEF = (double *)Calloc((size_t)(p), double);  
    SE   = (double *)Calloc((size_t)(p), double);
    P = p;
  }
  if (nobs > NOBS) {
    if (MU)      Free(MU);
    MU   = (double *)Calloc((int)nobs, double);
    NOBS = nobs;
  }
}

// Called by R_unload_BAS. 
void freeStaticStructs() {
  if (DEV)     Free(DEV);
  if (COEF)    Free(COEF);
  if (SE)      Free(SE);
  if (MU)      Free(MU);
  if (RANKP)   Free(RANKP);
  if (REGSS)   Free(REGSS);
  if (G)       Free(G);
  if (SHRINKP) Free(SHRINKP);
  if (LOGMLIK) Free(LOGMLIK);
}

void R_unload_BAS(DllInfo *info) { freeStaticStructs(); }

// Core GLM code. Called by bayesglm_fit().
glmfitptr *glm_fit_C(double *X, double *Y, int n, int p, glmstptr *glmfamily, 
		     double *offset, double *weights, const char *prior_family, 
		     const char *prior_class, double *prior_hyper, int maxit, 
		     double epsilon)
{
  int    inc = 1, nmodels = 1, it=0;
  int    i, j, l, rank = 1, conv = 0;
  double one = 1.0, devold, devnew;
  double tol = fmin(1e-07, epsilon / 1000);
  char   trans[]="N";
  coefdistptr coefprior;

  // glmresult & its components
  glmfitptr *glmresult = (glmfitptr *)Calloc(1, glmfitptr);

  double *dev  = DEV;
  double *coef = COEF;
  double *se   = SE;
  double *mu   = MU;

  int    *rankp   = RANKP;
  double *g       = G;
  double *shrinkp = SHRINKP;
  double *regSS = REGSS;
  double *log_marg_lik = LOGMLIK;

  double mu_eta[n];
  double residuals[n];
  double effects[n];
  double eta[n];
  double w[n];
  double variance[n];
  int    pivot[p];
  double qraux[p];
  double coefwork[p];
  double Xwork[n * p];
  double Ywork[n];
  double dqwork[2 * p];
  double R[p * p];
  double cov[p * p];

  if  (strcmp(prior_class, "gprior") == 0) {
    coefprior.shrinkage = shrinkage_gprior;
    coefprior.log_marginal_likelihood = log_marginal_likelihood_gprior;
    if  (strcmp(prior_family, "fixed-g-prior") == 0) 
      coefprior.g = g_gprior;
    else 
      coefprior.g = g_EB_local;
  }	

  if (strcmp(prior_class, "IC") == 0) {
    coefprior.shrinkage = no_shrinkage;
    coefprior.log_marginal_likelihood = log_marginal_likelihood_IC;
    coefprior.g = no_g;
  }

  if  (strcmp(glmfamily->family, "binomial") == 0) {
    glmfamily->dev_resids = binomial_dev_resids;
    glmfamily->dispersion = binomial_dispersion;
    glmfamily->initialize = binomial_initialize;
    if (strcmp(glmfamily->link, "logit") == 0) {
       glmfamily->linkfun = logit_link;	
       glmfamily->mu_eta_fn = logit_mu_eta;
       glmfamily->variance = logit_variance; 
       glmfamily->linkinv =  logit_linkinv;
    }	
   else  Rprintf("no other links implemented yet\n");
  }
  else  Rprintf("no other families implemented yet\n");

  for (int m = 0; m < nmodels; m++) {

    glmfamily->initialize(Y, mu, weights, n);
    glmfamily->linkfun(mu, eta, n);
    glmfamily->linkinv(eta, mu, n);
    glmfamily->dev_resids(Y, mu, weights, residuals, n);
    devold = deviance(residuals, n);
    devnew = devold;
    conv = 0.0;
    it = 0;

    while ( conv < 1 && it < maxit) {
      glmfamily->mu_eta_fn(eta, mu_eta, n);
      glmfamily->variance(mu, variance, n);

      for (i = 0, l = 0; i < n; i++) {
	w[i] = sqrt(weights[i] * mu_eta[i] * mu_eta[i] / variance[i]);
	Ywork[i] = w[i] * (eta[i] - offset[i] + (Y[i] - mu[i])/mu_eta[i]);
	residuals[i] = (Y[i] - mu[i]) / mu_eta[i];
      }

      for (j = 0, l = 0; j < p; j++) {
	pivot[j] = j + 1;
	for (i = 0; i < n; i++, l++) {
	  Xwork[l] = X[l] * w[i];
	}
      }

      rank = 1;
      for (j = 0; j < p; j++) {
	pivot[j] = j + 1;
      }

      F77_NAME(dqrls)(&Xwork[0], &n, &p, &Ywork[0], &inc, &tol,  &coefwork[0],
		      &residuals[0], &effects[0], &rank, &pivot[0], &qraux[0], 
		      &dqwork[0]);

      if (n < rank) {
	Rprintf("X has rank %ld but there are only %ld observations");
	conv = 1;
      }

      for (j=0; j < p; j++) { 
	coef[pivot[j] - 1] = coefwork[j];
      }


      F77_NAME(dcopy)(&n, &offset[0], &inc, &eta[0], &inc);
      F77_NAME(dgemv)(trans, &n, &p, &one, &X[0], &n, &coef[0], &inc, &one, 
		      &eta[0],&inc);
      
      glmfamily->linkinv(eta, mu, n);
      glmfamily->dev_resids(Y, mu, weights, residuals, n);
      devnew = deviance(residuals, n);
      glmfamily->mu_eta_fn(eta, mu_eta, n);
      glmfamily->variance(mu, variance, n);

      devnew = deviance(residuals, n);

      if (fabs(devnew - devold) / (0.1 + fabs(devnew)) < epsilon) {
	conv = 1;
      }
      else { 
	devold=devnew;
      }

      it += 1;
    }	

    dev[m] = devnew;

    if (rank == p)   chol2se(&Xwork[0], &se[0], &R[0], &cov[0], p, n);
    else	{  
      QR2cov(&Xwork[0], &R[0], &cov[0], rank, n);
      for (j = 0; j < rank; j++)  se[pivot[j] - 1] = sqrt(cov[j * rank + j]);
    }
  
    regSS[m] = quadform(coefwork, R, rank);
    g[m] = coefprior.g(dev[m], regSS[m], n, p, rank, prior_hyper);
    shrinkp[m] = 
      coefprior.shrinkage(dev[m],  regSS[m],  n,  p, rank, g[m], prior_hyper);
    log_marg_lik[m] = 
      coefprior.log_marginal_likelihood(dev[m], regSS[m], n,  p, rank, g[m], 
					 prior_hyper);
    rankp[m] = rank;
  }

  glmresult->coef = coef;
  glmresult->se = se;
  glmresult->mu = mu;
  glmresult->deviance = dev;
  glmresult->rank = rankp;
  glmresult->g = g;
  glmresult->regSS = regSS;
  glmresult->shrinkage = shrinkp;
  glmresult->log_marg_lik = log_marg_lik;

  return(glmresult);
} // end glm_fit_C


/* R interface to glm_fit_C(). A plug-in replacement for
   stats::glm.fit. Also returns marginal likelihoods for Bayesian
   model comparison. Called from converge.EMC, fit.EMC, and post.prob
   in the MISA package. Example:    
     fit <- .Call("bayesglm_fit", RX=as.matrix(X), RY=data[[1]],
                  Rfamily=binomial(), Roffset=NULL, Rweights=NULL,
                  Rpriorcoef=bic.prior(nobs), Rcontrol=glm.control(),
                  PACKAGE="BAS") */

SEXP bayesglm_fit(SEXP RX, SEXP RY, SEXP Rfamily, SEXP Roffset, SEXP Rweights, 
		  SEXP Rpriorcoef, SEXP Rcontrol)
{
  double *X = REAL(RX);
  double *Y = REAL(RY);
  int n = nrows(RX), 
      p = ncols(RX);
  glmstptr *glmfamily = (struct glmfamilystruc *)Calloc(1, glmstptr);
  glmfamily->family = CHAR(STRING_ELT(getListElement(Rfamily, "family"), 0));
  glmfamily->link   = CHAR(STRING_ELT(getListElement(Rfamily, "link"), 0));
  int maxit;
  double *weights, *offset;
  double tol, *epsilon;
  const char *prior_family, *prior_class;
  double *prior_hyper;  
  glmfitptr *fit;
  int gc_weights = FALSE, gc_offset = FALSE;
  int nProtected = 0;
  int nmodels = 1;

  allocStaticStructs(1, n, p);

  if (Rweights == R_NilValue) {
    weights = (double *)Calloc(n, double);
    for (int i = 0; i < n; i++) weights[i] = 1.0;
    gc_weights = TRUE;
  }
  else {
    weights = REAL(Rweights);
  }

  if (Roffset == R_NilValue) {
    offset = (double *)Calloc(n, double);
    for (int i = 0; i < n; i++) offset[i] = 0.0;
    gc_offset = TRUE;
  }
  else {
    offset = REAL(Roffset);
  } 	

  epsilon = REAL(getListElement(Rcontrol,"epsilon"));
  prior_family = CHAR(STRING_ELT(getListElement(Rpriorcoef, "family"), 0));
  prior_class  = CHAR(STRING_ELT(getListElement(Rpriorcoef, "class") , 0));
  prior_hyper  = REAL(getListElement(Rpriorcoef, "hyper"));	
  tol = fmin(1e-07, REAL(getListElement(Rcontrol,"epsilon"))[0]/1000);
  maxit = REAL(getListElement(Rcontrol, "maxit"))[0];

  if (!strcmp("BIC", prior_family)) *prior_hyper = (double)n;

  fit = glm_fit_C(X, Y, n, p, glmfamily, offset, weights, prior_family, 
		  prior_class, prior_hyper, maxit, *epsilon);

  if (Roffset  == R_NilValue) Free(offset);
  if (Rweights == R_NilValue) Free(weights);
  
  SEXP ANS = PROTECT(allocVector(VECSXP, 9)); ++nProtected;
  SEXP ANS_names = PROTECT(allocVector(STRSXP, 9)); ++nProtected;

  SEXP Rse   = PROTECT(allocVector(REALSXP, p)); ++nProtected;  
  double *se = REAL(Rse);
  memcpy(se, fit->se, p * sizeof(double));

  SEXP Rcoef   = PROTECT(allocVector(REALSXP, p)); ++nProtected;
  double *coef = REAL(Rcoef);
  memcpy(coef, fit->coef, p * sizeof(double));

  SEXP Rmu   = PROTECT(duplicate(RY)); ++nProtected;
  double *mu = REAL(Rmu);
  memcpy(mu, fit->mu, n * sizeof(double));

  SEXP Rdeviance   = PROTECT(allocVector(REALSXP,nmodels)); ++nProtected; 
  double *dev = REAL(Rdeviance);
  memcpy(dev, fit->deviance, nmodels * sizeof(double));

  SEXP Rrank   = PROTECT(allocVector(INTSXP,1)); ++nProtected;
  int *rank    = INTEGER(Rrank);
  memcpy(rank, fit->rank, nmodels * sizeof(int));

  SEXP Rg   = PROTECT(allocVector(REALSXP, nmodels)); ++nProtected; 
  double *g = REAL(Rg);
  memcpy(g, fit->g, nmodels * sizeof(double));

  SEXP Rshrinkage = PROTECT(allocVector(REALSXP, nmodels)); ++nProtected; 
  double *shrinkage = REAL(Rshrinkage);
  memcpy(shrinkage, fit->shrinkage, nmodels * sizeof(double));

  SEXP RregSS = PROTECT(allocVector(REALSXP,nmodels)); ++nProtected; 
  double *regSS = REAL(RregSS);
  memcpy(regSS, fit->regSS, nmodels * sizeof(double));

  SEXP Rlog_marg_lik = PROTECT(allocVector(REALSXP,nmodels)); ++nProtected; 
  double *log_marg_lik = REAL(Rlog_marg_lik);
  memcpy(log_marg_lik, fit->log_marg_lik, nmodels * sizeof(double));

  Free(fit);

  SET_VECTOR_ELT(ANS, 0, Rcoef);
  SET_VECTOR_ELT(ANS, 1, Rse);
  SET_VECTOR_ELT(ANS, 2, Rmu);
  SET_VECTOR_ELT(ANS, 3, Rdeviance);
  SET_VECTOR_ELT(ANS, 4, Rrank);
  SET_VECTOR_ELT(ANS, 5, Rg);
  SET_VECTOR_ELT(ANS, 6, Rshrinkage);
  SET_VECTOR_ELT(ANS, 7, RregSS);
  SET_VECTOR_ELT(ANS, 8, Rlog_marg_lik);
  
  SET_STRING_ELT(ANS_names, 0, mkChar("coefficients"));
  SET_STRING_ELT(ANS_names, 1, mkChar("se"));
  SET_STRING_ELT(ANS_names, 2, mkChar("mu"));	
  SET_STRING_ELT(ANS_names, 3, mkChar("deviance"));
  SET_STRING_ELT(ANS_names, 4, mkChar("rank"));
  SET_STRING_ELT(ANS_names, 5, mkChar("g"));
  SET_STRING_ELT(ANS_names, 6, mkChar("shrinkage"));
  SET_STRING_ELT(ANS_names, 7, mkChar("RegSS"));
  SET_STRING_ELT(ANS_names, 8, mkChar("logmarglik"));
 
  setAttrib(ANS, R_NamesSymbol, ANS_names);

  Free(glmfamily);
  UNPROTECT(nProtected);
  return(ANS);
   
} // end bayesglm_fit()

