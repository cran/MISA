/*********************************************/
/* bf4assoc.c       Last Modified  03/16/10  */
/* MISA Version                              */
/* Use Laplace Approx to calculate Bayes     */
/* factors in favor of 3 genetic models of   */
/* association.                              */
/* 'Features:'                               */
/*   (1) Laplace approximation based ests.   */
/*   (2) Output log-ORs + SE(logOR)s.        */
/*   (3) Option for Cauchy prior.            */
/*********************************************/

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include<R_ext/Rdynload.h>
#include<R_ext/BLAS.h>
#include<R_ext/Linpack.h>
#include<R_ext/Lapack.h>

#include <stdlib.h>
#include <math.h>
#include <errno.h>
#define SGN(x) ((x)>0.0) ? 1.0 : (((x)<0) ? -1.0 : 0.0)
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) > (b) ? (b) : (a))

/* mem alloc */
double *vecalloc(int nr);
double **matalloc(int das,int dbs);
void   matfree(double **x,int nr, int nc);
int    *ivecalloc(int nr);
int    **imatalloc(int das,int dbs);
void   imatfree(int **x,int nr, int nc);
void   veczero(int nr, double *x);
void   iveczero(int nr, int *x);
void   matzero(int nr, int nc, double **x);
/* lapack */
/* void   F77_NAME(dpotrf)(char *ul, int *nd, double *a, int *lda, int *info);
   void   F77_NAME(dpotri)(char *ul, int *nd, double *a, int *lda, int *info); */


void bf4assoc(double *D, double *x, double *xs, double *output, int *nobs, 
	      int *ns, int *nx, double *ssd, int *prior, double *mc, 
	      int *mfe, double *reltol, double *nullstartval, int *scoring)
{
  int    NumFit=4;
  double *theta; 
  double *mu, *OR, *seOR;
  double **S, ratio=1.0;
  double logINull=0.0, logIRec=0.0, logIAdd=0.0, logIDom=0.0;
  double *thetaNull; 
  double *muNull, **mulast;
  double **SNull;
  double nc=0.0, ladet=1.0;
  double *mean, cflag=0.0;
  double *meanNull;
  int    i, j;
  void minItLS(double *x, double **S, double epsilon, int maxiter, 
	       double ssd, double *lad, int prior, double *d, double **xmat,
	       double **xs, int genmodel, int nobs, int nx, int snp,
	       double **zz, double **zztzz, double *zzty, double *AA, int alg,
	       double *lgpost, double *cf);
  double tol;       /*relative error, governs convergence of optimizer*/
  int    maxnfe;    /*maximum number of function evaluations*/
  double *var;      /*variance estimates*/
  double *varNull;  /*variance estimates*/
  double fx;        /*value of f at optimum*/
  double **LaMargBF;
  int    ndim;
  int     *numcat;
  double  **table;
  double  *bf, *labf, sumbf=0.0;
  double  **X=NULL, **XS, MinCount, snpsd;
  int    geno, Nobs, Ns, Nx, Prior;
  int    GenModel, SNP, d, optalg;
  double **z, **ztz, *zty, *A;
  double **zNull, **ztzNull, *ztyNull, *ANull;

  tol=(*reltol);
  maxnfe=(*mfe);
  Nobs=(*nobs);
  Ns=(*ns);
  Nx=(*nx);
  snpsd=(*ssd);
  Prior=(*prior);
  MinCount=(*mc);
  optalg=(*scoring);

  printf("Number of observations=%d \n",Nobs);
  printf("Number of SNPs=%d \n",Ns);
  printf("Number of confounder variables=%d \n",Nx);
  if (Prior==0) printf("Prior=Normal \n");
  if (Prior==1) printf("Prior=Cauchy \n");
  printf("Prior standard deviation on genetic log-odds parameter=%5.3f \n",snpsd);
  printf("Maximum number of function evaluations per optimization=%d \n",maxnfe);
  printf("Relative tolerance=%4.3e\n",tol);
  if (optalg==0) printf("Newton-like algorithm used for optimization.\n");
  if (optalg==1) printf("Scoring algorithm used for optimization.\n");
  for (i=0; i<(Nx+1); i++){
    printf("Null model initial value [%d]=%5.3f\n",i,nullstartval[i]);
  }

  z=matalloc(Nobs,Nx+2);
  matzero(Nobs,Nx+2,z);
  zty=vecalloc(Nx+2);
  veczero(Nx+2,zty);
  ztz=matalloc(Nx+2,Nx+2);
  matzero(Nx+2,Nx+2,ztz);
  A=vecalloc((Nx+2)*(Nx+2));
  veczero((Nx+2)*(Nx+2),A);

  zNull=matalloc(Nobs,Nx+1);
  matzero(Nobs,Nx+1,zNull);
  ztyNull=vecalloc(Nx+1);
  veczero(Nx+1,ztyNull);
  ztzNull=matalloc(Nx+1,Nx+1);
  matzero(Nx+1,Nx+1,ztzNull);
  ANull=vecalloc((Nx+1)*(Nx+1));
  veczero((Nx+1)*(Nx+1),ANull);

  numcat=ivecalloc(Ns);  /* number of unique genotype categories */
  iveczero(Ns,numcat);
  table=matalloc(Ns,8);  /* table of genotypes by SNP and Disease outcome.*/
  matzero(Ns,8,table);

  XS=matalloc(Nobs,Ns);
  matzero(Nobs,Ns,XS);
  for (i=0;i<Ns;i++){
    for (j=0;j<Nobs;j++){
      geno=((int)xs[i*Nobs + j]);
      if ((geno!=0)&&(geno!=1)&&(geno!=2)&&(geno!=3)){
	printf("ERROR: Illegal genotype value = %d\n",geno);
	exit(-1);
      }
      /* first 4 columns for controls, last four for cases */
      d=((int)D[j]);
      table[i][4*d + geno]+=1.0;
      XS[j][i]=xs[i*Nobs + j];
    }
  }
  if (Nx>0){
    X=matalloc(Nobs,Nx);
    matzero(Nobs,Nx,X);
    for (i=0;i<Nx;i++){
      for (j=0;j<Nobs;j++){
	X[j][i]=x[i*Nobs + j];
      }
    }
  }
  
  LaMargBF=matalloc(Ns,3);
  matzero(Ns,3,LaMargBF);
  
  for (i=0;i<Ns;i++){
    numcat[i]=3; /* default: all models fit */
    /* genotypes with very small counts reduce the number of genetic models fit */
    if ((table[i][0]<=MinCount)||(table[i][1]<=MinCount)||(table[i][2]<=MinCount)||
	(table[i][4]<=MinCount)||(table[i][5]<=MinCount)||(table[i][6]<=MinCount)) numcat[i]=2;
    if (((min(table[i][0],table[i][4])<=MinCount)&&(min(table[i][1],table[i][5])<=MinCount))||
	((min(table[i][0],table[i][4])<=MinCount)&&(min(table[i][2],table[i][6])<=MinCount))||
	((min(table[i][1],table[i][5])<=MinCount)&&(min(table[i][2],table[i][6])<=MinCount))) numcat[i]=1;
    /* printf("numcat[%d]=%d\n",i+1,numcat[i]); */
  }
  
  bf=vecalloc(3);
  OR=vecalloc(3);
  seOR=vecalloc(3);
  labf=vecalloc(3);
  /* Genotype model structures */
  ndim=2+Nx; /*intercept + fixed variables + SNP */
  mu=vecalloc(ndim);
  mulast=matalloc(3,ndim);
  var=vecalloc(ndim);
  S=matalloc(ndim,ndim);
  theta=vecalloc(ndim);
  mean=vecalloc(ndim);
  
  /*Null model structures */
  ndim=1+Nx; /*intercept + fixed variables */
  muNull=vecalloc(ndim);
  varNull=vecalloc(ndim);
  SNull=matalloc(ndim,ndim);
  thetaNull=vecalloc(ndim);
  meanNull=vecalloc(ndim);
  
  /***************************/
  /* Begin cycle over SNPs:  */
  for (SNP=1; SNP<=Ns;SNP++){
    if (SNP%100 == 0) printf("Processing SNP %d\n",SNP);
    
    /************************************************************/
    /*********** Genetic model = 0;  Null model *****************/
    /************************************************************/
    /* Note Null model calcs SNP-specific b/c sample sizes may  */
    /*   vary if missing observations are dropped (not imputed) */
    /************************************************************/
    GenModel=0;
    ndim=1+Nx; /*intercept+fixed variables*/
    
    veczero(3,bf);
    veczero(3,labf);
    veczero(ndim,muNull);
    veczero(ndim,varNull);
    matzero(ndim,ndim,SNull);
    veczero(ndim,thetaNull);
    veczero(ndim,meanNull);
    for (i=0; i<ndim; i++){
      muNull[i]=nullstartval[i];
    }
    cflag=0.0;
    minItLS(muNull,SNull,tol,maxnfe,snpsd,&ladet,Prior,D,X,XS,GenModel,Nobs,Nx,SNP,
	    zNull,ztzNull,ztyNull,ANull,optalg,&fx,&cflag);
    /* printf("(1) NegLL=%6.4e \n", fx); */
    logINull=(-fx)+log(ladet)+(((double)ndim)/2)*log(2.0*PI);
    
    output[(SNP-1)*18 + 14 + GenModel]=cflag;  /* null model convergence flag */
    output[(SNP-1)*18]=((double)SNP);
    output[(SNP-1)*18+1]=((double)numcat[SNP-1]);
    
    for (i=0;i<3;i++){
      OR[i]=0.0; /* initialize log-ORs as missing */
      seOR[i]=0.0; /* initialize se(log-OR)s as missing */
    }
    
    /*************************************************************/
    /* Begin cycle over Genetic Models:                          */
    /* GenModel: 0=null, 1=log additive, 2=dominant, 3=recessive */
    /*************************************************************/
    NumFit=4;
    if (numcat[SNP-1]==2) NumFit=2; /*only 2 of 3 possible genotypes observed in adequate numbers: fit additive only*/
    if (numcat[SNP-1]==1) NumFit=1; /*evidently monoallelic, or nearly so, dont fit anything!*/
    for (GenModel=1; GenModel<NumFit; GenModel++){
      ratio=1.0;
      ndim=2+Nx; /*intercept+fixed variables+SNP*/
      nc=0.0;
      
      veczero(ndim,mu);
      veczero(ndim,var);
      matzero(ndim,ndim,S);
      veczero(ndim,theta);
      veczero(ndim,mean);
      /* Starting Values */
      for (i=0; i<ndim; i++){
	mu[i]=0.0;
	/* null model fits as starting values*/
	if ((i<(ndim-1))&&(GenModel==1))  mu[i]=muNull[i];
	if ((i==(ndim-1))&&(GenModel==1)) mu[i]=0.0;
	/* Use gen effect estimate from previous gen model if */
	/* dominant or recessive:  use log-additive estimates to start: */
	if (GenModel>1) mu[i]=mulast[0][i]; 
      }
      
      cflag=0.0;
      minItLS(mu,S,tol,maxnfe,snpsd,&ladet,Prior,D,X,XS,GenModel,Nobs,Nx,SNP,
	      z,ztz,zty,A,optalg,&fx,&cflag);
      for (i=0; i<ndim; i++){
	/* save as start vals for next imputation */
	mulast[GenModel-1][i]=mu[i];
      }
      /* OR[GenModel-1]=exp(mu[ndim-1]); */ /* SNP OR */
      output[(SNP-1)*18 + 14 + GenModel]=cflag;  /* convergence flags, models of association */
      
      OR[GenModel-1]=mu[ndim-1];  /* SNP log-OR */
      seOR[GenModel-1]=sqrt(S[ndim-1][ndim-1]);  /* se of SNP log-OR */
      
      if (GenModel==1){
	logIAdd=(-fx)+log(ladet)+(((double)ndim)/2)*log(2.0*PI);
      }
      if (GenModel==2){
	logIDom=(-fx)+log(ladet)+(((double)ndim)/2)*log(2.0*PI);
      }
      if (GenModel==3){
	logIRec=(-fx)+log(ladet)+(((double)ndim)/2)*log(2.0*PI);
      }
      
      if (GenModel==1){
	labf[0]=ratio*exp(logIAdd-logINull);
      }
      if (GenModel==2){
	labf[1]=ratio*exp(logIDom-logINull);
      }
      if (GenModel==3){
	labf[2]=ratio*exp(logIRec-logINull);
      }
    }  /* end of loop over Genetic Models */
    
    /* Laplace approximation BFs: */
    LaMargBF[SNP-1][0]=labf[0];
    LaMargBF[SNP-1][1]=labf[1];
    LaMargBF[SNP-1][2]=labf[2];      
    
    if (NumFit==4){
      /* Laplace approximation BFs: */
      output[(SNP-1)*18 + 2]=LaMargBF[SNP-1][0];
      output[(SNP-1)*18 + 3]=LaMargBF[SNP-1][1];
      output[(SNP-1)*18 + 4]=LaMargBF[SNP-1][2];
      sumbf=LaMargBF[SNP-1][0]+LaMargBF[SNP-1][1]+LaMargBF[SNP-1][2];
      output[(SNP-1)*18 + 5]=LaMargBF[SNP-1][0]/sumbf;
      output[(SNP-1)*18 + 6]=LaMargBF[SNP-1][1]/sumbf;
      output[(SNP-1)*18 + 7]=LaMargBF[SNP-1][2]/sumbf;
      output[(SNP-1)*18 + 8]=OR[0];
      output[(SNP-1)*18 + 9]=OR[1];
      output[(SNP-1)*18 + 10]=OR[2];
      output[(SNP-1)*18 + 11]=seOR[0];
      output[(SNP-1)*18 + 12]=seOR[1];
      output[(SNP-1)*18 + 13]=seOR[2];
    }
    else if (NumFit==2){
      /* Laplace approximation BFs: */
      output[(SNP-1)*18 + 2]=LaMargBF[SNP-1][0];
      output[(SNP-1)*18 + 3]=0.0;
      output[(SNP-1)*18 + 4]=0.0;
      output[(SNP-1)*18 + 5]=1.0;
      output[(SNP-1)*18 + 6]=0.0;
      output[(SNP-1)*18 + 7]=0.0;
      output[(SNP-1)*18 + 8]=OR[0];
      output[(SNP-1)*18 + 9]=0.0;
      output[(SNP-1)*18 + 10]=0.0;
      output[(SNP-1)*18 + 11]=seOR[0];
      output[(SNP-1)*18 + 12]=0.0;
      output[(SNP-1)*18 + 13]=0.0;
    }
    else {
      /* Laplace approximation BFs: */
      output[(SNP-1)*18 + 2]=0.0;
      output[(SNP-1)*18 + 3]=0.0;
      output[(SNP-1)*18 + 4]=0.0;
      output[(SNP-1)*18 + 5]=0.0;
      output[(SNP-1)*18 + 6]=0.0;
      output[(SNP-1)*18 + 7]=0.0;
      output[(SNP-1)*18 + 8]=0.0;
      output[(SNP-1)*18 + 9]=0.0;
      output[(SNP-1)*18 + 10]=0.0;
      output[(SNP-1)*18 + 11]=0.0;
      output[(SNP-1)*18 + 12]=0.0;
      output[(SNP-1)*18 + 13]=0.0;
    }
  }   /* end of cycle over SNPs */
  
  /* Free(mu);
  matfree(mulast,3,ndim);
  Free(var);
  matfree(S,ndim,ndim);
  Free(theta);
  Free(mean);
  Free(muNull);
  Free(varNull);
  matfree(SNull,ndim,ndim);
  Free(thetaNull);
  Free(meanNull);
  Free(OR);
  Free(seOR);
  Free(bf);
  Free(labf);
  Free(numcat);
  matfree(table,Ns,8);
  Free(D);
  if (Nx>0) matfree(X,Nobs,Nx);
  matfree(XS,Nobs,Ns);
  matfree(LaMargBF,Ns,3);
  */
} /* Th, th, th, that's all folks! */


/*************************************************************/
/* minItLS.c      Last Modified 02/10/10                     */
/* Minimize log posterior using iteratively reweighted least */
/*  squares (either scoring or a Newton-like method).        */
/*************************************************************/
void minItLS(double *x, double **S, double epsilon, int maxiter, 
	     double sdsnp, double *ladet, int Prior, double *D,
	     double **X, double **XS, int GenModel, int Nobs, int Nx, int SNP,
	     double **z, double **ztz, double *zty, double *A, int score, 
	     double *lg, double *convflag)
{
  int    i, j, k, ndim, one=1, iter=0;
  double linpred, denom, pid, geno, det=1.0;
  /*  void    F77_NAME(dgesv)(int *n, int *nrhs, double *A, int *lda, int *ipiv, 
      double *B, int *ldb, int *info);
      void    F77_NAME(dposv)(char *uplo,int *n, int *nrhs, double *AP, int *lda,  
      double *B, int *ldb, int *info); */
  int     info=0,switched=0;
  double  delta=10e10, wi=0.0, yhati=0.0;
  double  deltal=10e10, lpost=0.0, lpostlast=10e10;
  char    uplo[]="U";
  
  if (GenModel==0) ndim=Nx+1; /* fixed vars + intercept */
  if (GenModel>0) ndim=Nx+2;  /* fixed vars + intercept + SNP*/
  while ((delta>epsilon)&&(deltal>epsilon)&&(iter<maxiter)){
    iter++;
    matzero(Nobs,ndim,z);
    veczero(ndim,zty);
    veczero(ndim*ndim,A);
    matzero(ndim,ndim,ztz);
    lpost=0.0;
    for (i=0;i<Nobs;i++) { /* cycle over samples */
      geno=XS[i][SNP-1];  /* snp number SNP, additive model */
      if (!(geno==3)){  /* ignore observation if geno=3  */
	linpred=x[0]; /* intercept */
	if (Nx>0){
	  for (j=0;j<Nx;j++)  linpred+=X[i][j]*x[j+1]; /* always-in variables */
	}
	/* observed genotypes */
	if ((GenModel==2)&&(geno==2)) geno=1;  /* dominant model */
	if ((GenModel==3)&&(geno==1)) geno=0;  /* recessive model */
	if ((GenModel==3)&&(geno==2)) geno=1;  /* recessive model */
	if (GenModel>0)  linpred+=x[Nx+1]*geno;  /* snp */
	pid=exp(linpred)/(1.0+exp(linpred));
	if (score==0){ /* if newton-like algorithm */
	  for (j=0;j<ndim;j++){
	    if      (j==0)               z[i][j]=(D[i] - pid); /* intercept */
	    else if ((j<(Nx+1))&&(Nx>0)) z[i][j]=X[i][j-1]*(D[i] - pid); /*fixed variables*/
	    else                         z[i][j]=geno*(D[i] - pid); /* geno variable */
	    zty[j]+=z[i][j];
	  }
	}  /* end of if newton-like algorithm */
	if (score==1){ /* if scoring algorithm */
	  wi=sqrt(pid*(1-pid)); /* sqrt of weights in W, McCullough & Nelder p 116 */
	  yhati=(D[i]-pid);
	  for (j=0;j<ndim;j++){
	    if      (j==0)               z[i][j]=wi;           /* intercept */
	    else if ((j<(Nx+1))&&(Nx>0)) z[i][j]=wi*X[i][j-1]; /* fixed variables*/
	    else                         z[i][j]=wi*geno;      /* geno variable */
	    zty[j]+=(z[i][j]/wi)*yhati;
	  }
	}  /* end of is scoring algorithm */
	for (j=0;j<ndim;j++){ /* common to both algorithms */
	  for (k=0;k<ndim;k++){ 
	    ztz[j][k]+=z[i][j]*z[i][k];
	  }
	}
	lpost+=(D[i]*log(pid/(1-pid))+log(1-pid));
      } /* end of not ignore observation */
    }
    /* Prior: */
    if (GenModel>0){
      if (Prior==0){ 
	/* normal prior */
	zty[ndim-1] += (-x[Nx+1]/(sdsnp*sdsnp));
	ztz[ndim-1][ndim-1] += 1.0/(sdsnp*sdsnp);
	lpost+=(-0.5*(x[Nx+1]*x[Nx+1])/(sdsnp*sdsnp) - 0.5*log(2.0*PI*sdsnp*sdsnp));
      }
      else{
	/* Cauchy prior  -- see 6/15/09 notes.        */
	denom=(1.0+((x[Nx+1]*x[Nx+1])/(sdsnp*sdsnp)));
	zty[ndim-1] += (-2.0*x[Nx+1]/(sdsnp*sdsnp))/denom;
	ztz[ndim-1][ndim-1] += (2.0/(sdsnp*sdsnp))*(1.0-((x[Nx+1]*x[Nx+1])/(sdsnp*sdsnp)))/(denom*denom);
	lpost+=(-log(denom)-log(PI*sdsnp));
      }
    }
    /* solve system: */
    for (j=0;j<ndim;j++){
      for (k=0;k<ndim;k++){ 
	/* printf("%4.3e ",ztz[j][k]);*/
	A[j*ndim + k]=ztz[j][k];
      }
      /* printf("\n"); */
    }
    /* for (j=0;j<ndim;j++){
      printf("%4.3e ",zty[j]);
    }
    printf("\n"); */
    /* on output:  A is the Cholesky square root of ztz */
    /* on output:  zty is solution to Ax=zty */
    F77_NAME(dposv)(uplo,&ndim,&one,&A[0],&ndim,&zty[0],&ndim,&info);
    det=1.0; /* Determinant of the square-root of the inverse covariance matrix.   */
    for (j=0;j<ndim;j++) det*=A[ndim*j+j];
    det=1.0/det;  /* Determinant of the square-root of the (normal approx) covariance matrix. */
    if (info!=0){
      /* printf("******************************************************\n");
	 printf("* ERROR:  Solution to linear system not found,       *\n");
	 printf("*  dposv_() exited w/ flag=%i                        *\n",info);
	 printf("******************************************************\n"); */
      printf(" NOTE: switching to Newton-like iteration\n");
      /* switch to newton-like iteration if apparently diverging */
      score=0;
      switched=1;
      x[Nx+1]=0.0; /* reset */
    }
    /* for (j=0;j<ndim;j++){
       for (k=0;k<ndim;k++){ 
       printf("%4.3e ",ztz[j][k]);
       }
       printf("\n");
       } */
    delta=0.0;
    deltal=fabs(lpost-lpostlast)/(fabs(lpostlast)+epsilon);
    for (j=0;j<ndim;j++){
      delta+=fabs(zty[j])/(fabs(x[j])+epsilon);
      x[j]=x[j]+zty[j];
    }
    /* polish if switched */
    if ((switched==1)&&(delta<epsilon)){
      printf(" NOTE: switching back to scoring for final iterations.\n");
      score=1;
      delta=10.0*epsilon;
      switched=0;
    }

    /********************************************************************************/
    /* At convergence inverse(ztz) is an estimate of the covariance matrix.         */
    /* The structure A currently contains the cholesky square root of inverse(ztz). */
    /********************************************************************************/
  } /* end of while loop */

  /* set convergence code */
  if (delta<=epsilon)       (*convflag)=1.0;
  else if (deltal<=epsilon) (*convflag)=2.0;
  else if (iter>=maxiter)   (*convflag)=3.0;
  else                      (*convflag)=-9.0;

  /* Inverse of ztz using its Cholesky square root*/
  F77_NAME(dpotri)(uplo,&ndim,&A[0],&ndim,&info);
  if (info!=0){
    printf("ERROR:  inversion routine exited with error %i\n",info);
    printf("  see man page for dpotri\n");
    exit(-1);
  }
 
  /* The upper triangle of A holds the inverse of ztz. */
  /* Note that Fortran is column major:  if the matrix is 2 by 2 its four
     elements have index in A:
          | 0  2 |
          | 1  3 |    */
  /* The following puts the covariance matrix into the matrix S */
  for (j=0;j<ndim;j++){
    for (k=0;k<ndim;k++){ 
      if (k>=j){
	S[j][k]=A[k*ndim + j];
	S[k][j]=A[k*ndim + j];
      }
    }
  }

  /* for (j=0;j<ndim;j++){
     for (k=0;k<ndim;k++){ 
     printf("%4.3e ",S[j][k]);
     }
     printf("\n");
     }
     printf("\n");
     for (j=0;j<ndim;j++) printf("%4.3e ",x[j]);
     printf("\n"); */

  (*ladet)=det;
  (*lg)=-lpost; /*minimization, not maximization*/

}  /* end of minItLS */




#if 0
/******************************************/
/* mem.c   LAST MODIFIED:  10/07/98       */
/* Memory allocation subroutines from the */
/* test Numerical Computation Using C by  */
/* R. Glassey and elsewhere...            */
/******************************************/
/* Modified to be compatible with R memory allocation by MAC */
/* double *vecalloc(int nr)
{
  double *x;
  x=Calloc( nr, double);
  return x;
}
*/
/* int *ivecalloc(int nr)
{
  int *x;
  x=Calloc( nr, int);
  return x;
}
*/

int  **imatalloc(int nr, int nc)
{
int k;
int **x;
x = (int **) R_alloc(nr, sizeof(int *));
for (k=0;k<nr;k++){
  x[k] = (int *) R_alloc(nc, sizeof(int));
}
return x;
}

unsigned char  **cmatalloc(int nr, int nc)
{
int k;
unsigned char  **x;
x = (unsigned char **) R_alloc(nr, sizeof(unsigned char *));
for (k=0;k<nr;k++){
  x[k] = (unsigned char  *) R_alloc(nc, sizeof(unsigned char));
}
return x;
}

double *vecalloc(int nr)
{
  double *x;
  x=(double *) R_alloc(nr,sizeof(double));
  return x;
}

long *lvecalloc(int nr)
{
  long *x;
  x=(long *) R_alloc((unsigned) nr,sizeof(long));
  return x;
}



float *fvecalloc(int nr)
{
  float *x;
  x=(float *) R_alloc((unsigned) nr,sizeof(float));
  return x;
}

int *ivecalloc(int nr)
{
  int *x;
  x= (int *) R_alloc( nr, sizeof(int));
  return x;
}

double **matalloc(int nr, int nc)
{
int k;
double **x;
x=(double **) R_alloc((unsigned) nr, sizeof(double *));

for (k=0;k<nr;k++){
  x[k]=(double *) R_alloc((unsigned) nc, sizeof(double));
}
return x;
}

double **matalloc_C(int nr, int nc)
{
  int k;
  double **x;
  x= Calloc( nr, double *);
  for (k=0;k<nr;k++){
    x[k]= Calloc(nc, double);
  }
  return x;
}

int **imatalloc_C(int nr, int nc)
{
  int k;
  int **x;
  x= Calloc( nr, int *);
  for (k=0;k<nr;k++){
    x[k]= Calloc(nc, int);
  }
  return x;
}


void  matfree(double **mat, int  nr, int nc)
{
int k;

for (k=0;k<nr;k++){
  Free(mat[k]);
  }
Free(mat);

}

void  imatfree(int **mat, int  nr, int nc)
{
int k;

for (k=0;k<nr;k++){
  Free(mat[k]);
  }
Free(mat);

}

#endif

void veczero(int nr, double *x)
{
  int i;

  for (i=0; i<nr; i++) x[i]=0.0;
}

void iveczero(int nr, int *x)
{
  int i;

  for (i=0; i<nr; i++) x[i]=0;
}

void matzero(int nr, int nc, double **x)
{
  int i, j;

  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++) x[i][j]=0.0;
}








