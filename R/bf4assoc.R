bf4assoc <- function(D, X=NULL, XS, Ns, Nx, snpsd, Prior, MinCount, MaxIt=100,
                     RelTol=1e-4, scoring=1){
  
  ## Performs a marginal screen on the SNPs of interest by using Laplace
  ## approximations to estimate the marginal Bayes Factors (BFs) of each SNP.
  ##   D:        response vector coded (0=control, 1=case) of length N.
  ##   X:        numeric matrix of counfounder/design variables of dimension
  ##             N by Nx; is NULL if Nx==0.
  ##   XS:       numeric matrix of SNP variables of dimension N by Ns
  ##             coded (0=common homozygote, 1=heterozygote, 2=rare
  ##             homozygote, 3=missing.
  ##   Ns:       Number of SNPs in the input data set.
  ##   Nx:       Number of design/confounder variables included in all models.
  ##   snpsd:    Standard deviation of mean zero prior on the genetic effect
  ##             parameter when Prior==0 and scale when Prior==1.
  ##   Prior:    Set to 0 to chose Normal prior and 1 to chose Cauchy prior.
  ##   MinCount: Count below which a genotype is treated as absent for a given
  ##             SNP. Effect is to reduce the number of unique genetic models
  ##             that are discernible for that SNP.

  MaxIt <- max(MaxIt,10)
  RelTol <- min(RelTol,0.001)
  
  N <- length(D)
  if (nrow(XS) != N) stop ("Genetic design matrix, XS, has different number of rows than the response variable has elements.")
  #convert design matrices to vector format
  xs<-as.numeric(XS)
  x<-0
  if (Nx > 0){
    if (nrow(X) != N) stop ("Confounder design matrix, X, has different number of rows than the response variable, D, has elements.")
    if (nrow(X) != nrow(XS)) stop("Genetic and confounder design matrices (XS and X) have a different number of rows.\n")
    x <- as.numeric(X)
  }

  if (any((D != 0) & (D != 1))) stop("Response vector D improperly coded.\n")
  if (any((XS != 0) & (XS != 1) & (XS != 2) & (XS != 3)))
    stop("SNP genotypes are improperly coded.\n")
  if ((Prior != 0) && (Prior != 1))
    stop("Prior must be 0 (Normal) or 1 (Cauchy).\n")
  if ((scoring != 0) && (scoring != 1))
    stop("Optimization algorithm must be either scoring (scoring=1) or Newton-like (scoring=0).\n")

  #starting values for the null model:
  if (Nx == 0) startval <-
    as.numeric(coefficients(glm(D~1, family=binomial(link=logit))))
  if (Nx > 0)
    startval <- as.numeric(coefficients(glm(D~X,family=binomial(link=logit))))
  if (any(is.na(startval))) stop("Null regression possibly ill-determined.")
  if (length(startval) != (Nx + 1))
    stop("Null model starting values have wrong length -- possible multicollinearity.")
  print(round(startval, 5))
  bfout<-.C("bf4assoc",
            as.double(D),
            as.double(x),
            as.double(xs),
            output=as.double(rep(0, 18*Ns)),
            as.integer(N),
            as.integer(Ns),
            as.integer(Nx),
            as.double(snpsd),
            as.integer(Prior),
            as.double(MinCount),
            as.integer(MaxIt),
            as.double(RelTol),
            as.double(startval),
            as.integer(scoring),
            PACKAGE="MISA")

  temp<-matrix(bfout$output,nrow=Ns,ncol=18,byrow=TRUE)
  colnames(temp)<-c("SNP","N.geno ","bfAtoN","bfDtoN","bfRtoN",
                    "PrAgvnAssoc","PrDgvnAssoc","PrRgvnAssoc",
                    "logOR.LogAdd","logOR.Dom","logOR.Rec",
                    "SE.lOR.LogAdd","SE.lOR.Dom","SE.lOR.Rec",
                    "ConvergeCode.Null","Convergecode.LogAdd",
                    "ConvergeCode.Dom","ConvergeCode.Rec")
  return(temp)
}
