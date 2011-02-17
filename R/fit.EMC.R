## For backward-compatibility
`fit` <- function(...) { fit.EMC(...) }

`fit.EMC` <-
function(samp, force= NULL, data, fitness="AIC", impute=impute, b=NULL, a=1,
         homozyg.rare) {
  
  ## Calculate fitness of a model.
  ##   samp:    Vector of dim p(number of snps) where for each snp
  ##              0 = snp is not in the model,
  ##              1 = log-additive snp is in the model,
  ##              2 = dom snp is in the model,
  ##              3 = rec snp is in the model.
  ##   force:   The variables that are to be forced into the model.
  ##   data:    Data frame of type .snp with first column the response 
  ##            variable, next columns the forced variables and last columns
  ##            the sampled predictors (where we first put the snp.la,
  ##            then the snp.dom, and last the snp.rec).  
  ##   fitness: the fitness function to be used in the algorithm
  ##   impute:  The number of imputed data sets in the calculations in which
  ##            case the AIC returned is the average across all AIC's for 
  ##            each data set

  null.mdl <- as.logical(sum(samp) == 0)

  p <- length(samp)
  n.force <- length(force)
  start <- n.force + 2
  
  if (impute > 1) {
    data.1 <- data[[1]]
  }
  if (impute == 1) {
    data.1 <- data
  }

  cov <- c(force)
  if (!null.mdl) {  # At least one SNP in model.
    for (k in 1:3) {
      covmdl <-
        names(data.1[, (start + (k - 1) * p) : (start + k * p - 1)]
             ) [samp == k]
      cov <- c(cov, covmdl)
    } 
  }
    
  k <- length(cov) 
  n <- dim(data.1)[1]


  # If any of the SNPs in the current model have a dominant or recessive
  # model (samp > 1) AND there are no homozygous rare genotypes for this SNP
  # (homozyg.rare == FALSE), then dev is zero.
  dev <- 0
  if (sum(homozyg.rare == 0 & samp > 1) == 0) {
    if (impute == 1) {
      X <- matrix(c(rep(1, n),
                    unlist(data.1[cov], recursive=FALSE, use.names=FALSE)),
                  nrow=n)
      fit <- .Call("bayesglm_fit", RX=X, RY=data.1[[1]], 
                   Rfamily=binomial(), Roffset=NULL,
                   Rweights=NULL, Rpriorcoef=bic.prior(n),
                   Rcontrol=glm.control(), PACKAGE="MISA");
      
      dev <- fit$deviance
    }
    else {
      BASglm <- function(data, cov) {
        ## A wrapper for bayesglm.fit() to be called via lapply().
        X <- data[cov]
        X1 <- cbind(1, as.matrix(X))
        fit1 <- .Call("bayesglm_fit", RX=X1, RY=data.1[[1]], 
                      Rfamily=binomial(), Roffset=NULL,
                      Rweights=NULL, Rpriorcoef=bic.prior(n),
                      Rcontrol=glm.control(), PACKAGE="MISA");
      }
      fit <- lapply(data, BASglm, cov=cov)
      
      `get.dev` <- function(fit) {
        return(fit$dev)
      }
      
      dev   <- unlist(lapply(fit, get.dev), recursive=FALSE, use.names=FALSE)
      dev.m <- mean(dev)
      dev   <- mean(exp(-.5 * (dev - dev.m)))
      dev   <- -2 * (log(dev)) + dev.m
    }
  }
  
  if (fitness == "AIC.BB") {
    kstar <- k - n.force
    
    gp.prior <- sum(samp > 0 & homozyg.rare) * log(1/3)
    l.prior <- lgamma(a + b) - lgamma(a) - lgamma(b) - lgamma(a + b + p) +
               gp.prior + lgamma(a + kstar) + lgamma(b + p - kstar)
    fit <- dev + k * 2 - l.prior * 2
  }
    
  if (fitness == "AIC") {
    fit <- dev + k * 2
  }
  
  if (fitness == "BIC") {
    fit <- dev + k * log(n)
  }
  
  if (fitness == "DEV") {
    fit <- dev
  }

  return(.5 * fit)
}
