
## For backward-compatibility.
`emc.converge` <- function(...) { converge.EMC(...) }


`converge.EMC` <- function(emc.out, plot.type, bandwidth=1000, a=1,
                           b=NULL, ...)
{
  ## Takes the output of two independent runs of the Genetic EMC sampling
  ## algorithm and creates plots to assess the convergence of the algorithm.
  ##   emc.out:   output from Gene.EMC().
  ##   plot.type: indicates the type of plot; options are "iter",
  ##              "gelman.rubin","bayes.factor", "snp.inc", or "all".
  ##   bandwidth: If the "bayes.factor" plot.type is chosen, user must
  ##              indicate the bandwidth of iterations over which to calculate
  ##              the global Bayes Factor over.
  ##   a:         If "bayes.factor" or "snp.inc" plot.type is chosen, and the
  ##              fitness function used in the Genetic EMC algorithm is
  ##              "AIC.BB" the user must specify the value of a for the beta
  ##              hyper-parameter.
  ##   b:         If "bayes.factor" or "snp.inc" plot.type is chosen, and the
  ##              fitness function used in the Genetic EMC algorithm is
  ##              "AIC.BB" the user must specify the value of b for the beta
  ##              hyper-parameter.

  if (!(plot.type %in%
        c("iter", "gelman.rubin", "bayes.factor", "snp.inc", "all"))) {
    warning("Unknown plot type \"", plot.type, "\"")
    return()
  }

  library(coda)
  num.chains <- dim(summary(emc.out))[1]
  data <- emc.out[[1]]$data
  if(length(dim(data)) < 2) { data <- data[[1]]}

  ## Combine the iterations from the two unique chains
  iter.aic <- NULL
  iter <- NULL
  for(i in 1:num.chains) {
    sum.i <- length(emc.out[[i]]$iter.aic)
    iter.aic.i <- cbind(rep(i,sum.i), emc.out[[i]]$iter.aic)
    iter.aic <- rbind(iter.aic, iter.aic.i)
    iter <- c(iter, emc.out[[i]]$iter)
    fitness <- emc.out[[i]]$fitness
    force <- emc.out[[i]]$force
  }
  
  if(plot.type == "all") {
    par(mfrow=c(2,2))
  }
  
  ## Make iteration plots
  if(plot.type == "iter" || plot.type == "all") {
    red.plot <- round(iter.aic[,1]/2) == iter.aic[,1] / 2
    blue.plot <- red.plot == FALSE
    aic <- iter.aic[, 2]
    aic.red.ind  <- c(1:(sum(red.plot)  / bandwidth)) * bandwidth
    aic.blue.ind <- c(1:(sum(blue.plot) / bandwidth)) * bandwidth
    total.iter <- dim(iter.aic)[1]
    ylim <- c(min(aic), max(aic))
    plot(y=iter.aic[red.plot, 2][aic.red.ind],
         x=c(1:total.iter)[red.plot][aic.red.ind], type="l",
         xlim=c(0, total.iter), col="red",
         main="Cost vs Iteration", ylab="Cost",
         xlab="Iteration", ylim=ylim)
    lines(y=iter.aic[blue.plot, 2][aic.blue.ind],
          x=c(1:total.iter)[blue.plot][aic.blue.ind], type="l", col="blue")
  }

  ## Make Gelman Rubin Plot
  if (plot.type == "gelman.rubin" || plot.type == "all") {
    aic.list <- list(NA)
    for (i in 1:max(iter.aic[,1])) {
      aic.list[[i]] <- iter.aic[iter.aic[,1] == i, 2]
    }
    iter.min <- min(as.numeric(summary(aic.list)[, 1]))
    sub.aic <- function(aic.list) {
      aic.list <- aic.list[1:iter.min]
      return(aic.list)
    }
    aic.list <- lapply(aic.list, sub.aic)
    mcmc.list <- lapply(aic.list, mcmc, start=1, end=iter.min)
    
    gelman.plot(mcmc.list, xlab="Iteration", ylab="Shrink Factor",
                auto.layout=FALSE, main="Gelman-Rubin Shrink Factor",
                ylim=c(0, 4))}
  
  ## Function to get the BF or SNP Inclusion Prob
  calc.bf <- function(emc.out, bandwidth=NULL, calc, b=b, a=a) {
    which <- emc.out$which
    cost <- which[, dim(which)[2]]
    which <- which[, -dim(which)[2]]
    snps <- colnames(which)
    snp.info <- emc.out$snp.info
    fitness <- emc.out$fitness
    force <- emc.out$force
    p <- dim(which)[2]

    ## If the null model has been sampled, find at which iteration it was
    ## sampled.
    ind.null <- NULL
    for (i in 1:dim(which)[1]) {
      num.pred <- sum(which[i,] != 0)
      if (num.pred == 0) {
        ind.null <- i
        break
      }
    }

    ## Move sampled null models to the top of "which", and remove the entries
    ## in emc.out$iter.unique that corresponds to the sampled null models.
    if (!is.null(ind.null)) {
      which <- rbind(which[ind.null,],which[-ind.null,])
      cost <- c(cost[ind.null],cost[-ind.null])
      emc.out$iter.unique = emc.out$iter.unique[-ind.null]
    }
    
    ## Compute penalties and null.priors for all fitness functions
    if (fitness=="AIC.BB") {
      null.prior <- exp(lgamma(a+b) + lgamma(b+p) - lgamma(b) - lgamma(a+b+p))
      l.nc <- lgamma(a+b) - lgamma(a) - lgamma(b) - lgamma(a+b+p)
      penalty <- 2*length(force) - 2*(lgamma(a) + lgamma(b+p)) - 2*l.nc
    }
    
    n <- dim(data)[1]

    ## If null model was not sampled, fill it in!! 
    if (is.null(ind.null)) {
      ##Compute the Null estimated cost
      which <- rbind(rep(0, length(snps)), which)
      X <- cbind(1, as.matrix(data[-1]))
      ##fit <- bayesglm.fit(x=X, y=data[[1]], family=binomial())
      fit <- .Call("bayesglm_fit", x=X, y=data[[1]], family=binomial(),
                   Roffset=NULL, Rweights=NULL, Rpriorcoef=bic.prior(n),
                   Rcontrol=glm.control(), PACKAGE="MISA")
      dev <- fit$deviance

      cost  <- c(dev + penalty, cost)
    }
    
    ## Calculate Global BF and Model posteriors
    prior.odds <- null.prior / (1 - null.prior)
    bf.inner   <- exp(-(1/2) * (cost[-1] - cost[1])) * prior.odds
    num.prob   <- exp(-(1/2) * (cost - cost[1]))
    postprob   <- num.prob / sum(num.prob)
    
    if (calc=="snp") {
      probne0 <- as.numeric((postprob * 10000) %*% (which > 0)) / 10000
      probe0 <- as.numeric((postprob * 10000) %*% (which == 0)) / 10000
      snp.bf <- (probne0 / probe0)/(a / b)
      return(log(snp.bf))
    }
    
    if (calc=="bayes") {
      bf.subset <- function(cutpoint, bf.inner, iter.unique) {
        bf.i <- bf.inner[iter.unique <= cutpoint]
        return(sum(bf.i))
      }
      cutpoints <- c(1:floor((max(emc.out$iter.unique) / bandwidth))) *
        bandwidth
      bf <- unlist(lapply(cutpoints, bf.subset, bf.inner=bf.inner,
                          iter.unique=emc.out$iter.unique))
      return(list(bf,cutpoints))
    }
  }

  ## Make Bayes Factor Plot
  if (plot.type == "bayes.factor" || plot.type =="all") {
    bf <- lapply(emc.out, calc.bf, bandwidth=bandwidth, calc="bayes", b=b, a=a)
    min.cut <- min(length(bf[[1]][[1]]), length(bf[[2]][[1]]))
    if (sum(c(bf[[1]][[1]], bf[[2]][[1]]) == 0) > 0) {
      cat("Not enough models have been sampled to calculate Global Bayes Factor \n")
    }
    
    if (sum(c(bf[[1]][[1]], bf[[2]][[1]]) == 0) == 0) { 
      main <- "Convergence of Global Bayes Factors"
      plot(bf[[1]][[2]][1:min.cut], log(bf[[1]][[1]][1:min.cut]), type="l",
           xlab="Iteration", ylab="log(Global Bayes Factor)", main=main,
           col="blue")
      lines(bf[[2]][[2]][1:min.cut], log(bf[[2]][[1]][1:min.cut]),type="l",
            col="red")
    }
  }

  ## SNP Inclusion Probability Plot
  if (plot.type=="snp.inc" || plot.type== "all") {
    snp.inc <- lapply(emc.out,calc.bf, calc="snp", b=b, a=a)
    if (sum(c(snp.inc[[1]], snp.inc[[2]]) == -Inf) > 0) {
      cat("Not enough models have been sampled to calculate Marginal Bayes Factors \n")
    }
    if (sum(c(snp.inc[[1]], snp.inc[[2]]) == -Inf) == 0) {
      min.bf <- min(c(snp.inc[[1]], snp.inc[[2]]))
      max.bf <- max(c(snp.inc[[1]], snp.inc[[2]]))
      plot(snp.inc[[1]], snp.inc[[2]], xlab="Log Marg. BF from 1st Run",
           ylab="Log Marg. BF from 2nd Run", main="Convergence of Marg. BF",
           ylim=c(min.bf, max.bf), xlim=c(min.bf, max.bf))
      lines(c(min.bf, max.bf), c(min.bf, max.bf), col="purple")
      lines(c(min.bf, max.bf), c(log(3.2), log(3.2)), col="red")
      lines(c(log(3.2), log(3.2)), c(min.bf, max.bf), col="red")
    }
  }
}

