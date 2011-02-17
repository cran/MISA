`post.prob` <-
function(emc.out, sim.info, b=NULL, a=1){
  
  ## Calculates the global and marginal posterior probabilities and Bayes
  ## Factors that give the evidence of there being an association in the
  ## overall set of SNPs of interest, the individual genes of interest and
  ## the individual SNPs of interest.
  ##   emc.out:  Output from Gene.EMC().
  ##   sim.info: Vector of character strings giving the names of the genes for
  ##             each of the SNPs appearing in the same order that they appear
  ##             in the data frame for Gene.EMC() OR a data frame giving the
  ##             gene for each SNP, with one column named "Gene" and the
  ##             other named "SNP".
  ##   b:        If the fitness function used in the Genetic EMC algorithm is
  ##             "AIC.BB" the user must specify the value of b for the beta
  ##             hyper-parameter.
  ##   a:        If the fitness function used in the Genetic EMC algorithm is
  ##             "AIC.BB" the user must specify the value of a for the beta
  ##             hyper-parameter.

  ## Extract pertinent info from emc.out.
  which    <- emc.out$which
  fitness  <- emc.out$fitness
  force    <- emc.out$force
  n.accept <- emc.out$n.accept 
  data     <- emc.out$data
  if(length(dim(data)) < 2) {
    data <- data[[1]]
  }

  if (is.vector(sim.info)) {
    snp.info <- sim.info
  } else if (is.data.frame(sim.info)) {
    ## Select rows from sim.info that contain the SNPs of interest.
    snp.info <- as.character(sim.info$Gene[sim.info$SNP %in%
                                           colnames(emc.out$which)])
  } else {
    stop("post.prob: sim.info unrecognized data type.")
  }

  # Separate fitness column (last column) from SNPs.
  cost  <- which[,dim(which)[2]]
  which <- which[,-dim(which)[2]] # matrix of indicators; row=model, col=SNP
  
  snps <- colnames(which)
  p <- dim(which)[2]  

  ##### Find the Posterior Probability of the models sampled #####
  ## Find the null model and rearrange so that null model comes first
  ind.null <- NULL
  for(i in 1:dim(which)[1]) {
    num.pred <- sum(which[i,] != 0)
    if(num.pred == 0) {
      ind.null <- c(ind.null, i)
    }
  }

  ## Move null model to the top row.
  if (length(ind.null) == 1) {
    which    <- rbind(which[ind.null,], which[-ind.null,])
    n.accept <- c(n.accept[ind.null], n.accept[-ind.null])
    cost     <- c(cost[ind.null], cost[-ind.null])
  } else if (length(ind.null) > 1) {
    stop(paste("More than one null model (", length(ind.null),
               "in emc.out.", sep=""))
  }
  
  ## Compute penalties and null.priors for all fitness functions
  if(fitness=="AIC") {
    penalty <- 2 * length(force)
  }
  if(fitness=="BIC") {
    penalty <- (length(force) + length(snps)) * log(dim(data)[1])
  }
  if(fitness=="AIC.BB") {
    null.prior <- exp(lgamma(a+b) + lgamma(b+p) - lgamma(b) - lgamma(a+b+p))
    l.nc <- lgamma(a+b) - lgamma(a) - lgamma(b) - lgamma(a+b+p)
    penalty <- 2 * length(force) - 2 * (lgamma(a) + lgamma(b+p)) - 2 * l.nc
  }

  ## If null model was not sampled, fill it in!! 
  if(length(ind.null) == 0) {
    ## Add row of zeroes to top of which.
    which <- rbind(rep(0, length(snps)), which)    
    X <- cbind(1, data[-1])
    ##fit <- bayesglm.fit(x=X, y=data[[1]], family=binomial())
    nobs <- dim(data)[1]
    fit <- .Call("bayesglm_fit", RX=as.matrix(X), RY=data[[1]],
               Rfamily=binomial(), Roffset=NULL, Rweights=NULL,
               Rpriorcoef=bic.prior(nobs), Rcontrol=glm.control(), 
               PACKAGE="MISA");
    cost <- c(fit$deviance + penalty, cost)
    n.accept <- c(0, n.accept)
  }
  
  model.size <- apply(which > 0, 1, sum)
  
  ## Calculate Global BF and Model posteriors 
  prior.odds <- null.prior / (1 - null.prior)
  post.odds <- sum(exp(-(1/2)*(cost[-1] - cost[1])))
  bf <- post.odds * prior.odds
  num.prob <- exp(-(1/2) * (cost - cost[1]))
  postprob <- num.prob / sum(num.prob)
  
  ## Calculate Prior and posterior on model size
  model.pr <- function(i, b, p, prob="prior") {
    i.ind <- c(1:length(model.size))[model.size == i]
    
    if(prob=="prior") {
      if(fitness == "AIC.BB"){
        return(exp(lgamma(a+b) - lgamma(a) - lgamma(b) - lgamma(a+b+p) +
                   lchoose(p,i) + lgamma(a+i) + lgamma(b+p-i)))
      }
    }
    if(prob=="post") {
      return(sum(postprob[i.ind]))
    }
    if(prob=="freq") {
      return(length(i.ind))
    }
  }
  
  prior.size <- unlist(lapply(0:p, model.pr, b=b, p=p, prob="prior"))
  post.size  <- unlist(lapply(0:p, model.pr, b=b, p=p, prob="post"))
  freq.size  <- unlist(lapply(0:p, model.pr, b=b, p=p, prob="freq"))

  size <- matrix(NA, nrow=length(prior.size), ncol=3)
  colnames(size) <- c("Prior", "Post", "Freq")
  size[,1] <- prior.size
  size[,2] <- post.size
  size[,3] <- freq.size
  
  post.model <- matrix(NA, nrow=length(postprob), ncol=1)
  colnames(post.model) <- "Postprob"
  post.model[,1] <- postprob
  
  ## SNP inclusion probabilities
  snps <- colnames(which)
  probne0     <- as.numeric((postprob * 10000) %*% (which > 0)) / 10000
  probe0      <- as.numeric((postprob * 10000) %*% (which == 0)) / 10000
  probne0.la  <- (as.numeric(postprob %*% (which == 1))) / probne0
  probne0.dom <- (as.numeric(postprob %*% (which == 2))) / probne0
  probne0.rec <- (as.numeric(postprob %*% (which == 3))) / probne0
  
  prior.snp <- sum(((0:p) / p) * prior.size)
  bf.snp <- (probne0 / (probe0)) * ((1 - prior.snp) / prior.snp)
  
  ## make a matrix of inclusion probabilities
  snp.inc <- matrix(NA, nrow=length(snps), ncol=6)
  rownames(snp.inc) <- snps
  colnames(snp.inc) <- c("Gene", "Inc.Pr", "Pr.LA", "Pr.Dom", "Pr.Rec", "BF")
  snp.inc[,1] <- snp.info
  snp.inc[,2] <- round(probne0,5)
  snp.inc[,3] <- round(probne0.la, 3)
  snp.inc[,4] <- round(probne0.dom, 3)
  snp.inc[,5] <- round(probne0.rec, 3)
  snp.inc[,6] <- bf.snp
  
  ## Find gene inclusion prob
  genes <- unique(as.character(snp.info))
  which.genes <- matrix(0, nrow <- length(snps), ncol <- length(genes))
  for(j in 1:length(genes)) {
    which.genes[snp.info == genes[j], j] <- 1
  }
  which.g <- which %*% which.genes
  probne0.g <- as.numeric((postprob * 10000) %*% (which.g > 0)) / 10000
  probe0.g <- as.numeric((postprob * 10000) %*% (which.g == 0)) / 10000
  
  size.gene <- apply(which.genes, 2, sum)
  
  gene.pr <- function(size.gene, p, prior.size) {
    pr.g <- 0
    for(j in 1:size.gene) {
      i <- j:p
      inner <- sum((choose((p - size.gene), (i - j)) / choose(p, i)) *
                   prior.size[i+1])
      pr.g <- pr.g + choose(size.gene, j) * inner
    }
    return(pr.g)
  }
  
  prior.gene <- unlist(lapply(size.gene, gene.pr, p=p, prior.size=prior.size))
  j.g.prior <- prior.gene[size.gene == 2][1]
  j.g.prior <- 2 * prior.snp - j.g.prior
  bf.gene <- (probne0.g / (probe0.g)) * ((1 - prior.gene) / prior.gene)
  
  gene.results <- matrix(NA, nrow <- length(genes), ncol=2)
  rownames(gene.results) <- genes
  colnames(gene.results) <- c("Inc.Pr", "BF")
  gene.results[,1] <- probne0.g
  gene.results[,2] <- bf.gene
  
  ## Find bf and other post prob. of alternative
  bf.results <- matrix(NA, nrow <- 3, ncol=1)
  colnames(bf.results) <- c(paste("Prior=", round(prior.odds, 3)))
  rownames(bf.results) <- c("Post.Prob", "Bayes.Factor", "Prior.Odds")
  bf.results[1,] <-  sum(postprob[-1])
  bf.results[2,1] <- bf
  bf.results[3,1] <- prior.odds 
  
  post.prob <- list(post.model, snp.inc, bf.results, gene.results, size, which,
                    cost, j.g.prior,  n.accept)
  names(post.prob) <- c("Post.Model", "Post.SNP", "BF.Assoc", "Post.Gene",
                        "Size", "which", "cost", "j.g.prior", "n.accept")
  
  return(post.prob)
}

