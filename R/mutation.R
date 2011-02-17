`mutation` <-
function(chain, pop, pop.fit, force, data, fitness, t, impute=impute, b=NULL,
         a=1, homozyg.rare, burnin=TRUE){
  

  ## Perform pointwise mutation on selected chain.
  ##         chain: Indicator of which chain to perform the mutation on.
  ##           pop: Matrix containing one row of indicators for each chain.
  ##       pop.fit: Vector of fitness values for each of the models specified
  ##                in the current population
  ##         force: Character vector specify the variables to force in the
  ##                models
  ##          data: Data frame of the same form as in 'Gene.EMC'. Columns are:
  ##                  case
  ##                  age
  ##                  log-additive indicators for all SNPs
  ##                  dominant predictors for all SNPs
  ##                  recessive predictors for all SNPs
  ##       fitness: Character string specifying the fitness function to use in
  ##                the algorithm
  ##             t: Temperature vector specifying the temperature value for
  ##                each chain in the population
  ##        impute: Number of imputed data sets
  ##             b: If the fitness function is "AIC.BB", the user must specify
  ##                the value for the beta hyper-parameter b
  ##             a: If the fitness function is "AIC.BB", the user must specify

  ## Select model, fit, and temperature for current chain.
  samp <- pop[chain,]
  samp.fit <- pop.fit[chain]
  tmpr <- t[chain]
     
  p <- length(samp)
  ## Randomly select a SNP.
  snpID <- sample(c(1:p), 1)

  snp.model <- samp[snpID] ## Current SNP model.

  ## Propose a move to a new snp genetic model or to remove the snp from
  ## the model. Codes are:
  ##   0  not in model
  ##   1  log-additive
  ##   2  dominant
  ##   3  recessive

  no.rcssv.mdl <- !homozyg.rare[snpID]

  ## If there are no homozygous rare genotypes for this SNP, just choose
  ## from models 0 and 1. Otherwise, choose from all 4 models.
  if (!homozyg.rare[snpID]){
    choose <- 0:1
  } else {
    choose <- 0:3
  }

  ## Create a 2 or 4-column matrix of indicators. Each row is a SNP; each
  ## column is a model with different values at each SNP. Each column is 
  ## identical, except at the selected SNP, which has values 0:1 or 0:3.
  samp.star <- matrix(NA, nrow=length(samp), ncol=length(choose))
  samp.star[, 1:(length(choose))] <- rep(samp, length(choose))
  samp.star[snpID,] <- choose

  ## Apply fit() to each column of samp.star; i.e. calculate a fitness for
  ## each of the 2 or 4 SNP models, excluding the current model (whose fitness
  ## we already know). Samp.fit.star is a vector of 4 fitness values.
  if (length(choose) > 2) {
    samp.fit.star <- apply(samp.star[,-(snp.model + 1)], 2, fit.EMC,
                           force=force, data=data, fitness=fitness,
                           impute=impute, b=b, a=a, homozyg.rare=homozyg.rare)
  } else {
    samp.fit.star <- fit.EMC(samp.star[,-(snp.model + 1)], force=force,
                             data=data, fitness=fitness, impute=impute, b=b,
                             a=a, homozyg.rare=homozyg.rare)
  }

  samp.fit.star <- append(samp.fit.star, samp.fit, after=snp.model)
  
  ## Normalize the fitnesses by subtracting the min fitness of the model.
  samp.fit.min  <- min(samp.fit.star)
  samp.fit.star <- samp.fit.star - samp.fit.min

  ## Remove the fitness for the current model from the vector of fitnesses,
  ## and apply Boltzmann distr to get probability for each of three new models.
  samp.prob <- exp(-samp.fit.star[-(snp.model +1)] / tmpr)
  samp.prob <- samp.prob / sum(samp.prob)
  
  ## Select from among 3 possible new models, weighted by samp.prob.new.
  new.model <- sample(choose[-(snp.model + 1)], 1, prob=samp.prob)
     
  rm.top    <- sum(exp(-(samp.fit.star[-(snp.model + 1)] / tmpr)))
  rm.bottom <- sum(exp(-(samp.fit.star[-(new.model + 1)] / tmpr)))
  rm <- rm.top/rm.bottom

  # Decide whether to accept this model.
  rm <- min(1,rm)
  mut.a <- 0
  if (runif(1) <= rm) {
    samp[snpID] <- new.model
    samp.fit <- samp.fit.star[(new.model + 1)] + samp.fit.min
    mut.a <- 1 
  }
     
  x <- list(samp, samp.fit, mut.a)
  names(x) <- c("samp", "samp.fit", "mut.a")
  return(x)
}
