`crossover` <-
function(pop, pop.fit, cross.a, force, data, fitness, t, impute=impute, b=NULL,
         a=1, homozyg.rare) {
#, cores=cores) {

  ## Performs uniform crossover on two randomly selected individuals.
  ##   pop:     Population of models.
  ##   pop.fit: Fitness value for each model.
  ##   cross.a: Number of accepted crossovers.

  N <- dim(pop)[1]
  p <- dim(pop)[2]

  ## selection of (x_k, x_j)
  prob <- pop.fit
  prob.exp <- exp(-(prob - mean(prob)))
  prob.exp <- prob.exp / sum(prob.exp)  # vector of acceptance probs
  ## Select a chromosome weighted by the probabilities.
  k <- sample(c(1:N), 1, prob=prob.exp)
  ## Select any other chromosome.
  j <- sample(c(1:N)[-k], 1)

  fit.k <- pop.fit[k]
  fit.j <- pop.fit[j]
  
  ## Make new population with uniform crossover.
  ## Randomly select a subset of the SNPs. For the SNPs not in the subset,
  ## swap the indicators between j and k. These will comprise our new
  ## chromosomes: j.star and k.star.
  ind <- sample(c(0, 1), p, replace=TRUE)
  k.star <- rep(0, p)
  j.star <- rep(0, p)
  k.star[ind==1] <- pop[k, ind==1]
  k.star[ind==0] <- pop[j, ind==0]
  j.star[ind==1] <- pop[j, ind==1]
  j.star[ind==0] <- pop[k, ind==0]

  #if (cores == 1) {
    fit.star <- lapply(list(k.star, j.star), fit.EMC, force=force, data=data,
                       fitness=fitness, impute=impute, b=b, a=a,
                       homozyg.rare=homozyg.rare)
  #} else {
    #fit.star <- multicore::mclapply(list(k.star, j.star), fit.EMC, force=force,
    #                                data=data, fitness=fitness, impute=impute,
    #                                b=b, a=a, homozyg.rare=homozyg.rare,
    #                                mc.cores=cores)
  #}
  
  fit.k.star <- fit.star[[1]]
  fit.j.star <- fit.star[[2]]

  pop.star  <- pop
  prob.star <- pop.fit
  
  ## If (fit.k, fit.j) and (fit.k.star, fit.j.star) have the same order,
  ## replace old k,j rows with new k,j rows; otherwise, swap & replace.
  if ((fit.k <= fit.j) == (fit.k.star <= fit.j.star)) {
    pop.star[k,] <- k.star
    pop.star[j,] <- j.star
    prob.star[k] <- fit.k.star
    prob.star[j] <- fit.j.star
  } else {
    pop.star[k,] <- j.star
    pop.star[j,] <- k.star
    prob.star[k] <- fit.j.star
    prob.star[j] <- fit.k.star
  }
  
  prob.exp.star <- exp(-(prob.star - mean(prob.star)))
  prob.exp.star <- prob.exp.star / sum(prob.exp.star)
  
  ## selection probabilites
  select.x <- (prob.exp[k] + prob.exp[j]) / (N-1)
  select.y <- (prob.exp.star[k] + prob.exp.star[j]) / (N-1)
  
  rc <- min(1, exp(-((prob.star[k] - prob[k]) / t[k]) -
                    ((prob.star[j] - prob[j]) / t[j])) * (select.y / select.x))
  
  if(runif(1) <= rc) {
    pop <- pop.star
    pop.fit <- prob.star
    cross.a <- cross.a + 1
  }

  return(list(pop, pop.fit, cross.a, k, j))
}

