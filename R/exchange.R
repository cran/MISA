`exchange` <-
function(pop, pop.fit, exc.a, temp) {
  
  ## Pick a row and exchange it with a neighbor, with a probability determined
  ## by the difference of the fitness values.
  ##   pop:     matrix of indicators; one row for each chain.
  ##   pop.fit: corresponding fitness values.
  ##   exc.a:   vector of the number of exchanges at each chain.
  ##   temp:    vector of N temperatures

  N <- dim(pop)[1]
  p <- dim(pop)[2]
  
  ## Randomly pick a row, and then randomly pick a neighbor (one choice at
  ## endpoints; two otherwise).
  ex.num <- sample(c(1:N),1)
  samp.ex.num <- pop[ex.num,]
  if (ex.num == 1) {
    ex.pair <- ex.num+1
  } else if (ex.num == N) {
    ex.pair <- ex.num - 1
  }
  else {
    ex.pair <- sample(c(-1,1),1) + ex.num
  }
  samp.pair <- pop[ex.pair,]
  
  re <- min(1, exp((pop.fit[ex.num] - pop.fit[ex.pair]) *
                   ((1 / temp[ex.num]) - (1 / temp[ex.pair]))))

  if (runif(1) <= re) {
    ## Swap row ex.num with row ex.pair.
    row.tmp <- pop[ex.num,]
    pop[ex.num,] <- pop[ex.pair,]
    pop[ex.pair,] <- row.tmp

    fit.tmp <- pop.fit[ex.num]
    pop.fit[ex.num]  <- pop.fit[ex.pair]
    pop.fit[ex.pair] <- fit.tmp
    
    exc.a[ex.num] <- exc.a[ex.num] + 1
  }

  return(list(pop, pop.fit, exc.a, ex.num))
}

