combine.EMC <- function(emc.list)
{
  ## Combine results from multiple runs.
  ##   emc.list: A list of output structures from Gene.EMC().
  ## Example: combine.EMC(list(emc.out.1, emc.out.2, emc.out.3))
  ## Author: Gary Lipton
  ## Date: 07/29/10

  ## Make sure all have the same number of columns.
  if (!all(unlist(lapply(lapply(emc.list, function(abc){abc$which}), ncol)) ==
           ncol(emc.list[[1]]$which))) {
    stop("Design matrices have differing numbers of columns.")
  }

  ## Initialize
  which.full       <- c()
  cost.full        <- c()
  iter.aic.full    <- c()
  iter.unique.full <- c()
  n.accept.full    <- c()
  iter.full        <- 0

  n.out = length(emc.list)
  
  ## Create design matrix containing all models from all runs.
  for (i in 1:n.out) {
    cost.col <- dim(emc.list[[i]]$which)[2]
    which.full <- rbind(which.full, emc.list[[i]]$which[, -cost.col])
    cost.full <- c(cost.full, emc.list[[i]]$which[, cost.col])
    iter.aic.full <- c(iter.aic.full, emc.list[[i]]$iter.aic)
    iter.unique.full <- c(iter.unique.full,
                          emc.list[[i]]$iter.unique + iter.full)
    iter.full <- iter.full + emc.list[[i]]$iter
    n.accept.full <- c(n.accept.full, emc.list[[i]]$n.accept)
  }
  
  rownames(which.full) <- c(1:dim(which.full)[1])

  ## Find unique models.
  which <- unique(which.full)
  unique.ind <- as.numeric(rownames(which))
  cost <- cost.full[unique.ind]
  iter.unique <- iter.unique.full[unique.ind]
  
  which <- cbind(which, cost)
  iter <- sum(iter.full)
  
  ## Recalculate the number of times each of the unique models were accepted
  n.accept <- n.accept.full[unique.ind]
  n.accept.redun <- n.accept.full[-unique.ind]
  n.accept[match(apply(which.full[-unique.ind, ], 1, paste, collapse=""),
                 apply(which.full[ unique.ind, ], 1, paste, collapse=""))] <-
    n.accept[match(apply(which.full[-unique.ind,], 1, paste, collapse=""),
                   apply(which.full[unique.ind,], 1, paste, collapse=""))] + 
                     n.accept.redun
 
  ## Create combined emc.out
  emc.out <- list(which, emc.list[[1]]$data, iter.aic.full, iter, iter.unique, 
                  emc.list[[1]]$fitness, emc.list[[1]]$force, n.accept)
  names(emc.out) <- names(emc.list[[1]])
  emc.out
}
