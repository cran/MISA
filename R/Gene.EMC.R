`Gene.EMC` <-
function(data, force=NULL, fitness="AIC.BB", b=NULL, a=1, impute=1,
         snp.subset=NULL, start.snps=NULL, n.iter=10000, N=5, tmax=1, tlone=0,
         qm=.25, display.acc=FALSE, display.acc.ex=FALSE, status.out=NULL,
         chkpt.rda="chkpt.rda", save.iter=10000, burnin=1, log=FALSE, pb=FALSE)
{

  ## Run the Evolutionary Monte Carlo algorithm on SNP association data.
  ##   data:   a data frame (or a list of data frames for multiple imputed
  ##           data sets) where the first column corresponds to the response
  ##           variable of interest (disease status), the next columns are
  ##           the forced variables of interest, and the last columns are
  ##           the SNPs of interest with a log-additive parameterization.
  ##   force:  a character vector of variable names to be forced into the
  ##           sampled models.  These variable names should correspond to
  ##           the column headers for these variables in data.
  ##   fitness: a character string that specifies the fitness function to be
  ##           used.  Options are "AIC","BIC","AIC.BB".
  ##   b:      the Beta hyper-parameter b, when the fitness function is
  ##           "AIC.BB".
  ##   a:      the Beta hyper-parameter a, when the fitness function is
  ##           "AIC.BB".
  ##   impute: the number of data sets on which to run the algorithm when
  ##           data is a list of data frames for multiple imputed data sets.
  ##   snp.subset: a logical vector that indicates which SNPs (in the same
  ##           order as the data frame) to run the algorithm on. This
  ##           variable can be useful for searching through only a number of
  ##           SNPs that have passed a set screen.
  ##   start.snps: a logical vector of the SNPs that will be in the initial
  ##           model as a log-additive parameterization. If NULL the
  ##           algorithm initializes with 5 random SNPs with the
  ##           log-additive parametrization.
  ##   n.iter: the number of iterations to run the algorithm.
  ##   N:      the number of parallel chains.
  ##   tmax:   the maximum temperature value of the chains. Default is to
  ##           set tmax=1 so that there is a constant temperature ladder.
  ##   tlone:  the number of chains that the user wants to run at a
  ##           temperature value less than zero.
  ##   qm:     the probability of the mutation update.
  ##   display.acc: logical flag indicating whether to output the acceptance
  ##           rates of the mutation and crossover steps at each iteration.
  ##   display.acc.ex: logical flag indicating whether to output the
  ##           acceptance rates of the exchanges between each parallel chain
  ##           for each iteration.  Helps identify whether the chains are
  ##           too far apart and the temperature scheme needs to be changed.
  ##   status.out: character string giving the pathname of the file to write
  ##           the status of the algorithm and the acceptance rates to,
  ##           instead of stdout.
  ##   chkpt.rda: character string giving the pathname of the checkpoint file
  ##           to save the output of the algorithm to.
  ##   save.iter: the number of iterations between each checkpoint. A
  ##           checkpoint file is written every save.iter iterations.
  ##   burnin: integer indicating the length of the burnin.
  ##   log:    logical flag indicating whether to write warnings and errors
  ##           to a time-stamped log file.
  ##   pb:     display X11 progress bar.


  ######################## Supporting Functions ########################

  ## Return a string giving the version of a library.
  ##   libnm: name of library (string)
  `libVersion` <-
    function(libnm)
    {
      pkginfo <- installed.packages()
      nrow <- dim(pkginfo)[1]
      row <- (1:nrow)[pkginfo[, "Package"] == libnm]
      pkginfo[row, "Version"][1]
    }

  
  ## Create a time-stamped log file if log ==TRUE.
  ## Write to stderr if log == FALSE.
  `init.log` <-
    function(log=log)
      {
        fname <- format(Sys.time(), "log.%d%m%y.%H:%M:%S")
        if (log) { conn = fname }
        else     { conn = stdout() }
        
        cat(paste("MISA v", libVersion("MISA"), "\nDate: ", date(),
                  "\n", sep=""), file=conn)
        cat(paste("User: ", Sys.info()["user"], "\n", sep=""), file=conn,
            append=TRUE)
        cat(paste("Node: ", Sys.info()["nodename"], "\n", sep=""), file=conn,
            append=TRUE)
        #cat(paste("Libraries: BAS v", libVersion("BAS"), sep=""), file=conn,
        #    append=TRUE)
        cat("\n\n", file=conn, append=TRUE)
        cat(paste("Iter = ", n.iter, "\n", sep=""))
        conn
      }


  ## Print a log message.
  `logmsg` <-
  function(conn, msg)
  {
    for (i in 1:length(msg)) {
      cat(msg[i], file=conn, append=TRUE)
    }
  }


  ## Perform sanity checks on input data.
  `check.data` <-
  function(conn, data, force)
  {
    logmsg(conn, "Checking data ... ")
    emsg <- c()
    wmsg <- c()
    
    case   <- data[1]
    force  <- data[2:(n.force + 1)]
    snps   <- data[, (n.force + 2):dim(data)[2]]
    status <- 0
    
    snperr <- FALSE
    ## All SNP genotypes must be integers.
    if (!all(as.integer(as.matrix(as.vector(snps))) ==
             as.vector(as.matrix(as.vector(snps))))) {
      snperr <- TRUE
    }
    
    ## All SNP genotypes must be {0, 1, 2}.
    if (!all(as.vector(as.matrix(as.vector(snps))) >= 0) ||
        !all(as.vector(as.matrix(as.vector(snps))) <= 2)) {
      snperr <- TRUE
    }
    
    if (snperr) {
      emsg <- c(emsg,
                paste("All SNP genotypes must have integer values 0, 1, or 2"))
    }
    
    ## Case values must be two-valued.
    if (any(is.na(case))) {
      na.cases <- (1:length(case[[1]]))[is.na(case)]
      str <- ""
      for (i in 1:length(na.cases)) str <- paste(str, na.cases[i])
      wmsg <- c(wmsg, paste("Case variable \"", names(case),
                            "\" has value NA at ", str, sep=""))
    }
    
    nlevels <- length(attributes(as.factor(case[[1]]))$levels)
    if (nlevels != 2) {
      emsg <- c(emsg, paste("Case variable \"", names(case),
                            "\" has more than 2 values.", sep=""))
    }    

    ## Report conversion of factor to numeric. 
    if (is.factor(case[[1]])) {
      lvls <- unique(case[[1]])
      if (as.numeric(lvls)[2] == 2) {
        caseval <- lvls[2]
        contrval <- lvls[1]
      } else {
        caseval <- lvls[1]
        contrval <- lvls[2]
      }
      wmsg <- c(wmsg,
                paste("Converting levels into numeric case/control status: ",
                      "case=", caseval, " control=", contrval, ".\n",
                      sep=""))
    }    
    
    ## Eliminate monomorphic columns.
    snp <- 1
    snps.mono <- c()
    while (snp <= length(snps)) {
      if (all(snps[,snp] == snps[1, snp])) {
        snps.mono <- c(snps.mono, names(snps[snp]))
        snps[[snp]] <- NULL
        snp = snp - 1
      }
      snp <- snp + 1
    }
    if (length(snps.mono) > 0) {
      str = ""
      for (i in 1:length(snps.mono)) str <- paste(str, snps.mono[i])
      wmsg <- c(wmsg, paste("Removing monomorphic SNPs", str))
    }
    
    ## If the designated minor allele frequency (MAF) is greater than the
    ## designated major allele frequency, swap definitions.
    f1 <- apply(data==1, 2, sum)
    f2 <- apply(data==2, 2, sum)
    
    MAF <- (f1 + 2 * f2) / (2 * nrow(data))
    snps.recode <- c()
    if (sum(MAF > 0.5) > 0) {
      data[, MAF > 0.5] <- (2 - data[, MAF > 0.5])
      snps.recode <- c(snps.recode, colnames(data[, MAF > 0.5]))
    }
    if (length(snps.recode) > 0) {
      str <- ""
      for (i in 1:length(snps.recode)) str <- paste(str, snps.recode[i])
      wmsg <- c(wmsg, paste("Minor/major allele recode at ", str))
    }  
    
    if (is.null(wmsg) & is.null(emsg)) logmsg(conn, "OK\n")
    
    ## Print warnings.
    if (!is.null(wmsg)) {
      logmsg(conn, "\nWarnings:\n")
      for (msg in wmsg) {
        logmsg(conn, paste(indent, msg, "\n", sep=""))
      }
    }
    
    ## Print errors.
    if (!is.null(emsg)) {
      logmsg(conn, "\nErrors:\n")
      for (msg in emsg) {
        logmsg(conn, paste(indent, msg, "\n", sep=""))
      }
    }
    
    if (!is.null(emsg)) logmsg(conn, "Fatal errors in data.\n")
    return(!is.null(emsg))  
    
  } # end check.data()
  
  
  ## Check all args to Gene.EMC for type, range, and consistency.
  `check.args` <-
  function(log, data, force, fitness, b, a, impute,
             snp.subset, start.snps, n.iter, N, tmax, tlone,
             qm, display.acc, display.acc.ex, status.out,
             chkpt.rda, save.iter, burnin, pb)  
  {
    logmsg(conn, "\nChecking parameters ... ")
    emsg <- c()
    wmsg <- c()
    fit.fns <- c("AIC.BB", "AIC", "BIC", "DEV")
    
    n.snps <- dim(data)[2] - n.force - 1
    d.names = names(data)
    
    ## force
    force.err = FALSE
    if (n.force > 0) {
      if (!is.vector(force)) {
        emsg <- c(emsg, "Parameter \"force\" is not vector mode.")
        force.err <- TRUE
      }
      if (!is.character(force)) {
        emsg <- c(emsg, "Parameter \"force\" is not character type.")
        force.err <- TRUE
      }
      if (is.vector(force) && is.character(force)) {
        if (!all(force %in% names(data))) {
          emsg <- c(emsg, "Some \"force\" variables not in data.")
          force.err <- TRUE
        }
      } 
    }
    
    ## fitness
    if (! fitness %in% fit.fns) {
      emsg <- c(emsg, paste("Unknown fitness function \"", fitness, "\".",
                            sep=""));
    }
    
    ## beta hyper-parameters a & b
    if (!is.numeric(a) || (a <= 0)) {
      emsg <- c(emsg, "Beta hyper-parameter \"a\" must be a positive number.")
    }
    if (!is.numeric(b) || (b <= 0)) {
      emsg <- c(emsg, "Beta hyper-parameter \"b\" must be a positive number.")
    }
    
    if  ((as.integer(impute) != impute) || (impute < 1)) {
      emsg <- c(emsg, "Parameter \"impute\" must be a positive integer.")
    }
    
    ## snp.subset
    if (!is.null(snp.subset)) {
      if (!is.logical(snp.subset)) {
        emsg <- c(emsg, "Parameter \"snp.subset\" must be a logical vector.")
      }
      if (!force.err & length(snp.subset) !=
          (dim(data)[2] - length(force) - 1)) {
        emsg <- c(emsg, "Length of vector \"snp.subset\" must be equal to the number of SNPs.")
      }
    }
    
    ## start.snps
    if (!is.null(start.snps)) {
      if (!is.logical(start.snps)) {
        emsg <- c(emsg, "Parameter \"start.snps\" must be a logical vector.")
        start.snps <- NULL  # To skip consistency check.
      }
      if (!force.err & length(start.snps) != (dim(data)[2] - length(force) - 1)) {
        emsg <- c(emsg, "Length of vector \"start.snps\" must be equal to the number of SNPs.")
        start.snps <- NULL  # To skip consistency check.
      }
    }
    
    ## Consistency between start.snps and snp.subset
    if (!is.null(snp.subset) && !is.null(start.snps)) {
      if (any(start.snps[!snp.subset])) {
        wmsg <- c(wmsg, "Vector \"snps.subset\" excludes some SNPs in \"start.snps\".")
      }
    }
        
    ## n.iter
    if ((as.integer(n.iter) != n.iter) || n.iter < 1) {
      emsg <- c(emsg, "Parameter \"n.iter\" must be a positive integer.")
    }
    
    ## N: number of chains
    if (!is.numeric(N) || (as.integer(N) != N) || (N < 1)) {
      emsg <- c(emsg, "Parameter \"N\" must be a positive integer.")
    }
    
    ## tmax
    if (!is.numeric(tmax) || tmax <= 0) {
      emsg <- c(emsg, "Parameter \"tmax\" must be a positive number.")
    }
    
    ## tlone
    if (!is.numeric(tlone) || as.integer(tlone) != tlone || tlone < 0) {
      emsg <- c(emsg, "Parameter \"tlone\" must be a nonnegative integer.")
    }
    
    ## qm
    if (!is.numeric(qm) || qm > 1 || qm < 0) {
      emsg <- c(emsg, "Parameter \"qm\" must be within the interval [0, 1].")
    }
    
    ## display.acc
    if (!is.logical(display.acc)) {
      emsg <- c(emsg, "Parameter \"display.acc\" must be a logical constant.")
    }
    
    ## display.acc.ex
    if (!is.logical(display.acc.ex)) {
      emsg <- c(emsg,
                "Parameter \"display.acc.ex\" must be a logical constant.")
    }
    
    ## status.out
    if (!is.null(status.out)) {
      if (!is.character(status.out) || nchar(status.out) == 0) {
        emsg <- c(emsg, "Parameter \"status.out\" must be a non-null string.");
      } else {
        fexists <- file.exists(status.out)
        if (!file.create(status.out)) {
          emsg <- c(emsg, paste("Cannot open output file \"", status.out,
                                "\" for writing.", sep=""))
        } else {
          if (!fexists) unlink(status.out)
        }
      }
    }
    
    ## chkpt.rda
    if (!is.null(chkpt.rda)) {
      if (!is.character(chkpt.rda) || nchar(chkpt.rda) == 0) {
        emsg <- c(emsg, "Parameter \"chkpt.rda\" must be a non-null string.");
      } else {
        fexists <- file.exists(chkpt.rda)
        if (!file.create(chkpt.rda)) {
          emsg <- c(emsg, paste("Cannot open output file \"", chkpt.rda,
                                "\" for writing.", sep=""))
        } else {
          if (!fexists) unlink(chkpt.rda)
        }
      }
    }
    
    ## save.iter
    if ((as.integer(save.iter) != save.iter) || (save.iter < 1)) {
      emsg <- c(emsg, "Parameter \"save.iter\" must be a positive integer.")
    }
    
    ## burnin
    if ((as.integer(burnin) != burnin) || (burnin < 1)) {
      emsg <- c(emsg, "Parameter \"burnin\" must be a positive integer.")
    }
    
    ## cores
    #if ((as.integer(cores) != cores) || (cores < 1)) {
    #  emsg <- c(emsg, "Parameter \"cores\" must be a positive integer.")
    #}
    #if ((cores > 1) & pb) {
      #wmsg <- c(wmsg, paste("The progress bar is not permitted with multicore processing.\n", indent, ">> Turning off progress bar.\n", sep=""))
      #pb = FALSE
    #}
    
    if (is.null(wmsg) & is.null(emsg)) logmsg(conn, "OK\n")
    
    ## Print warnings.
    if (!is.null(wmsg)) {
      logmsg(conn, "\nWarnings:\n")
      for (msg in wmsg) {
        logmsg(conn, paste(indent, msg, sep=""))
      }
    }
    
    ## Print errors.
    if (!is.null(emsg)) {
      logmsg(conn, "\nErrors:\n")
      for (msg in emsg) {
        logmsg(conn, paste(indent, msg, "\n", sep=""))
      }
      logmsg(conn, "Fatal errors in parameters.\n\n")
      return(TRUE)
    }
    return(FALSE)
    
  }  # end check.args()


  ## Initialize progress bar.
  `pb.init` <-
   function() {
    if (require(tcltk) & capabilities()["X11"] & capabilities()["tcltk"]) {
      pbar <- tkProgressBar(title="MISA", min=0, max=100, width=300)
      setTkProgressBar(pbar, 0, label=paste("0%  Cost=---------",
                                            sep=""))
      pbscale <- 100 / n.iter
      logmsg(conn, "\n")
    } else {
      return(NULL)
    }
    pblist <-  list(on = pb, widget=pbar, last.stat=0, scale=pbscale,
                    start.time=proc.time()[3], start.iter=0)
    return(pblist)
  }

  ## Update progress bar.
  `pb.update` <-
  function(pb, iter, popfit) {
    if (pb$on) {
      pbstat <- floor(iter * pb$scale)
      if (pbstat > pb$last.stat) {
        if (capabilities()["X11"] & capabilities()["tcltk"]) {
          setTkProgressBar(pb$widget, pbstat,
                           label=paste(pbstat, "%   Cost=",
                                       format(popfit, digits=6), "  ", sep=""))
        } else {
          pb$on = FALSE
        }
        
        pb$last.stat = pbstat
      }
    }
    pb
  }

  ## Extract columns selected by logical vector snp.subset.
  `sub.snp` <- function(data, snp.subset){
    snp.subset <- c(rep(TRUE, (n.force + 1)), snp.subset)
    data <- data[, snp.subset]
    return(data)
  }
  
  `expand.data` <- function(data, force) {
    p <- dim(data)[2] - length(force) - 1
    ## Translate log-additive parametrization (i.e. genotype) to possible
    ## log-additive, dominant, and recessive parametrizations.
    data <- expand.data.snp(data.snp=data[, (n.force + 2):dim(data)[2]],
                            ind.info=data[, 1:(n.force + 1)], force=force)
    return(data)
  }

  `transcribe.case` <- function(data) {
     if (is.factor(data[[1]])) {
       data[[1]] <- as.numeric(data[[1]]) - 1
     }
     if (is.logical(data[[1]])) {
       data[[1]] <- as.numeric(data[[1]])
     }
     return(data)
   }

  ## For each SNP, assign a Boolean value designating whether it can
  ## have a recessive model; i.e. whether any of the observations have a
  ## copy number of 2.
  `has.homozyg.rare` <- function(data, p) {
    homozyg.rare <- rep(FALSE, p)
    for (i in 1:p) {
      homozyg.rare[i] <- (sum(data[, i + 2 * p + n.force + 1]) != 0)
    }
    return(homozyg.rare)
  }

  `report.accept` <-
   function(m, c, e) {
     if(display.acc == TRUE) {
       cat(paste("  Accept rates of mutation and crossover are:",
                 format(m, nsmall=3, digits=3), format(c, nsmall=3, digits=3),
                 "\n"), file=status.out, append=TRUE)
     }    
     if(display.acc.ex == TRUE) {
       cat(paste("  Accept rates exchange are:", format(e, digits=4), "\n"),
           file=status.out, append=TRUE)
     }
  }
 
  ###################### End Supporting Functions ######################

  
  ######################################################################
  ########################### Main Function ############################
  ######################################################################

  #require(BAS)
  #if (cores > 1) {
   # require(multicore)
  #}
  indent <- "    "

  if (impute == 1) { data.1 <- data }
  else             { data.1 <- data[[1]] }

  n.force <- length(force)
  conn <- init.log(log)
  logmsg(conn, paste(format(match.call()), "\n", sep=""))
  status1 <- 
    check.args(conn, data=data.1, force=force, fitness=fitness, b=b, a=a,
               impute=impute, snp.subset=snp.subset, start.snps=start.snps,
               n.iter=n.iter, N=N, tmax=tmax, tlone=tlone, qm=qm,
               display.acc=display.acc, display.acc.ex=display.acc.ex,
               status.out=status.out, chkpt.rda=chkpt.rda,
               save.iter=save.iter, burnin=burnin, pb=pb)
  status2 <- check.data(conn, data.1, force) 
  if (status1 | status2) {
    logmsg(conn, "\nStopping due to fatal errors.\n")
    stop("Fatal errors.")
  }
	
  if (impute == 1) {
    data <- transcribe.case(data)
  }
  else {
    data <- lapply(data, transcribe.case)
  }

  emc.out.names <- c("which", "data", "iter.aic", "iter", "iter.unique",
                     "fitness", "force", "n.accept")

  checkpt <- !is.null(chkpt.rda) && (save.iter > 0)

  ## Set output connection.
  if (length(status.out) == 0) {
    status.out <- stdout()
  }

  if (length(snp.subset) > 0) {
    start.snps <- start.snps[snp.subset]
    if (impute == 1) {
      data <- sub.snp(data, snp.subset)
    } else {
      data <- lapply(data, sub.snp, snp.subset=snp.subset)
    }
  }
  
  gc()

  if(impute == 1) {
    data <- expand.data(data=data, force=force)
    data.expand <- data
  }
  
  if(impute > 1) {
    data <- lapply(data, expand.data, force=force)
    data.expand <- data[[1]]
  }
  
  start <- 1 + length(force)
  p <- (dim(data.expand)[2] - 1 - length(force))/3
  cov <- names(data.expand[, (start + 1):(p+start)])
  snps <- unlist(strsplit(cov, split=".", fixed=TRUE),
                 use.names=FALSE)[c(1:p) * 2 - 1]
  nobs <- dim(data.expand)[1]

  ## Start the algorithm with the list of start snps in the model with la
  ## genetic effect.
  samp <- rep(0, length(snps))     
  if(length(start.snps) == 0) {
    start.snps <- sample(c(1:length(snps)), min(5, length(snps)))
  }     
  samp[start.snps] <- 1 

  pop <- matrix(NA, nrow=N, ncol=p)
  for(i in 1:N) pop[i, ] <- samp
       
  if(impute == 1) {
    homozyg.rare <- has.homozyg.rare(data.expand, p)
  } else {
    homozyg.rare.i <- lapply(data.expand, has.homozyg.rare, p)
    homozyg.rare <- rep(1, p)
    for(i in 1:impute) {
      homozyg.rare <- homozyg.rare * homozyg.rare.i[[i]]
    }
  }
  
  ## Initialize the temperature based on N, tmax, and tlone.
  temp <- exp(log(tmax)*(c((0-tlone):(N-1-tlone))/(N-1-tlone)))
  
  ## Initialize the values that will calculate acceptance rates.
  mut <- 0
  mut.a <- 0
  cross <- 0
  cross.a <- 0
  exc <- rep(0, N)
  exc.a <- rep(0, N)
  
  ## Initialize Matrix of the unique visited models, their fitness value, and 
  ## the first generations on which they were visited.
  pop.emc.nrow <- N * n.iter
  pop.emc.ncol <- p + 3
  pop.emc.matrix <- matrix(NA, nrow=pop.emc.nrow, ncol=pop.emc.ncol)
  pop.emc.rowIDs <- rep("", pop.emc.nrow)
  colnames(pop.emc.matrix) <- c(snps, "Fitness", "Iter", "n.accept")

  pop.emc.rptr <- 1  ## current row in pop.emc.matrix

  samp.fit <- fit.EMC(samp=samp, force=force, data=data.expand, fitness=fitness,
                  impute=impute, b=b, a=a, homozyg.rare=homozyg.rare)

  pop.fit <- rep(samp.fit, N)
  pop.emc.matrix[1, ] <- c(pop[(tlone+1), ], 2 * pop.fit[(tlone+1)], 0, 1)

  ## This gives a result identical to paste(), but 37x faster.
  samp.str <- rawToChar(as.raw(pop[(tlone + 1),] + 48))

  pop.emc.rowIDs[1] <- "INIT"
  pop.emc.rptr <- 2  ## current row in pop.emc.matrix
  cost.iter <- NULL

  last.pbstat <- 0
  if (pb) { pblist <- pb.init() }

  
  ## Begin iterations
  for (iter in 1:n.iter) {
    popfit2 <- 2 * pop.fit[tlone + 1]
    cat("Iter ", format(iter, width=5),": Cost = ", format(popfit2, nsmall=3),
        "\n", file=status.out, append=TRUE, sep="")
    if (pb) {
      pblist <- pb.update(pblist, iter, popfit2)
      pb <- pblist$on
    }
    
    ## STEP 1) apply either MUTATION or CROSSOVER operator with prob qm and
    ## (1 - qm) respectively    
    flip <- runif(1)
    if (flip <= qm) {
      ## Run mutation step on all chains of the population.
      mut <- mut + N

      x <- lapply(c(1:N), mutation, pop=pop, pop.fit=pop.fit, force=force, 
                  data=data.expand, fitness=fitness, t=temp, impute=impute,
                  b=b, a=a, homozyg.rare=homozyg.rare)
      
      xv <- unlist(x, use.names=FALSE)
      pop <- matrix(xv[-c(c(1:N) * (p + 2) - 1, c(1:N) * (p + 2))],
                    nrow=N, ncol=p, byrow=TRUE)
      pop.fit <- as.numeric(xv[c(1:N) * (p + 2) - 1])
      mut.a <- sum(mut.a, as.numeric(xv[c(1:N) * (p + 2)]))
    } else {
      ## Run crossover step
      cross <- cross + 1 
      x <- crossover(pop, pop.fit, cross.a, force=force, data=data.expand,
                    fitness=fitness, t=temp, impute=impute, b=b, a=a,
                    homozyg.rare=homozyg.rare)
      pop     <- x[[1]]
      pop.fit <- x[[2]]
      cross.a <- x[[3]]
    }
    
    ## STEP 2) EXCHANGE xi with xj for N pairs (i, j) where i is samp
    ## uniformly 
    for(i in 1:N) {
      x <- exchange(pop, pop.fit, exc.a, t=temp)
      ex.num <- x[[4]]
      exc[ex.num] <- exc[ex.num] + 1
      pop <- x[[1]]
      pop.fit <- x[[2]]
      exc.a <- x[[3]]
    }
    
    ## Calculate ACCEPTANCE rates    
    if (mut != 0)   { m <- mut.a/mut }     else { m <- 0 }
    if (cross != 0) { c <- cross.a/cross } else { c <- 0 }
    
    e <- rep(0, N)
    for(k in 1:N) {
      if (exc[k] != 0) { e[k] <- exc.a[k]/exc[k] } else { e[k] <- 0 }
    }
    
    report.accept(m, c, e)
  
    ## If the current model is new, or we are at burn-in, add the model
    ## to pop.emc.matrix.
    if (iter >= burnin) {
      samp.model <- pop[(tlone + 1), ]
      ## This gives a result identical to paste(), but 37x faster.
      samp.str <- rawToChar(as.raw(samp.model + 48))

      ## Find rows (if any) that match samp.str.
      row.match <- (pop.emc.rowIDs == samp.str)

      ## If no rows match, add this model to the table.
      if (!any(row.match)) {
        if (pop.emc.rptr <= pop.emc.nrow) {
          new.row <- c(samp.model, 2 * pop.fit[(tlone + 1)], iter, 1)
          pop.emc.matrix[pop.emc.rptr,] <- new.row
          pop.emc.rowIDs[pop.emc.rptr] <- samp.str
          pop.emc.rptr <- pop.emc.rptr + 1
        } else {
          warning("Model results table full. Cannot add results.")
        }
      }

      ## If there was a match, increment the match counter for this row.
      if (sum(as.integer(row.match)) > 0) {
        pop.emc.matrix[row.match, p + 3] <-
          pop.emc.matrix[row.match, p + 3] + 1
      }
    }
   
    cost.iter <- c(cost.iter, 2*pop.fit[(tlone+1)])
       
    ## Save emc.out results after every save.iter iterations.
    if(checkpt && (iter %% save.iter == 0)) {
      pop.emc.global <- pop.emc.matrix
      emc.out <- list(pop.emc.global[-1, 1:(p + 1)], data.expand, cost.iter, iter,
                     pop.emc.global[-1, p + 2], fitness, force,
                     pop.emc.global[-1, p + 3])
      names(emc.out) <- emc.out.names
      save(emc.out, file=chkpt.rda)
      gc()
    }
  } ## end for (iter in 1:n.iter)

  pop.emc.matrix <- pop.emc.matrix[pop.emc.rowIDs != "",]
  emc.out <- list(pop.emc.matrix[-1, 1:(p + 1)], data.expand, cost.iter, iter,
                  pop.emc.matrix[-1, p + 2], fitness, force,
                  pop.emc.matrix[-1, p + 3])

  names(emc.out) <- emc.out.names
  if (pb) { close(pblist$widget) }
  
  emc.out
} # end Gene.EMC()
