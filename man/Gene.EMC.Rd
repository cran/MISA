\name{Gene.EMC}
\alias{Gene.EMC}
\title{Bayesian Model Search Algorithm for Case-Control Genetic Association
  Studies}
\description{
  This function performs a sampling algorithm that is based on
  Evolutionary Monte Carlo (EMC) (Liang and Wong, 2000) which is
  a combination of parallel tempering (Geyer, 1991) and the genetic
  algorithm (Holland, 1975). The basic idea behind the genetic algorithm
  is taken from evolution in that individuals in a population compete
  and mate in order to produce increasingly stronger individuals.  Here,
  the individuals correspond to models of SNP data and the population
  corresponds to a set of models of interest.  The strength of each
  model is determined by the user specified fitness function and we are
  therefore interested in sampling the strongest models on basis of the
  given fitness function.
}

\usage{
Gene.EMC(data, force = NULL, fitness = "AIC.BB", b = NULL, a = 1,
         impute = 1, snp.subset = NULL, start.snps = NULL, n.iter = 10000, 
         N = 5, tmax = 1, tlone = 0, qm = 0.25, display.acc = FALSE,
         display.acc.ex = FALSE, status.out = NULL, 
         chkpt.rda = "chkpt.rda", save.iter = 10000, burnin = 1, 
         # cores = 1,
         log = FALSE, pb = FALSE)
}

\arguments{
  
  \item{data}{ 
    a data frame (or a list of data frames for multiple imputed data sets)
    where the first column corresponds to the response variable of interest
    (disease status), the next columns are the forced variables of interest,
    and the last columns are the SNPs of interest with a log-additive
    parameterization.}

  \item{force}{ 
    a character vector of variable names to be forced into the sampled
    models.  These variable names should correspond to the column headers
    for these variables in data.}
  
  \item{fitness}{ 
    a character string that specifies the fitness function to
    be used.  Options are "AIC","BIC","AIC.BB".}
  
  \item{b}{
    the Beta hyper-parameter b, when the fitness function is "AIC.BB".}
  
  \item{a}{
    the Beta hyper-parameter a, when the fitness function is "AIC.BB".}
  
  \item{impute}{
    the number of data sets on which to run the algorithm when data is a
    list of data frames for multiple imputed data sets.}
  
  \item{snp.subset}{
    a logical vector that indicates which SNPs (in the same order as the data
    frame) to run the algorithm on. This variable can be useful for searching
    through only a number of SNPs that have passed a set screen.}
  
  \item{start.snps}{
    a logical vector of the SNPs that will be in the initial model as a
    log-additive parameterization. If NULL the algorithm initializes with 5
    random SNPs with the log-additive parametrization.}
  
  \item{n.iter}{
    the number of iterations to run the algorithm.}
  
  \item{N}{
    the number of parallel chains.}
  
  \item{tmax}{
    the maximum temperature value of the chains. Default is to set tmax=1 so
    that there is a constant temperature ladder.}
  
  \item{tlone}{
    the number of chains that the user wants to run at a temperature value
    less than zero.}
  
  \item{qm}{
    the probability of the mutation update.}
  
  \item{display.acc}{
    logical flag indicating whether to output the acceptance rates of the
    mutation and crossover steps at each iteration.}
  
  \item{display.acc.ex}{
    logical flag indicating whether to output the acceptance rates of the
    exchanges between each parallel chain for each iteration.  Helps identify
    whether the chains are too far apart and the temperature scheme needs to
    be changed.}
  
  \item{status.out}{
    character string giving the pathname of the file to write the status of
    the algorithm and the acceptance rates to, instead of stdout.}
  
  \item{chkpt.rda}{
    character string giving the pathname of the checkpoint file 
    to save the output of the algorithm to.}

  \item{save.iter}{
    the number of iterations between each checkpoint. A checkpoint file is
    written every save.iter iterations.}
  
  \item{burnin}{
    integer indicating the length of the burnin.}
  
  #\item{cores}{
   # the number of cores to use; i.e. the maximum number of processes to spawn.}

  \item{log}{
    logical flag indicating whether to write warnings and errors to a 
    time-stamped log file.}

  \item{pb}{
    display X11 progress bar.}
}

\details{
  The algorithm is run for a chosen number of iterations where we update
  the population via the mutation, crossover, and exchange steps of the
  genetic algorithm. For each iteration the algorithm runs the mutation
  operator with probability qm and the crossover operator with probability
  (1-qm) and then performs the exchange step with all adjacent chains.
  The user must specify a the population size, N, the maximum and minimum
  temperature for the ladder T. (We make the assumption that the jumps
  in the temperature ladder have the form tj-ti = exp(tj/ti) where ti =
  1 for some value i in the ladder). These parameters are chosen so that
  the algorithm converges at a fast rate and can normally be determined
  in test runs of the algorithm by examining overall acceptance rates
  of the mutation and crossover rates and acceptance rates of exchanges
  between adjacent chains.  Low exchange rates between chains indicates
  that the temperature values are too far apart for adjacent chains and
  either the max temp. should be decreased or the number of parallel
  chains should be increased.  Also, convergence of the algorithm can
  be determined by running two independent runs of the algorithm with
  different starting values and examining the convergence plots produced
  in \code{\link{emc.converge}}.
}

\value{
  This function outputs a list of the following values to the file
  write.out if this file is specified for every save.iter number of
  iterations:
  
  \item{which}{
    Matrix where each row corresponds to a model specification vector for
    the unique models visited and the value of the fitness function for
    the model.}
  
  \item{data}{
    Data frame used to run the algorithm of the same form as the imputed
    data frame.}
  
  \item{iter.aic}{The value of the cost function for each of the models
    visited at each of the iterations.}
  
  \item{iter}{
    The total number of iterations run.}
  
  \item{iter.unique}{
    Vector indicating the number of the iteration in which each of the unique
    models were found.}
  
  \item{force}{
    Character vector indicating the names of the forced variables.}
  
  \item{fitness}{
    Character string indicating the fitness function used for the algorithm.}
}

\references{ 
  Geyer C (1991). \emph{Markov chain Monte Carlo maximum likelihood}. Computing
  Science and Statistics:156.

  Holland J (1975). \emph{Adaptation in Natural and Artificial Systems}.
  University of Michigan Press.

  Liang F, Wong W (2000). \emph{Evolutionary Monte Carol:Applications to Cp
    Model Sampling and Change Point Problem}. Statistica Sinica 10:317.
}

\author{ Melanie Wilson <maw27@stat.duke.edu>}

\examples{
## Load the data
data(dna.snp)
## Find the number of snps in the data set
p = (dim(dna.snp)[2] - 2)

## Set the algorithm to start with SNPS 3,7, and 21.
start.snp = rep(FALSE,p)

## Run algorithm for 100 iterations and save output to emc.out.
emc.out = Gene.EMC(data=dna.snp, force=c("age"), fitness="AIC.BB", b=p, a=1,
                   start.snps=start.snp, n.iter=100, N=5, tmax=5, tlone=1,
                   qm=.25, display.acc=TRUE, display.acc.ex=TRUE) 
}

\keyword{ models }

