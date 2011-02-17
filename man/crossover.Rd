\name{crossover}
\alias{crossover}

\title{ EMC: Crossover Step }
\description{
  This function takes the current state of the population in the Genetic
  EMC algorithm and performs the crossover step.
}

\usage{
crossover(pop, pop.fit, cross.a, force, data, fitness, t,
          impute = impute, b = NULL, a = 1, homozyg.rare)
          # cores = cores)
}

\arguments{
  \item{pop}{matrix specifying the current status of each chain (or
    model) of the population.}

  \item{pop.fit}{vector of fitness values for each of the models
    specified in the current population.}

  \item{cross.a}{current number of accepted crossover steps.}

  \item{force}{character vector specify the variables to force in the models }

  \item{data}{data frame of the same form in \code{Gene.EMC}  }

  \item{fitness}{character string specifying the fitness function to use
    in the algorithm.}

  \item{t}{temperature vector specifying the temperature value for each
    chain in the population.}

  \item{impute}{number of imputed data sets.}

  \item{b}{If the fitness function is "AIC.BB", the user must specify
    the value for the beta hyper-parameter b.}

  \item{a}{If the fitness function is "AIC.BB", the user must specify
    the value for the beta hyper-parameter a.}

  \item{homozyg.rare}{indicator vector that indicates which SNPs should
    not have a recessive parameter since the power is too weak.}

  \item{cores}{number of cores to use; i.e. the maximum number of
    processes to spawn.}
}

\details{
  In the crossover step, one of the top current models from the
  population is chosen to mate with another random model and two new
  models are formed by some composition of the two parental models.
}

\value{
  This function outputs a list of the following values:
  
  \item{pop}{current status of the population.}

  \item{pop.fit}{vector of fitness values for each of the models
    specified in the current population.}

  \item{cross.a}{current number of accepted crossover steps.}
}

\references{ 
  Geyer C (1991). \emph{Markov chain Monte Carlo maximum likelihood}.
  Computing Science and Statistics:156.

  Holland J (1975). \emph{Adaptation in Natural and Artificial Systems}.
  University of Michigan Press.

  Liang F, Wong W (2000). \emph{Evolutionary Monte Carlo:Applications to
    Cp Model Sampling and Change Point Problem}. Statistica Sinica 10:317.
}

\author{ Melanie Wilson <maw27@stat.duke.edu> }

\keyword{ methods }
