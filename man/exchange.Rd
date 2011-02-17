\name{exchange}
\alias{exchange}

\title{ EMC: Exchange Step }
\description{
  This function takes the current state of the population in the Genetic
  EMC algorithm and performs the exchange step between each neighboring chain.
}

\usage{
exchange(pop, pop.fit, exc.a, temp)
}

\arguments{
  \item{pop}{matrix specifying the current status of each chain (or
    model) of the population.}
  \item{pop.fit}{vector of fitness values for each of the models
    specified in the current population.}
  \item{exc.a}{current number of accepted crossover steps.}
  \item{temp}{temperature vector specifying the temperature value for
    each chain in the population.}
}

\details{
  This step corresponds to the normal parallel tempering exchange step,
  where we allow the models, or current states of the chains to move up
  and down the temperature ladder.
}

\value{This function outputs a list of the following values:
  \item{pop}{current status of the population.}
  \item{pop.fit}{vector of fitness values for each of the models
    specified in the current population.}
  \item{exc.a}{current number of accepted exchange steps.}
}

\references{ 
  Geyer C (1991). \emph{Markov chain Monte Carlo maximum likelihood}.
  Computing Science and Statistics:156.

  Holland J (1975). \emph{Adaptation in Natural and Artificial Systems}.
  University of Michigan Press.

  Liang F, Wong W (2000). \emph{Evolutionary Monte Carol:Applications to
    Cp Model Sampling and Change Point Problem}. Statistica Sinica 10:317.
}

\author{ Melanie Wilson <maw27@stat.duke.edu> }

\keyword{ methods }

