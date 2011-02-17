\name{mutation}
\alias{mutation}

\title{ EMC: Mutation Step }
\description{
  This function takes a current state of one of the chains of the
  population in the Genetic EMC algorithm and performs the mutation step.
}

\usage{
mutation(chain, pop, pop.fit, force, data, fitness, t , impute = impute,
         b = NULL, a = 1, homozyg.rare, burnin = TRUE)
}

\arguments{
  \item{chain}{integer specifying which chain (or model) of the
    population to perform the mutation on. }
  \item{pop}{matrix specifying the current status of each chain of the
    population. }
  \item{pop.fit}{vector of fitness values for each of the models
    specified in the current population. }
  \item{force}{character vector specify the variables to force in the models.}
  \item{data}{data frame of the same form in \code{Gene.EMC}  }
  \item{fitness}{character string specifying the fitness function to use
    in the algorithm.}
  \item{t}{temperature vector specifying the temperature value for each
    chain in the population. }
  \item{impute}{number of imputed data sets. }
  \item{b}{If the fitness function is "AIC.BB", the user must specify
    the value for the beta hyper-parameter b.}
  \item{a}{If the fitness function is "AIC.BB", the user must specify
    the value for the beta hyper-parameter a.}
  \item{homozyg.rare}{indicator vector that indicates which SNPs should
    not have a recessive parameter since the power is too weak.}
  \item{burnin}{integer indicating the length of the burnin.}
}

\details{
  In the mutation step we are performing a Metropolis update on the
  population by choosing a model, or current value of one of the chains
  and taking one of the SNP indicators and mutating its state in the
  chosen model. 
}

\value{
  This function outputs a list of the following values:
  \item{samp}{current status of the chain of the current population.}
  \item{samp.fit}{fitness value for the model specified in the chain of
    the current population }
  \item{mut.a}{indicator if the mutation proposal was accpeted.}
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



