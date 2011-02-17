\name{fit.EMC}
\alias{fit.EMC}
\alias{fit}

\title{ EMC: Fitness Function }
\description{
This function takes one of the models in the population (or a current state of a chain) and calculates the fitness/cost value of the model.  
}

\usage{
fit.EMC(samp, force = NULL, data, fitness = "AIC", impute = impute, b =
        NULL, a = 1, homozyg.rare)
}

\arguments{
  \item{samp}{vector specifying the current value of one of the chains
    of the population: or a model specification vector. }
  \item{force}{character vector specify the variables to force in the
    model.}
  \item{data}{data frame of the same form in \code{Gene.EMC}.}
  \item{fitness}{character string specifying the fitness function to
    use in the algorithm.}
  \item{impute}{number of imputed data sets.}
  \item{b}{If the fitness function is "AIC.BB", this specifies the
    hyper-parameter in the beta distribution on the model prior. }
  \item{a}{If the fitness function is "AIC.BB", this specifies the
    hyper-parameter in the beta distribution on the model prior. }
  \item{homozyg.rare}{Indicator vector that indicates which SNPs should
    not have a recessive parameter since the power is too weak.}  }

\details{
  The different options for the fitness function are "AIC", "BIC", or
  "AIC.BB" where we add in a penalty to each of the models that
  corresponds to a beta-binomial prior on the model size with a and b
  chosen by the user.  }

\value{
  The function returns the cost value (-1/2*fitness value) of the model
  corresponding to the model specification vector given in samp.
}

\references{ 
  Geyer C (1991). \emph{Markov chain Monte Carlo maximum likelihood}.
  Computing Science and Statistics:156.

  Holland J (1975). \emph{Adaptation in Natural and Artificial Systems}.
  University of Michigan Press.

  Liang F, Wong W (2000). \emph{Evolutionary Monte Carol:Applications to
    Cp Model Sampling and Change Point Problem}. Statistica Sinica 10:317.
}

\author{ Melanie Wilson <maw7@stat.duke.edu>}



\keyword{ methods }

