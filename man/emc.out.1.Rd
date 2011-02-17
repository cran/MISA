\name{emc.out.1}
\alias{emc.out.1}
\docType{data}

\title{Example Output I from the Genetic EMC algorithm}

\description{
  Output from the first independent run of the Genetic EMC algorithm for the 
  dna.snp data set.  The algorithm was run for approximately 500,000
  iterations.  The assumed fitness function was "AIC.BB" where we
  assumed a Beta-Binomial prior on the size of the models sampled with
  hyper-parameters on the beta distribution, a=1 and b=p (p=82).
}

\usage{data(emc.out.1)}

\format{
  Output from the function \code{\link{Gene.EMC}} that is a list of the
  following values:
  \describe{
    \item{which:}{Matrix where each row corresponds to a model
      specification vector for the unique models visited and the value
      of the fitness function for the model.}
    \item{data:}{Data frame of the same form as the imputed data frame.}
    \item{iter.aic:}{The value of the fitness function for each of the
      models visited at each of the iterations.}
    \item{iter:}{The total number of iterations run.}
    \item{iter.unique:}{Vector indicating the number of the iteration in
      which each of the unique models was found.}
    \item{force:}{Character vector indicating the names of the forced
      variables.}
    \item{fitness:}{Character string indicating the fitness function
      used for the algorithm.}
  }
}

\examples{
  data(emc.out.1)
}

\keyword{datasets}
