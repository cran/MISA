\name{converge.EMC}
\alias{converge.EMC}
\alias{emc.converge}
\title{ Convergence Diagnostic Plots for Genetic EMC }

\description{
  Takes the output of two independent runs of the Genetic EMC sampling
  algorithm and creates plots to assess the convergence of the algorithm.   
}

\usage{
converge.EMC(emc.out, plot.type, bandwidth = 1000, a = 1, b = NULL,...)
}


\arguments{
  \item{emc.out}{output from \code{\link{Gene.EMC}} }

  \item{plot.type}{indicates the type of plot; options are "iter",
    "gelman.rubin","bayes.factor", "snp.inc", or "all". }

  \item{bandwidth}{If the "bayes.factor" plot.type is chosen, user must
    indicate the bandwidth of iterations over which to calculate the
    global Bayes Factor over.}

  \item{a}{If "bayes.factor" or "snp.inc" plot.type is chosen, and the
    fitness function used in the Genetic EMC algorithm is "AIC.BB" the
    user must specify the value of a for the beta hyper-parameter.}

  \item{b}{If "bayes.factor" or "snp.inc" plot.type is chosen, and the
    fitness function used in the Genetic EMC algorithm is "AIC.BB" the
    user must specify the value of b for the beta hyper-parameter. }

  \item{...}{general parameters for plotting functions.}
}

\details{
  The four plot types are described as follows:

  iter: Plot of the cost values over each iteration for each of the
  independent runs.  These trace plots help to examine if a balanced has
  been reached between exploring the model space and convergence rates.
  Make sure that the cost values are not sticking too much in one area
  and are moving around freely to explore the space properly.  If the
  trace plots do show that the algorithm is tending to get stuck in on
  area this may be a sign that you need to increase the max temperature
  or increase the number of parallel chains so that adjacent chains can
  communicate better.
  
  gelman.rubin: Plot of the gelman rubin convergance diagnostic (see
  Gleman, Rubin (1992)) of the cost values of the two independent
  chains.
  
  bayes.factor: Plot of the global Bayes factor computed across
  iterations for each independent chain.  Since our global Bayes factor is
  a lower bound for the global Bayes factor computed if we were able to
  enumerate all models the Bayes factor will increase for every new unique
  model that we find at any of the iterations.  Therefore, we are
  interested in seeing if the Bayes factor begins to converge after a
  given number iterations and is no longer making large jumps.

  snp.inc: Plot of the Marginal Bayes factors for one independent
  run vs. another independent run.  This plot enables us to 
  determine if the values of the marginal SNP inclusion probabilities
  are consistent across two independent runs of the algorithm. If the
  plot does not following the line y = x then this is an indication that
  the algorithm as not yet converged and is still exploring the space.  
}

\value{
  This function plots four convergence diagnostics for the Genetic EMC
  algorithm.  
}

\references{ 
  Gelman A, Rubin D (1992). \emph{Inference from iterative simulation
  using multiple sequences}. Statistical Science 7:457.
}

\author{ Melanie Wilson <maw27@stat.duke.edu> }

\examples{
##Load the emc.out files for two independent runs of the dna.snp data 
data(emc.out.1)
p <- dim(emc.out.1$which)[2] - 1
data(emc.out.2)
emc.out <- list(emc.out.1, emc.out.2)

##Look at all of the convergence plots for the dna.snp data
converge.EMC(emc.out,plot.type="all", bandwidth=100 ,b=p, a=1)
}

\keyword{ hplot }
