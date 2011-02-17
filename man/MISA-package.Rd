\name{MISA-package}
\alias{MISA-package}
\alias{MISA}
\docType{package}
\title{Multilevel Inference for SNP Association Studies}
\description{
The functions in this package focus on intermediate
throughput case-control association studies, where the outcome
of interest is often a binary disease state and where the
genetic markers have been chosen to capture variation in a set
of related genes, such as those involved in a specific
biochemical pathway. Given this data, we are interested in
addressing two questions: "To what extent does the data support
an overall association between the pathway and outcome of
interest?" and "Which markers or genes are most likely to be
driving this association?" To address both of these
questions, this package performs a Bayesian model search
technique that utilizes Evolutionary Monte Carlo and searches
over models including main effects of all genetic markers and
marker-specific genetic effects in a computationally efficient
manner.  The package incorporates functions that perform a marginal
screen on the genetic markers, summarize the
output of the model search algorithm, including image plots of
the models with the highest posterior probability, marginal
summaries of SNP and gene inclusion probabilities and Bayes
Factors, and global summaries of the posterior probability and
Bayes Factor giving evidence of an association in the set of
SNPs of interest.
}

\details{
\tabular{ll}{
Package: \tab MISA\cr
Version: \tab 2.6.5\cr
Date: \tab 2010-07-08\cr
Depends: \tab R (>= 2.10.0)\cr
License: \tab GPL-2\cr
URL: \tab http://www.stat.duke.edu/gbye/MISA.html\cr
}
}

\author{Melanie Wilson, 
  \cr
  Maintainer: Gary Lipton <gl37@stat.duke.edu>
}
\references{
  Geyer C (1991). \emph{Markov chain Monte Carlo maximum likelihood}. Computing
  Science and Statistics:156.

  Holland J (1975). \emph{Adaptation in Natural and Artificial Systems}.
  University of Michigan Press.

  Liang F, Wong W (2000). \emph{Evolutionary Monte Carol:Applications to Cp
    Model Sampling and Change Point Problem}. Statistica Sinica 10:317.
}

\keyword{package}
\keyword{snp, SNP}
