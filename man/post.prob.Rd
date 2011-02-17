\name{post.prob}
\alias{post.prob}

\title{ Calculates Global and Marginal Summaries }

\description{
  This function calculates the global and marginal posterior
  probabilities and Bayes Factors that give the evidence of there being an
  association in the overall set of SNPs of interest, the individual genes
  of interest and the individual SNPs of interest.
}

\usage{
post.prob(emc.out, sim.info, b = NULL, a = 1) }

\arguments{
  \item{emc.out}{Output from \code{\link{Gene.EMC}} }
  \item{sim.info}{Vector of character strings giving the names of the
    genes for each of the SNPs in the model, or a data frame that
    a column labeled "SNP" and a column labeled "Gene". This
    data frame may contain all SNPs in the data, not just the ones in
    the model. SNPs must appear in the same order that they appear in 
    the input data frame for Gene.EMC.}
  \item{b}{If the fitness function used in the Genetic EMC algorithm is
    "AIC.BB" the user must specify the value of b for the beta
    hyper-parameter. }
  \item{a}{If the fitness function used in the Genetic EMC algorithm is
    "AIC.BB" the user must specify the value of a for the beta
    hyper-parameter. }
}
\details{
  Global and marginal summaries are computed based on calculating the
  posterior probabilities of each of the unique models that were visited
  in the Genetic EMC algorithm.  The global summaries included a posterior
  probability of association in the overall set of SNPs and Bayes Factor
  for the hypothesis that there is an association in the overall set.  The
  marginal summaries are calculated at the gene and the SNP level.  At the
  gene level, posterior probabilities and Bayes Factors are computed for
  the overall evidence of at least one of the SNPs within the gene of
  interest being associated.  At the SNP level, posterior probabilities
  and Bayes Factors are computed for the evidence of an association within
  the given SNP and posterior probabilities of the most likely genetic
  mode of inheritance of the SNP given that it is associated is computed
  (for the log-additive, dominant, and recessive models).
}

\value{
  The output of the function is a list of the following values:
  \item{Post.Model}{the posterior probability of each of the unique
    models visited by the Genetic EMC algorithm}
  \item{BF.Assoc}{matrix with the global posterior probability and Bayes
    Factor of an overall association and prior odds of H0:Ha}
  \item{Post.SNP}{matrix where each row gives the gene name, inclusion
    probability, probability of log-additive, dominate, and recessive
    genetic mode of inheritance respectively, and the Bayes Factor for
    association for a given SNP. } 
  \item{Post.Gene}{matrix where each row gives the inclusion probability
    and Bayes Factor for the evidence of an association for a given gene}
}

\author{ Melanie Wilson <maw27@stat.duke.edu> }

\examples{
## Combine emc.out.1 and emc.out.2 to get results from both runs.
data(emc.out.1, emc.out.2)
emc.out <- combine.EMC(list(emc.out.1, emc.out.2))

data(sim.info)
p <- dim(emc.out.1$which)[2] - 1

##Calculate Posterior Quantities
post.prob.out <- post.prob(emc.out, sim.info, b=p)
}


\keyword{ methods }
