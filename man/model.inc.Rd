\name{model.inc}
\alias{model.inc}

\title{ Image Plots for top SNP and Gene Inclusions }
\description{
  This function allows the user to create image plots of the top SNPs
  and top Genes included in the top models.  For the SNP inclusion plots,
  the color of the inclusion block signifies genetic mode of inheritance
  for the specified SNP in each of the modes with: Purple = Log-Add.,
  Red= Dom., Blue = Rec. SNPs and Genes are ordered based on marginal
  SNP inclusion probabilities which are plotted on the right axis.  The
  width of the inclusion blocks are proportional to the posterior model
  probability that the SNP or Gene is included in.
}

\usage{
model.inc(prob.out, num.models = 100, num.snps = 20, num.genes = 20, 
          inc.typ = "s", hide.name = FALSE, ...)
}

\arguments{
  \item{prob.out}{output list from \code{\link{post.prob}}.}
  \item{num.models}{the number of the top models to place on the x-axis.}
  \item{num.snps}{ If inc.type="s", the number of the top SNPs to place
    on the y-axis.}
  \item{num.genes}{ If inc.type="g",the number of the top genes to place
    on the y-axis. }
  \item{inc.typ}{specifies if we want to plot the SNP inclusion ("s") or
    gene inclusion ("g") }
  \item{hide.name}{logical indicator determining if we should hide the
    specific gene and SNP names. }
  \item{...}{General parameters for plotting functions}
}


\author{ Melanie Wilson <maw27@stat.duke.edu> }

\examples{
  ## Load the emc.out files for one run of the dna.snp data (100,000
  ## iterations)
  data(emc.out.1)
  data(sim.info)
  p = dim(emc.out.1$which)[2]-1

  ## Calculate posterior summaries for the output of the Genetic EMC algorithm
  post.prob.out = post.prob(emc.out.1, sim.info, b=p)

  ## Plot the SNP Inclusions in the top 100 models for the top 20 SNPs
  model.inc(post.prob.out, num.models=100, num.snps=20, inc.typ="s")

  ## Plot the Gene Inclusions in the top 100 models for the top 20 Genes
  model.inc(post.prob.out,num.models=100,num.genes=6,inc.typ="g")
}


\keyword{ hplot }