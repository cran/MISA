\name{dna.snp}
\alias{dna.snp}
\docType{data}
\title{Simulated Example SNP data set.}

\description{
  Simulated data set of 26 SNP variables in the DNA repair pathway.  The
  first column represents the disease status of the individual, the
  second column the forced variable of age and in the remaining columns
  we have the log-additive parameterizations of the SNP variable.  This
  is the subset data set of dna.snp.full once the marginal Bayes factor
  screen was performed.}

\usage{data(dna.snp)}

\format{
  A data frame with 1197 observations on the following 19 variables
(case, age, snp.la).
}

\details{
  The simulated samples comprised a binary outcome and genetic data on 399
  cases and 798 controls where the genetic data was simulated at the same
  60 tag SNPs as genotyped in a particular study of interest on the DNA
  repair pathway.  We simulated the genotypes in two stages. First, for
  each of the 6 genes represented in the data set, we phased the NCOCS
  control SNP genotype data and estimated recombination rates using PHASE
  (Stephens et al., 2001). Second, given a model of association, we
  generated case-control data at these tags using HAPGEN (Marchini and Su,
  2006).  For this simulation, we assumed that a randomly chosen subset of
  2 genes were associated and that, within the associated genes, a single,
  randomly chosen tag was the source of the association. One of the
  associated tags SNPs were accorded an odds ratio (OR) of 1.75 and an
  assumed genetic log-additive genetic parametrization and one an OR of
  1.25 and a dominate genetic parametrization.  To facilitate the mixing
  of the Genetic EMC algorithm and ease computation time we screened the
  SNPs in each full simulation to the final subset of SNPs seen in this
  data set with a marginal Bayes factor estimated to be 1.0 or above by a
  procedure described in \code{\link{bf4assoc}}.  Information on the SNPs
  in this simulation can be found in the data set \code{\link{sim.info}}.
}

\references{
  Stephens M, Smith N, Donnelly P (2001). A New Statistical Method for
  Haplotype Reconstruction from Population Data. The American Journal of
  Human Genetics 68:978-989.

  Marchini J, Su Z (2006). HAPGEN, a C++ program for simulating case and
  control SNP haplotypes.  
}

\examples{
data(dna.snp)
}

\keyword{datasets}
