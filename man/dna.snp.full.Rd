\name{dna.snp.full}
\alias{dna.snp.full}
\docType{data}
\title{Simulated Example SNP data set. }

\description{
  Simulated data set of 60 SNP variables in the DNA repair pathway.  The
  first column represents the disease status of the individual, the second
  column the forced variable of age, and the remaining columns the
  log-additive parametrizations of each SNP variable.  This is the
  full simulated data set before the marginal screen was performed.
}

\usage{data(dna.snp)}
\format{
  A data frame with 1197 observations on the following 60 variables
  (case, age, snp.la).
 }

\details{
  The simulated samples comprised a binary outcome and genetic data on 399
  cases and 798 controls where the genetic data  was simulated at the same
  60 tag SNPs as genotyped in a particular study of interest on the DNA
  repair pathway. We simulated the genotypes in two stages. First, for
  each of the 6 genes represented in the data set, we phased the NCOCS
  control SNP genotype data and estimated recombination rates using PHASE
  (Stephens et al., 2001). Second, given a model of association, we
  generated case-control data at these tags using HAPGEN (Marchini and Su,
  2006).  For this simulation, we assumed that a randomly chosen subset of
  2 genes were associated and that, within the associated genes, a single,
  randomly chosen tag was the source of the association. One of the
  associated tags SNPs were accorded an odds ratio (OR) of 1.75 and an
  assumed genetic log-additive genetic parametrization and one an OR of
  1.25 and a dominate genetic parametrization. Information on the SNPs in
  the simulation can found in the data set \code{\link{sim.info}}.  
}

\references{
  Stephens M, Smith N, Donnelly P (2001). A New Statistical Method for
  Haplotype Reconstruction from Population Data. The American Journal of
  Human Genetics 68:978-989.

  Marchini J, Su Z (2006). HAPGEN, a C++ program for simulating case and
  control SNP haplo- types.  }

\examples{
  data(dna.snp.full)
}
\keyword{datasets}
