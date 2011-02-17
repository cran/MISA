\name{bf4assoc}
\alias{bf4assoc}
\title{Bayes Factors for 3 Association Models: Marginal SNP Screen}

\description{
  This function performs a marginal screen on the SNPs of interest by
  using Laplace approximations to estimate the marginal Bayes Factors
  (BFs) of each SNP. In particular, we estimated the marginal likelihood
  of each of the three genetic models of association (log--additive,
  dominant and recessive) and under the null model (model of no genetic
  association).  The BF for a model of association is defined as the ratio
  of the marginal likelihood of that model of association to the marginal
  likelihood of the null model.  We can then use a decision rule such as
  including only the SNPs with a maximum marginal Bayes factor for each
  genetic model greater than 1.
}

\usage{
bf4assoc(D, X = NULL, XS, Ns, Nx, snpsd, Prior, MinCount, MaxIt = 100,
         RelTol = 1e-4, scoring = 1)
}

\arguments{
  \item{D}{
    response vector coded (0=control, 1=case) of length N.}

  \item{X}{
    numeric matrix of confounder/design variables of dimension N by Nx
    is NULL if Nx==0.}
  
  \item{XS}{
    numeric matrix of SNP variables of dimension N by Ns coded (0 = common
    homozygote, 1 = heterozygote, 2 = rare homozygote, 3 = missing).}
  
  \item{Ns}{
    number of SNPs in the input data set.}

  \item{Nx}{
    number of design/confounder variables included in all models.}

  \item{snpsd}{standard deviation of mean zero prior on the genetic
    effect parameter when Prior==0 and scale when Prior==1.}

  \item{Prior}{
    set to 0 to chose Normal prior and 1 to chose Cauchy prior.}

  \item{MinCount}{
    count below which a genotype is treated as absent for a given
    SNP. Effect is to reduce the number of unique genetic models that
    are discernable for that SNP.}
  
  \item{MaxIt}{
    maximum number of function evaluations per optimization.}

  \item{RelTol}{
    relative tolerance.}

  \item{scoring}{optimization algorithm; must be either scoring
    (scoring = 1) or Newton-like (scoring = 0).}
}

\details{
  Use Laplace Approx to calculate Bayes factors in favor of 3 genetic
  models of association. Features:(1) Laplace approximation based
  estimates. (2) Output log-ORs + SE(logOR)s.(3) Option for Cauchy
  prior.
}

\value{
  This function outputs a matrix of the following values:
  \item{SNP}{
    SNP ID number.}
  
  \item{N.geno}{
    number of genotypes with counts exceeding MinCount
    among both cases and controls. If 1, the snp is treated as if it
    were monomorphic -- no models are fit; if 2, only the log-additive
    model is fit.}
  
  \item{bfAtoN}{
    marginal likelihood ratio: Pr(Data|logadditive model)/Pr(Data|null
    model) ie BF for log-additive model.}
  
  \item{bfDtoN}{
    marginal likelihood ratio: Pr(Data|dominant model)/Pr(Data|null
    model) or BF for dominant model.}
  
  \item{bfRtoN}{
    marginal likelihood ratio: Pr(Data|recessive model)/Pr(Data|null
    model) or BF for recessive model.}
  
  \item{PrAgvnAssoc}{
    posterior probability of a log-additive genetic model given an
    association.}
  
  \item{PrDgvnAssoc}{
    posterior probability of a dominant genetic model given an
    association.}
  
  \item{PrRgvnAssoc}{
    posterior probability of a recessive genetic model given an association.}

  \item{logOR.LogAdd}{
    modal estimate of log odds ratio under the log-additive model.}
  
  \item{logOR.Dom}{
    modal estimate of log odds ratio under the dominant model.}
  
  \item{logOR.Rec}{modal estimate of log odds ratio under the recessive
    model.}
  
  \item{SE.lOR.LogAdd}{
    estimate of the standard error for the modal estimate of the log
    odds ratio for the SNP variable under the log-additive model.}
  
  \item{SE.lOR.Dom}{
    estimate of the standard error for the modal estimate of the log
    odds ratio for the SNP variable under the dominant model.}

  \item{SE.lOR.Rec}{
    estimate of the standard error for the modal estimate of the log
    odds ratio for the SNP variable under the recessive model.}
}

\examples{
## Load the data
data(dna.snp.full)
## Find the number of snps in the data set
p <- (dim(dna.snp.full)[2] - 2)

## Calculate the Marginal BF's for the SNPs
marg.bf <- bf4assoc(D=dna.snp.full$case,
                    X=as.matrix(dna.snp.full$age),
                    XS=as.matrix(dna.snp.full[,-c(1,2)]),
                    Ns=p, Nx=1, snpsd=0.25, Prior=0, MinCount=1.9,
                    MaxIt=1000, RelTol=1e-7)
## Calculate the Maximum BF for each SNP (LA, Dom, Rec)
max.bf <- apply(marg.bf[,c(3:5)], 1, max)

## Screen the data based on max.bf > 1
dna.snp <- dna.snp.full[, c(TRUE, TRUE, max.bf > 1)]
}

\author{ Ed Iversen <iversen@stat.duke.edu> }

\keyword{ methods }
