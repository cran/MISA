\name{expand.data.snp}
\alias{expand.data.snp}
\alias{data.snp}

\title{Data: Genetic Parametrization}
\description{
  This function takes the data set with log-additive parameterization
  for each SNP and makes a data set with log-additive, dominant, and
  recessive parameterization.  
}
\usage{
expand.data.snp(data.snp, ind.info, force, subset)
}

\arguments{
  \item{data.snp}{(n x p) matrix of log-additive SNP genotypes for each
    individual where 0 is the hom. common, 1 het. and 2 hom. rare genotype.}
  \item{ind.info}{(n x (q+1)) matrix of information on each individual
    where q is the number of forced variables in the study. The first
    column must represent the case/control status followed by a column
    for each of the forced variables.}
  \item{force}{vector of variable names that you wish to force into the
    final model. These variables must be found in the ind.info matrix. }
  \item{subset}{vector of TRUE/FALSE indicating the individuals you
    want to include in the study }
% \item{mdldef}{3 x m matrix, where each column is a model and each row
%   is a genotype (0, 1, 2). Each entry in the matrix determines the value 
%   that goes into the design matrix for the given genotype and model.}
}

\value{
This function outputs a data frame where the first column is the
case/control variable, the next columns are the forced variables and the
last columns are the SNP.la, SNP.dom and SNP.rec variables, respectively.  
}

\author{ Melanie Wilson <maw27@stat.duke.edu> }

\keyword{ methods }

