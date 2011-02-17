## For backward-compatibility
`data.snp` <- function(...) { expand.data.snp(...) }

`expand.data.snp` <- 
function(data.snp, ind.info, force, subset=rep(TRUE, dim(data.snp)[1])){
  ## Makes a data frame of type data.snp for use with EMC 
  ##   data.snp: (n x p) matrix of log-additive snp genotypes for each
  ##             individual where p is the number of snps under investigation
  ##   ind.info: (n x c) matrix of information on each individual where
  ##             c is the number of covariates in the study. The first column
  ##             must be case/control status, the next length(force) columns
  ##             must be forced variables, and the remaining columns must
  ##             be SNP copy numbers.
  ##   force:    Vector of variable names that you wish you force into the
  ##             final model.  These variables must be found in the ind.info
  ##             matrix.
  ##   subset:   Vector of Trues and Falses based on if we want to consider
  ##             that case or not.
  ##   ### NOT YET IMPLEMENTED ##
  ##   mdldef:   A 3 x m matrix, where each column is a model and each row
  ##             is a genotype (0, 1, 2). Each entry in the matrix determines
  ##             the value that goes into the design matrix for the given
  ##             genotype and model. Example:
  ##
  ##                Genotype        Code  LogAdd  Dom Rec Overdom Underdom
  ##             homozyg wild         0      0     0   0     0       1
  ##             heterozygous         1      1     1   0     1       0
  ##             homozyg mutant       2      2     1   1     0       1
  ##
  ## Returns an expanded n x (m * p) data frame, which may be thought of as
  ## a set of m n x p matrices, each coded according to the corresponding
  ## model.

  mdldef <- cbind(c(0, 1, 2), c(0, 1, 1), c(0, 0, 1))
  colnames(mdldef) <- c("la", "dom", "rec")

  ## get forced variables
  n.force <- length(force)
  f <- NULL
  if (n.force > 0) {
    for(i in 1:n.force) {
      f <- cbind(f,ind.info[, names(ind.info) == force[i]])
    }
    f <- data.frame(f)
    names(f) <- force
  }
  
  ## get case variable
  if (n.force == 0) {
    case <- data.frame(ind.info)
  } else {
    case <- data.frame(ind.info[[1]])
  }
  names(case) <- "case"

  # Generate a command to cbind() all the design matrices for each model.
  snpmdl <- list()
  cbstr <- "cbind(case"
  if (n.force > 0) {
    cbstr <- paste(cbstr, ",f", sep="")
  }
  for (m in 1:dim(mdldef)[2]) {
    snpmdl[[m]] = matrix(mdldef[as.matrix(data.snp) + 1, m],
                         nrow=dim(data.snp)[1])
    cbstr <- paste(cbstr, ",snpmdl[[", m, "]]", sep="")
    colnames(snpmdl[[m]]) <- paste(names(data.snp), colnames(mdldef)[m],
                                   sep=".")
  }
  cbstr <- paste(cbstr, ", deparse.level=0)", sep="")
  data = eval(parse(text=cbstr))
  
  data <- data[subset,]
  return(data)
}
		
	
