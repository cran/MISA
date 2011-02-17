`model.inc` <-
function(prob.out, num.models=100, num.snps=20, num.genes=20, inc.typ="s",
         hide.name=FALSE, ...) {
  
  ## Create image plots of the top SNPs and top Genes included in the top
  ## models.
  ##   prob.out:   output list from post.prob(). 
  ##   num.models: the number of the top models to place on the x-axis.
  ##   num.snps:   If inc.type="s", the number of the top SNPs to place on the
  ##               y-axis.
  ##   num.genes:  If inc.type="g",the number of the top genes to place on the
  ##               y-axis.
  ##   inc.typ:    specifies if we want to plot the SNP inclusion ("s") or gene
  ##               inclusion ("g")
  ##   hide.name:  logical indicator determining if we should hide the specific
  ##               gene and SNP names.

  # Fixed for single-gene case. gl37 5/24/10

  which <- prob.out$which
  cost <- prob.out$cost
  snp.info <- as.character(prob.out$Post.SNP[, 1])
 
  snps <- colnames(which)
  postprob <- as.numeric(prob.out$Post.Model[, 1])
  snp.inc <- prob.out$Post.SNP
  gene.inc <- prob.out$Post.Gene
  bf <- as.numeric(snp.inc[, 6])
  probne0 <- as.numeric(snp.inc[, 2])
  bf.g <- as.numeric(gene.inc[, 2])
  probne0.g <- as.numeric(gene.inc[, 1])
 

  genes <- unique(as.character(snp.info))
  which.genes <- matrix(0, nrow=length(snps), ncol=length(genes))
   for(j in 1:length(genes)) {
     which.genes[(snp.info == genes[j]), j] <- 1
   }
  which.g <- which %*% which.genes

  null.post <- postprob[1]
  model.ord <- order(-postprob[-1])
  which <- which[-1,][model.ord, ]  
  which.g <- as.matrix(which.g[-1,])[model.ord,]
  cost <- cost[-1][model.ord]
  postprob <- postprob[-1][model.ord]

  ## Graphic Parameters
  nmodel <- num.models

  if (inc.typ == "s") {
    nvar <- num.snps
    ordr <- order(-bf)
    rownms <- rownames(snp.inc)[ordr][1:nvar]
    genes.s <- snp.inc[ordr, 1][1:nvar]
    rownms <- paste(genes.s, rownms, sep=":  ")
  
  
    if (hide.name == TRUE) {
      snp.hide <- paste("S", c(1:num.snps), sep="")
      gene.hide <- paste("G", c(1:num.genes), sep="")
      genes.o <- genes[order(-bf.g)][1:num.genes]
      rownms.hide <- NULL
      for (s in 1:num.snps) {
        snp.gene <- c(1:num.genes)[genes.o == snp.inc[ordr,][s,1]]
        snp.gene <- paste("G", snp.gene, sep="")
        rownms.hide <- c(rownms.hide, paste(snp.hide[s], snp.gene, sep="."))
      }
      rownms <- rownms.hide[1:nvar]
    }

    clr <- c("#FFFFFF", "#A020F0", "#FF0000", "#0000CD")
    ordr <- ordr[1:nvar]
    color.matrix <- which[1:(nmodel), ordr[1:nvar]] + 1
    prob.labels <- as.numeric(round(bf[ordr[1:nvar]], 2))
  } 

  if (inc.typ == "g") {
    nvar <- min(length(genes), num.genes)
    rownms <- genes[order(-bf.g)][1:nvar]

    if (hide.name==TRUE) {
      rownms <- paste("G", c(1:nvar), sep="")
    }
    
    clr <- c("#FFFFFF", "#6495ED")
    ordr <- order(-bf.g)[1:nvar]
    which.g[which.g > 0] <- 1
    color.matrix <- as.matrix(which.g)[1:(nmodel), ordr[1:nvar]] + 1
    prob.labels <- as.numeric(round(bf.g[ordr[1:nvar]], 2))
  }
    
  ## matrix of colors white, purple(la), blue(dom), red(rec)    
  keep.mar <- par(mar=c(5, 6, 4, 2) + 0.1)
  par(las=1, mar=c(8, 10, 5, 10), ps=12, font=2)
  
  if (inc.typ == "s") {
    maintitle=paste("SNP Inclusions of Top Models")
    if (hide.name==TRUE) {
      maintitle=""
    }
  }
  
  if (inc.typ == "g") {
    maintitle=paste("Gene Inclusions of Top Models")
    if (hide.name==TRUE) {
    	maintitle=""
      }
  }

  if (hide.name == TRUE) {
    maintitle <- NULL
  }
   
  prob.axis <- postprob[1:(nmodel)]
  prob.axis <- prob.axis/sum(prob.axis) 
  image(c(0, cumsum(prob.axis)), 1:nvar, as.matrix(color.matrix), col=clr, 
        xlab="Model", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1),
        main=maintitle)

  xat <- (cumsum(prob.axis) + c(0, cumsum(prob.axis[-nmodel]))) / 2
  axis(1, at=xat, labels=1:nmodel)
  axis(2, at=1:nvar, labels=rownms)
  axis(4, at=1:nvar, labels=prob.labels)
  par(mar=keep.mar)
}

