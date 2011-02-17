bic.prior = function(n=NULL) {
  structure(list(family="BIC", class="IC", hyper=log(n)), class="prior")
}

aic.prior = function() {
  structure(list(family="AIC", class="IC", hyper=2.0), class="prior")
}

IC.prior = function(penalty) {
  structure(list(family="IC", class="IC", hyper=as.numeric(penalty)), 
  class="prior")
}
