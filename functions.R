# Functions needed for INLA analyses
# Part of mtDNA project
# Update: May 2020


SummarizeFun = function(x, Quantiles = c(0.025, 0.975)) {
  c(mean(x), sd(x), quantile(x, probs = Quantiles))
}


SummarizeInlaVars = function(x, nSamples = 10000) {
  # Summarize INLA effects "precisions" in form of Standard deviations, Variances, and Proportions (var / sum(all vars))
  Terms = names(x$marginals.hyperpar)
  Terms = Terms[grepl(pattern = "Precision for ", x = Terms)]
  TermsShort = sub(pattern = "Precision for ", replacement = "", x = Terms)
  TermsShort = sub(pattern = "the ",           replacement = "", x = TermsShort)
  nTerms = length(Terms)
  Samples = matrix(data = numeric(), nrow = nSamples, ncol = nTerms + 1)
  Out = vector(mode = "list", length = 3)
  names(Out) = c("Sd", "Var", "Proportion")
  Out[[1]] = Out[[2]] = Out[[3]] = matrix(data = numeric(), nrow = nTerms + 1, ncol = 4)
  dimnames(Out[[1]]) = dimnames(Out[[2]]) = dimnames(Out[[3]]) = list(c(TermsShort, "Total"), c("Mean", "Sd", "Q0.025", "Q0.975"))
  for (Term in 1:nTerms) {
    # Term = 3
    Samples[, Term] = 1 / inla.rmarginal(n = nSamples, marginal = x$marginals.hyperpar[[Terms[Term]]])
  }
  Samples[, Term + 1] = rowSums(x = Samples[, 1:nTerms, drop = FALSE])
  Out$Var[]        = t(apply(X = Samples,                         MARGIN = 2, FUN = SummarizeFun))
  Out$Sd[]         = t(apply(X = sqrt(Samples),                   MARGIN = 2, FUN = SummarizeFun))
  Out$Proportion[] = t(apply(X = Samples / Samples[, nTerms + 1], MARGIN = 2, FUN = SummarizeFun))
  return(Out)
}


SummarizeInlaSpdeVars = function(x, nSamples = 10000, name = 'spatial', spde = SpdeStat) {
  # Summarize INLA SPDE hyper-parameters
  SpdeParam = inla.spde2.result(inla = x, name = name, spde = spde)
  Samples = matrix(data = numeric(), nrow = nSamples, ncol = 2)
  Out = vector(mode = "list", length = 3)
  names(Out) = c("Sd", "Var", "Range")
  Out[[1]] = Out[[2]] = Out[[3]] = matrix(data = numeric(), nrow = 1, ncol = 4)
  colnames(Out[[1]]) = colnames(Out[[2]]) = colnames(Out[[3]]) = c("Mean", "Sd", "Q0.025", "Q0.975")
  Samples[, 1] = inla.rmarginal(n = nSamples, marginal = SpdeParam$marginals.variance.nominal[[1]])
  Samples[, 2] = inla.rmarginal(n = nSamples, marginal = SpdeParam$marginals.range.nominal[[1]])
  Out$Sd[]    = SummarizeFun(x = sqrt(Samples[, 1]))
  Out$Var[]   = SummarizeFun(x =      Samples[, 1])
  Out$Range[] = SummarizeFun(x =      Samples[, 2])
  return(Out)
}


# Summarize from inla.posterior.sample
SummarizeINLApostsample = function(x, nSamples = 10000) {
  # Summarize INLA effects "precisions" in form of Standard deviations, Variances, and Proportions (var / sum(all vars))
  Terms = names(x[[nSamples]]$hyperpar)
  Terms = Terms[grepl(pattern = "Precision for ", x = Terms)]
  TermsShort = sub(pattern = "Precision for ", replacement = "", x = Terms)
  TermsShort = sub(pattern = "the ",           replacement = "", x = TermsShort)
  nTerms = length(Terms)
  Samples = matrix(data = numeric(), nrow = nSamples, ncol = nTerms + 1)
  Out = vector(mode = "list", length = 3)
  names(Out) = c("Sd", "Var", "Proportion")
  Out[[1]] = Out[[2]] = Out[[3]] = matrix(data = numeric(), nrow = nTerms + 1, ncol = 4)
  dimnames(Out[[1]]) = dimnames(Out[[2]]) = dimnames(Out[[3]]) = list(c(TermsShort, "Total"), c("Mean", "Sd", "Q0.025", "Q0.975"))
  for (Term in 1:nTerms) {
    # Term = 3
    for (i in 1:nSamples) {
    Samples[i, Term] = 1 / x[[i]]$hyperpar[[Terms[Term]]]
    }
  }
  Samples[, Term + 1] = rowSums(x = Samples[, 1:nTerms, drop = FALSE])
  Out$Var[]        = t(apply(X = Samples,                         MARGIN = 2, FUN = SummarizeFun))
  Out$Sd[]         = t(apply(X = sqrt(Samples),                   MARGIN = 2, FUN = SummarizeFun))
  Out$Proportion[] = t(apply(X = Samples / Samples[, nTerms + 1], MARGIN = 2, FUN = SummarizeFun))
  return(Out)
}




#### Calculate if effect different that zero ####

ProbDiffFromZero = function(x) {
  # Probability of effect being different from zero
  # x - numeric vector
  Sel = !is.na(x)
  x = x[Sel]
  n = sum(Sel)
  ProbG = sum((0 + .Machine$double.eps) < x) / n
  ProbS = sum(x < (0 - .Machine$double.eps)) / n
  # Prob0 = sum((0 - .Machine$double.eps) < x & x < (0 + .Machine$double.eps)) / n
  # Prob0 = 1 - ProbG - ProbS
  return(max(c(ProbG, ProbS)))
}

# FATRidgeMarEffP = apply(X = RidgeMarEffList$FAT, MARGIN = 2, FUN = ProbDiffFromZero)


MinusLogProb = function(x, Min = NULL) {
  # x - numeric vector of p-values
  # Min - numeric value below which we cut x (to avoid log(0) = Inf)
  if (!is.null(Min)) {
    x[x < Min] = Min
  }
  return(-log(x))
}

# FATRidgeMarEffMLP = MinusLogProb(x = 1 - FATRidgeMarEffP, Min = (1 / nIter) / 2)

# Samples are from one posterior distribution and we are here calculating one p-value
# I remember now. nIter  in MinusLogProb is to avoid doing log(0) .
# I asummed that if prob is smaller than (1 / nIter) / 2 then I just set i to that value.
# nIter from MCMC is the same as nSamp  you take from INLA posteriors



##### FUNKCIJE ZA TREE MODEL #####

expandHaplotypeDf = function(haplotypeDf, noMut, phantomStart = 300){
  # haplotypeDf should be a dataframe with 2 columns (parent haplotype in column 1 and child haplotype in column 2)
  # noMut should be a vector with no of mutations between each haplotype pair in haplotypeDf
  # phantomStart should be higher than highest haplotype code
  
  myHap = as.matrix(haplotypeDf)
  
  for(i in 1:nrow(haplotypeDf)){
    numMut = noMut[i]
    
    if(numMut>1){
      
      phantom = c(phantomStart:(phantomStart+numMut-2))
      nPhantom = length(phantom)
      
      
      if(length(phantom)==1){
        myHap = rbind(myHap,  c(myHap[i,1] ,phantom[1]),  c(phantom[1], myHap[i,2]) )
      }else{
        middlePhantom = matrix( c(phantom[1:(nPhantom-1)], phantom[2:nPhantom]  ) ,ncol = 2)
        
        myHap = rbind(myHap, c(myHap[i,1] ,phantom[1]), middlePhantom,c(phantom[nPhantom], myHap[i,2]) )
        
      }
      phantomStart = phantomStart + numMut - 1
    }
  }
  
  myHap = myHap[  -c(which(noMut>1)) ,]
  
  
  return(myHap) 
}  



# Construct Q so can sample from model
build.Q.haplotype.network = function(hapNet, p = list(prec = 1, rho = 0.6)){
  
  nHap = length(unique(c(hapNet[,1], hapNet[,2])))
  
  A = matrix(0, ncol =nHap,nrow = nHap)
  diag(A) =1
  for(i in 1:nrow(hapNet)){
    theChild = hapNet[i,2]
    nrParentsOfChild = sum(hapNet[,2] ==theChild)
    A[theChild,hapNet[i,1]] = -p$rho/nrParentsOfChild # The parent
    
  }
  
  # Most haplotypes have 1/sigma_h^2 on diagonal,which is 1(sigma_0^2*(1-rho^2)). Here we use the inverse variance
  P = rep(p$prec/(1-p$rho^2),nHap) 
  # Ancestor haplotypes have 1/sigma_0^2 on the diagonal. Here we use the inverse variance
  myAnc = unique(hapNet[ which( !(hapNet[,1] %in% hapNet[,2] )),1])
  P[myAnc]=p$prec
  # Haplotypes with multiple parents have, lower variance. Here, we use the inverse variance, so we multiply rather than dividing.
  nParents = sapply(1:nHap, function(X) length(which( hapNet[,2] == X )   ))
  nParents[nParents==0] = 1
  P = P*nParents
  
  Q =( (t(A)%*%diag(P))  %*%A) 
  return(inla.as.sparse(Q))
}

# regeneric function implementing the haplotype network model
inla.rgeneric.haplotype.network = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),  theta = NULL){
  # for reference and potential storage for objects to
  # cache, this is the environment of this function
  # which holds arguments passed as `...` in
  # `inla.rgeneric.define()`.
  envir = environment(sys.call()[[1]])
  
  library(plyr)
  
  interpret.theta = function(){
    return(list(prec = exp(theta[1L]),rho = 2 * exp(theta[2L])/(1 + exp(theta[2L])) - 1))
  }
  
  graph = function(){
    return(Q())
  }
  
  Q = function(){
    p = interpret.theta()
    
    nHap = length(unique(c(hapNet[,1], hapNet[,2])))
    # Make matrix that defines the tree
    A = matrix(0, ncol =nHap,nrow = nHap)
    diag(A) =1
    for(i in 1:nrow(hapNet)){
      theChild = hapNet[i,2]
      nrParentsOfChild = sum(hapNet[,2] ==theChild)
      A[theChild,hapNet[i,1]] = -p$rho/nrParentsOfChild # The parent
      
    }
    
    # Most haplotypes have 1/sigma_h^2 on diagonal,which is 1(sigma_0^2*(1-rho^2)). Here we use the inverse variance
    P = rep(p$prec/(1-p$rho^2),nHap) 
    # Ancestor haplotypes have 1/sigma_0^2 on the diagonal. Here we use the inverse variance
    myAnc = unique(hapNet[ which( !(hapNet[,1] %in% hapNet[,2] )),1])
    P[myAnc]=p$prec
    # Haplotypes with multiple parents have, lower variance. Here, we use the inverse variance, so we multiply rather than dividing.
    nParents = sapply(1:nHap, function(X) length(which( hapNet[,2] == X )   ))
    nParents[nParents==0] = 1
    P = P*nParents
    
    Q =( (t(A)%*%diag(P))  %*%A) 
    
    
    return(inla.as.sparse(Q))
  }
  
  mu = function(){
    return(numeric(0))
  }
  
  log.norm.const = function(){
    # Since the normalizing constant depends on the phylogeny, we let INLA compute it. 
    return(numeric(0))
  }
  
  log.prior = function(){ # This function must return the log-prior of the prior density for theta
    # Log PC priors for the precision and auto-correlation
    
    p = interpret.theta()
    
    val =  inla.pc.dprec(p$prec, uPrec, alphaPrec, log = TRUE) +  theta[1L]  + inla.pc.dcor1(p$rho,uRho, alphaRho, log = T) + log(2) - abs(theta[2L]) - 2*log(1 + exp(-abs(theta[2L])))       
    # The term, theta[1] is the Jacobian for the change of variable from precision to log precision, and the final term is the Jacobian from the change of variable from rho to theta_2
    
    return (val)
  }
  
  initial = function(){
    # The initial values of theta
    return (c(4,1))
  }
  
  quit = function(){
    return (invisible())
  }
  
  # Sometimes this is useful, as argument 'graph' and 'quit'
  # will pass theta=NULL as the values of theta are not
  # required for defining the graph. However, this statement
  # will ensure that theta is always defined.
  if (is.null(theta)) {theta = initial() }
  val = do.call(match.arg(cmd), args = list())
  return (val)
}


