#F2R, t2R, r.CI
# made by joe hilgard in early 2014 *flexes biceps*
F2R = function(Fstat, N, width=.95, neg=F) {
  r.equiv = sqrt(Fstat/(Fstat + N - 2)) # is this where the loss of fidelity happens?
  if(neg==T) {r.equiv=-(r.equiv)}
  zScore = 1/2 * log((1+r.equiv)/(1-r.equiv))
  z.se = 1/sqrt(N-3)
  margin = -qnorm((1-width)/2)
  z.low = r.equiv-margin*z.se
  z.hi = r.equiv+margin*z.se
  r.equiv.low = (exp(2*z.low)-1)/(exp(2*z.low)+1)
  r.equiv.hi = (exp(2*z.hi)-1)/(exp(2*z.hi)+1)
  return(c(r.equiv.low, r.equiv, r.equiv.hi))
}

t2R = function(tstat, N, digits=2) {
  neg = tstat<0
  Fstat = tstat^2
  r.equiv = sqrt(Fstat/(Fstat + N - 2))
  if(neg==T) {r.equiv=-(r.equiv)}
  zScore = 1/2 * log((1+r.equiv)/(1-r.equiv))
  z.se = 1/sqrt(N-3)
  z.low = r.equiv-1.96*z.se
  z.hi = r.equiv+1.96*z.se
  r.equiv.low = (exp(2*z.low)-1)/(exp(2*z.low)+1)
  r.equiv.hi = (exp(2*z.hi)-1)/(exp(2*z.hi)+1)
  print(paste("Point estimate:", r.equiv))
  print(paste("95% CI: [", r.equiv.low, ", ", r.equiv.hi, "]", sep=""))
  return(c(r.equiv.low, r.equiv, r.equiv.hi))
}

r.CI=function(r.equiv, N) {
  zScore = 1/2 * log((1+r.equiv)/(1-r.equiv))
  z.se = 1/sqrt(N-3)
  z.low = r.equiv-1.96*z.se
  z.hi = r.equiv+1.96*z.se
  r.equiv.low = (exp(2*z.low)-1)/(exp(2*z.low)+1)
  r.equiv.hi = (exp(2*z.hi)-1)/(exp(2*z.hi)+1)
  return(c(r.equiv.low, r.equiv, r.equiv.hi))
}

pool.sd = function (sds, ns) {
  SSlist = sds^2 %*% (ns-1)
  pool.var = sum(SSlist) / sum(ns-1)
  return(sqrt(pool.var))
}

Cramer = function(chisq, n, k=2) {
  v = sqrt(chisq/(n*(k-1)))
  return(v)
}

#Computation of r according to CMA

stderr.d = function(d, n1, n2) {
  term1 = (n1+n2)/(n1*n2) + d^2 / (2*(n1+n2-2))
  term2 = (n1+n2)/(n1+n2-2)
  return(term1*term2)
}

# Odds Ratio into Cohen's d
# Hasselblad & Hedges (1995) technique
OR.to.d = function(OR=NULL, b=NULL) {
  if (!is.null(OR)) b1 = log(OR)
  if (!is.null(b) & !is.null(OR)) if (b1 != b) print("Nonmatching OR and b! One or the other, please.")
  if (is.null(b)) b = b1
  return(b * sqrt(3)/pi)
}


d2r = function(d, n1, n2, width=.95) {
  r = d / sqrt(d ^ 2 + 4)
  term1 = (n1+n2)/(n1*n2) + d^2/(2*(n1+n2-2))
  term2 = (n1+n2)/(n1+n2-2)
  StdErr.d = sqrt(term1*term2)
  a = ((n1+n2)^2)/(n1*n2)
  StdErr.r = sqrt(a^2 * StdErr.d ^ 2 / ((d ^ 2 + a) ^ 3))
  return(list("r"=r, "StdErr.r"=StdErr.r))
}

d2r2z = function(d, n1, n2) {
  r = d2r(d, n1, n2)[[1]]
  StdErr.r = d2r(d, n1, n2)[[2]]
  FisherZ = 0.5 * log((1 + r) / (1 - r))
  StdErr.z = StdErr.r / (1 - r ^ 2) 
  return(list("Z" = FisherZ, "StdErr.z" = StdErr.z))
}

fucker = function(r, n) {
  se.r = sqrt((1-r^2)/(n-2))
  se.r.2 = sqrt(1-r^2) / sqrt(n-2)
  print(paste(round(se.r,2), round(se.r.2, 2)))
  stopifnot(se.r - se.r.2 < .001)
  #
  se.z = se.r / (1 - r^2)
  se.z.2 = 1 / (sqrt(1-r^2) * sqrt(n-2))
  stopifnot(se.z - se.z.2 < .001)
  return(list("se.r"=se.r, "se.z"=se.z))
}
# 
# r = 3,000 / Sqr(3,000 ^ 2 + 4) = 0,832 # assuming n1 = n2 
# StdErr(r) = Sqr(16 * 0,149 ^ 2 / ((3,000 ^ 2 + 4) ^ 3)) = 0,013
# 
# Computation of Z
# FisherZ = 0.5 * Log((1 + r) / (1 - r))
# StdErr(FisherZ) = StdErr(r) / (1 - r ^ 2)
# 
# FisherZ = 0.5 * Log((1 + 0,832) / (1 - 0,832)) = 1,195
# StdErr(FisherZ) = (0,013) / (1 - 0,832 ^ 2) = 0,041

# ask_joe, for moral and practical guidance in hard times.
# Hilgard fortune cookie
ask_joe = function() {
  advice = c(
    "Are you cleaning your data with dplyr?",
    "Did you look at histograms of your DV?",
    "Maybe the DV is Poisson-distributed.",
    "Are you looking at accuracy rates? You should probably be using a logistic HLM.",
    "I hope you've got this gathered into an R Project.",
    "I hope you're backing everything up on GitHub.",
    "I hope you're going to share the data and code on OSF.",
    "Did you know? You can link your GitHub account to your OSF account.",
    "Did you know? If you need a private GitHub repo, go to collaborate.missouri.edu.",
    "Should those columns be rows? Use tidyr.",
    "Follow me on Twitter at @JoeHilgard.",
    "Follow Hadley Wickham on Twitter at @HadleyWickham",
    "Publish your null results, in PLOS or Frontiers if you need to.",
    "If your project's worth doing, it's worth taking the time now to preregister it.",
    "Starting a new project? Figure out what statistical test you'll do. Decide now, not later.",
    "What would it take to convince you that your hypothesis is right? What would it take to convince you that your hypothesis is wrong?",
    "You should probably be using ggplot2.",
    "You should probably be using a heirarchical linear model.",
    "Don't turn a continuous variable categorical unless you have a damn good reason."
    
  )
  cat(sample(advice, 1))
}
