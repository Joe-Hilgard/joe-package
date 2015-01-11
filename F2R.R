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
  SSlist = sds %*% (ns-1)
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
  if (!is.null(b)) if (b1 != b) print("Nonmatching OR and b! One or the other, please.")
  if (is.null(b)) b = b1
  return(b * sqrt(3)/pi)
}


# d2r = function(d)
# 
# r = d / Sqr(d ^ 2 + 4)
# StdErr(r) = Sqr(16 * StdErr(d) ^ 2 / ((d ^ 2 + 4) ^ 3))
# 
# r = 3,000 / Sqr(3,000 ^ 2 + 4) = 0,832
# StdErr(r) = Sqr(16 * 0,149 ^ 2 / ((3,000 ^ 2 + 4) ^ 3)) = 0,013
# 
# Computation of Z
# FisherZ = 0.5 * Log((1 + r) / (1 - r))
# StdErr(FisherZ) = StdErr(r) / (1 - r ^ 2)
# 
# FisherZ = 0.5 * Log((1 + 0,832) / (1 - 0,832)) = 1,195
# StdErr(FisherZ) = (0,013) / (1 - 0,832 ^ 2) = 0,041
