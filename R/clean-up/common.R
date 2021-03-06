# Functions shared in different models
library(flexsurv)

inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

# Simpson's method of integration
# Far better performance that "life-table" aka trapezoidal method
alt_simp_coef <- function(i) c(17, 59, 43, 49, rep(48, i-8), 49, 43, 59, 17) / 48
alt_simp      <- function(x,h) h*sum(alt_simp_coef(length(x)) * x)


# secular death rate over time interval
gompertz_ratio <- function(t0, t1, shape, rate)
{
  r <- (pgompertz(t1, shape, rate) - pgompertz(t0, shape, rate)) / (1 - pgompertz(t0, shape, rate))
  r[is.na(r)] <- 1
  r
}

# Secular death from t-1 to t (scaled by interval)
# A modified version of the gompertz_ratio function to return probability of secular death in each time cycle accommodating different cycle lengths.   
# Default interval value 1 means 1-year cycle length, 12 means monthly, 365 means daily. 
gompertz_ratio2 <- function(t, interval, shape, rate) # p.sd for each step
{
  t1 <- t/interval
  t0 <- (t-1)/interval
  r  <- (pgompertz(t1, shape, rate) - pgompertz(t0, shape, rate)) /
        (1 - pgompertz(t0, shape, rate))
  r[is.na(r)] <- 1  
  
  r
}

# Sample states for multiple individuals simultaneously, this function was developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Citations:
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 
samplev <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}


# from "heemod" package code
# Citation: Filipović-Pierucci A, Zarca K, Durand-Zaleski I (2017). “Markov Models for Health Economic Evaluation: The R Package heemod.” ArXiv e-prints. R package version 0.11.0, 1702.03252.
rate_to_prob <- function(r, to=1, per=1)
{
  1 - exp(-(r/per) * to)
}

# Once secular death mortablity reaches 1, all other probs need to put 0, no effect for default 40-year time horizon
# NOTE: This is used a lot in cost effectiveness research and is the result
# of failing to treat it as competing risks when embedding a continuous rate into a probability
cap_max <- function(value, sd) ifelse(value+sd>1, 1-sd, value)

#trying to replicate life-table method, and other numerical integration methods
integrator <- function(x, method="beginning")
{
  if(method=="life-table")
  {
    n0 <- x[- nrow(x), ]
    n1 <- x[-1, ]
    (n0 + n1) / 2
  } else if(method=="end")
  {
    x[-1, ]
  } else if(method=="beginning")
  {
    x[- nrow(x), ]
  } else if(method=="simpson")
  {
    if(length(x) < 8) stop("Too few time steps for Alternate Extended Simpson method")
    apply(x, 2, function(y) alt_simp_coef(length(y)) * y)
  } else stop("undefined method") 
}

# Run a discounted integrator on differences in an accululator
disc_acc <- function(m, cols, disc, method)
{
  sum(diag(disc) %*% diff(m[,cols, drop=FALSE]))
}