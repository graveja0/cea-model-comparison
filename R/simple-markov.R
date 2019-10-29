library(here)
source(here("R/common.R"))

#### 02 Functions ####
gompertz_ratio2 <- function(t, interval, shape, rate) # p.sd for each step
{
  t1 <- t/interval
  t0 <- (t-1)/interval
  r  <- (pgompertz(t1, shape, rate) - pgompertz(t0, shape, rate)) /
    (1 - pgompertz(t0, shape, rate))
  r[is.na(r)] <- 1  
  
  r
}

# Once secular death mortablity reaches 1, all other probs need to put 0, no effect for default 40-year time horizon
# NOTE: This is used a lot in cost effectiveness research and is the result
# of failing to treat it as competing risks when embedding a continuous rate into a probability
cap_max <- function(value, sd) ifelse(value+sd>1, 1-sd, value)


#### 03 Main Simulation
# Unweighted model
markov0 <- function(params, N=NULL, gene=0, test=0, method="beginning")
{
  if (!is.null(N)) params$n <- N
  
  with(params, {
    
    #### 1. Model Setting ####
    n.t       <- horizon*interval           # number of cycles
    v.n       <- c("H",paste0("A",1:(d_at*interval+1)),"BS1","BS2","BD1","BD2","D")    # state names
    n.s       <- length(v.n)          # number of states
    
    #### 2. Inputs ####
    #secular death risk
    p.HD      <- gompertz_ratio2(1:n.t, interval, shape, rate)
    pA        <- cap_max(value=rate_to_prob(r_a,1/interval), sd=p.HD)
    
    rr        <- ifelse(gene==1 & test==1, rr_b, 1)
    pB        <- cap_max(value=rate_to_prob(r_b,1/interval), sd=p.HD)
    pBS       <- rr*pB*(1-p_bd)
    pBD       <- rr*pB*p_bd
    
    #### 3. Transition Matrix ####
    # cohort counts
    m.M <- matrix(0, 
                  nrow = n.t + 1 ,                # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
                  ncol = n.s,                  
                  dimnames = list(0:n.t, v.n)) 
    
    ####
    # The cohort starts from state 1
    m.M[1,] <- c(1, rep(0, n.s-1))      # initialize first cycle of Markov trace accounting for the tunnels
    ####
    
    # transition probability array
    tp <-  matrix(0,nrow = n.s, ncol = n.s, dimnames = list(v.n, v.n))
    l.P <- lapply(1:n.t, function(x)
    {
      tp["H", "A1"] <- pA[x]
      tp["H", "D"]  <- p.HD[x]
      tp["H", "H"]  <- 1-pA[x]-p.HD[x]
      
      for(y in 1:(d_at*interval+1))
      {
        tp[paste0("A",y),"BS1"] <- pBS[x]  
        tp[paste0("A",y),"BD1"] <- pBD[x]
        tp[paste0("A",y),"D"]   <- p.HD[x]
      }
      for(y in 1:(d_at*interval))
      {
        tp[paste0("A",y), paste0("A",y+1)] <- 1-pBS[x]-pBD[x]-p.HD[x]
      }
      
      tp[paste0("A", d_at*interval+1), paste0("A", d_at*interval+1)] <- 1-pBS[x]-pBD[x]-p.HD[x]
      
      tp["BS1","D"]   <- tp["BS2","D"]   <- p.HD[x]
      tp["BS1","BS2"] <- tp["BS2","BS2"] <- 1-p.HD[x]
      tp["BD1","BD2"] <- tp["BD2","BD2"] <- 1
      tp["D","D"]     <- 1
      
      return(tp)
    })
    
    #### 4. Run ####
    for (t in 1:n.t)                        # loop through the number of cycles
      m.M[t+1, ] <- m.M[t, ] %*% l.P[[t]] # estimate the Markov trace for cycle t + 1 using the t-th matrix from the probability array
    
    #### 5. Computation ####
    d.r   <- (1 + disc)^(1/interval)-1
    v.dwc <- 1 / (1 + d.r) ^ (0:(n.t-1)) # calculate discount weights for costs for each cycle based on discount rate d.r
    v.dwe <- 1 / (1 + d.r) ^ (0:(n.t-1))
    
    #adjust counts based on method
    mm        <- integrator(m.M, method=method)
    dd        <- ifelse(gene==1 & test==1, c_alt*365/interval,c_tx*365/interval) #conditional drug cost
    tt        <- ifelse(test==1, c_t,0)
    
    cc        <- mm %*% c(0,c_a+dd+tt,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0)
    ee        <- mm %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)
    cc        <- as.vector(cc) * v.dwc
    ee        <- as.vector(ee) * v.dwe
    c.test    <- mm[,"A1"]*tt*v.dwc
    c.drug    <- mm %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0)
    c.drug    <- as.vector(c.drug) * v.dwc
    
    possible  <- mm %*% c(1,rep(1,d_at*interval+1),1,1,0,0,0)
    possible  <- as.vector(possible) * v.dwe
    fatal_b   <- sum(mm[n.t,c("BD1","BD2")])
    living    <- 1 - sum(mm[n.t,c("BD1","BD2","D")])
    disutil_a <- mm %*% c(0,rep(d_a,d_at*interval),0,0,0,0,0,0)
    disutil_a <- as.vector(disutil_a) * v.dwe
    disutil_b <- mm %*% c(0,rep(0,d_at*interval+1),d_b,d_b,0,0,0)
    disutil_b <- as.vector(disutil_b) * v.dwe
    c.treat   <- mm %*% c(0,c_a,rep(0,d_at*interval),c_bs,0,c_bd,0,0)
    c.treat   <-  as.vector(c.treat) * v.dwc
    
    list(
      m.M     = m.M,
      l.P     = l.P,
      mm      = mm,
      results = c(dCOST       = sum(unlist(cc)),
                  dQALY       = sum(unlist(ee)),
                  possible    = sum(unlist(possible)),
                  fatal_b     = unname(fatal_b),
                  living      = unname(living),
                  disutil_a   = sum(unlist(disutil_a)),
                  disutil_b   = sum(unlist(disutil_b)),
                  dCOST.test  = sum(unlist(c.test)),
                  dCOST.drug  = sum(unlist(c.drug)),
                  dCOST.treat = sum(unlist(c.treat)))
    )
  }) # Closes "with(params, {
}

# Combined model
markov_comb <- function(params, method="life-table")
{
  ref1        <- markov0(params, gene=0, test=0, method=method)
  test1       <- markov0(params, gene=1, test=1, method=method)
  test2       <- markov0(params, gene=0, test=1, method=method)
  
  o           <- params$p_o
  g           <- params$p_g
  wt          <- c(1-o,g*o,(1-g)*o)
  dCOST.ref   <- unname(ref1$results["dCOST"])
  dQALY.ref   <- unname(ref1$results["dQALY"])
  dCOST.test  <- sum(wt*c(ref1$results["dCOST"], test1$results["dCOST"], test2$results["dCOST"]))
  dQALY.test  <- sum(wt*c(ref1$results["dQALY"], test1$results["dQALY"], test2$results["dQALY"]))
  
  c(ICER       = (dCOST.test-dCOST.ref)/(dQALY.test-dQALY.ref),
    NMB        = (dCOST.ref-dCOST.test)+params$wtp*(dQALY.test-dQALY.ref),
    dCOST.ref  = unname(dCOST.ref),
    dCOST.test = unname(dCOST.test),
    dQALY.ref  = unname(dQALY.ref),
    dQALY.test = unname(dQALY.test)
  )
}
