
library(here)
library(Matrix)
source(here("R/common.R"))

markov_corr_tp <- function(params, v.n, t)
{
  # To Debug
  # sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
  with(params, {
    
    #### 2. Inputs ####
    #secular death risk -- All same as original Markov computation
    r.death   <- inst_rate(gompertz_ratio2(t, interval, shape, rate), 1)
    rr        <- ifelse(p_g==1 & p_o==1, rr_b, 1)

    # Always start on zeros
    tr <-  matrix(0, nrow = length(v.n), ncol = length(v.n), dimnames = list(v.n, v.n))

    # Health to Indication or secular death
    tr["H", "A"]  <- r_a
    tr["H", "D"]  <- r.death
      
    # Indication to adverse event or secular death
    tr["A", "BS"] <- rr*r_b*(1-p_bd)
    tr["A", "BD"] <- rr*r_b*p_bd
    tr["A", "D"]  <- r.death
      
    # Post adverse event to secular death
    tr["BS","BSD"]  <- r.death
    
    diag(tr)      <- -rowSums(tr)        # Markovian rate specification
      
    # Add accumulators
    tr["H", "CUM_A"]   <- tr["H", "A"]        # Accumulator for A is non-Markovian
    #tr["A", "CUM_BS"]  <- rr*r_b*(1-p_bd)     # Accumulator for B is non-Markovian
    tr["H", "CUM_T"]   <- p_o*tr["H", "A"]    # Accumulator for testing is non-Markovian
    
    # Tunnel accumulator for A
    tr["H", "TUN"]     <- tr["H", "A"]
    tr["TUN", "PTUN"]  <- tr["A", "BS"]
    tr["TUN", "PTUN"]  <- tr["A", "BD"]
    tr["TUN", "PTUN"]  <- tr["A", "D"]
    tr["TUN", "TUN"]   <- -sum(tr["A", c("BS", "BD", "D")])
    
    x <- as.matrix(expm(tr))

    # Attempted modification for accurate discounting integration
    eff_disc             <- (exp(-disc_rate*(t-1)) - exp(-disc_rate*t))/disc_rate
    x["A", "CUM_BS"]    <- x["A", "BSD"]+x["A", "BS"]
    x["H", "CUM_BS"]    <- x["H", "BSD"]+x["H", "BS"]
    x["A", "DCUM_BS"]   <- eff_disc*x["A", "CUM_BS"]
    x["H", "DCUM_BS"]   <- eff_disc*x["H", "CUM_BS"]
    
    # Leave tunnel state on yearly cycle (zero is null transfer back)
    x["TUN", "TUN"]    <- 0
    
    # The alpha / beta was handled by the matrix exponent

    x
  })
}

markov_corr_sim <- function(params)
{
  # To Debug
  # sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
  with(params, {
        
    #### 1. Model Setting ####
    n.t       <- horizon*interval          # number of cycles
    v.n       <- c("H", "A", "BS", "BD", "D", "BSD", "CUM_A", "CUM_BS", "CUM_T", "DCUM_BS", "TUN", "PTUN")  # state names, and accumulators
    n.s       <- length(v.n)               # number of states
        
    # transition probability matrix list
    l.P <- lapply(1:n.t, function(x) markov_corr_tp(params, v.n, x))
        
    #### 4. Run Probabilistic Matrix ####
    m.M <- matrix(0, 
                  nrow = n.t + 1 ,                # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
                  ncol = n.s,                  
                  dimnames = list(0:n.t, v.n)) 

    # The cohort starts from state 1, aka 'Healthy'
    m.M[1,] <- c(1, rep(0, n.s-1))      # initialize first cycle of Markov trace accounting for the tunnels
    for (t in 1:n.t) m.M[t+1, ] <- m.M[t, ] %*% l.P[[t]]
    
    m.M
  })
}

#### 03 Main Simulation
# Unweighted model
markov_corr <- function(params, N=NULL, gene=0, test=0, method="beginning")
{
  if (!is.null(N)) params$n <- N
  
  params$disc_rate <- inst_rate(1-1/(1 + params$disc), 1)
  params$p_g <- gene
  params$p_o <- test
  
  m.M <- markov_corr_sim(params)
  
  # To Debug
  # sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
  with(params, {

    n.t <- dim(m.M)[1]
    #### 5. Computation ####
    d.r  <- inst_rate(1-1/(1 + params$disc), 1)
    v.dw <- (exp(-d.r*(0:(n.t-2))) - exp(-d.r*(1:(n.t-1))))/d.r
    #v.dw <- exp( 0:(n.t-1) * -d.r)

    # adjust counts using integration methods
    mm        <- integrator(m.M,                method=method) # Total counts
    dmm       <- integrator(diag(exp( 0:(n.t-1) * -d.r)) %*% m.M, method=method) # Discounted counts

    # conditional drug costs
    dd        <- ifelse(gene==1 & test==1, c_alt*365/interval, c_tx*365/interval)
    tt        <- ifelse(test==1, c_t,0)

    living    <- sum(m.M[n.t, c("H", "A", "BS")])         # Proportion alive at end of sim
    fatal_b   <- sum(m.M[n.t, c("BD")])
    possible  <- sum(dmm[,c("H", "A", "BS")]) / interval  # Sum of possible dQALY

    # Discounted Costs    
    c.test    <- tt   * sum(disc_acc(m.M, "CUM_T", v.dw, method)) 
    c.drug    <- dd   * sum(dmm[,c("A", "BS")]) 
    c.treat   <- c_a  * disc_acc(m.M, "CUM_A",  v.dw, method) +
                 c_bs * disc_acc(m.M, "CUM_BS", v.dw, method) + 
                 c_bd * disc_acc(m.M, "BD",     v.dw, method)

    
    disutil_a <- d_a*sum(dmm[,"TUN"])

    # This is correct
    disutil_b <- d_b*sum(dmm[,c("BS")]) / interval 

    list(
      m.M     = m.M,
      mm      = mm,
      results = c(dCOST       = c.test+c.drug+c.treat,
                  dQALY       = possible - disutil_a - disutil_b,
                  possible    = possible,
                  fatal_b     = fatal_b,
                  living      = living,
                  disutil_a   = disutil_a,
                  disutil_b   = disutil_b,
                  dCOST.test  = c.test,
                  dCOST.drug  = c.drug,
                  dCOST.treat = c.treat,
                  dCOST.a     = c_a  * disc_acc(m.M, "CUM_A",  v.dw, method),
                  dCOST.bs    = c_bs * disc_acc(m.M, "CUM_BS", v.dw, method),
                  dCOST.bd    = c_bd * disc_acc(m.M, "BD",     v.dw, method))
    )
  }) # Closes "with(params, {
}

# Combined model
markov_comb_corr <- function(params, method="life-table")
{
  ref1        <- markov_corr(params, gene=0, test=0, method=method)
  test1       <- markov_corr(params, gene=1, test=1, method=method)
  test2       <- markov_corr(params, gene=0, test=1, method=method)
     
  o           <- params$p_o
  g           <- params$p_g
  wt          <- c(1-o,g*o,(1-g)*o)
  dCOST.ref   <- ref1$results["dCOST"]
  dQALY.ref   <- ref1$results["dQALY"]
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


