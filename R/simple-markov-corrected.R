
library(here)
library(Matrix)
source(here("R/common.R"))

markov_corr_sim <- function(params)
{
  # To Debug
  # sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
  with(params, {
        
    #### 1. Model Setting ####
    n.t       <- horizon*interval          # number of cycles
    v.n       <- c("H", "A", "BS", "BD", "D", "CUM_A", "CUM_B", "CUM_T")  # state names, and accumulators
    n.s       <- length(v.n)               # number of states
        
    #### 2. Inputs ####
    #secular death risk -- All same as original Markov computation
    r.HD      <- inst_rate(gompertz_ratio2(1:n.t, interval, shape, rate), 1)
    rr        <- ifelse(p_g==1 & p_o==1, rr_b, 1)
   
        
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
    l.P <- lapply(1:n.t, function(x)
    {
      # Always start on zeros
      tr <-  matrix(0,nrow = n.s, ncol = n.s, dimnames = list(v.n, v.n))

      # Health to Indication or secular death
      tr["H", "A"]  <- r_a
      tr["H", "D"]  <- r.HD[x]
      
      # Indication to adverse event or secular death
      tr["A", "BS"] <- rr*r_b*(1-p_bd)
      tr["A", "BD"] <- rr*r_b*p_bd
      tr["A", "D"]  <- r.HD[x]
      
      # Post adverse event to secular death
      tr["BS","D"]  <- r.HD[x]
      
      diag(tr) <- -rowSums(tr)        # Markovian rate specification
      
      # Add accumulators
      tr["H", "CUM_A"] <- tr["H", "A"]      # Accumulator for A is non-Markovian
      tr["A", "CUM_B"] <- rr*(1-p_bd)*r_b   # Accumulator for B is non-Markovian
      tr["H", "CUM_T"] <- p_o*tr["H", "A"]  # Accumulator for testing is non-Markovian
      
      as.matrix(expm(tr))
    })
        
    #### 4. Run Probabilistic Matrix ####
    for (t in 1:n.t) m.M[t+1, ] <- m.M[t, ] %*% l.P[[t]]
    
    m.M
  })
}

#### 03 Main Simulation
# Unweighted model
markov_corr <- function(params, N=NULL, gene=0, test=0, method="beginning")
{
  if (!is.null(N)) params$n <- N
  
  params$p_g <- gene
  params$p_o <- test
  
  m.M <- markov_corr_sim(params)
  
  # To Debug
  # sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
  with(params, {
        
    n.t <- dim(m.M)[1]
    #### 5. Computation ####
    d.r  <- (1 + disc)^(1/interval)-1
    v.dw <- 1 / (1 + d.r) ^ (0:(n.t-1)) # calculate discount weights for costs for each cycle based on discount rate d.r
    
    # adjust counts using integration methods
    mm        <- integrator(m.M,                method=method) # Total counts
    dmm       <- integrator(diag(v.dw) %*% m.M, method=method) # Discounted counts

    # conditional drug costs
    dd        <- ifelse(gene==1 & test==1, c_alt*365/interval,c_tx*365/interval)
    tt        <- ifelse(test==1, c_t,0)

    living    <- sum(m.M[n.t, c("H", "A", "BS")])         # Proportion alive at end of sim
    fatal_b   <- sum(m.M[n.t, "BD"])
    possible  <- sum(dmm[,c("H", "A", "BS")]) / interval  # Sum of possible dQALY

    # Discounted Costs    
    c.test    <- tt   * sum(disc_acc(m.M, "CUM_T", v.dw, method)) 
    c.drug    <- dd*sum(dmm[,c("A", "BS")]) / interval
    c.treat   <- c_a  * sum(disc_acc(m.M, "CUM_A", v.dw, method)) +
                 c_bs * sum(disc_acc(m.M, "CUM_B", v.dw, method)) + 
                 c_bd * sum(disc_acc(m.M, "BD",    v.dw, method))

    # FIXME: Temp hack, not correct
    disutil_a <- d_a*sum(disc_acc(m.M, "CUM_A", v.dw, method))

    # This is correct
    disutil_b <- d_b*sum(dmm[,c("BS")]) / interval 

    list(
      m.M     = m.M,
      l.P     = l.P,
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
                  dCOST.treat = c.treat)
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


