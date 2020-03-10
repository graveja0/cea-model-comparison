# Load the following packages and scripts to run this model independently
# library(here)
# library(Matrix)
# source(here("common.R"))
# source(here("simple-params.R"))

# NOTE: Due to nature of temporary disutility in our case, 
# different cycle lengths require different construction of tunnel accumulators (non-Markovian parts).
# It might be feasible to build a more flexible version of the model, 
# but please note this script only works for a yearly cycle.
params$interval <- 1

#### Embed Transition Probability Matrices ####
markov_corr_tp <- function(params, v.n, t)
{
  # To Debug
  # sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
  with(params, {
    
    #secular death risk -- All same as original Markov computation
    r.death   <- inst_rate(gompertz_ratio2(t, interval, shape, rate), 1)
    rr        <- ifelse(p_g==1 & p_o==1, rr_b, 1)
    r_a <- r_a/interval # Rates are in terms of interval units
    r_b <- r_b/interval # Rates are in terms of interval units

    # Always start on zeros, continuous transition rates, i.e. tr
    tr <-  matrix(0, nrow = length(v.n), ncol = length(v.n), dimnames = list(v.n, v.n))

    #########################################
    ## transition _RATES_ into matrix form
    
    # Health to Indication or secular death
    tr["H", "A"]  <- r_a
    tr["H", "D"]  <- r.death
      
    # Indication to adverse event or secular death
    tr["A", "BS"] <- rr*r_b*(1-p_bd)
    tr["A", "BD"] <- rr*r_b*p_bd
    tr["A", "D"]  <- r.death
      
    # Post adverse event to secular death
    tr["BS","BSD"]  <- r.death
    
    # Diagnolize to balance Markovian transitions
    # Makes sure no one is created or destroyed
    diag(tr)      <- -rowSums(tr)        # Markovian rate specification for rates
    
    ###########################################
    ## Non-Markovian Accumulators (as rates) Added to matrix
    ## I.e., do not include in inverse diagonal
    
    # Add accumulators (Non-Markovian)
    tr["H", "CUM_A"]   <- tr["H", "A"]        # Accumulator for A is non-Markovian
    tr["H", "CUM_T"]   <- p_o*tr["H", "A"]    # Accumulator for testing is non-Markovian
    
    # Tunnel accumulator for A
    tr["H", "TUN"]     <- tr["H", "A"]
    # TUN Now duplicates transitions in from H, just like A
    # PTUN is a "waste" bucket necessary for bookkeeping
    # Transitions into PTUN keep track of things that exited tunnel state
    # for reasons other than expiration of the tunnel, i.e. external risks
    tr["TUN", "PTUN"]  <- tr["A", "BD"] + tr["A", "D"] + tr["A", "BS"]
    tr["TUN", "TUN"]   <- -sum(tr["A", c("BD", "D", "BS")]) # semi-Markovian
    
    ###########################################
    ## Embed into unit (interval) timestep the matrix
    ## x is now probabilities
    x <- as.matrix(expm(tr))

    ###########################################
    # Adding skip overs back to probabilities for the accumulators
    # Attempted modification for accurate discounting integration
    x["A", "CUM_BS"]    <- x["A", "BSD"]+x["A", "BS"]
    x["H", "CUM_BS"]    <- x["H", "BSD"]+x["H", "BS"]

    # Leave tunnel state on yearly cycle (zero is null transfer back)
    x["TUN", "TUN"]    <- 0 # This patches the semi-Markovian 
    # PROBLEM: This doesn't work if interval is not yearly
    # E.G. one needs a TUN1 through TUN12 for monthly, etc to work on other intervals
    
    # The alpha / beta was handled by the matrix exponent

    x
  })
}

#### Main Simulation ####
# Unweighted model: Markov trace  
markov_corr_sim <- function(params)
{
  # To Debug
  # sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
  with(params, {
        
    #### 1. Model Setting ####
    n.t       <- horizon*interval          # number of cycles
    v.n       <- c("H", "A", "BS", "BD", "D", "BSD", "CUM_A", "CUM_BS", "CUM_T", "TUN", "PTUN")  # state names, and accumulators
    n.s       <- length(v.n)               # number of states
        
    #### 2. List of Transition Probability Matrices #### 
    l.P <- lapply(1:n.t, function(x) markov_corr_tp(params, v.n, x))
        
    #### 3. Initiate Markov Trace ####
    m.M <- matrix(0, 
                  nrow = n.t + 1 ,                # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
                  ncol = n.s,                  
                  dimnames = list(0:n.t, v.n)) 

    # The cohort starts from state 1, aka 'Healthy'
    m.M[1,] <- c(1, rep(0, n.s-1))      # initialize first cycle of Markov trace accounting for the tunnels
    #### 4. Run the Model ####
    for (t in 1:n.t) m.M[t+1, ] <- m.M[t, ] %*% l.P[[t]]
    
    m.M
  })
}

# Discounting & Summation
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
    # discount rate transformation
    d.r  <- inst_rate(1-1/(1 + params$disc), interval)
    v.dw <- (exp(-d.r*(0:(n.t-2))) - exp(-d.r*(1:(n.t-1))))/d.r
    
    # adjust counts using integration methods
    mm        <- integrator(m.M,                method=method) # Total counts
    dmm       <- integrator(diag(exp( 0:(n.t-1) * -d.r)) %*% m.M, method=method) # Discounted counts
    
    # drug/testing costs conditional on gene/testing decision
    dd        <- ifelse(gene==1 & test==1, c_alt*365/interval, c_tx*365/interval) #conditional drug cost
    tt        <- ifelse(test==1, c_t,0) #conditional testing cost
    
    # Survival Metrics
    living    <- sum(m.M[n.t, c("H", "A", "BS")])         # Proportion alive at end of sim
    fatal_b   <- sum(m.M[n.t, c("BD")])
    possible  <- sum(dmm[,c("H", "A", "BS")]) / interval  # Discounted Life Years

    # Discounted Costs    
    c.test    <- tt   * sum(disc_acc(m.M, "CUM_T", v.dw, method))  #discounted testing cost
    c.drug    <- dd   * sum(dmm[,c("A", "BS")]) #discounted drug cost
    c.treat   <- c_a  * disc_acc(m.M, "CUM_A",  v.dw, method) + 
                 c_bs * disc_acc(m.M, "CUM_BS", v.dw, method) + 
                 c_bd * disc_acc(m.M, "BD",     v.dw, method) #discounted adverse event costs

    # Discounted Disutilities
    disutil_a <- d_a*sum(dmm[,"TUN"]) / interval #discounted disutility from A
    disutil_b <- d_b*sum(dmm[,c("BS")]) / interval #discounted disutility from B

    list(
      m.M     = m.M,
      mm      = mm,
      dmm     = dmm,
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
# Under testing strategy, cost/QALY implications are different conditional on gene type/testing decision
# 1.no testing because of physician not ordering test (p_o); 2. no-variant gets tested not prescribed alternative drug;
# 3.variant-carrier gets tested and prescribed alternative drug;
markov_corr_icer <- function(params, method="life-table")
{
  ref1        <- markov_corr(params, gene=0, test=0, method=method)
  test1       <- markov_corr(params, gene=1, test=1, method=method)
  test2       <- markov_corr(params, gene=0, test=1, method=method)
     
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

# markov_corr_icer(params) #sample code to run the model
