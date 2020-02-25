# This code built on a model example developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Citations:
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 

# Load the following packages and scripts to run this model independently
# library(here)
# source(here("common.R")) #load shared functions
# source(here("simple-params.R")) #load inputs

#### Main Simulation ####
# Unweighted model
markov0 <- function(params, N=NULL, gene=0, test=0, method="beginning")
{
  if (!is.null(N)) params$n <- N
  
  with(params, {
    
    #### 1. Model Setting ####
    n.t       <- horizon*interval           # number of cycles
    v.n       <- c("H",paste0("A",1:(d_at*interval+1)),"BS1","BS2","BD1","BD2","D")    # state names
    n.s       <- length(v.n)          # number of states
    
    #### 2. Probabilities ####
    p.HD      <- gompertz_ratio2(1:n.t, interval, shape, rate)
    pA        <- cap_max(value=rate_to_prob(r_a,1/interval), sd=p.HD)
    
    rr        <- ifelse(gene==1 & test==1, rr_b, 1)
    pB        <- cap_max(value=rate_to_prob(r_b,1/interval), sd=p.HD)
    pBS       <- rr*pB*(1-p_bd)
    pBD       <- rr*pB*p_bd
    
    #### 3. Markov Trace ####
    m.M <- matrix(0, 
                  nrow = n.t + 1 ,       
                  ncol = n.s,                  
                  dimnames = list(0:n.t, v.n)) 
    
    ####
    # The cohort starts from everyone in the healthy state 
    m.M[1,] <- c(1, rep(0, n.s-1))
    ####
    
    #### 3. Transition Probability Matrices #### 
    tp <-  matrix(0,nrow = n.s, ncol = n.s, dimnames = list(v.n, v.n))
    
    # fill transition probability matrix for each cycle to avoid storing all matrices in one giant object
    get_transition_prob <- function(x) 
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
    }
    
    #### 4. Run ####
    for (t in 1:n.t)   {
      # loop through the number of cycles
      l.P.tmp <- get_transition_prob(t)
      m.M[t+1, ] <- m.M[t, ] %*% l.P.tmp # estimate the Markov trace for cycle t + 1 using the t-th transition probability matrix
    }
    
    #### 5. Discounting & Summation ####
    d.r  <- inst_rate(1-1/(1 + disc), interval) # annual discount rate adjusted to cycle length and converted to instaneous rate
    
    #adjust counts based on integration method (need to dicount before this step if applicable)
    mm        <- integrator(m.M, method=method)
    dmm       <- integrator(diag(exp( 0:n.t * -d.r)) %*% m.M, method=method)
    
    # drug/testing costs conditional on gene/testing decision
    dd        <- ifelse(gene==1 & test==1, c_alt*365/interval,c_tx*365/interval) #conditional drug cost
    tt        <- ifelse(test==1, c_t,0) #conditional testing cost
    
    # discounted cost
    cc        <- as.vector(dmm %*% c(0,c_a+dd+tt,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0)) #discounted total cost by cycle
    c.test    <- as.vector(dmm[,"A1"]*tt) #discounted testing cost
    c.drug    <- as.vector(dmm %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0)) #discounted drug cost
    c.treat   <- as.vector(dmm %*% c(0,c_a,rep(0,d_at*interval),c_bs,0,c_bd,0,0)) #discounted adverse event costs
    
    # discounted effectiveness
    ee        <- as.vector(dmm%*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)) #discounted QALYs by cycle
    possible  <- as.vector(dmm %*% c(1,rep(1,d_at*interval+1),1,1,0,0,0))/interval #discounted life years 
    disutil_a <- as.vector(dmm %*% c(0,rep(d_a,d_at*interval),0,0,0,0,0,0))/interval #discounted disutility from A
    disutil_b <- as.vector(dmm %*% c(0,rep(0,d_at*interval+1),d_b,d_b,0,0,0))/interval #discounted disutility from B
    
    # other metrics
    living    <- 1 - sum(mm[n.t,c("BD1","BD2","D")])/interval #proportion alive at end of sim
    fatal_b   <- sum(mm[n.t,c("BD1","BD2")]) #fatal B count
    
    list(
      # m.M     = m.M,
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
  }) 
}

# Combined model
# Under testing strategy, cost/QALY implications are different conditional on gene type/testing decision
# 1.no testing because of physician not ordering test (p_o); 2. no-variant gets tested not prescribed alternative drug;
# 3.variant-carrier gets tested and prescribed alternative drug;
markov_icer <- function(params, method="life-table")
{
  ref1        <- markov0(params, gene=0, test=0, method=method)
  test1       <- markov0(params, gene=1, test=1, method=method)
  test2       <- markov0(params, gene=0, test=1, method=method)
  
  o           <- params$p_o
  g           <- params$p_g
  wt          <- c(1-o,g*o,(1-g)*o) #Under testing strategy, results are weighted by gene prevalence and physician behavior 
  dCOST.ref   <- unname(ref1$results["dCOST"])
  dQALY.ref   <- unname(ref1$results["dQALY"])
  dCOST.test  <- sum(wt*c(ref1$results["dCOST"], test1$results["dCOST"], test2$results["dCOST"]))
  dQALY.test  <- sum(wt*c(ref1$results["dQALY"], test1$results["dQALY"], test2$results["dQALY"]))
  
  c(ICER       = (dCOST.test-dCOST.ref)/(dQALY.test-dQALY.ref),
    NMB        = params$wtp*(dQALY.test-dQALY.ref)-(dCOST.test-dCOST.ref),
    dCOST.ref  = unname(dCOST.ref),
    dCOST.test = unname(dCOST.test),
    dQALY.ref  = unname(dQALY.ref),
    dQALY.test = unname(dQALY.test)
  )
}

markov_icer(params) #sample code to run the model
