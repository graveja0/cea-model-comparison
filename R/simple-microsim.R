library(data.table)
library(tidyverse)
library(flexsurv)
library(heemod)

#### 01 inputs ####
source("./R/simple-params.R")

#### 02 Sim Functions ####
#---------------------------------------------------------------------------#
#### R function to sample states for multiple individuals simultaneously ####
#---------------------------------------------------------------------------#
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
### from workshop code ###


gompertz_ratio2 <- function(t, interval, shape, rate) # p.sd for each step
{
  t1 <- t/interval
  t0 <- (t-1)/interval
  r <- (pgompertz(t1, shape, rate) - pgompertz(t0, shape, rate)) / (1 - pgompertz(t0, shape, rate))
  if(is.na(r)) r=1
  return(r)
}

# Event draw function:
# looks up probability of transition
ProbsH <- function(M_t,t,p.a,p.sd,interval){
  # M_t    :    current health state 
  # t    :  cycle
  p.die.base <- p.sd[t]
  p <- matrix(0,nrow=length(M_t),ncol=2,dimnames=list(c(),c("A1","D")))
  pp.a <- rate_to_prob(p.a,1/interval)
  
  #cap probs to match markov version: is it appropriate????
  #given time horizon of 40, no impact
  if(p.die.base+pp.a>1) {pp.a <- 1-p.die.base}
  
  p[,"A1"] <- pp.a
  p[,"D"] <- p.die.base
  return(p)
}

ProbsA <- function(M_t,v.drug,t,p.b,rr.b,p.bd,p.sd,interval) {
  # v.drug :    treatment
  p.die.base <- p.sd[t]
  p <- matrix(0,nrow=length(M_t),ncol=3,dimnames=list(c(),c("BS1","BD1","D")))
  
  pp.b <- rate_to_prob(p.b,1/interval)
  if(p.die.base+pp.b>1) {pp.b <- 1-p.die.base}
  
  p[,"BS1"] <- ifelse(v.drug=="Alternate", pp.b*rr.b*(1-p.bd), pp.b*(1-p.bd))
  p[,"BD1"] <- ifelse(v.drug=="Alternate", pp.b*rr.b*p.bd, pp.b*p.bd)
  p[,"D"] <- p.die.base
  return(p)
}


ProbsBS <- function(M_t,t,p.sd){
  rep(p.sd[t],length(M_t))
}


#draw event for competing risks
draw.events <- function(current,p) {
  pp <- rbinom(n=length(current),size=1,prob=rowSums(p))
  if(sum(pp)>0) { # if at least one person has a event
    # randomly sample the event type for those who experience progression events this cycle
    # uses samplev() from the "Functions Folder"
    event.type <- samplev(p[pp==1,,drop=FALSE]/rowSums(p[pp==1,,drop=FALSE]),m=1)   
    # update state based on event type
    current[pp==1] = event.type
  }
  return(unname(current))
}

#trying to replicate life-table method
ltm <- function(x,method="beginning") {
  if(method=="life-table") {
    n0 <- x[- nrow(x), ]
    n1 <- x[-1, ]
    (n0 + n1) / 2
  } else {
    x[-1, ]
  }
}

microsim_run <- function(params, N = NULL, method="beginning") {
  if (!is.null(N)) params$n <- N
  with(params, {
    
    #### 1. Model Setting ####
    n.i       <- n              # number of individuals
    n.t       <- horizon*interval           # number of cycles
    v.n       <- c("H",paste0("A",1:(d_at*interval+1)),"BS1","BS2","BD1","BD2","D")    # state names
    n.s       <- length(v.n)          # number of states
    
    #secular death risk
    p.HD <- 1:n.t %>% map_dbl(~gompertz_ratio2(.x,interval,shape,rate))
    
    #### 2. Sample individual level characteristics ####
    #Static characteristics
    df.Pop <-  data.frame(name=1:n.i,variant=sample(c(FALSE,TRUE),n.i,prob=c(1-p_g,p_g),replace = TRUE)) %>% 
      mutate(tested = FALSE,treat = "")
    
    #Dynamic characteristics 
    v.M_Init <- rep("H", n.i)       # everyone begins in the healthy state 
    
    #### 3. Start Simultation ####
    # m.M: health state for each patient at each cycle
    m.M  <- matrix(nrow = n.i, ncol = n.t + 1, 
                   dimnames = list(paste(1:n.i),     # could name the rows ind1, ind2, ind3, etc.
                                   paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.
    
    m.M[, 1] <- v.M_Init         # initial health state for all individuals
    
    # Time for-loop from cycle 1 through n.t
    for (t in 1:n.t) {
      
      ## Simulate Events ##
      # Default Update: if no events occur, health state remains unchanged
      m.M[,t+1] = m.M[,t]
      
      # From H
      fromH <- (m.M[,t]=="H")
      if(sum(fromH)>0) {
        pH <- ProbsH(m.M[,t][fromH],t,r_a,p.HD,interval)
        m.M[,t+1][fromH] <- draw.events(m.M[,t][fromH],pH)
        
        # assign treatment upon indication
        elig <- (m.M[,t]=="H" & m.M[,t+1]=="A1")
        df.Pop$tested[elig] <- sample(c(TRUE,FALSE),sum(elig),prob=c(p_o,1-p_o),replace = TRUE)
        df.Pop$treat[elig] <- ifelse(df.Pop$tested[elig] & df.Pop$variant[elig],"Alternate","Primary")
      }
      
      # From A
      fromA <- m.M[,t] %in% paste0("A",1:(d_at*interval+1))
      
      if(sum(fromA)>0) {
        pA <- ProbsA(m.M[,t][fromA],df.Pop$treat[fromA],t,r_b,rr_b,p_bd,p.HD,interval)
        
        # move forward in A tunnel as needed
        drawA <- draw.events(m.M[,t][fromA],pA)
        intunnel <- (substr(drawA,1,1)=="A")
        ttunnel <- as.numeric(substr(drawA[intunnel],2,2))
        if(any(ttunnel<=d_at*interval)) {
          drawA[intunnel][ttunnel<=d_at*interval] <- paste0("A",ttunnel[ttunnel<=d_at*interval]+1)
        }
        m.M[,t+1][fromA] <- drawA
      }
      
      # From BS
      fromBS <- m.M[,t] %in% c("BS1","BS2") 
      if(sum(fromBS)>0) {
        pBS <- ProbsBS(m.M[,t][fromBS],t,p.HD)
        drawBS <- rbinom(n=length(pBS),size=1,prob=pBS)
        m.M[,t+1][fromBS][drawBS==1] <- "D"
        m.M[,t+1][fromBS][drawBS==0] <- "BS2"
      }
      
      # From BD
      fromBD1 <- (m.M[,t] == "BD1")
      if(sum(fromBD1)>0) {
        m.M[,t+1][fromBD1] <- "BD2"
      }
      
    }        
    
    
    
    
    ### compute C/E
    zerocol <- as.data.frame(setNames(replicate(length(v.n),numeric(0), simplify = F),v.n))
    
    mm <- data.frame(m.M) %>% mutate(name=as.integer(rownames(.))) %>%
      gather("cycle","state",starts_with("cycle")) %>% mutate(cycle=as.integer(str_replace(cycle,"cycle.",""))) %>%
      arrange(name,cycle)
    
    pop.test.alt <- df.Pop %>% filter(treat=="Alternate" | is.na(treat))
    pop.test.prim <- df.Pop %>% filter(tested==TRUE, treat=="Primary")
    pop.none <- df.Pop %>% filter(tested==FALSE)
    cc <- list()
    ee <- list()
    c.test <- list()
    c.drug <- list()
    d.r <- (1 + disc)^(1/interval)-1
    v.dwc <- 1 / (1 + d.r) ^ (0:(n.t-1)) # calculate discount weights for costs for each cycle based on discount rate d.r
    v.dwe <- 1 / (1 + d.r) ^ (0:(n.t-1))
    
    #no testing folks
    if(nrow(pop.none)>0) {
      m.test.none <- mm %>% filter(name %in% pop.none$name) %>% group_by(cycle) %>% count(state) %>%
        spread(state,n) %>% merge(zerocol,all=TRUE) %>% replace(is.na(.),0) %>% 
        select_at(c("cycle",v.n)) %>% ltm(method=method)
      
      m1 <- m.test.none[,-1] %>% as.matrix()
      dd <- c_tx*365/interval
      cc[[1]] <- m1 %*% c(0,c_a+dd,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0)
      ee[[1]] <- m1 %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)
      cc[[1]] <- as.vector(cc[[1]]) * v.dwc
      ee[[1]] <- as.vector(ee[[1]]) * v.dwe
      c.test[[1]] <- 0
      c.drug[[1]] <- m1 %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0)
      c.drug[[1]] <- as.vector(c.drug[[1]]) * v.dwc
    }
    
    
    if(nrow(pop.test.prim)>0) {
      m.test.prim <- mm %>% filter(name %in% pop.test.prim$name) %>% group_by(cycle) %>% count(state) %>%
        spread(state,n) %>% merge(zerocol,all=TRUE) %>% replace(is.na(.),0) %>% 
        select_at(c("cycle",v.n)) %>% ltm(method=method)
      m2 <- m.test.prim[,-1] %>% as.matrix()
      dd <- c_tx*365/interval
      cc[[2]] <- m2 %*% c(0,c_a+dd+c_t,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0)
      ee[[2]] <- m2 %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)
      cc[[2]] <- as.vector(cc[[2]]) * v.dwc
      ee[[2]] <- as.vector(ee[[2]]) * v.dwe
      c.test[[2]] <- m2[,"A1"]*c_t*v.dwc
      c.drug[[2]] <- m2 %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0)
      c.drug[[2]] <- as.vector(c.drug[[2]]) * v.dwc
    }
    
    if(nrow(pop.test.alt)>0) {
      m.test.alt <- mm %>% filter(name %in% pop.test.alt$name) %>% group_by(cycle) %>% count(state) %>%
        spread(state,n) %>% merge(zerocol,all=TRUE) %>% 
        replace(is.na(.),0) %>% select_at(c("cycle",v.n)) %>% ltm(method=method)
      
      m3 <- m.test.alt[,-1] %>% as.matrix()
      dd <- c_alt*365/interval
      cc[[3]] <- m3 %*% c(0,c_a+dd+c_t,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0)
      ee[[3]] <- m3 %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)
      cc[[3]] <- as.vector(cc[[3]]) * v.dwc
      ee[[3]] <- as.vector(ee[[3]]) * v.dwe
      c.test[[3]] <- m3[,"A1"]*c_t*v.dwc
      c.drug[[3]] <- m3 %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0)
      c.drug[[3]] <- as.vector(c.drug[[3]]) * v.dwc
      
    }
    
    #other
    all <- mm %>% group_by(cycle) %>% count(state) %>%
      spread(state,n) %>% merge(zerocol,all=TRUE) %>% arrange(cycle) %>% 
      replace(is.na(.),0) %>% select_at(c("cycle",v.n)) %>% ltm(method=method)
    mmm <- all[,-1] %>% as.matrix()
    possible <- mmm %*% c(1,rep(1,d_at*interval+1),1,1,0,0,0)
    possible <-  as.vector(possible) * v.dwe
    fatal_b <- all[n.t,c("BD1","BD2")] %>% sum()
    living <- n.i - all[n.t,c("BD1","BD2","D")] %>% sum()
    disutil_a <- mmm  %*% c(0,rep(d_a,d_at*interval),0,0,0,0,0,0)
    disutil_a <- as.vector(disutil_a) * v.dwe
    disutil_b <- mmm  %*% c(0,rep(0,d_at*interval+1),d_b,d_b,0,0,0)
    disutil_b <- as.vector(disutil_b) * v.dwe
    c.treat <- mmm  %*% c(0,c_a,rep(0,d_at*interval),c_bs,0,c_bd,0,0)
    c.treat <-  as.vector(c.treat) * v.dwc
    
    return(list(raw=m.M,pop=df.Pop,count=all,
                results=                c(dCOST       = sum(unlist(cc)),
                                          dQALY       = sum(unlist(ee)),
                                          possible    = sum(unlist(possible)),
                                          fatal_b     = unname(fatal_b),
                                          living      = unname(living),
                                          disutil_a   = sum(unlist(disutil_a)),
                                          disutil_b   = sum(unlist(disutil_b)),
                                          dCOST.test  = sum(unlist(c.test)),
                                          dCOST.drug  = sum(unlist(c.drug)),
                                          dCOST.treat = sum(unlist(c.treat)))/n.i
    ))
  })   
}

microsim_icer <- function(params, reference=NULL, genotype=NULL, N = NULL, seed=123, method="beginning")
{
  set.seed(seed)
  params$p_o <- 0.0 # No testing, reference
  reference <- microsim_run(params,method=method)
  reference <- reference$results
  
  set.seed(seed)
  params$p_o <- 1.0 # Genotype testing upon indication
  genotype   <- microsim_run(params,method=method)
  genotype <- genotype$results
  
  c( ICER       = unname((genotype['dCOST'] - reference['dCOST']) / (genotype['dQALY'] - reference['dQALY'])),
     NMB        = unname((reference['dCOST'] - genotype['dCOST']) + params$wtp*(genotype['dQALY'] - reference['dQALY'])),
     dCOST.ref  = unname(reference['dCOST']),
     dCOST.test = unname(genotype['dCOST']),
     dQALY.ref  = unname(reference['dQALY']),
     dQALY.test = unname(genotype['dQALY'])
  )
}


