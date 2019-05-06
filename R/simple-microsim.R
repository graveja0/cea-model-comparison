# library(data.table)
# library(tidyverse)
# library(flexsurv)

#### 01 inputs ####
source("./simple-params.R")

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
ProbsE <- function(M_t,v.drug,t,p.a,p.b,rr.b,p.bd,p.sd,interval,dur_a){
  # st : all states name
  # M_t    :    current health state
  # t    :  cycle
  # v.drug :    treatment
  p.die.base <- p.sd[t]
  p <- matrix(0,nrow=length(M_t),ncol=6+dur_a,dimnames=list(c(),c(paste0("A",1:(dur_a+1)),"BS1","BS2","BD1","BD2","D")))
  
  #not sure about using this heemod func
  library(heemod)
  pp.a <- rate_to_prob(p.a,1/interval)
  pp.b <- rate_to_prob(p.b,1/interval)
  
  #cap probs to match markov version: is it appropriate????
  if(p.die.base+p.a>1) {pp.a <- 1-p.die.base}
  if(p.die.base+p.b>1) {pp.b <- 1-p.die.base}
  
  
  p[M_t == "H","A1"] <- pp.a
  p[M_t == "A1","BS1"] <- ifelse(v.drug[M_t == "A1"]=="Alternate", pp.b*rr.b*(1-p.bd),pp.b*(1-p.bd))
  p[M_t == "A2","BS1"] <- ifelse(v.drug[M_t == "A2"]=="Alternate", pp.b*rr.b*(1-p.bd),pp.b*(1-p.bd))
  p[M_t == "A1","BD1"] <- ifelse(v.drug[M_t == "A1"]=="Alternate", pp.b*rr.b*p.bd,pp.b*p.bd)
  p[M_t == "A2","BD1"] <- ifelse(v.drug[M_t == "A2"]=="Alternate", pp.b*rr.b*p.bd,pp.b*p.bd)
  
  p[M_t == "H","D"] <- p[M_t == "A1","D"] <- p[M_t == "A2","D"] <- p[M_t == "BS1","D"] <- p[M_t == "BS2","D"] <- p.die.base
  p[M_t == "A1","A2"] <- 1-p.die.base-pp.b
  p[M_t == "BS1","BS2"] <- 1-p.die.base
  p[M_t == "BD1","BD2"] <- 1
  
  if(dur_a>1) {
    i <- 2
    while(i<=dur_a) {
      p[M_t == paste0("A",i+1),"BS1"] <- ifelse(v.drug[M_t == paste0("A",i+1)]=="Alternate", pp.b*rr.b*(1-p.bd),pp.b*(1-p.bd))
      p[M_t == paste0("A",i+1),"BD1"] <- ifelse(v.drug[M_t == paste0("A",i+1)]=="Alternate", pp.b*rr.b*p.bd,pp.b*p.bd)
      p[M_t == paste0("A",i),paste0("A",i+1)] <- 1-p.die.base-pp.b
      i <- i+1
    }
  }
  
  return(p)
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
    v.n       <- c("H",paste0("A",1:(d_at+1)),"BS1","BS2","BD1","BD2","D")    # state names
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
    
    a.P <- array(0,                                          # Create 3-D array
                 dim = c(n.i, length(v.n)-1, n.t),
                 dimnames = list(c(), v.n[-1], 1:n.t))
    
    m.M[, 1] <- v.M_Init         # initial health state for all individuals
    
    # Time for-loop from cycle 1 through n.t
    for (t in 1:n.t) {
      
      ## Simulate Events ##
      # Default Update: if no events occur, health state remains unchanged
      m.M[,t+1] = m.M[,t]
      # Progression of illness 
      p <- ProbsE(m.M[,t],df.Pop$treat,t,r_a,r_b,rr_b,p_bd,p.HD,interval,d_at)
      a.P[,,t] <- p
      # note: if someone died this cycle, they can't experience other events
      draw.event <- rbinom(n=n.i,size=1,prob=rowSums(p))
      if(sum(draw.event)>0) { # if at least one person has a event
        # randomly sample the event type for those who experience progression events this cycle
        # uses samplev() from the "Functions Folder"
        event.type <- samplev(p[draw.event==1,,drop=FALSE]/rowSums(p[draw.event==1,,drop=FALSE]),m=1)   
        # update state based on event type
        m.M[draw.event==1,t+1] = event.type
        
        # assign treatment upon indication
        elig <- (draw.event==1 & m.M[,t+1]=="A1")
        if(sum(elig)>0) {
          df.Pop$tested[elig] <- sample(c(TRUE,FALSE),sum(elig),prob=c(p_o,1-p_o),replace = TRUE)
          df.Pop$treat[elig] <- ifelse(df.Pop$tested[elig] & df.Pop$variant[elig],"Alternate","Primary")
        }
      }
    }
    
    ### depending on p_o
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
    if(p_o==1) {
      zerocol <- matrix(rep(0,(n.t+1)*(length(v.n)-2)),nrow=n.t+1,ncol=length(v.n)-2) %>% 
        data.frame() %>% set_names(v.n[!(v.n %in% c("H","D"))])
      m.test.none <- mm %>% filter(name %in% pop.none$name) %>% group_by(cycle) %>% count(state) %>%
        spread(state,n) %>% bind_cols(zerocol) %>% replace(is.na(.),0) %>% select_at(c("cycle",v.n)) %>% ltm(method=method)
    } else {
      m.test.none <- mm %>% filter(name %in% pop.none$name) %>% group_by(cycle) %>% count(state) %>%
        spread(state,n) %>% replace(is.na(.),0) %>% select_at(c("cycle",v.n)) %>% ltm(method=method)
    }
    
    m1 <- m.test.none[,-1] %>% as.matrix()
    dd <- c_tx*365/interval
    cc[[1]] <- m1 %*% c(0,c_a+dd,rep(dd,d_at),c_bs+dd,dd,c_bd,0,0)
    ee[[1]] <- m1 %*% c(1/interval,rep((1-d_a)/interval,d_at),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)
    cc[[1]] <- as.vector(cc[[1]]) * v.dwc
    ee[[1]] <- as.vector(ee[[1]]) * v.dwe
    c.test[[1]] <- 0
    c.drug[[1]] <- m1 %*% c(0,dd,rep(dd,d_at),dd,dd,0,0,0)
    c.drug[[1]] <- as.vector(c.drug[[1]]) * v.dwc
    
    
    if(nrow(pop.test.prim)>0) {
      m.test.prim <- mm %>% filter(name %in% pop.test.prim$name) %>% group_by(cycle) %>% count(state) %>%
        spread(state,n) %>% replace(is.na(.),0) %>% select_at(c("cycle",v.n)) %>% ltm(method=method)
      m2 <- m.test.prim[,-1] %>% as.matrix()
      dd <- c_tx*365/interval
      cc[[2]] <- m2 %*% c(0,c_a+dd+c_t,rep(dd,d_at),c_bs+dd,dd,c_bd,0,0)
      ee[[2]] <- m2 %*% c(1/interval,rep((1-d_a)/interval,d_at),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)
      cc[[2]] <- as.vector(cc[[2]]) * v.dwc
      ee[[2]] <- as.vector(ee[[2]]) * v.dwe
      c.test[[2]] <- m2[,"A1"]*c_t*v.dwc
      c.drug[[2]] <- m2 %*% c(0,dd,rep(dd,d_at),dd,dd,0,0,0)
      c.drug[[2]] <- as.vector(c.drug[[2]]) * v.dwc
    }
    
    if(nrow(pop.test.alt)>0) {
      m.test.alt <- mm %>% filter(name %in% pop.test.alt$name) %>% group_by(cycle) %>% count(state) %>%
        spread(state,n) %>% replace(is.na(.),0) %>% select_at(c("cycle",v.n)) %>% ltm(method=method)
      
      m3 <- m.test.alt[,-1] %>% as.matrix()
      dd <- c_alt*365/interval
      cc[[3]] <- m3 %*% c(0,c_a+dd+c_t,rep(dd,d_at),c_bs+dd,dd,c_bd,0,0)
      ee[[3]] <- m3 %*% c(1/interval,rep((1-d_a)/interval,d_at),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)
      cc[[3]] <- as.vector(cc[[3]]) * v.dwc
      ee[[3]] <- as.vector(ee[[3]]) * v.dwe
      c.test[[3]] <- m3[,"A1"]*c_t*v.dwc
      c.drug[[3]] <- m3 %*% c(0,dd,rep(dd,d_at),dd,dd,0,0,0)
      c.drug[[3]] <- as.vector(c.drug[[3]]) * v.dwc
      
    }
    
    #other
    all <- mm %>% group_by(cycle) %>% count(state) %>%
      spread(state,n) %>% replace(is.na(.),0) %>% select_at(c("cycle",v.n)) %>% ltm(method=method)
    mmm <- all[,-1] %>% as.matrix()
    possible <- mmm %*% c(1,rep(1,d_at+1),1,1,0,0,0)
    possible <-  as.vector(possible) * v.dwe
    fatal_b <- all[n.t,c("BD1","BD2")] %>% sum()
    living <- n.i - all[n.t,c("BD1","BD2","D")] %>% sum()
    disutil_a <- mmm  %*% c(0,rep(d_a,d_at),0,0,0,0,0,0)
    disutil_a <- as.vector(disutil_a) * v.dwe
    disutil_b <- mmm  %*% c(0,rep(0,d_at+1),d_b,d_b,0,0,0)
    disutil_b <- as.vector(disutil_b) * v.dwe
    c.treat <- mmm  %*% c(0,c_a,rep(0,d_at),c_bs,0,c_bd,0,0)
    c.treat <-  as.vector(c.treat) * v.dwc
    
    return(list(raw=m.M,pop=df.Pop,probs=a.P,count=all,
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


