rm(list =ls()) 

library(data.table)
library(tidyverse)
library(flexsurv)

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
ProbsE <- function(M_t,v.drug,t,p.a,p.b,rr.b,p.bd,p.sd){
  # M_t    :    current health state
  # t    :  cycle
  # v.drug :    treatment
  p.die.base <- p.sd[t]
  p <- matrix(0,nrow=length(M_t),ncol=4,dimnames=list(c(),c("A","BS","BD","D")))
  
  #cap probs to match markov version: is it appropriate????
  if(p.die.base+p.a>1) {p.a <- 1-p.die.base}
  if(p.die.base+p.b>1) {p.b <- 1-p.die.base}
  
  p[M_t == "H","A"] <- p.a
  p[M_t == "A","BS"]  <- ifelse(v.drug[M_t == "A"]=="Alternate", 
                                p.b*rr.b*(1-p.bd),
                                p.b*(1-p.bd))
  p[M_t == "A","BD"]  <- ifelse(v.drug[M_t == "A"]=="Alternate", 
                                p.b*rr.b*p.bd,
                                p.b*p.bd)
  p[M_t == "H","D"] <- p[M_t == "A","D"] <- p[M_t == "BS","D"] <- p.die.base
  return(p)
}


microsim_run <- function(params, N = NULL) {
  if (!is.null(N)) params$n <- N
  with(params, {
    
    #### 1. Model Setting ####
    n.i       <- n              # number of individuals
    n.t       <- horizon*interval           # number of cycles
    v.n       <- c("H","A","BS","BD","D")    # state names
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
                   dimnames = list(paste(1:n.i),     # could ame the rows ind1, ind2, ind3, etc.
                                   paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.
    
    # Intialize all individuals
    m.M[, 1] <- v.M_Init         # initial health state for all individuals
    
    # Time for-loop from cycle 1 through n.t
    for (t in 1:n.t) {
      
      ## Simulate Events ##
      # Default Update: if no events occur, health state remains unchanged
      m.M[,t+1] = m.M[,t]
      # Progression of illness 
      p <- ProbsE(m.M[,t],df.Pop$treat,t,r_a,r_b,rr_b,p_bd,p.HD)
      # note: if someone died this cycle, they can't experience other events
      draw.event <- rbinom(n=n.i,size=1,prob=rowSums(p))
      if(sum(draw.event)>0) { # if at least one person has a event
        # randomly sample the event type for those who experience progression events this cycle
        # uses samplev() from the "Functions Folder"
        event.type <- samplev(p[draw.event==1,,drop=FALSE]/rowSums(p[draw.event==1,,drop=FALSE]),m=1)   
        # update state based on event type
        m.M[draw.event==1,t+1] = event.type
        
        # assign treatment upon indication
        elig <- (draw.event==1 & m.M[,t+1]=="A")
        df.Pop$tested[elig] <- sample(c(TRUE,FALSE),sum(elig),prob=c(p_o,1-p_o),replace = TRUE)
        df.Pop$treat[elig] <- ifelse(df.Pop$tested[elig] & df.Pop$variant[elig],"Alternate","Primary")
      }}
    
    
    ### reframe results to feed des summary functions
    re <- df.Pop %>% arrange(name)
    mm <- data.frame(m.M) %>% mutate(name=as.integer(rownames(.))) %>%
      gather("cycle","state",starts_with("cycle")) %>% mutate(cycle=as.integer(str_replace(cycle,"cycle.",""))) %>%
      arrange(name,cycle)

    re$secular_death <- find_time(mm,"D",n.i)
    re$indication <- find_time(mm,"A",n.i)
    re$adverse <- find_time(mm,c("BS","BD"),n.i)
    re$adverse_death <- find_time(mm,"BD",n.i)

    re$death <- re$secular_death
    re$death[!is.na(re$adverse_death)] <- re$adverse_death[!is.na(re$adverse_death)]

    re$end_of_sim <- re$death
    re$end_of_sim[is.na(re$end_of_sim)] <- n.t+1
    re$cutoff <- re$indication + d_at
    re$cutoff <- ifelse(!is.na(re$cutoff) & !is.na(re$adverse) & re$adverse < re$cutoff, re$adverse, re$cutoff)
    re$cutoff <- ifelse(!is.na(re$cutoff) & !is.na(re$secular_death) & re$secular_death < re$cutoff, re$secular_death, re$cutoff)
    re$cutoff[re$cutoff > n.t] <- n.t+1
    
    return(list(raw=m.M,results=re))
  })   
}


####### end cycle value ###### 

find_time <- function(dt,st,n.i) {
  ss <- dt %>% filter(state %in% st) %>% group_by(name) %>% summarise(when=min(cycle)) %>% arrange(name)
  dd <- vector("integer",n.i)
  dd[ss$name] <- ss$when
  dd[dd==0] <- NA
  return(dd)
}


####03 Post-sim summary####


# Discounted integral from A to B at annual rate ar of discounting
discount_int2 <- function(dr, A, B)
{
  q <- 1/(1+dr)
  (q^(A-1)-q^(B-1))/(1-q)
}

# Discounted Integral over interval
disum2 <- function(dr, A, B) {
  sum(discount_int2(dr, A, B), na.rm=TRUE)
}

# Discounted Discrete sum (aka single event discounting)
# Each value is a time of event in units of days (365 per year)
ddsum2 <- function(dr,A) {
  sum(1 / ((1 + dr)^(A-1)), na.rm=TRUE)
}

microsim_summary <- function(df, params, N = NULL)
{
  if (!is.null(N)) params$n = N
  with(c(df[["results"]], params), {
    
    disc <- (1 + disc)^(1/interval)-1
    
    # Testing costs
    test.cost  <- c_t  * ddsum2(disc,indication[tested])
    treat.cost <- c_a  * ddsum2(disc,indication) + 
      c_bs * ddsum2(disc,adverse[is.na(adverse_death)]) +
      c_bd * ddsum2(disc,adverse_death) 
    
    pri <- !is.na(treat) & treat=="Primary"
    alt <- !is.na(treat) & treat=="Alternate"
    
    drug.cost <- sum(c_tx*365*disum2(disc, indication[pri], end_of_sim[pri]),
                     c_alt*365*disum2(disc, indication[alt], end_of_sim[alt]))
    
    # Total living in model
    life <- (n - sum(!is.na(death)))
    
    # Total possible discounted life units is integral of discounted time
    pQALY <- disum2(disc, 1, end_of_sim) / interval
    
    # Temp disutility of Indication
    disA   <- d_a*disum2(disc, indication, cutoff) / interval
    
    # Permanent disutility for Event
    disB   <- d_b*disum2(disc, adverse, end_of_sim) / interval
    
    c(dCOST       = unname(treat.cost+test.cost+drug.cost),
      dQALY       = unname(pQALY-disA-disB),
      possible    = unname(pQALY),
      fatal_b     = sum(!is.na(adverse_death)),
      living      = unname(life),
      disutil_a   = unname(disA),
      disutil_b   = unname(disB),
      dCOST.test  = unname(test.cost),
      dCOST.drug  = unname(drug.cost),
      dCOST.treat = unname(treat.cost)
    )/n
  })
}


microsim_icer <- function(params, reference=NULL, genotype=NULL, N = NULL, seed=123)
{
  set.seed(seed)
  params$p_o <- 0.0 # No testing, reference
  reference <- microsim_summary(if(is.null(reference)) microsim_run(params,N=N) else reference, params,N=N)
  
  set.seed(seed)
  params$p_o <- 1.0 # Genotype testing upon indication
  genotype   <- microsim_summary(if(is.null(genotype)) microsim_run(params,N=N) else genotype, params,N=N)
  
  c( ICER       = unname((genotype['dCOST'] - reference['dCOST']) / (genotype['dQALY'] - reference['dQALY'])),
     NMB        = unname((reference['dCOST'] - genotype['dCOST']) + params$wtp*(genotype['dQALY'] - reference['dQALY'])),
     dCOST.ref  = unname(reference['dCOST']),
     dCOST.test = unname(genotype['dCOST']),
     dQALY.ref  = unname(reference['dQALY']),
     dQALY.test = unname(genotype['dQALY'])
  )
}

# ### compare to Markov ###
# params$n = 1000
library(heemod)
source("./R/simple-markov.R")
# ### icer
# icer_micro <- microsim_icer(params)
# icer_mark <- markov_icer(params)

# micro_comp <- vector("list",5)
# mark_comp <- vector("list",5)
# ss <- c(32,532,13,3267,8098)
# for(s in 1:5) {
#   set.seed(ss[s])
#   micro_comp[[s]] <- microsim_icer(params)
#   mark_comp[[s]] <- markov_icer(params)
# }


### event
#compare testing strategy
# params$n = 1000000
# params$p_o = 1
# micro <- microsim_run(params)
# mak <- markov_simulation(params)
# comp <-  mak$combined_model$eval_strategy_list$reference$counts %>% mutate(cycle=as.integer(rownames(.)))
# c2 <-  micro$raw %>% data.frame() %>% gather() %>% count(key,value)
# c2 <- c2 %>% spread(value,n) %>% mutate(key=as.integer(str_replace(key,"cycle.",""))) %>% arrange(key) %>% rename(cycle=key)
# c3 <- c2[-1,] %>% replace_na(list(A=0,BS=0,BD=0,D=0))
# 
# dt <- comp %>% gather("state","ct",-cycle) %>% mutate(model="mark") %>% bind_rows(
#   c3 %>% gather("state","ct",-cycle) %>% mutate(model="micro",ct=ct/params$n*1000)
# )
# ggplot(dt,aes(x=cycle)) + geom_line(aes(y=ct,group=state,color=state)) + facet_grid(~model) + theme_bw()
# 
# ### life-table method: I think microsim current version is like beginning methods ###
# ck <- dt %>% spread(model,ct) %>% mutate(dd=mark-micro)
# ggplot(ck,aes(x=cycle)) + geom_line(aes(y=dd,group=state,color=state)) + theme_bw()
