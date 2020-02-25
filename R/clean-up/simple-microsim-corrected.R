# Load the following packages and scripts to run this model independently
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(Matrix)
source(here("common.R"))

### 01 inputs ####
source(here("simple-params.R"))
source(here("simple-markov-corrected.R"))

#Get embedded transition probability matrices defined in Markov model (only run yearly cycle, i.e. interval=1)
markov_corr_get <- function(params, gene=0, test=0)
{
  params$p_g <- gene
  params$p_o <- test

  with(params, {

    #### 1. Model Setting ####
    n.t       <- horizon*interval          # number of cycles
    v.n       <- c("H", "A", "BS", "BD", "D", "BSD", "CUM_A", "CUM_BS", "CUM_T", "TUN", "PTUN")  # state names, and accumulators
    n.s       <- length(v.n)               # number of states

    # transition probability matrix list
    l.P <- lapply(1:n.t, function(x) markov_corr_tp(params, v.n, x))

    return(l.P)
  })
}

m00 <- markov_corr_get(params,0,0) #no testing
m11 <- markov_corr_get(params,1,1) #tested carrier


#### Model-Specific Functions ####
# Transition Probability Functions:
# Feed in current Markov trace, time cycle, and other parameters to look up probability of transition from one state
# Plug in embedded transition probability matrices based on individual gene type
ProbsH_corr <- function(M_t, t, gene, test)
{
    p <- matrix(0,nrow=length(M_t),ncol=5,dimnames=list(c(),c("A1","HBS","HBD","D","HBSD"))) #additional jump-over states added compared to non-embedded microsim
    #assume behavior parameter only takes 0 or 1, otherwise do not know which matrix to use
    if(test==1) {
        b0 <- !gene
        b1 <- gene
        
        p[b0,"A1"] <- m00[[t]]["H","A"]
        p[b0,"HBS"] <- m00[[t]]["H","BS"]
        p[b0,"HBD"] <- m00[[t]]["H","BD"]
        p[b0,"D"] <- m00[[t]]["H","D"]
        p[b0,"HBSD"] <- m00[[t]]["H","BSD"]
        
        p[b1,"A1"] <- m11[[t]]["H","A"]
        p[b1,"HBS"] <- m11[[t]]["H","BS"]
        p[b1,"HBD"] <- m11[[t]]["H","BD"]
        p[b1,"D"] <- m11[[t]]["H","D"]
        p[b1,"HBSD"] <- m11[[t]]["H","BSD"]
    } else if(test==0) {
        p[,"A1"] <- m00[[t]]["H","A"]
        p[,"HBS"] <- m00[[t]]["H","BS"]
        p[,"HBD"] <- m00[[t]]["H","BD"]
        p[,"D"] <- m00[[t]]["H","D"]
        p[,"HBSD"] <- m00[[t]]["H","BSD"]
    } else stop("invalid p_o")
    return(p)
}

ProbsA_corr <- function(M_t, t, gene, test) {
    p <- matrix(0,nrow=length(M_t),ncol=4,dimnames=list(c(),c("ABS","ABD","D","ABSD")))
    if(test==1) {
        b0 <- !gene
        b1 <- gene
        
        p[b0,"ABS"] <- m00[[t]]["A","BS"]
        p[b0,"ABD"] <- m00[[t]]["A","BD"]
        p[b0,"D"] <- m00[[t]]["A","D"]
        p[b0,"ABSD"] <- m00[[t]]["A","BSD"]
        
        p[b1,"ABS"] <- m11[[t]]["A","BS"]
        p[b1,"ABD"] <- m11[[t]]["A","BD"]
        p[b1,"D"] <- m11[[t]]["A","D"]
        p[b1,"ABSD"] <- m11[[t]]["A","BSD"]
    } else if(test==0) {
        p[,"ABS"] <- m00[[t]]["A","BS"]
        p[,"ABD"] <- m00[[t]]["A","BD"]
        p[,"D"] <- m00[[t]]["A","D"]
        p[,"ABSD"] <- m00[[t]]["A","BSD"]
    } else stop("invalid p_o")
    return(p)
}


ProbsBS_corr <- function(M_t, t, gene, test){
    p <- matrix(0,nrow=length(M_t),ncol=1,dimnames=list(c(),c("D")))
    if(test==1) {
        b0 <- !gene
        b1 <- gene
        p[b0,"D"] <- m00[[t]]["BS","BSD"]
        p[b1,"D"] <- m11[[t]]["BS","BSD"]
    } else if(test==0) {
        p[,"D"] <- m00[[t]]["BS","BSD"]
    } else stop("invalid p_o")
    return(p)
}

# Function to Draw Events for Competing Risks
draw.events <- function(current,p) { #feed in current state and transition probabilities
    pp <- rbinom(n=length(current),size=1,prob=rowSums(p))
    if(sum(pp)>0) { # if at least one person has a event then randomly sample the event type for those who experience any event
        event.type <- samplev(p[pp==1,,drop=FALSE]/rowSums(p[pp==1,,drop=FALSE]),m=1)   
        # update state based on event type
        current[pp==1] = event.type
    }
    return(unname(current))
}


#### Main Simulation ####
microsim_run_corr <- function(params, N = NULL, method="beginning")
{
    if (!is.null(N)) params$n <- N
    
    # TO DEBUG USE
    # sapply(names(params), function(p) assign(p, params[[p]], 1))
    with(params, {
        
        #### 1. Model Setting ####
        n.i       <- n                          # number of individuals
        n.t       <- horizon*interval           # number of cycles
        v.n       <- c("H",paste0("A",1:(d_at*interval+1)),"HBS","ABS","BS2","HBD","ABD","BD2","D","HBSD","ABSD")    # state names
        n.s       <- length(v.n)                # number of states
        
        #secular death risk
        p.HD      <- gompertz_ratio2(1:n.t, interval, shape, rate)
        
        #### 2. Sample individual level characteristics ####
        df.Pop    <- data.frame(name    = 1:n.i,
                                variant = sample(c(FALSE,TRUE), n.i, prob=c(1-p_g,p_g), replace = TRUE),
                                tested  = FALSE,
                                treat   = factor("",levels=c("","Alternate","Primary")))
        # "tested" and "treat" are dynamic attributes and will be updated along with the simulation
        v.M_Init <- rep("H", n.i)       # everyone begins in the healthy state 
        
        #### 3. Initiate Simulation ####
        # m.M: state for each patient, only capturing the current cycle and the one following to avoid storing giant matrices
        m.M  <- matrix(nrow = n.i, ncol = 2, 
                       dimnames = list(paste(1:n.i),     # could name the rows ind1, ind2, ind3, etc.
                                       c("current","forward")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.
        
        m.M[, 1] <- v.M_Init         # initial health state for all individuals
        
        # create a placeholder list to record state counts in each cycle
        ct <- vector("list",n.t)
      
        #### 4. Run Simultation ####
        for (t in 1:n.t)
        {
          ## Simulate Events ##
          # Default Update: if no events occur, current state is copied to the next cycle
          m.M[,2] = m.M[,1]
          
          # From H
          fromH <- m.M[,1]=="H"
            if(sum(fromH) > 0)
            {
                pH <- ProbsH_corr(m.M[,1][fromH],t,df.Pop$variant[fromH],p_o)
                m.M[,2][fromH] <- draw.events(m.M[,1][fromH],pH)
                
                # assign treatment upon indication (including jumping thru A)
                elig <- (m.M[,2] %in% c("A1","HBS","HBD","HSBD"))
                df.Pop$tested[elig] <- sample(c(TRUE,FALSE),sum(elig),prob=c(p_o,1-p_o),replace = TRUE)
                df.Pop$treat[elig]  <- ifelse(df.Pop$tested[elig] & df.Pop$variant[elig],"Alternate","Primary")
            }

            # From A
            fromA <- m.M[,1] %in% paste0("A",1:(d_at*interval+1))
            
            if(sum(fromA)>0) 
            {
                pA <- ProbsA_corr(m.M[,1][fromA],t,df.Pop$variant[fromA],p_o)
            
                # move forward in A tunnel as needed
                drawA <- draw.events(m.M[,1][fromA],pA)
                intunnel <- (drawA %in% paste0("A",1:(d_at*interval+1)))
                ttunnel <- as.numeric(substring(drawA[intunnel],2)) # which tunnel cycle
                if(any(ttunnel<=d_at*interval)) {
                    drawA[intunnel][ttunnel<=d_at*interval] <- paste0("A",ttunnel[ttunnel<=d_at*interval]+1)
                }
                m.M[,2][fromA] <- drawA
            }
            
            # From BS 
            fromBS <- m.M[,1] %in% c("HBS","ABS","BS2") 
            if(sum(fromBS)>0)
            {
                pBS <- ProbsBS_corr(m.M[,1][fromBS],t,df.Pop$variant[fromBS],p_o)
                drawBS <- rbinom(n=length(pBS),size=1,prob=pBS)
                m.M[,2][fromBS][drawBS==1] <- "D"
                m.M[,2][fromBS][drawBS==0] <- "BS2"

            }
            
            # From BD
            m.M[,2][m.M[,1] %in% c("HBD","ABD")] <- "BD2"
            m.M[,2][m.M[,1] %in% c("HBSD","ABSD")] <- "D"
            
            # Depending on testing decision/treatment, individuals incur different costs/disutilities,
            # therefore they are grouped into three tracks: 
            # 1.not tested; 2. tested & on primary drug; 3. tested & on alternative drug.
            pop.none      <- df.Pop[df.Pop$tested==FALSE,                            'name']
            pop.test.prim <- df.Pop[df.Pop$tested==TRUE & df.Pop$treat=="Primary",   'name']
            pop.test.alt  <- df.Pop[df.Pop$treat=="Alternate" | is.na(df.Pop$treat), 'name']
            
            ct[[t]] <- list(pop.none, pop.test.prim, pop.test.alt) %>% 
              map(function(x) {
                if(length(x)==0) {
                  return(as.data.frame(as.list(rep(0,length(v.n))),row.names = "1",col.names = v.n))
                  } else {
                  m.M[x,"forward"] %>% as.data.frame() %>% set_names("state") %>%
                    mutate(state=factor(state,levels=v.n)) %>% 
                    table() %>% as.data.frame() %>% spread_(".","Freq")
                }
              })
            
            # move one cycle forward to reset the matrix
            m.M[,1] <- m.M[,2]
        }        
        
        
        # condense the counter list to get state counts for all cylces in each track
        mm <- map(1:3, function(x) map_df(1:n.t,~ct[[.x]][[x]]))
        
        ####5. Discounting & Summation ####
        cc            <- numeric(3)
        ee            <- numeric(3)
        c.test        <- numeric(3)
        c.drug        <- numeric(3)
        d.r  <- inst_rate(1-1/(1 + disc), interval)
        
        for(i in 1:3)
        {
          if(nrow(mm[[i]]) > 0)
          {
                # Jump-over states accrues additional costs/disutilities compared to non-embedded microsim
                dd        <- (if(i==3) c_alt else c_tx)*365/interval
                m0        <- rbind(c(sum(mm[[i]][1,]),rep(0,ncol(mm[[i]]))),mm[[i]]) %>% as.matrix() #add back t=0 row to ensure all groups sum up to n.i in healthy state at t=0 cycle
                m1        <- integrator(diag(exp( 0:n.t * -d.r)) %*% m0, method=method) #first discount then integrate
                cc[i]     <- sum(as.vector(m1 %*% c(0,if(i==1) c_a+dd else c_a+dd+c_t,rep(dd,d_at*interval),c_t+c_a+c_bs+dd,c_bs+dd,dd,c_t+c_a+c_bd,c_bd,0,0,c_t+c_a+c_bs,c_bs)))
                ee[i]     <- sum(as.vector(m1 %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,rep((1-d_b)/interval,2),(1-d_b)/interval,rep(0,6))))
                c.test[i] <- if(i==1) 0 else sum(as.vector(m1 %*% c(0,c_t,rep(0,d_at*interval),c_t,0,0,c_t,0,0,0,c_t,0)))
                c.drug[i] <- sum(as.vector(m1 %*% c(0,dd,rep(dd,d_at*interval),rep(dd,3),rep(0,6))))
            }
        }

        tout  <- map(mm, function(x) mutate(x,cycle=as.integer(row.names(x)))) %>% do.call("bind_rows",.) %>% 
          group_by(cycle) %>% summarise_all(sum) %>% select(-cycle) #combine counts from all three tracks
        mmm <- rbind(c(sum(tout[1,]),rep(0,ncol(tout[1,])-1)),tout) %>% as.matrix() #add back t=0 when all inviduals are in healthy state
        mmm2 <- integrator(mmm,method = method) #non-discounted integration
        dmm <- integrator(diag(exp( 0:n.t * -d.r)) %*% mmm, method=method) #discounted integration
        
        # other metrics 
        possible  <- as.vector(dmm %*% c(1,rep(1,d_at*interval+1),rep(1,3),rep(0,6)))/interval
        fatal_b   <- sum(mmm2[n.t,c("HBD","ABD","BD2")])
        living    <- n.i - sum(mmm2[n.t,c("HBD","ABD","BD2","D","HBSD","ABSD")])
        disutil_a <- as.vector(dmm  %*% c(0,rep(d_a,d_at*interval),rep(0,10)))/ interval
        disutil_b <- as.vector(dmm  %*% c(0,rep(0,d_at*interval+1),d_b,d_b,d_b,rep(0,6)))/ interval
        c.treat   <- as.vector(dmm  %*% c(0,c_a,rep(0,d_at*interval),c_a+c_bs,c_bs,0,c_a+c_bd,c_bd,0,0,c_a+c_bs,c_bs))
        
        list(
            #raw_ct = mmm,
            count = mmm2,
            pop     = df.Pop,
            results = c(dCOST       = sum(cc),
                        dQALY       = sum(ee),
                        possible    = sum(possible),
                        fatal_b     = unname(fatal_b),
                        living      = n.i*unname(living),
                        disutil_a   = sum(disutil_a),
                        disutil_b   = sum(disutil_b),
                        dCOST.test  = sum(c.test),
                        dCOST.drug  = sum(c.drug),
                        dCOST.treat = sum(c.treat))/n.i
        )
        
    })  
}

# Combined model: Run Both Strategies
microsim_corr_icer <- function(params, reference=NULL, genotype=NULL, method="life-table",seed=NULL,...)
{
    if(!is.null(seed)) set.seed(seed)
    
    seed         <- .Random.seed  
    
    params$p_o   <- 0.0 # No testing, reference
    reference    <- microsim_run_corr(params,method=method,...)
    reference    <- reference$results
    
    .Random.seed <- seed
    params$p_o   <- 1.0 # Genotype testing upon indication
    genotype     <- microsim_run_corr(params,method=method,...)
    genotype     <- genotype$results
    
    c( ICER       = unname((genotype['dCOST'] - reference['dCOST']) / (genotype['dQALY'] - reference['dQALY'])),
       NMB        = unname((reference['dCOST'] - genotype['dCOST']) + params$wtp*(genotype['dQALY'] - reference['dQALY'])),
       dCOST.ref  = unname(reference['dCOST']),
       dCOST.test = unname(genotype['dCOST']),
       dQALY.ref  = unname(reference['dQALY']),
       dQALY.test = unname(genotype['dQALY'])
    )
}

# microsim_corr_icer(params,N=1000,seed=123) #sample code to run the model
