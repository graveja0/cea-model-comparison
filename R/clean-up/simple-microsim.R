# This code built on a model example developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Citations:
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 

# Load the following packages and scripts to run this model independently
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(here)
# source(here("common.R")) #load shared functions
# source(here("simple-params.R")) #load inputs

#### Model-Specific Functions ####
# Transition Probability Functions:
# feed in current Markov trace, time cycle, and other parameters to look up probability of transition from one state
ProbsH <- function(M_t, t, p.a, p.sd, interval)
{
    # M_t  :    current health state 
    # t    :  cycle
    p.die.base <- p.sd[t]
    p <- matrix(0,nrow=length(M_t),ncol=2,dimnames=list(c(),c("A1","D")))
    pp.a <- rate_to_prob(p.a,1/interval)
    
    # Once secular death mortablity reaches 1, all other probs need to put 0, no effect for default 40-year time horizon
    # NOTE: This is used a lot in cost effectiveness research and is the result
    # of failing to treat it as competing risks when embedding a continuous rate into a probability
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
microsim_run <- function(params, N = NULL, method="beginning")
{
    if (!is.null(N)) params$n <- N
    
    # TO DEBUG USE
    # sapply(names(params), function(p) assign(p, params[[p]], 1))
    with(params, {
        
        #### 1. Model Setting ####
        n.i       <- n                          # number of individuals
        n.t       <- horizon*interval           # number of cycles
        v.n       <- c("H",paste0("A",1:(d_at*interval+1)),"BS1","BS2","BD1","BD2","D")    # state names
        n.s       <- length(v.n)                # number of states
        
        # secular death risk
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
                       dimnames = list(paste(1:n.i),     
                                       c("current","forward"))) 
        
        m.M[, 1] <- v.M_Init         # initial state for all individuals 
        
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
                pH <- ProbsH(m.M[,1][fromH],t,r_a,p.HD,interval)
                m.M[,2][fromH] <- draw.events(m.M[,1][fromH],pH)
                
                # assign treatment upon indication
                elig <- (m.M[,1]=="H" & m.M[,2]=="A1")
                df.Pop$tested[elig] <- sample(c(TRUE,FALSE),sum(elig),prob=c(p_o,1-p_o),replace = TRUE)
                df.Pop$treat[elig]  <- ifelse(df.Pop$tested[elig] & df.Pop$variant[elig],"Alternate","Primary")
            }
            
            # From A
            fromA <- m.M[,1] %in% paste0("A",1:(d_at*interval+1))
            
            if(sum(fromA)>0) 
            {
                pA <- ProbsA(m.M[,1][fromA],df.Pop$treat[fromA],t,r_b,rr_b,p_bd,p.HD,interval)
                
                # move forward in A tunnel as needed
                drawA <- draw.events(m.M[,1][fromA],pA)
                intunnel <- (substr(drawA,1,1)=="A")
                ttunnel <- as.numeric(substring(drawA[intunnel],2)) # which tunnel cycle
                if(any(ttunnel<=d_at*interval)) {
                    drawA[intunnel][ttunnel<=d_at*interval] <- paste0("A",ttunnel[ttunnel<=d_at*interval]+1)
                }
                m.M[,2][fromA] <- drawA
            }
            
            # From BS
            fromBS <- m.M[,1] %in% c("BS1","BS2") 
            if(sum(fromBS)>0)
            {
                pBS <- ProbsBS(m.M[,1][fromBS],t,p.HD)
                drawBS <- rbinom(n=length(pBS),size=1,prob=pBS)
                m.M[,2][fromBS][drawBS==1] <- "D"
                m.M[,2][fromBS][drawBS==0] <- "BS2"
                
            }
            
            # From BD
            m.M[,2][m.M[,1] == "BD1"] <- "BD2"
            
            # Depending on testing decision/treatment, individuals incur different costs/disutilities,
            # therefore they are grouped into three tracks: 
            # 1.not tested; 2. tested & on primary drug; 3. tested & on alternative drug.
            pop.none      <- df.Pop[df.Pop$tested==FALSE,                            'name']
            pop.test.prim <- df.Pop[df.Pop$tested==TRUE & df.Pop$treat=="Primary",   'name']
            pop.test.alt  <- df.Pop[df.Pop$treat=="Alternate" | is.na(df.Pop$treat), 'name']
            
            # summarize counts by track for the next cycle
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
        
        for(i in 1:3) #depending on track, different drug/testing costs are accrued 
        {
            if(nrow(mm[[i]]) > 0)
            {
                dd        <- (if(i==3) c_alt else c_tx)*365/interval #drug cost
                m0        <- rbind(c(sum(mm[[i]][1,]),rep(0,ncol(mm[[i]]))),mm[[i]]) %>% as.matrix() #add back t=0 row to ensure all groups sum up to n.i in healthy state at t=0 cycle
                m1        <- as.matrix(integrator(diag(exp( 0:n.t * -d.r)) %*% m0, method=method)) #first discount then integrate
                cc[i]     <- sum(as.vector(m1 %*% c(0,if(i==1) c_a+dd else c_a+dd+c_t,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0)))
                ee[i]     <- sum(as.vector(m1 %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)))
                c.test[i] <- if(i==1) 0 else sum(m1[,"A1"]*c_t)
                c.drug[i] <- sum(as.vector(m1 %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0)))
            }
        }
        
        tout  <- map(mm, function(x) mutate(x,cycle=as.integer(row.names(x)))) %>% do.call("bind_rows",.) %>% 
            group_by(cycle) %>% summarise_all(sum) %>% select(-cycle) #combine counts from all three tracks
        mmm <- rbind(c(sum(tout[1,]),rep(0,ncol(tout[1,])-1)),tout) %>% as.matrix() #add back t=0 when all inviduals are in healthy state
        mmm2 <- integrator(mmm,method = method) #non-discounted integration
        dmm <- integrator(diag(exp( 0:n.t * -d.r)) %*% mmm, method=method) #discounted integration
        
        # other metrics 
        possible  <- as.vector(dmm %*% c(1,rep(1,d_at*interval+1),1,1,0,0,0))/interval #discounted life years
        fatal_b   <- sum(mmm2[n.t,c("BD1","BD2")]) #fatal B count
        living    <- n.i - sum(mmm2[n.t,c("BD1","BD2","D")]) #living individuals by the end of simulation
        disutil_a <- as.vector(dmm  %*% c(0,rep(d_a,d_at*interval),0,0,0,0,0,0))/ interval #discounted disutility from A
        disutil_b <- as.vector(dmm  %*% c(0,rep(0,d_at*interval+1),d_b,d_b,0,0,0))/ interval #discounted disutility from B
        c.treat   <- as.vector(dmm  %*% c(0,c_a,rep(0,d_at*interval),c_bs,0,c_bd,0,0)) #discounted adverse event costs

        list(
            # raw_ct = mmm,
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
microsim_icer <- function(params, reference=NULL, genotype=NULL, seed=NULL, method="life-table",...)
{
    if(!is.null(seed)) set.seed(seed)
    
    seed         <- .Random.seed  
    
    params$p_o   <- 0.0 # No testing, reference
    reference    <- microsim_run(params,method=method,...)
    reference    <- reference$results
    
    .Random.seed <- seed
    
    params$p_o   <- 1.0 # Genotype testing upon indication
    genotype     <- microsim_run(params,method=method,...)
    genotype     <- genotype$results
    
    c( ICER       = unname((genotype['dCOST'] - reference['dCOST']) / (genotype['dQALY'] - reference['dQALY'])),
       NMB        = unname((reference['dCOST'] - genotype['dCOST']) + params$wtp*(genotype['dQALY'] - reference['dQALY'])),
       dCOST.ref  = unname(reference['dCOST']),
       dCOST.test = unname(genotype['dCOST']),
       dQALY.ref  = unname(reference['dQALY']),
       dQALY.test = unname(genotype['dQALY'])
    )
}

# microsim_icer(params,N=1000,seed=123) #sample code to run the model
