# Load the following packages and scripts to run this model independently
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(here)
# source(here("R/common.R")) #load shared functions
# source(here("R/simple-params.R")) #load inputs

#### 01 inputs ####
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


# Event draw function:
# looks up probability of transition
ProbsH <- function(M_t, t, p.a, p.sd, interval)
{
    # M_t  :    current health state 
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
        
        #secular death risk
        p.HD      <- gompertz_ratio2(1:n.t, interval, shape, rate)
        
        #### 2. Sample individual level characteristics ####
        #Static characteristics
        df.Pop    <- data.frame(name    = 1:n.i,
                                variant = sample(c(FALSE,TRUE), n.i, prob=c(1-p_g,p_g), replace = TRUE),
                                tested  = FALSE,
                                treat   = factor("",levels=c("","Alternate","Primary")))
        
        #Dynamic characteristics 
        v.M_Init <- rep("H", n.i)       # everyone begins in the healthy state 
        
        #### 3. Start Simultation ####
        # m.M: health state for each patient at each cycle
        m.M  <- matrix(nrow = n.i, ncol = 2, 
                       dimnames = list(paste(1:n.i),     # could name the rows ind1, ind2, ind3, etc.
                                       c("current","forward")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.
        
        m.M[, 1] <- v.M_Init         # initial health state for all individuals 
        
        zerocol <- as.data.frame(setNames(replicate(length(v.n),numeric(0), simplify = F),v.n))
        zerocol[1,] <- 0
        
        ct <- vector("list",n.t)
        
        # Time for-loop from cycle 1 through n.t
        for (t in 1:n.t)
        {
            ## Simulate Events ##
            # Default Update: if no events occur, health state remains unchanged
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
            
            #summarize counts at each cycle
            pop.none      <- df.Pop[df.Pop$tested==FALSE,                            'name']
            pop.test.prim <- df.Pop[df.Pop$tested==TRUE & df.Pop$treat=="Primary",   'name']
            pop.test.alt  <- df.Pop[df.Pop$treat=="Alternate" | is.na(df.Pop$treat), 'name']
            
            ct[[t]] <- list(pop.none, pop.test.prim, pop.test.alt) %>% 
                map(function(x) {
                    if(length(x)==0) {return(zerocol)} else {
                        m.M[x,"forward"] %>% as.data.frame() %>% set_names("state") %>%
                            mutate(state=factor(state,levels=v.n)) %>% 
                            table() %>% as.data.frame() %>% spread_(".","Freq")
                    }
                })
            
            ### reset matrix
            m.M[,1] <- m.M[,2]
            
        }        
        
        mm <- map(1:3, function(x) map_df(1:n.t,~ct[[.x]][[x]]))
        
        cc            <- numeric(3)
        ee            <- numeric(3)
        c.test        <- numeric(3)
        c.drug        <- numeric(3)
        d.r  <- inst_rate(1-1/(1 + disc), interval)
        
        for(i in 1:3)
        {
            if(nrow(mm[[i]]) > 0)
            {
                dd        <- (if(i==3) c_alt else c_tx)*365/interval
                m0        <- rbind(c(sum(mm[[i]][1,]),rep(0,ncol(mm[[i]]))),mm[[i]]) %>% as.matrix() #add back t=0 row (initial states)
                m1        <- as.matrix(integrator(diag(exp( 0:n.t * -d.r)) %*% m0, method=method)) #first discount then integrate
                cc[i]     <- sum(as.vector(m1 %*% c(0,if(i==1) c_a+dd else c_a+dd+c_t,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0)))
                ee[i]     <- sum(as.vector(m1 %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0)))
                c.test[i] <- if(i==1) 0 else sum(m1[,"A1"]*c_t)
                c.drug[i] <- sum(as.vector(m1 %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0)))
            }
        }
        
        tout  <- map(mm, function(x) mutate(x,cycle=as.integer(row.names(x)))) %>% do.call("bind_rows",.) %>% 
            group_by(cycle) %>% summarise_all(sum) %>% select(-cycle)
        mmm <- rbind(c(sum(tout[1,]),rep(0,ncol(tout[1,])-1)),tout) %>% as.matrix()
        mmm2 <- integrator(mmm,method = method)
        dmm <- integrator(diag(exp( 0:n.t * -d.r)) %*% mmm, method=method)
            
        possible  <- as.vector(dmm %*% c(1,rep(1,d_at*interval+1),1,1,0,0,0))/interval
        fatal_b   <- sum(mmm2[n.t,c("BD1","BD2")])
        living    <- n.i - sum(mmm2[n.t,c("BD1","BD2","D")])
        disutil_a <- as.vector(dmm  %*% c(0,rep(d_a,d_at*interval),0,0,0,0,0,0))/ interval
        disutil_b <- as.vector(dmm  %*% c(0,rep(0,d_at*interval+1),d_b,d_b,0,0,0))/ interval
        c.treat   <- as.vector(dmm  %*% c(0,c_a,rep(0,d_at*interval),c_bs,0,c_bd,0,0))

        list(
            raw_ct = mm,
            count = mmm,
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

# microsim_icer(params,N=1000) #sample code to run the model
