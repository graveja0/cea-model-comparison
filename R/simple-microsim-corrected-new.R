# library(here)
# library(magrittr)
# library(dplyr)
# library(stringr)
# library(tidyr)
# source(here("R/common.R"))

#### 01 inputs ####
# source(here("R/simple-params.R"))
# source(here("R/simple-markov-corrected.R"))
#Option 1: get corrected transition matrices by running corrected markov
markov_corr_get <- function(params, gene=0, test=0)
{
  params$disc_rate <- inst_rate(1-1/(1 + params$disc), 1)
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

# Option 2: load matrices
# load(here("output/markov_matrices.rda"))
# m00 <- corrected_tp$g0t0
# m11 <- corrected_tp$g1t1


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
ProbsH_corr <- function(M_t, t, gene, test)
{
    p <- matrix(0,nrow=length(M_t),ncol=5,dimnames=list(c(),c("A1","BS1","BD1","D","BSD")))
    #assume behavior parameter only takes 0 or 1, otherwise do not know which matrix to use
    if(test==1) {
        b0 <- !gene
        b1 <- gene
        
        p[b0,"A1"] <- m00[[t]]["H","A"]
        p[b0,"BS1"] <- m00[[t]]["H","BS"]
        p[b0,"BD1"] <- m00[[t]]["H","BD"]
        p[b0,"D"] <- m00[[t]]["H","D"]
        p[b0,"BSD"] <- m00[[t]]["H","BSD"]
        
        p[b1,"A1"] <- m11[[t]]["H","A"]
        p[b1,"BS1"] <- m11[[t]]["H","BS"]
        p[b1,"BD1"] <- m11[[t]]["H","BD"]
        p[b1,"D"] <- m11[[t]]["H","D"]
        p[b1,"BSD"] <- m11[[t]]["H","BSD"]
    } else if(test==0) {
        p[,"A1"] <- m00[[t]]["H","A"]
        p[,"BS1"] <- m00[[t]]["H","BS"]
        p[,"BD1"] <- m00[[t]]["H","BD"]
        p[,"D"] <- m00[[t]]["H","D"]
        p[,"BSD"] <- m00[[t]]["H","BSD"]
    } else stop("invalid p_o")
    return(p)
}

ProbsA_corr <- function(M_t, t, gene, test) {
    p <- matrix(0,nrow=length(M_t),ncol=4,dimnames=list(c(),c("BS1","BD1","D","BSD")))
    if(test==1) {
        b0 <- !gene
        b1 <- gene
        
        p[b0,"BS1"] <- m00[[t]]["A","BS"]
        p[b0,"BD1"] <- m00[[t]]["A","BD"]
        p[b0,"D"] <- m00[[t]]["A","D"]
        p[b0,"BSD"] <- m00[[t]]["A","BSD"]
        
        p[b1,"BS1"] <- m11[[t]]["A","BS"]
        p[b1,"BD1"] <- m11[[t]]["A","BD"]
        p[b1,"D"] <- m11[[t]]["A","D"]
        p[b1,"BSD"] <- m11[[t]]["A","BSD"]
    } else if(test==0) {
        p[,"BS1"] <- m00[[t]]["A","BS"]
        p[,"BD1"] <- m00[[t]]["A","BD"]
        p[,"D"] <- m00[[t]]["A","D"]
        p[,"BSD"] <- m00[[t]]["A","BSD"]
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


microsim_run_corr <- function(params, N = NULL, method="beginning")
{
    if (!is.null(N)) params$n <- N
    
    # TO DEBUG USE
    # sapply(names(params), function(p) assign(p, params[[p]], 1))
    with(params, {
        
        #### 1. Model Setting ####
        n.i       <- n                          # number of individuals
        n.t       <- horizon*interval           # number of cycles
        v.n       <- c("H",paste0("A",1:(d_at*interval+1)),"BS1","BS2","BD1","BD2","D","BSD")    # state names
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
                pH <- ProbsH_corr(m.M[,1][fromH],t,df.Pop$variant[fromH],p_o)
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
                pA <- ProbsA_corr(m.M[,1][fromA],t,df.Pop$variant[fromA],p_o)
                
                # move forward in A tunnel as needed
                drawA <- draw.events(m.M[,1][fromA],pA)
                intunnel <- (substr(drawA,1,1)=="A")
                ttunnel <- as.numeric(substring(drawA[intunnel],2)) # which tunnel cycle
                if(any(ttunnel<=d_at*interval)) {
                    drawA[intunnel][ttunnel<=d_at*interval] <- paste0("A",ttunnel[ttunnel<=d_at*interval]+1)
                }
                m.M[,2][fromA] <- drawA
                #a.P[fromA,c("BS1","BD1","D"),t] <- pA
            }
            
            # From BS 
            fromBS <- m.M[,1] %in% c("BS1","BS2") 
            if(sum(fromBS)>0)
            {
                pBS <- ProbsBS_corr(m.M[,1][fromBS],t,df.Pop$variant[fromBS],p_o)
                drawBS <- rbinom(n=length(pBS),size=1,prob=pBS)
                m.M[,2][fromBS][drawBS==1] <- "D"
                m.M[,2][fromBS][drawBS==0] <- "BS2"
                
                #a.P[fromBS,"D",t] <- pBS
            }
            
            # From BD
            m.M[,2][m.M[,1] == "BD1"] <- "BD2"
            m.M[,2][m.M[,1] == "BSD"] <- "D"
            
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
        
        
        ### compute C/E
        mm <- map(1:3, function(x) map_df(1:n.t,~ct[[.x]][[x]]))
        
        cc            <- numeric(3)
        ee            <- numeric(3)
        c.test        <- numeric(3)
        c.drug        <- numeric(3)
        d.r           <- (1 + disc)^(1/interval)-1
        v.dwc         <- 1 / (1 + d.r) ^ (0:(n.t-1)) # calculate discount weights for costs for each cycle based on discount rate d.r
        v.dwe         <- v.dwc  # Use same discount for qaly

        
        for(i in 1:3)
        {
          if(nrow(mm[[i]]) > 0)
          {
                dd        <- (if(i==3) c_alt else c_tx)*365/interval
                m0 <- rbind(c(sum(mm[[i]][1,]),rep(0,ncol(mm[[i]]))),mm[[i]])
                m1        <- as.matrix(integrator(m0,method=method))
                cc[i]     <- sum(as.vector(m1 %*% c(0,if(i==1) c_a+dd else c_a+dd+c_t,rep(dd,d_at*interval),c_bs+dd,dd,c_bd,0,0,c_bs)) * v.dwc)
                ee[i]     <- sum(as.vector(m1 %*% c(1/interval,rep((1-d_a)/interval,d_at*interval),1/interval,(1-d_b)/interval,(1-d_b)/interval,0,0,0,0)) * v.dwe)
                c.test[i] <- if(i==1) 0 else sum(m1[,"A1"]*c_t*v.dwc)
                c.drug[i] <- sum(as.vector(m1 %*% c(0,dd,rep(dd,d_at*interval),dd,dd,0,0,0,0)) * v.dwc)
            }
        }
        
        tout  <- map(mm, function(x) mutate(x,cycle=as.integer(row.names(x)))) %>% do.call("bind_rows",.) %>% 
          group_by(cycle) %>% summarise_all(sum) %>% select(-cycle)
        mmm <- rbind(c(sum(tout[1,]),rep(0,ncol(tout[1,])-1)),tout) %>% 
          integrator(method=method) %>% as.matrix()
        
        #other
possible  <- mmm %*% c(1,rep(1,d_at*interval+1),1,1,0,0,0,0)
        possible  <- as.vector(possible) * v.dwe
        fatal_b   <- sum(mmm[n.t,c("BD1","BD2")])
        living    <- n.i - sum(mmm[n.t,c("BD1","BD2","D","BSD")])
        disutil_a <- mmm  %*% c(0,rep(d_a,d_at*interval),0,0,0,0,0,0,0)
        disutil_a <- as.vector(disutil_a) * v.dwe
        disutil_b <- mmm  %*% c(0,rep(0,d_at*interval+1),d_b,d_b,0,0,0,0)
        disutil_b <- as.vector(disutil_b) * v.dwe
        c.treat   <- mmm  %*% c(0,c_a,rep(0,d_at*interval),c_bs,0,c_bd,0,0,c_bs)
        c.treat   <- as.vector(c.treat) * v.dwc

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


microsim_icer_corr <- function(params, reference=NULL, genotype=NULL, method="life-table",seed=123,...)
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

