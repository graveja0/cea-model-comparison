# This script can run DES and microsimulation in parallel to generate results for the convergence plot.
# In our manuscript, we ran each model was run 1k times for each N.
# For the microsimulation model, we ran a 1000-person model 1,000,000 times to get a pool of 1b and sample from it to get 1k estimates for each N.
library(parallel)
library(doParallel)

####Iterate DES runs####
# How many batches/sessions to run in parallel
## Grab batch number for random seed (makes job repeatible)
x <- 1:5

nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
clusterExport(cl, c("n"))

icer_des_N <-
    foreach(i = icount(length(x)), .combine=bind_rows,            
            .packages = c("tcltk",#for progress bar (optional)
                          "here","dplyr","purrr")
    ) %dopar% {
        #set up each sessionl
        source(here("common.R")) #load shared functions
        source(here("simple-params.R")) #load inputs
        source(here("simple-des.R")) 
        
        #progress bar by session (optional)
        if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=length(x))
        setTkProgressBar(pb, i)
        
       10^seq(3,6,by=0.25) %>% 
                 map_df(function(s) des_icer(params,seed=x[i],N=s) %>% t() %>% data.frame() %>% mutate(N = s))
    }
stopCluster(cl)

####Sample from microsim results####
# run microsimulation models a few times to get a pool of 1000-person cohort estimates to sample from
if(!exists("icer_micro_y") | !exists("icer_micro_corr_y")) source(here("parallel_microsim"))

# sample (with replacement) from the pool and get NMB estimates for target sample size
get_microsim_nmb <- function(df, n, seed = 123, wtp = params$wtp) {
    set.seed(seed)
    sampled_rows <- sample(1:nrow(df),n/1e3, replace = TRUE)
    tmp <- apply(df[sampled_rows,] ,2, mean)
    nmb <- unname((tmp['dCOST.ref'] - tmp['dCOST.test']) + wtp*(tmp['dQALY.test'] - tmp['dQALY.ref']))
    nmb.test <- wtp*(tmp['dQALY.test']) -  (tmp['dCOST.test'])
    nmb.ref <- wtp*(tmp['dQALY.ref']) -  (tmp['dCOST.ref'])
    icer <- unname((tmp['dCOST.test'] - tmp['dCOST.ref']) / (tmp['dQALY.test'] - tmp['dQALY.ref']))
    return(c(nmb))
}

# How many batches/sessions to run in parallel
## Grab batch number for random seed (makes job repeatible)
x <- 1:5
micro_sample_sizes = 10^seq(4,6,0.25)

nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
clusterExport(cl, c("n"))

icer_micro_N <-
    foreach(i = icount(length(x)), .combine=bind_rows,            
            .packages = c("tcltk",#for progress bar (optional)
                          "here","dplyr","purrr")
    ) %dopar% {
        
        if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=length(x))
        setTkProgressBar(pb, i)
        
        data.frame(N=micro_sample_sizes,
            NMB=map_dbl(micro_sample_sizes,~get_microsim_nmb(df=icer_micro_y,n = .x)))
    }
stopCluster(cl)

cl <- makeCluster(nc)
registerDoParallel(cl)
clusterExport(cl, c("n"))
icer_micro_corr_N <-
    foreach(i = icount(length(x)), .combine=bind_rows,            
            .packages = c("tcltk",#for progress bar (optional)
                          "here","dplyr","purrr")
    ) %dopar% {
        
        if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=length(x))
        setTkProgressBar(pb, i)
        
        data.frame(N=micro_sample_sizes,
                   NMB=map_dbl(micro_sample_sizes,~get_microsim_nmb(df=icer_micro_corr_y,n = .x)))
    }
stopCluster(cl)
