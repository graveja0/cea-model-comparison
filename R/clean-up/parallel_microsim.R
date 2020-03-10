# This script can run micro simulation model in parallel and then average results to form an estimate.
# In our manuscript, we ran a 1000-person model 100,000 times to get an estimate for 10 million individuals on ACCERE.
library(parallel)
library(doParallel)

# How many batches/sessions to run in parallel
## Grab batch number for random seed (makes job repeatible)
x <- 1:5

nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
clusterExport(cl, c("n"))

icer_micro_y <-
    foreach(i = icount(length(x)), .combine=bind_rows,
            .packages = c("tcltk",#for progress bar (optional)
                          "here","dplyr","tidyr","purrr")
            ) %dopar% {
        #set up each session
        source(here("run_microsim.R"))
        
        #progress bar by session (optional)
        if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=length(x))
        setTkProgressBar(pb, i)
        
        #set random seed and run simulation multiple times in sequence    
        set.seed(x[i])
        simulation(10) %>% #here specify how many times to run in each session
            as.data.frame() #return one df from each session then to be combined using bind_rows()
    }
stopCluster(cl)

icer_micro_y %>% summarise_all(mean) #scale up small samples by averaging results


nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
clusterExport(cl, c("n"))

icer_micro_corr_y <-
    foreach(i = icount(length(x)), .combine=bind_rows,
            .packages = c("tcltk",#for progress bar (optional)
                          "here","dplyr","tidyr","purrr","Matrix")
    ) %dopar% {
        #set up each session
        source(here("run_microsim_corr.R"))
        
        #progress bar by session (optional)
        if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=length(x))
        setTkProgressBar(pb, i)
        
        #set random seed and run simulation multiple times in sequence    
        set.seed(x[i])
        simulation(10) %>% #here specify how many times to run in each session
            as.data.frame() #return one df from each session then to be combined using bind_rows()
    }
stopCluster(cl)

icer_micro_corr_y %>% summarise_all(mean) #scale up small samples by averaging results
