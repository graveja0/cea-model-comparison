# Run microsimulation model results

#library(tidyverse)
library(tidyverse)
library(survival)
library(here)
library(tictoc)
library(Matrix)
library(deSolve)
#library(eha)
library(flexsurv)
source(here("R/common.R"))
source(here("R/simple-params.R"))
source(here("R/simple-des.R"))
source(here("R/simple-deq.R"))
source(here("R/simple-markov.R"))
source(here("R/simple-markov-corrected.R"))
source(here("R/simple-microsim-corrected-new.R"))
source(here("R/simple-microsim-new.R"))

params_orig <- params

# Differntial Equations
icer_deq <- deq_icer(params)

# Embedded Markov
icer_mar_corr <- markov_comb_corr(params)

# Uncorrected Markov 
params$interval = 1
icer_mar_1y <- markov_comb(params) 
params$interval = 12
icer_mar_1m <- markov_comb(params) 
params$interval = 365
icer_mar_1d <- markov_comb(params) 

# Discrete Event Simulation
icer_des <- des_icer(params,N = 1e7, seed = 123)

# Microsimulation 
load("~/box/microsim_runs/microsim_1b_daily.rda")
icer_micro_d <- runs[1:(1e7/1e3),] %>% apply(2,mean)

load("~/box/microsim_runs/microsim_1b_yearly.rda")
icer_micro_y <- runs[1:(1e7/1e3),] %>% apply(2,mean)

load("~/box/microsim_runs/microsim_corr_yr_nov_1b.rda")
icer_micro_corr_y <- runs[1:(1e7/1e3),] %>% apply(2,mean)


# 1-Day Intervals for state transition models
params$interval = 365
icer_mar_1d <- markov_comb(params) 
params$interval = 12
icer_mar_1m <- markov_comb(params)


tic()
icer_microsim_1m <- microsim_icer(params, N = 1e7, seed=123)
toc()
icer_mar_corr_1d <- markov_comb_corr(params) 
#icer_microsim_corr_1d <- microsim_icer(params, N = 2e6, seed=123)

names(icer_mar) <- names(icer_mar_corr) <- names(icer_des) <- names(icer_microsim) <- names(icer_deq)
icer_results <- 
  list("DEQ" = icer_deq,
       "MRKCHT_UADJ" = icer_mar,
       "MRKCHT_EMB" = icer_mar_corr,
       "DES" = icer_des,
       "MICRO" = icer_microsim)
icer_results %>% 
  write_rds(here("output/icer-results.R"))


