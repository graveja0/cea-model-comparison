# This script obtains runtime estimates for each model.
# Please note that some model settings take a long time to run.
# Due to computational burden, we ran everything on ACCRE.
library(microbenchmark)

rtimes <- 2L
rn <- 1000

print("deq")
b_deq <- microbenchmark("deq" = { icer_deq <- deq_icer(params) },times=rtimes)

print("des")
b_des <- microbenchmark("des" = { icer_des <- des_icer(params,N = rn) },times=rtimes)

print("mark_y")
params$interval <- 1
b_mark_y <- microbenchmark("markov" = { icer_mar1 <- markov_icer(params) },times=rtimes)

print("mark_m")
params$interval = 12
b_mark_m <- microbenchmark("markov1m" = { icer_mar2 <- markov_icer(params) },times=rtimes)

print("mark_d")
params$interval = 365
b_mark_d <- microbenchmark("markov1d" = { icer_mar3 <- markov_icer(params) },times=rtimes)

print("mark_corr_y")
params$interval <- 1
b_mark_corr_y <- microbenchmark("markov_corr" = { icer_mar_corr <- markov_corr_icer(params) },times=rtimes)

print("mark_corr_m")
params$interval <- 12
b_mark_corr_m <- microbenchmark("markov_corr_m" = { icer_mar_corr_m <- markov_corr_m_icer(params) },times=rtimes)

print("micro_y")
b_micro_y <- microbenchmark("microsim" = {icer_microsim1 <- microsim_icer(params, N = rn)},times=rtimes)

print("micro_corr_y")
b_micro_corr_y <- microbenchmark("microsim_corr" = {icer_microsim2 <- microsim_corr_icer(params, N = rn)},times=rtimes)

rt <- b_deq %>% mutate(model="DEQ") %>% 
    bind_rows(b_des %>% mutate(model="DES (10mil)")) %>% 
    bind_rows(b_mark_y %>% mutate(model="MARKOV-Y")) %>% 
    bind_rows(b_mark_m %>% mutate(model="MARKOV-M")) %>% 
    bind_rows(b_mark_d %>% mutate(model="MARKOV-D")) %>% 
    bind_rows(b_mark_corr_y %>% mutate(model="MARKOV-EMB-Y")) %>% 
    bind_rows(b_mark_corr_m %>% mutate(model="MARKOV-EMB-M")) %>%     
    bind_rows(b_micro_y %>% mutate(model="MICROSIM-Y (10 mil)")) %>% 
    bind_rows(b_micro_corr_y %>% mutate(model="MICROSIM-EMB-Y (10 mil)")) %>% 
    mutate(times=time/10^9)
