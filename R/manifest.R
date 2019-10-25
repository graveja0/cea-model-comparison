# 
# library(tidyverse)
# library(directlabels)
# library(patchwork)
# library(rlang)
# library(randtoolbox)
# library(microbenchmark)
# library(progress)
# library(here)
# library(fs)
# library(deSolve)
# library(heemod)
# library(dplyr)
# library(purrr)
# library(flexsurv)
# library(furrr)
# library(conflicted)
# library(tangram)
# library(furrr)
# library(microbenchmark)
# library(kableExtra)
# library(gt)
# source(here("R/plot-theme.R"))
# # conflicted::conflict_scout()
# # heemood and purrr: modify
# conflict_prefer("filter", "dplyr")
# conflict_prefer("set_names", "magrittr")
# 

library(here)
source(here("R/common.R"))
source(here("R/simple-params.R"))
source(here("R/simple-des.R"))
source(here("R/simple-deq.R"))
#source(here("R/simple-markov-heemod.R")) 
source(here("R/simple-markov.R"))
source(here("R/simple-markov-corrected.R"))
source(here("R/simple-microsim.R"))


get_params <-  function(x,params) {
  map2(params,x,~(
    .x(.y)
  ))
} 

#this function add calculated parameters to param list from halton draw
get_params2 <-  function(x,params) {
  c(x,map(params,~(
    .x(x)
  )))
} 

dist_lut <- c("qnorm" = "Normal",
              "qunif" = "Uniform",
              "qnorm01" = "Normal",
              "qbeta" = "Beta",
              "qpois" = "Poisson")

get_psa_values <- function(x,params) {
  
  foo <- deparse(x, width  = 200)
  parameter_name <- str_split(gsub("list\\(","",foo[1])," = ")[[1]][1]
  dist = dist_lut[str_split(foo[2],"\\(")[[1]][1]]
  
  gsub(" ","",str_split(str_split(foo[2],"\\(|\\)")[[1]][2],",")[[1]]) %>% data.frame() %>% 
    set_names("parameter") %>% 
    separate(parameter,"=",into = c("parameter","value")) %>% 
    mutate(from_params_list = as.integer(grepl("params",value))) %>% 
    mutate(value = gsub("params\\$","",value))  %>% 
    mutate(rowname = parameter_name) %>% 
    filter(value != "x")  %>% 
    mutate(value = ifelse(from_params_list==1, gsub(rowname,params[rowname],value),value)) %>% 
    select(rowname,parameter,value) %>% 
    spread(parameter,value) %>% 
    mutate(distribution = dist) %>% 
    select(-rowname) %>% 
    select(distribution,everything()) %>% 
    tbl_df() %>% 
    mutate_at(vars(-distribution),list(~paste0(unlist(.)))) %>% 
    mutate_at(vars(-distribution),list(~eval(parse(text=.))))
}
