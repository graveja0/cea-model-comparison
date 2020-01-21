library(tidyverse)
library(here)
library(deSolve)
library(dplyr)
library(purrr)
library(flexsurv)
library(furrr)
library(conflicted)
library(ggnewscale)
# conflicted::conflict_scout()
# heemood and purrr: modify
conflict_prefer("filter", "dplyr")
conflict_prefer("set_names", "magrittr")
source(here::here("R/manifest.R"))

#### Convergence ####
icer_deq <- deq_icer(params)
icer_mar <- markov_comb(params) 
icer_mar_corr <- markov_comb_corr(params)

icer_des <- read_rds("~/Box Sync/microsim_runs/convergence/icer_des_1k.rds")

load("~/Box Sync/microsim_runs/convergence/icer_micro_N_e6_9.rds")
icer_micro_d <- icer_micro1d_N
icer_micro_y <- icer_micro_N
load("~/Box Sync/microsim_runs/convergence/icer_micro_N_e9_10.rds")
icer_micro_d <- icer_micro_d %>% bind_rows(icer_micro1d_N)
icer_micro_y <- icer_micro_y %>% bind_rows(icer_micro_N)
load("~/Box Sync/microsim_runs/convergence/icer_micro_N_e10_11.rds")
icer_micro_d <- icer_micro_d %>% bind_rows(icer_micro1d_N)
icer_micro_y <- icer_micro_y %>% bind_rows(icer_micro_N)

icer_micro_corr_y <- read_rds("~/Box Sync/microsim_runs/convergence/icer_micro_corr_y_e6_9_nov.rds") %>% 
    bind_rows(read_rds("~/Box Sync/microsim_runs/convergence/icer_micro_corr_y_e9_10_nov.rds")) %>% 
    bind_rows(read_rds("~/Box Sync/microsim_runs/convergence/icer_micro_corr_y_e10_11_nov.rds"))


#### p1 ####
lines <- data.frame(model=c("DEQ","Markov","Markov Corrected"),
                    NMB=c(icer_deq["NMB"],icer_mar["NMB"],icer_mar_corr["NMB"]),stringsAsFactors = F) 

dots <-   icer_des %>% mutate(model="DES") %>% 
    bind_rows(icer_micro_y %>% mutate(model="Microsim Yearly")) %>% 
    bind_rows(icer_micro_d %>% mutate(model="Microsim Daily")) %>% 
    bind_rows(icer_micro_corr_y %>% mutate(model="Microsim Yearly Corrected")) %>% 
    group_by(model,N) %>% slice(1) %>% summarise(NMB=mean(NMB))


ggplot() +
    geom_hline(data=lines,aes(yintercept=NMB,lty=model,colour=model,size=model)) +
    scale_linetype_manual(name="",values=c("solid","dashed","dotdash")) + 
    scale_color_manual(name="",values = c("black","darkgreen","green")) + 
    scale_size_manual(name="",values = c(0.5,0.5,0.8)) +
    new_scale_color() +
    geom_point(data=dots,aes(x=N,y=NMB,shape=model,colour=model,alpha=model)) +
    scale_color_manual(name="",values = c("tan4","orange2","red3","blue")) +
    scale_shape_manual(name="",values =  c(20,10,4,18)) +
    scale_alpha_manual(name="",values = c(1,1,1,0.7)) +
    theme_bw() +
    ylab("NMB\n(Genotype vs. Reference)") + 
    xlab("Number of Simulated Patients\n(log10)") + 
    scale_x_continuous(breaks = 10^seq(2,12,1), labels = seq(2,12,1), trans="log10") +
    theme(panel.grid.minor.y = element_blank())


#### p2 ####
icer_micro_d$timestep <- "Microsim Daily"
icer_micro_y$timestep <- "Microsim Yearly"
icer_micro_corr_y$timestep <- "Microsim Yearly Corrected"
icer_des$timestep <- "DES"

micro <- rbind(icer_micro_d, icer_micro_y,icer_micro_corr_y)
micro <- rbind(micro, icer_des[,c("NMB", "N", "seed", "iteration", "timestep")])
micro$logN <- factor(round(log(micro$N)))

ll <- lines %>% mutate(timestep=model)

tiff("Micro Simulation Convergence.tiff",res=300,units="in",width=6.5, height=5)
# Micro Simulation Convergence
ggplot(micro, aes(x=logN,y=NMB,color=timestep,fill=timestep),alpha=0.7) +
    geom_boxplot(coef=1e30) +
    scale_color_manual(name="",values = c("grey40","orange2","red3","blue")) +
    scale_fill_manual(name="",values = c("grey80","gold1","tomato","skyblue")) +
    new_scale_color() +
    geom_hline(data=ll,aes(yintercept=NMB,color=timestep,linetype=timestep)) + 
    scale_linetype_manual(name="",values=c("solid","dashed","dotdash")) + 
    scale_color_manual(name="",values = c("black","darkgreen","green")) + 
    scale_size_manual(name="",values = c(0.5,0.5,0.8)) +
    xlab("Number of Simulated Patients\n(log)") + 
    theme_bw() +
    theme(panel.grid.minor.y=element_blank())
dev.off()

#### runtime ####
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_deq.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_des.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_y.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_m.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_d.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_corr_y.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_corr_m.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_micro_y.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_micro_corr_yr_nov.rda")

rt <- b_deq %>% mutate(model="DEQ") %>% 
    bind_rows(b_des %>% mutate(model="DES")) %>% 
    bind_rows(b_mark_y %>% mutate(model="Markov Yearly")) %>% 
    bind_rows(b_mark_m %>% mutate(model="Markov Monthly")) %>% 
    bind_rows(b_mark_d %>% mutate(model="Markov Daily")) %>% 
    bind_rows(b_mark_corr_y %>% mutate(model="Markov Yearly Corrected")) %>%
    bind_rows(b_mark_corr_m %>% mutate(model="Markov Yearly Corrected")) %>%
    bind_rows(b_micro_y %>% mutate(model="Microsim Yearly")) %>% 
    bind_rows(b_micro_corr_y %>% mutate(model="Microsim Yearly Corrected")) %>% 
    mutate(times=time/10^9)
    
tiff("Runtime Comparison.tiff",res=300,units="in",width=6, height=5)
ggplot(rt) +
    geom_boxplot(aes(x=model,y=times),outlier.shape = 1,outlier.alpha = 0.5) +
    scale_y_continuous(breaks = 10^seq(-3,4,0.5), labels = seq(-3,4,0.5), trans="log10") + ylab("log10(seconds)") + xlab("Method") +
    coord_flip() + theme_bw() +
    theme(panel.grid.major.y=element_blank())
dev.off()


