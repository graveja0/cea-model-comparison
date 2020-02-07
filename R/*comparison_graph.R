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

icer_des <- read_rds("~/box/microsim_runs/convergence/icer_des_1k.rds")

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

icer_micro_corr_m <- read_rds("~/Box Sync/microsim_runs/convergence/icer_micro_corr_m_e6_11.rds")

#### p1 ####
lines <- data.frame(model=c("DEQ","MARKOV","MARKOV-EMB"),
                    NMB=c(icer_deq["NMB"],icer_mar["NMB"],icer_mar_corr["NMB"]),stringsAsFactors = F) 

dots <-   icer_des %>% mutate(model="DES") %>% 
    bind_rows(icer_micro_y %>% mutate(model="Microsim Yearly")) %>% 
    bind_rows(icer_micro_d %>% mutate(model="Microsim Daily")) %>% 
    bind_rows(icer_micro_corr_y %>% mutate(model="Microsim Yearly Corrected")) %>%
    bind_rows(icer_micro_corr_m %>% mutate(model="Microsim Monthly Corrected")) %>% 
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
icer_micro_d$timestep <- "MICROSIM-D"
icer_micro_y$timestep <- "MICROSIM-Y"
icer_micro_corr_y$timestep <- "MICROSIM-EMB-Y"
icer_micro_corr_m$timestep <- "MICROSIM-EMB-M"
icer_des$timestep <- "DES"

micro <- rbind(icer_micro_d, icer_micro_y,icer_micro_corr_y)
micro <- rbind(micro, icer_des[,c("NMB", "N", "seed", "iteration", "timestep")])
micro$logN <- factor(round(log(micro$N,base = 10)))

ll <- lines %>% mutate(timestep=model)
micro$timestep <- factor(micro$timestep, 
                         levels = c('DES', 'MICROSIM-Y','MICROSIM-D','MICROSIM-EMB-Y','MICROSIM-EMB-M'))


tiff("Micro Simulation Convergence.tiff",res=300,units="in",width=6.5, height=5)
# Micro Simulation Convergence
library(directlabels)
p <- 
ggplot(data= micro %>% filter(N>1e3), aes(x=logN,y=NMB,color=timestep,fill=timestep),alpha=0.7) +
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(name="",values = c("grey40","red3","orange2","blue","skyblue")) +
    scale_fill_manual(name="",values = c("grey80","tomato","gold1","dodgerblue","paleturquoise1")) +
    new_scale_color() +
    geom_hline(data=ll,aes(yintercept=NMB,color=timestep,linetype=timestep)) + 
    scale_linetype_manual(name="",values=c("solid","dashed","dotdash")) + 
<<<<<<< HEAD
    scale_color_manual(name="",values = c("black","darkgreen","purple")) + 
    scale_size_manual(name="",values = c(0.25,0.25,2)) +
    xlab("Number of Simulated Patients\n(log10)") + 
    theme_bw() + 
    scale_x_discrete(breaks = 2:11, labels = c("2\n100","3\n1k","4\n10k","5\n100k","6\n1mil","7\n10mil","8\n100mil","9\n1bil","10\n10bil","11\n100bil")) + 
    ylim(c(-1100,1100))
p
=======
    scale_color_manual(name="",values = c("black","darkgreen","green")) + 
    scale_size_manual(name="",values = c(0.5,0.5,0.8)) +
    xlab("Number of Simulated Patients\n(log)") + 
    theme_bw() +
    theme(panel.grid.minor.y=element_blank())
dev.off()
>>>>>>> 6fddfa946fb04ebb5b5476b85f852596057c4cb7


#truncation convergence
icer_des_trunc <- read_rds("~/Box Sync/microsim_runs/convergence/icer_des_trunc_e2_8.rds")
icer_des_trunc$timestep <- "DES truncated"
icer_des$timestep <- "DES"

trunc <- bind_rows(icer_des,icer_des_trunc)
trunc$logN <- factor(round(log10(trunc$N)))

tiff("DES Truncation Simulation Convergence.tiff",res=300,units="in",width=7, height=5)
ggplot(trunc, aes(x=logN,y=NMB,color=timestep,fill=timestep),alpha=0.7) +
    geom_boxplot(coef=1e30,size=0.2) +
    scale_color_manual(name="",values = c("blue","red3")) +
    scale_fill_manual(name="",values = c("lightskyblue","pink")) +
    new_scale_color() +
    geom_hline(data=ll[ll$model=="DEQ",],aes(yintercept=NMB,color=timestep,linetype=timestep),size=0.2) + 
    scale_linetype_manual(name="",values=c("solid")) + 
    scale_color_manual(name="",values = c("black")) + 
    xlab("Number of Simulated Patients\n(log10)") + 
    theme_bw() +
    scale_x_discrete(breaks = 2:8, labels = c("2\n100","3\n1k","4\n10k","5\n100k","6\n1mil","7\n10mil","8\n100mil")) + 
    theme(panel.grid=element_blank())

dev.off()



#### runtime ####
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_deq.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_des.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_y.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_m.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_d.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_m.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_corr_y.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_mark_corr_m.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_micro_y.rda")
load("~/Box Sync/microsim_runs/run_time_1e7/run_t_micro_corr_yr_nov.rda")

rt <- b_deq %>% mutate(model="DEQ") %>% 
<<<<<<< HEAD
    bind_rows(b_des %>% mutate(model="DES (10mil)")) %>% 
    bind_rows(b_mark_y %>% mutate(model="MARKOV-Y")) %>% 
    bind_rows(b_mark_m %>% mutate(model="MARKOV-M")) %>% 
    bind_rows(b_mark_d %>% mutate(model="MARKOV-D")) %>% 
    bind_rows(b_mark_corr_y %>% mutate(model="MARKOV-EMB-Y")) %>% 
    bind_rows(b_mark_corr_m %>% mutate(model="MARKOV-EMB-M")) %>%     
    bind_rows(b_micro_y %>% mutate(model="MICROSIM-Y (10 mil)")) %>% 
    bind_rows(b_micro_corr_y %>% mutate(model="MICROSIM-EMB-Y (10 mil)")) %>% 
    mutate(times=time/10^9)
    
rt %>% group_by(model,expr) %>% 
    summarise(med = quantile(times,0.5,na.rm=TRUE)) %>% 
    arrange(desc(med))

value_labs <- paste(seq(-3,4,0.5),paste0("\n"),signif(10^seq(-3,4,0.5),2),"sec")

ggplot(rt) +
    geom_boxplot(aes(x=reorder(model,-times),y=times),outlier.shape = 1,outlier.alpha = 0.5) +
    scale_y_continuous(breaks = 10^seq(-3,4,0.5), 
                       labels = value_labs, trans="log10") + ylab("log10(seconds)") + xlab("Method") +
    coord_flip() + theme_bw()  + xlab("")
=======
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
>>>>>>> 6fddfa946fb04ebb5b5476b85f852596057c4cb7


