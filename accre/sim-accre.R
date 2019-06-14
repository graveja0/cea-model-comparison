  #############################################################################
 ##
## ACCRE batch run
##
## Assumes an "output" and "status" directory underneath the directory
## of this file

  #############################################################################
 ##
## Must define function simulation(x)
source("run.R")

  #############################################################################
 ##
## Append ACCRE local directory of R libraries
.libPaths(c("~/R/rlib-3.4.0", .libPaths()))

  #############################################################################
 ##
## Grab ACCRE batch number for random seed (makes job repeatible on failure)
args <- commandArgs(trailingOnly=TRUE)
x    <- as.numeric(args[1])
set.seed(x)

  #############################################################################
 ##
## Run the job and save to output directory
progress(paste0("Batch ", x, " starting."))
results <- simulation(10000)
save(results, file=paste0("output/run-", x, ".RData"))
progress(paste0("Batch ", x, " complete."))
