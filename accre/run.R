source("../R/simple-params.R")
source("../R/simple-microsim.R")

params$n <- 100

simulation <- function(n) do.call("rbind", lapply(1:n, function(y) microsim_icer(params, seed=NULL)))