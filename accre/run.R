source("../R/simple-params.R")
source("../R/simple-microsim.R")

params$n <- 100

progress <- function(...)
{
  cat(date(), ' ')
  invisible(lapply(list(...), function(x) cat(x,'\n')))
}


simulation <- function(n) do.call("rbind", lapply(1:n, function(y) {
  if(y %% 100 == 1) progress(paste0("Starting run ",y))
  microsim_icer(params, seed=NULL)
}))