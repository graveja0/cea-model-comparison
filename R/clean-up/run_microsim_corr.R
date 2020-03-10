# To run microsim in 1000/batch
source(here("common.R")) #load shared functions
source(here("simple-params.R")) #load inputs
source(here("simple-microsim-corrected.R"))

params$n <- 1000

progress <- function(...)
{
    cat(date(), ' ')
    invisible(lapply(list(...), function(x) cat(x,'\n')))
}

simulation <- function(n) do.call("rbind", lapply(1:n, function(y) {
    if(y %% 100 == 1) progress(paste0("Starting run ",y))
    microsim_corr_icer(params, seed=NULL)
}))
