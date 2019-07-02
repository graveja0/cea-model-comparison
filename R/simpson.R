# Simpson's method of integration
# Far better performance that "life-table" aka trapezoidal method
alt_simp_coef <- function(i) c(17, 59, 43, 49, rep(48, i-8), 49, 43, 59, 17) / 48
alt_simp      <- function(x,h) h*sum(alt_simp_coef(length(x)) * x)
