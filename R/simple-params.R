inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

##############################################################
##
## Standardized Parameter List for Simple Model Comparison
##
params_lut <- 
  list(
    # Controls for model execution
    n          = "DES simulations to perform (1000s)",
    resolution = "Diff Eq Time step for DEQ approach.",
    interval   = "Markov Interval.",
    horizon    = "Time horizon (years) of simulation.",
    wtp        = "Willingness to pay threshold used for NMB ($1000s)",
    
    # Gompertz model of secular death for 40yr female
    # fit from 2012 social security data
    shape   = "Shape Parameter from Gompertz model of secular death for 40yr female fit from 2012 Social Security data.",
    rate    = "Rate Parameter from Gompertz model of secular death for 40yr female fit from 2012 Social Security data.",
    
    # Probabilities and rates
    p_o  = "Probability of ordering the genetic test.",
    p_bd = "Probability of death from adverse drug event.",
    p_g  = "Population prevalence of genetic variant.",
    r_a_pct = "Condition indication percentage",
    r_a_dur = "Condition indication duration (in years)",
    r_a  = "Inst rate of developing condition.",
    r_b_pct = "Adverse drug event percentage",
    r_b_dur = "Advarse drug event duration (in years)",
    r_b  = "Inst rate of adverse drug event.",
    rr_b = "Relative risk of adverse drug event | Genetic variant present and alternative therapy prescribed.",
    
    # Costs
    c_a   = "Initial cost of condition ($1000s).",
    c_bs  = "Cost of adverse drug event survival ($1000s).",
    c_bd  = "Cost of adverse drug event case fatality ($1000s).",
    c_tx  = "Cost of standard drug therapy.",
    c_alt = "Cost of alternate drug therapy.",
    c_t   = "Cost of genetic test.",
    
    # Disutilities
    d_a   = "Disutility of developing condition.",
    d_at  = "Disutility duration (years) after developing condition.",
    d_b   = "Lifetime disutility of adverse drug event among survivors.",
    
    # Discounting
    disc  ="Discount Rate."
  )

params <- list(
  # Controls for model execution
  n          = 1e6,      # DES simulations to perform
  resolution = 7/365,        # Diff Eq Time step for DEQ approach
  interval   = 1,            # Markov Interval
  horizon    = 40,           # Time horizon of simulation
  wtp        = 100000,       # Willingness to pay threshold
  
  # Gompertz model of secular death for 40yr female
  # fit from 2012 social security data
  shape   = 0.1007511,
  rate    = 0.0008370717,
  
  # Probabilities and rates
  p_o  = 1.0,                # Probability of ordering test
  p_bd = 0.10,               # Probability of death from B
  p_g  = 0.2,                # Probability of genetic variant
  r_a_pct = 10,              # Percentage who develop condition.
  r_a_dur = 10,              # Condition indication rate duration (in years).
  r_a  = inst_rate(.9999, 1e-7), # 10% Rate of A over a 10 year period
  r_b_pct = 5,               # Adverse drug event percentage
  r_b_dur = 5,               # Adverse drug event duration (in years). 
  r_b  = inst_rate(0.05, 5), # 5% Rate of B over a 5 year period  
  rr_b = 0.7,                # Reduced relative risk of B
  
  # Costs
  c_a   = 10000,             # Cost of Event A
  c_bs  = 25000,             # Cost of Event B survival
  c_bd  = 15000,             # Cost of Event B death
  c_tx  = 0.5,               # Cost of normal treatment
  c_alt = 2.75,                 # Cost of alternate treatment
  c_t   = 100,               # Cost of test
  
  # Disutilities
  d_a   = 0.05,              # Disutility of A
  d_at  = 1,                 # Duration of A in years.
  d_b   = 0.15,              # Disutility of B
  
  # Discounting
  disc  = 0.03               # Annual Discount Rate
)

params_psa <- list(
  # Controls for model execution
  n          = function(x) qnorm(p = x, mean = params$n,sd =0),      # DES simulations to perform
  resolution = function(x) qnorm(p = x, mean = params$resolution ,sd = 0),        # Diff Eq Time step for DEQ approach
  interval   = function(x) qnorm(p = x, mean = params$interval, sd = 0),            # Markov Interval
  horizon    = function(x) qnorm(p = x, mean = params$horizon, sd = 0),           # Time horizon of simulation
  wtp        = function(x) qnorm(p = x, mean = params$wtp, sd = 0),       # Willingness to pay threshold
  
  # Gompertz model of secular death for 40yr female
  # fit from 2012 social security data
  shape   = function(x) qnorm(p = x, mean = params$shape, sd = 0), 
  rate    = function(x) qnorm(p = x, mean = params$rate,sd = 0),
  
  # Probabilities and rates
  p_o  = function(x) qunif(p = x, min = 0.0, max = 1),                # Probability of ordering test (overwritten by runs to 0 and 1)
  p_bd = function(x) qunif(p = x, min = 0.0, max = 1),                       # Probability of death from B
  p_g  = function(x) qunif(p = x, min = 0.0, max = 0.5),              # Probability of genetic variant
  r_a_pct = function(x) qunif(p = x, min = 0.0, max = 1),              # Percentage who develop condition.
  r_a_dur = function(x) qunif(p = x, min = 1e-7, max = params$horizon),              # Condition indication rate duration (in years).
  r_b_pct = function(x) qunif(p = x, min = 0.0, max = 1),               # Adverse drug event percentage
  r_b_dur = function(x) qunif(p = x, min = 1e-7, max = params$horizon),               # Adverse drug event duration (in years). 
  rr_b = function(x) qunif(p = x, min = 0, max = 1),                # Reduced relative risk of B
  
  # Costs
  c_a   = function(x) qunif(p = x, min = 1000, max = 100000),              # Cost of Event A
  c_bs  = function(x) qunif(p = x, min = 1000, max = 100000),              # Cost of Event B survival
  c_bd  = function(x) qunif(p = x, min = 1000, max = 100000),             # Cost of Event B death
  c_tx  = function(x) qunif(p = x, min = 0.25, max = 4),                 # Cost of normal treatment
  c_alt = function(x) qunif(p = x, min = 0.25, max = 4),                # Cost of alternate treatment
  c_t   = function(x) qunif(p = x, min = 10, max = 10000),               # Cost of test
  
  # Disutilities
  d_a   = function(x) qunif(p = x, min = 0.0, max = 0.2),              # Disutility of A
  d_at  = function(x) qunif(p = x, min = 0, max = params$horizon),                 # Duration of A in years.
  d_b   = function(x) qunif(p = x, min = 0, max = 0.3),              # Disutility of B
  
  # Discounting
  disc  = function(x) qnorm(p = x, mean = 0.03, sd = 0)               # Annual Discount Rate
)


params_psa_realistic <- list(
  # Controls for model execution
  n          = function(x) qnorm(p = x, mean = params$n,sd =0),      # DES simulations to perform
  resolution = function(x) qnorm(p = x, mean = params$resolution ,sd = 0),        # Diff Eq Time step for DEQ approach
  interval   = function(x) qnorm(p = x, mean = params$interval, sd = 0),            # Markov Interval
  horizon    = function(x) qnorm(p = x, mean = params$horizon, sd = 0),           # Time horizon of simulation
  wtp        = function(x) qnorm(p = x, mean = params$wtp, sd = 0),       # Willingness to pay threshold
  
  # Gompertz model of secular death for 40yr female
  # fit from 2012 social security data
  shape   = function(x) qnorm(p = x, mean = params$shape, sd = 0), 
  rate    = function(x) qnorm(p = x, mean = params$rate,sd = 0),
  
  # Probabilities and rates
  p_o  = function(x) qunif(p = x, min = 0.0, max = 1),                # Probability of ordering test (overwritten by runs to 0 and 1)
  p_bd = function(x) qbeta(p = x, shape1 = params$p_bd*100, shape2 = 100-100*params$p_bd),                       # Probability of death from B
  p_g  = function(x) qbeta(p = x, shape1 = params$p_g*100, shape2 = 100-100*params$p_g),              # Probability of genetic variant
  r_a_pct = function(x)  100* qbeta(p = x, shape1 = params$r_a_pct, shape2 = 100-params$r_a_pct),              # Percentage who develop condition.
  r_a_dur = function(x) qunif(p = x, min = 1e-7, max =20),              # Condition indication rate duration (in years).
  r_b_pct = function(x) 100* qbeta(p = x, shape1 = params$r_b_pct, shape2 = 100-params$r_b_pct),               # Adverse drug event percentage
  r_b_dur = function(x) qunif(p = x, min = 0, max = 10),               # Adverse drug event duration (in years). 
  rr_b = function(x) qunif(p = x, min = 0.4, max = 1),                # Reduced relative risk of B
  
  # Costs
  c_a   = function(x) qpois(p = x, lambda = params$c_a),              # Cost of Event A
  c_bs  = function(x) qpois(p = x, lambda = params$c_bs),              # Cost of Event B survival
  c_bd  = function(x) qpois(p = x, lambda = params$c_bd),             # Cost of Event B death
  c_tx  = function(x) qpois(p = x, lambda = params$c_tx*365)/365,                 # Cost of normal treatment
  c_alt = function(x) qpois(p = x, lambda = params$c_alt*365)/365,                # Cost of alternate treatment
  c_t   = function(x) qunif(p = x, min = 50,max=150),               # Cost of test
  
  # Disutilities
  d_a   = function(x) qbeta(p = x, shape1 = params$d_a*100, shape2 = 100-100*params$d_a),              # Disutility of A
  d_at  = function(x) qunif(p = x, min = 0, max = 2),                 # Duration of A in years.
  d_b   = function(x) qbeta(p = x, shape1 = params$d_b*100, shape2 = 100-100*params$d_b),              # Disutility of B
  
  # Discounting
  disc  = function(x) qnorm(p = x, mean = 0.03, sd = 0)               # Annual Discount Rate
)




params_psa_test <- list(
  # Controls for model execution
  n          = function(x) qnorm(p = x, mean = params$n,sd =0),      # DES simulations to perform
  resolution = function(x) qnorm(p = x, mean = params$resolution ,sd = 0),        # Diff Eq Time step for DEQ approach
  interval   = function(x) qnorm(p = x, mean = params$interval, sd = 0),            # Markov Interval
  horizon    = function(x) qnorm(p = x, mean = params$horizon, sd = 0),           # Time horizon of simulation
  wtp        = function(x) qnorm(p = x, mean = params$wtp, sd = 0),       # Willingness to pay threshold
  
  # Gompertz model of secular death for 40yr female
  # fit from 2012 social security data
  shape   = function(x) qnorm(p = x, mean = params$shape, sd = 0), 
  rate    = function(x) qnorm(p = x, mean = params$rate,sd = 0),
  
  # Probabilities and rates
  p_o  = function(x) qnorm(p = x, mean = params$p_o, sd = 0),                # Probability of ordering test (overwritten by runs to 0 and 1)
  p_bd = function(x) qnorm(p = x, mean = params$p_bd, sd = 0),  # Probability of death from B
  p_g  = function(x) qnorm(p = x, mean = params$p_g, sd = 0),              # Probability of genetic variant
  r_a = function(x) qnorm(p = x , mean = params$r_a, sd =0),
  r_b  = function(x) qnorm(p = x, mean = params$r_b, sd = 0),# Inst Rate of B
  #rr_b = function(x) qnorm(p = x, mean = params$rr_b, sd = 0),                # Reduced relative risk of B
  rr_b = function(x) qunif(p = x, min = 0.5, max = 1), 
  
  # Costs
  c_a   = function(x) qnorm(p = x, mean = params$c_a, sd = 0),              # Cost of Event A
  c_bs  = function(x) qnorm(p = x, mean = params$c_bs, sd = 0),              # Cost of Event B survival
  c_bd  = function(x) qnorm(p = x, mean = params$c_bd, sd = 0),             # Cost of Event B death
  c_tx  = function(x) qunif(p = x, min = params$c_tx, max = params$c_tx),                 # Cost of normal treatment
  #c_alt = function(x) qunif(p = x, min = params$c_alt, max = params$c_alt),                # Cost of alternate treatment
  c_alt = function(x) qunif(p = x, min = 1, max = 4),
  c_t   = function(x) qunif(p = x, min = params$c_t, max = params$c_t),               # Cost of test
  
  # Disutilities
  d_a   = function(x) qnorm(p = x, mean = params$d_a, sd = 0),              # Disutility of A
  d_at  = function(x) qnorm(p = x, mean = 1, sd = 0),                 # Duration of A in years.
  d_b   = function(x) qnorm(p = x, mean = params$d_b, sd = 0),              # Disutility of B
  
  # Discounting
  disc  = function(x) qnorm(p = x, mean = 0.03, sd = 0)               # Annual Discount Rate
)



