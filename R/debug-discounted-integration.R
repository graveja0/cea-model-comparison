source('R/simple-params.R')
source('R/simple-deq.R')
source('R/simple-markov-corrected.R')

params$disc_rate <- inst_rate(1-1/(1 + params$disc), 1)
params$p_o <- 0
params$p_g <- 0

mc  <- markov_corr(params, method="simpson")
m.M <- mc$m.M
ds  <- deq_simulation(params)


mc$results
deq_summary(ds, params)

round(100*(mc$results - deq_summary(ds, params)) / deq_summary(ds, params), 3)

sapply(names(params), function(n) assign(n, params[[n]], envir=baseenv()))
n.t       <- horizon*interval
d.r  <- inst_rate(1-1/(1 + params$disc), 1)
#v.dw <- c(1,(exp(-d.r*(0:(n.t-2))) - exp(-d.r*(1:(n.t-1))))/d.r)
v.dw <- -diff(exp(-d.r * (0:(n.t+1))))/d.r # calculate discount weights for costs for each cycle based on discount rate d.r

# Discounted accumulators
plot(ds[,"time"], disc_acc(ds, "b_c", ds[,"dsc"], "simpson"), typ="l", xlab="time", ylab="b_c disc diff")
points(0:40, disc_acc(m.M, "CUM_BS", v.dw, "simpson")*(ds[2,"time"]-ds[1,"time"] ))
points(1:40-0.5, disc_acc(m.M, "CUM_BS", v.dw, "life-table")*(ds[2,"time"]-ds[1,"time"] ), col="green")

#plot(ds[,"time"], ds[,"b_c"], typ="l", xlab="time", ylab="b_c")
#points(0:40, m.M[,"CUM_B"])

plot(ds[-1,"time"], diff(ds[,"b_c"])/(ds[2,"time"] - ds[1,"time"]), type="l", xlab="time", ylab="diff(b_c)")
points(1:40, diff(m.M[,"CUM_BS"]))


# Discounted accumulators
plot(ds[,"time"], disc_acc(ds, "a_c", ds[,"dsc"], "simpson"), typ="l", xlab="time", ylab="b_c disc diff")
points(0:40,      disc_acc(m.M, "CUM_A", v.dw, "simpson")*(ds[2,"time"]-ds[1,"time"] ))
points(1:40-0.5,  disc_acc(m.M, "CUM_A", v.dw, "life-table")*(ds[2,"time"]-ds[1,"time"] ), col="green")


plot(ds[,"time"], c(0,diff(ds[,"db_c"]))/(ds[2,"time"]-ds[1,"time"]), typ="l", xlab="time", ylab="discounted B survival velocity")
points(1:40-0.5,      disc_acc(m.M, "CUM_BS", v.dw[1:41], "life-table"))

eff <- function(t1, t2, eta, b) log((exp(-eta*exp(b*t2)))/(exp(-eta*exp(b*t1)))) / (t1-t2)

# Attempt to find effective discounting
eff.dw <- sapply(1:45, function(t) eff(t-1, t, -1/d.r, -d.r))
plot(ds[,"time"], c(0,diff(ds[,"db_c"]))/(ds[2,"time"]-ds[1,"time"]), typ="l", xlab="time", ylab="discounted B survival velocity")
points(1:40-0.5,      disc_acc(m.M, "CUM_BS", v.dw[1:41], "life-table"))

# Run a discounted integrator on differences in an accululator
take2 <- function(m, cols, disc, method) as.vector(diag(disc) %*% rbind(diff(m[,cols, drop=FALSE])))
points(1:40, take2(m.M, "CUM_BS", eff.dw[1:40]), col='green')
sum(take2(m.M, "CUM_BS", eff.dw[1:40]))

ds[nrow(ds),"db_c"]
sum(take2(m.M, "CUM_BS", eff.dw[1:40]))
m.M[nrow(m.M), "DCUM_BS"]

plot(ds[,"time"], c(0,diff(ds[,"db_c"]))/(ds[2,"time"]-ds[1,"time"]), typ="l", xlab="time", ylab="discounted B survival velocity")
points(1:40-0.5, diff(m.M[,"DCUM_BS"]))

plot(ds[,"time"], c(0,diff(ds[,"b_c"]))/(ds[2,"time"]-ds[1,"time"]), typ="l", xlab="time", ylab="discounted B survival velocity")
points(1:40, diff(m.M[,"CUM_BS"]))

plot(ds[,"time"], ds[,"b_c"], typ="l", xlab="time", ylab="discounted B survival velocity")
points(0:40, m.M[,"CUM_BS"])

plot(ds[,"time"], ds[,"b_c"], typ="l", xlab="time", ylab="discounted B survival velocity")
points(0:40, m.M[,"CUM_BS"])

plot(ds[,"time"], c(0,diff(ds[,"a_c"]))/(ds[2,"time"]-ds[1,"time"]), typ="l", xlab="time", ylab="discounted B survival velocity")
points(1:40-0.5, diff(m.M[,"CUM_A"]))

mc$results
deq_summary(ds, params)


cat("Rel Error", (mc$results[['living']] - deq_summary(ds, params)[['living']]) / deq_summary(ds, params)[['living']], "\n")

mc$results[['living']]
deq_summary(ds, params)[['living']]


plot(ds[,"time"], rowSums(ds[,c("h_u","a_p","a_a","b_p","b_a")])*ds[,c("dsc")], typ="l", xlab="time", ylab="discounted living")
d.r  <- inst_rate(1-1/(1 + params$disc), 1)
v.dw <- c(1,(exp(-d.r*(0:(n.t-2))) - exp(-d.r*(1:(n.t-1))))/d.r)
points(0:40, rowSums(m.M[, c("H", "A", "BS")])*v.dw)

# Works
plot(ds[,"time"], rowSums(ds[,c("h_u","a_p","a_a","b_p","b_a")]), typ="l", xlab="time", ylab="living")
points(0:40, rowSums(m.M[, c("H", "A", "BS")]))

# Now get discounting working
plot(ds[,"time"], rowSums(ds[,c("h_u","a_p","a_a","b_p","b_a")])*ds[,c("dsc")], typ="l", xlab="time", ylab="discounted living")
points(0:40, rowSums(m.M[, c("H", "A", "BS")])*exp(0:40 * -d.r) )

plot(ds[,"time"], ds[,c("a_c")], typ="l", xlab="time", ylab="Cum A")
points(0:40, m.M[, c("CUM_A")])

plot(ds[,"time"], c(0,diff(ds[,"a_c"]))/(ds[2,"time"]-ds[1,"time"])*ds[,"dsc"], typ="l", xlab="time", ylab="discounted A survival velocity")
#points(0:40, c(0,diff(m.M[,"CUM_A"]))*exp(0:40 * -d.r))
#points(0:40, c(0, diff(dmm[,"CUM_A"])), col='red')

plot(ds[,"time"], rowSums(ds[,c("a_da", "a_dp")]), type="l")
points(0:40, m.M[,"TUN"])

plot(ds[,"time"], ds[,c("dsc")]*rowSums(ds[,c("a_da", "a_dp")]), type="l", ylim=c(0, 0.0125))
points(0:40, mc$dmm[,"TUN"])
lines(0:40, mc$dmm[,"TUN"], col='red')


