

rm(list = ls())


# sourcing library and codes ----
sourcecode.path = ...    # specifies path ".../GPDR_sourcecode" to the GPDR_sourcecode folder
source(file.path(sourcecode.path, "GPDR_sourcecode", "GPDR_functions.R"))    # sources ".../GPDR_sourcecode/GPDR_functions.R"


# data generation specifics ----
alpha_DP = .1
k = 10
xseq = seq(0, 1, by = 0.01)
beta.truth = beta1(xseq)


# specifics for KDE ----
set.seed(0)
data = generate_data_dp(beta1, 100, 100)
K = MK(nu=Inf, sigma=1, l=0.25)
## x grid and final evaluate grid are always the same
## for all m and n, therefore we only need to calculate once
Kxx = K(data$kde[[1]]$x, data$kde[[1]]$x)
Kxc = K(xseq, data$kde[[1]]$x)
Kcc = K(xseq, xseq)


# combinations of n, m and replications ----
param.combinations = 
  as.data.frame(tidyr::expand_grid(n = c(50, 100, 200, 300, 400),
                     m = c(10, 50, 100, 250, 500, 1000, 2000),
                     it = 1:100
                     # it = 1:3
                     ))

# starting (parallel) simulation ----
## set cores=1 below to implement serially
doParallel::registerDoParallel(cores = 8)
simRes = foreach::foreach(l = 1:nrow(param.combinations), #.packages=c("rstan", "kernlab"),
                          .combine = 'rbind', .multicombine = T) %dopar% {
  
  set.seed(l)
  
  n = param.combinations$n[l]
  m = param.combinations$m[l]
  it = param.combinations$it[l]
  
  data = generate_data_dp(beta = beta1, n = n, m = m, alpha = alpha_DP)
  fit = GPfit(data$X, data$y, nu=Inf, l=0.25, k=k, sigma=0.1, verbose=FALSE)
  pred0 = predict(fit, xseq, mean.only=TRUE)$E
  
  fit = kde_fit(data, Kxx, Kxc, Kcc, sig=0.1, newx=xseq)
  pred1 = fit$E
  
  exp_l2 = sum((beta.truth - pred0)^2)
  kde_l2 = sum((beta.truth - pred1)^2)
  
  print(as.numeric(param.combinations[l,]))
  
  c(n, m, it, exp_l2, kde_l2)
  
                          }

# loading simulation output ----
save(simRes, file=paste0("simulation_results_alpha_", alpha_DP, ".RData"))
