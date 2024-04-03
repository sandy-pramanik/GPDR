

rm(list = ls())


# sourcing library and codes ----
sourcecode.path = ...    # specifies path ".../GPDR_sourcecode" to the GPDR_sourcecode folder
source(file.path(sourcecode.path, "GPDR_sourcecode", "GPDR_functions.R"))    # sources ".../GPDR_sourcecode/GPDR_functions.R"


# data generation specifics ----
k = 10
rhoX = 0.5
xseq = seq(0, 1, by = 0.01)
beta.truth = beta1(xseq)


# specifics for KDE ----
set.seed(0)
data = generate_data_cp_dep(beta_gen = beta1, n = 100, m = 100)
K = MK(nu=Inf, sigma=1, l=0.25)
## x grid and final evaluate grid are always the same
## for all m and n, therefore we only need to calculate once
Kxx = K(data$kde[[1]]$x, data$kde[[1]]$x)
Kxc = K(xseq, data$kde[[1]]$x)
Kcc = K(xseq, xseq)


# combinations of n, m and replications ----
param.combinations = 
  as.data.frame(tidyr::expand_grid(n = c(50, 100, 200#, 300, 400
  ),
  m = c(10, 50, 100
        #, 250, 500#, 1000, 2000
  ),
  # it = 1:3
  it = 1:100
  ))

# starting (parallel) simulation ----
## set cores=1 below to implement serially
doParallel::registerDoParallel(cores = 5)
simRes = foreach::foreach(l = 1:nrow(param.combinations), #.packages=c("rstan", "kernlab"),
                          .combine = 'rbind', .multicombine = T) %dopar% {
                            
                            set.seed(l)
                            
                            n = param.combinations$n[l]
                            m = param.combinations$m[l]
                            it = param.combinations$it[l]
                            
                            # cat(it, " ")
                            
                            set.seed(it)
                            data = generate_data_cp_dep(beta_gen = beta1, n = n, m = m, sig=0.1, 
                                                        rho_gen = rhoX)
                            fit = GPfit(data$X, data$y, nu=Inf, l=0.25, k=k, sigma=0.1, verbose=FALSE)
                            pred0 = predict(fit, xseq, mean.only=TRUE)$E
                            
                            fit = kde_fit(data, Kxx, Kxc, Kcc, sig=0.1, newx=xseq)
                            pred1 = fit$E
                            rm(fit)
                            
                            exp_l2 = sum((beta.truth - pred0)^2)
                            kde_l2 = sum((beta.truth - pred1)^2)
                            
                            print(param.combinations[l,])
                            
                            c(n, m, it, exp_l2, kde_l2)
                            
                          }

# loading simulation output ----
save(simRes, file=paste0("simulation_results_cp_dep_rho", rhoX, ".RData"))
