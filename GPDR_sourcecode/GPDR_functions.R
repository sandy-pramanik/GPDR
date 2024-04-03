

# sourcing necessary libraries and functions ====
library(tidyverse)
library(doParallel)
source(file.path(sourcecode.path, "GPDR_sourcecode", "GPDR_Matern.R"))


# user defined functions ----
## random samples from DP(alpha) ----
rdp = function(n, alpha){
  # sample 1 measure mu from DP(unif[0,1], alpha)
  # then sample n values from mu
  # return mu and n samples
  
  probs = c(); bk = 1
  while(sum(probs) < 1-1e-5){
    sk = rbeta(1, 1, alpha)  
    probs = c(probs, sk * bk)
    bk = bk * (1-sk)
  }
  probs = probs / sum(probs)
  locs = runif(length(probs))
  samples = sample(locs, n, replace=TRUE, prob=probs)
  list(samples=samples, mu=list(locs, probs))
}

## simulate data using DP(alpha) ----
generate_data_dp = function(beta, n, m, sig=0.1, alpha=25){
  X = list(); Y = c(); kde = locs_rdp = probs_rdp = list()
  
  for(i in 1:n){
    dpi = rdp(m, alpha)
    X[[i]] = dpi$samples
    Y[i] = sum(beta(dpi$mu[[1]]) * dpi$mu[[2]]) + rnorm(1, sd=sig)
    kde[[i]] = density(X[[i]], from=0, to=1, n=512)
    locs_rdp[[i]] = dpi$mu[[1]]
    probs_rdp[[i]] = dpi$mu[[2]]
  }
  list(X=X, y=Y, kde=kde, locs_rdp=locs_rdp, probs_rdp=probs_rdp)
}

## simulate data using independent continuous covariates ----
generate_data_cp_indep = function(beta_gen, n, m, sig=0.1, 
                            sd.mu_gen=2, sigma_gen=.3){
  
  muX = rnorm(n=n, sd=sd.mu_gen)
  X=kde=list()
  Y = c()
  for(i in 1:n){
    X[[i]] = 1/(1 + exp(-rnorm(n=m, mean = muX[i], sd=sigma_gen)))
    Y[i] = mean(beta_gen(1/(1 + exp(-rnorm(n=1e+5, mean = muX[i], sd=sigma_gen))))) +
      rnorm(n=1, sd=sig)
    kde[[i]] = density(X[[i]], from=0, to=1, n=512)
  }
  list(X=X, y=Y, kde=kde, mu=muX)
}

## simulate data using dependent continuous covariates ----
generate_data_cp_dep = function(beta_gen, n, m, sig=0.1, 
                                sd.mu_gen=2, sigma_gen=.3, rho_gen=.5){
  
  muX = rnorm(n=n, sd=sd.mu_gen)
  Sigma.mat = matrix(nrow = m, ncol = m)
  Sigma.mat = (sigma_gen^2)*(rho_gen^(abs(row(Sigma.mat)-col(Sigma.mat))))
  logitXallsamp.std = LaplacesDemon::rmvn(n=n, mu=rep(0,m), Sigma = Sigma.mat)
  # logitXallsamp.std = LaplacesDemon::rmvnp(n=n, mu=rep(0,m),
  #                                          Omega = as.matrix(ar.matrix::Q.AR1(M=m, sigma = sigma_gen, rho = rho_gen)))
  X=kde=list()
  Y = c()
  for(i in 1:n){
    # X[[i]] = 1/(1 + exp(-LaplacesDemon::rmvnp(n=1, mu=rep(muX[i],m),
    #                                           Omega = as.matrix(ar.matrix::Q.AR1(M=m, sigma = sigma_gen, rho = rho_gen)))))
    X[[i]] = 1/(1 + exp(-(muX[i] + as.numeric(logitXallsamp.std[i,]))))
    Y[i] = mean(beta_gen(1/(1 + exp(-rnorm(n=1e+5, mean = muX[i], sd=sigma_gen))))) +
      rnorm(n=1, sd=sig)
    kde[[i]] = density(X[[i]], from=0, to=1, n=512)
  }
  list(X=X, y=Y, kde=kde, mu=muX)
}

## true beta functions for simulations ----
beta = function(xvec){
  out = c()
  for(x in xvec){
    if(x < 1/3) out = c(out, 3 / 2 * x^2)
    else if(x < 2/3) out = c(out, 1/6 + x - 3 * (x-1/3)^2 - 1/3)
    else out = c(out, -x + 3/2 * (x-2/3)^2 + 5/6)
  }
  out
}

beta1 = function(xvec){
  10 * xvec * exp(-5 * xvec)
}

## fit model with kde method ----
kde_fit = function(data, Kxx, Kxc, Kcc, sig, newx=seq(0, 1, 0.01)){
  n = length(data$y)
  TKf = matrix(NA, length(newx), n); fTKf = matrix(NA, n, n)
  for(i in 1:n){
    TKf[, i] = Kxc %*% (data$kde[[i]]$y / sum(data$kde[[i]]$y))
    for(j in i:n){
      fKf = data$kde[[i]]$y / sum(data$kde[[i]]$y) * c(Kxx %*% data$kde[[j]]$y / sum(data$kde[[j]]$y))
      fTKf[i, j] = sum(fKf)
      fTKf[j, i] = fTKf[i, j]
    }
  }
  M = solve(fTKf + sig^2*diag(n))
  list(
    E=c(TKf %*% M %*% data$y), 
    COV=Kcc - TKf %*% M %*% t(TKf)
  )
}
