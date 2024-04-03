## Simulation code for BDR with landmark points = 10, 50

# suppressWarnings(suppressMessages(library(foreach)))
source(file.path(sourcecode.path, "GPDR_sourcecode", 'BDR', "bdr.R"))    # sources ".../GPDR_sourcecode/BDR/bdr.R"

## model ## BDR functions
bdr_landmarks = function(
    Xtrain, ytrain, Xtest=NULL, ytest=c(),
    model, # compiled stan model
    Rscale = 1, # this is the eta that scales R
    eta = 1, # this is the SD of the Gaussian measure
    sigmaK = 1, # sigma for kernel
    k = 30, # number of landmark points
    landmarks = NULL, # landmark points
    l = 1, # length scale
    it = 200, # iteration numbers
    chains = 4, # chain numbers
    adapt = 0.8, # adapt_delta
    grid = seq(0, 1, 0.01) # get the regression function values on the grid
) {
  ## dimensions
  ntrain = nrow(Xtrain); ptrain = ncol(Xtrain)
  if (!is.null(Xtest))
    ntest = nrow(Xtest); ptest = ncol(Xtest)
    
    ## landmarks u, evenly spaced in [0, 1]
    if (is.null(landmarks)) {
      u = seq(0, 1, length.out=k)
    } else {
      u = landmarks
      k = length(landmarks)
    }
    
    ## use non stationary R kernel  
    Rplus = exp(-.5 * outer(1:k, 1:k, Vectorize(function(i,j) .25 * (u[i]+u[j])^2)) / (.5 * l + eta^2))
    R = Rscale * sigmaK^2 * kernelMatrix(rbfdot(.5 / (sqrt(2)*l)^2), u) * Rplus  
    R = R + diag(1e-6, length(u))  # add some jitter to preserve positive definiteness
    
    ## empirical kernel mean embedding
    phi_train = array(
      t(apply(Xtrain, 1, function(s) sigmaK^2 * kernelMatrix(rbfdot(.5 / l^2), s, u))),
      dim = c(ntrain, ptrain, k)
    )
    mu_train = apply(phi_train, c(1, 3), mean)
    
    if (!is.null(Xtest)) {
      phi_test = array(
        t(apply(Xtest, 1, function(s) sigmaK^2 * kernelMatrix(rbfdot(.5 / l^2), s, u))),
        dim = c(ntest, ptest, k)
      )
      mu_test = apply(phi_test, c(1, 3), mean)
    } else {
      mu_test = matrix(NA, 0, k)
    }
    
    mu = rbind(mu_train, mu_test)
    
    ## build stan data list (part of)
    stan_data = list(
      d = ncol(mu),
      p = nrow(mu),
      ntrain = length(ytrain), 
      mu = mu,
      y = ytrain,
      ytrue = c(ytrain, ytest)
    )
    
    ## using global mean for kernel mean embedding
    mu0 = apply(phi_train, 3, mean)
    ## using empirical covariance for KME
    Sigma0 = cov(array(phi_train, dim=c(ntrain*ptrain, k)))
    
    ## build final mu and Sigma for stan model
    Sigma = array(0, dim=c(nrow(mu), ncol(mu), ncol(mu)))
    
    for(i in 1:nrow(mu)) {
      n = ptrain
      mu[i,] = R %*% solve(R+Sigma0/n,mu[i,] - mu0) + mu0
      Sigma[i,,] = R - R %*% solve(R+Sigma0/n,R)
    }
    stan_data$mu = mu
    stan_data$Sigma = Sigma
    
    ## Run the model
    ptm = proc.time()
    fit = sampling(
      model, data=stan_data, iter=it, warmup=round(it/2),
      chains=chains, control=list(adapt_delta=adapt),
      verbose=FALSE, show_messages=FALSE, refresh=0
    )
    elapsed = as.numeric((proc.time() - ptm)[3])
    
    ## get predicted value (posterior mean) for test data (if any)
    if (!is.null(Xtest)){
      y_pred = summary(fit, "yhat")$summary[(nrow(mu_train)+1):nrow(mu), 1]
    } else {
      y_pred = c()
    }
    
    ## get posterior regression curve
    beta = summary(fit, "beta")$summary[, 1]
    alpha = summary(fit, "alpha")$summary[, 1]
    curve = c(sigmaK^2 * kernelMatrix(rbfdot(.5 / l^2), grid, u) %*% beta) + alpha
    
    list(
      elapsed=elapsed, fit=fit, 
      u=u, grid=grid, curve=curve, 
      ypred=y_pred
    )
}


## run the simulation

# set.seed(54321)
set.seed(12345)

my.cluster <- parallel::makeCluster(
  12, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

results = c()
for(n in c(50, 100, 200, 300, 400)){
  for(m in c(50, 100, 250, 500, 1000, 2000)){
    print(paste("Processing n = ", n, "and m = ", m))

    # out_n_m <- foreach(
    #   it=1:100, .combine="rbind", .packages=c("rstan", "kernlab")
    # ) %dopar% {
    out_n_m <- foreach(
      it=1:25, .combine="rbind", .packages=c("rstan", "kernlab")
    ) %dopar% {
      load("bdr_stan_model.RData")
      
      data = generate_data_dp(beta1, n, m)
      
      suppressWarnings({
        fit1 <- bdr_landmarks(
          matrix(unlist(data$X), nrow=n, byrow=TRUE), data$y,
          model=bdr_model, k=10, it=200, l=0.25, chains=2
        )
        fit2 <- bdr_landmarks(
          matrix(unlist(data$X), nrow=n, byrow=TRUE), data$y,
          model=bdr_model, k=50, it=200, l=0.25, chains=2
        )
      })
      
      fit1_l2 = sum((beta1(seq(0,1,0.01)) - fit1$curve)^2)
      fit2_l2 = sum((beta1(seq(0,1,0.01)) - fit2$curve)^2)
      # results = rbind(results, c(n, m, it, fit1_l2, fit2_l2))
      c(n, m, it, fit1_l2, fit2_l2)
    }
    cat("\n")
    
    results = rbind(results, out_n_m)
    write.table(results, file="bdr_simulation_results", row.names=FALSE)
  }
}

parallel::stopCluster(cl = my.cluster)

BDR_simRes = results
save(BDR_simRes, file="BDR_simulation_results.RData")
