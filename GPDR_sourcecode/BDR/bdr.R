suppressWarnings(suppressMessages(library(rstan)))
suppressWarnings(suppressMessages(library(kernlab)))


suppressWarnings(suppressMessages(rstan_options(auto_write = TRUE)))
options(mc.cores = parallel::detectCores())
# suppressWarnings(suppressMessages(
#   options(mc.cores = 4)
# ))


suppressWarnings(suppressMessages(
  bdr_model <- stan_model("bdr-landmarks-conjugacy.stan")
))
save(bdr_model, file = "bdr_stan_model.RData")


## BDR functions
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
  grid = seq(0, 1, 0.01), # get the regression function values on the grid
  verbose = FALSE
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
  if(!verbose){
    fit = sampling(
      model, data=stan_data, iter=it, warmup=round(it/2),
      chains=chains, control=list(adapt_delta=adapt),
      verbose=FALSE, show_messages=FALSE, refresh=0
    )
  } else {
    fit = sampling(
      model, data=stan_data, iter=it, warmup=round(it/2),
      chains=chains, control=list(adapt_delta=adapt)
    )
  }
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


##
# Xtrain = matrix(rnorm(100 * 50), 100, 50)
# ytrain = rnorm(100)
# 
# res = bdr_landmarks(Xtrain, ytrain, Xtrain[1:10,], ytrain[1:10], model=bdr_model)
# res = bdr_landmarks(Xtrain, ytrain, model=bdr_model)

# 
# data = generate_data(beta1, 200, 200)
# res = bdr_landmarks(
#   matrix(unlist(data$X), nrow=200, byrow=TRUE),
#   data$y,
#   model=bdr_model,
#   k=50, it=200,
#   l=0.25, eta=0.5, Rscale=1
# )
# 
# plot(res$grid, res$curve, type="l", ylim=c(0, 1))
# lines(seq(0,1,0.01), beta1(seq(0,1,0.01)), col="red")

# xx = array(
#   t(apply(rbind(matrix(1, 5, 3), c(1,2,3)), 1, function(s) kernelMatrix(rbfdot(.5), s, c(1, 2)))), 
#   dim = c(6, 3, 2)      
# )
