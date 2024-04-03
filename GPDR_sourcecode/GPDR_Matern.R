###################################################################
## Simply  Y_i = 1/m \sum_j g(x_ij) + e_i
##         g ~ GP(0, K)
###################################################################


##################################################################
## Kernel function
matern_kernel <- function(d, nu = 1.5, sigma = 1, l = 1) {
  # d: distance (scalar, vector, matrix)
  # nu: smoothness (only accept 0.5, 1.5, 2.5 and Inf for now)
  # sigma: covariance scale
  # l: location scale
  p <- nu - 0.5
  if (p == 0) {
    sigma^2 * exp(- d / l)
  } else if (p == 1) {
    sigma^2 * (1 + sqrt(3)*d/l) * exp(- sqrt(3)*d/l)
  } else if (p == 2) {
    sigma^2 * (1 + sqrt(5)*d/l + 5*d^2 / (3*l^2)) * exp(-sqrt(5)*d/l)
  } else {
    sigma^2 * exp(-d^2 / (2 *l^2))
  }
}

MK = function(nu=1.5, sigma=1, l=1){
  function(x, y){
    d = outer(x, y, function(a, b) abs(a-b))
    matern_kernel(d, nu, sigma, l)
  }
}


########################################################################################
## Aux function:
##    return interval boundaries such that each one
##    contains roughly k data points
greedy_shrink = function(P, k=10){
  all_p = c()
  for(pi in P)
    all_p = c(all_p, pi)
  probs = seq(0, 1, length.out=round(length(unique(all_p)) / k))
  list(qs=quantile(all_p, probs), cdf=ecdf(all_p))
}

########################################################################################
## Aux function
##    return design matrix A
##    after shrinking
Amat = function(P, k=10, verbose=TRUE){
  if(inherits(P, "matrix")){
    Plist = list()
    for(i in 1:nrow(P)){
      Plist[[i]] = P[i,]
    }
    P = Plist
    rm(Plist)
  }
  if(verbose) cat("Start measurement shrinkage\n")
  gs = greedy_shrink(P, k)
  iv = gs$qs; cdf=gs$cdf
  l = length(iv) # will have l-1 intervals
  iv[l] = iv[l] + min(diff(iv))
  if(verbose) cat("Constructing design matrix\n")
  A = matrix(NA, length(P), l - 1)
  for(i in 1:length(P)){
    loci = ceiling(cdf(P[[i]]) * (l-1))
    loci[loci == 0] = 1
    loci = factor(loci, levels=1:(l-1))
    Ai = table(loci) / length(loci)
    names(Ai) = NULL
    A[i, ] = Ai
  }
  pts = (iv[1:(l-1)] + iv[2:l]) / 2
  list(A = A, xs = pts)
}


GPfit = function(
  P, y, kernel_gen=MK, sigma=0.1, nu=1.5, sigK=1, l=1, k=10, verbose=TRUE
){
  D = Amat(P, k=k, verbose=verbose); A = D$A
  kernel_func = kernel_gen(nu, sigK, l)
  K = kernel_func(D$xs, D$xs)
  if(verbose) cat("Begin solving\n")
  M = solve(A %*% K %*% t(A) + sigma^2 * diag(length(y)))
  My = c(M %*% y)
  
  if(verbose) cat("Finished\n")
  output = list(A=A, xs=D$xs, kf=kernel_func, M=M, My=My)
  attr(output, "class") = "MGPDR"
  output
}


predict.MGPDR = function(model, newx, mean.only=FALSE){
  Kp = model$kf(newx, model$xs)
  KpA = Kp %*% t(model$A)
  Egx = c(KpA %*% model$My)
  if(!mean.only){
    COVgx = model$kf(newx,newx) - KpA %*% model$M %*% t(KpA)
    list(E = Egx, COV=COVgx)
  } else {
    list(E = Egx)
  }
}


plot.MGPDR = function(model, scale=1.96, ...){
  xs = seq(min(model$xs), max(model$xs), length.out=1e3)
  res = predict(model, xs)
  yub = res$E + scale * sqrt(diag(res$COV))
  ylb = res$E - scale * sqrt(diag(res$COV))
  rmax = max(yub)
  rmin = min(ylb)
  param = list(...)
  param[["ylim"]] = c(rmin, rmax)
  param[["x"]] = xs
  param[["y"]] = res$E
  if(is.null(param[["xlab"]]))
    param[["xlab"]] = "x"
  if(is.null(param[["ylab"]]))
    param[["ylab"]] = "y"
  param[["type"]] = "l"
  do.call("plot", param)
  param[["lty"]] = 2
  param[["y"]] =yub
  do.call("lines", param)
  param[["y"]] = ylb
  do.call("lines", param)
}


generate_data = function(beta, n, m, sig=0.1){
  Xreal = list()
  X = list()
  for(i in 1:n){
    g1 = rgamma(1, 4, 1); g2 = rgamma(1, 4, 1);
    Xreal[[i]] = rbeta(1e5, rgamma(1, g1, g2), rgamma(1, g1, g2))
    X[[i]] = sample(Xreal[[i]], m, replace=TRUE)
  }
  Y = sapply(Xreal, function(xs){mean(beta(xs))})
  Y = Y + sig * rnorm(length(Y))
  list(X=X, y=Y)
}

beta = function(x){
  x + (1 + 1/2 * sin(2*pi * x^2))
}
# g1 = rgamma(1, 4, 1); g2 = rgamma(1, 4, 1)
# hist(rbeta(1e5, rgamma(1, g1, g2), rgamma(1, g1, g2)))


##########################################################################
## Some test associated with the model, not included in the paper
main = function(){

########################################################
## Simulation
data = generate_data(beta, 50, 3000, sig=0.5)

system.time(fit0 <- GPfit(data$X, data$y, sigma=1, nu=Inf, sigK=1, l=0.1, k=50))

plot(fit0, scale=1.96)
lines(seq(0,1,0.01), beta(seq(0,1,0.01)), col="red")

plot(seq(0,1,0.01), beta(seq(0,1,0.01)), type="l")

library(mgcv)

Xmat = c()
for(s in data$X){
  Xmat = rbind(Xmat, s)
}
df = list(X=Xmat, y=data$y, W=matrix(1/ncol(Xmat), nrow(Xmat), ncol(Xmat)))
dfi = list(X=Xmat[,2], y=data$y)

fit0 = gam(y~s(X, k=30), data=dfi)
plot(fit0)

fit1 = gam(y~s(X, k=30, by=W), data=df)
plot(fit1, rug=FALSE, ylab="f(X)")

#########################################################
## Canadian weather data
library(refund)

wdata = fda::CanadianWeather
start = 172; end = 265
temp = t(wdata$dailyAv[start:end,,1])
prec = wdata$dailyAv[start:end,,2]

temp_dens = c()
for(i in 1:nrow(temp)){
  di = density(temp[i,], n=512, from=-10, to=26)$y
  temp_dens = rbind(temp_dens, di)
}

CDF = ecdf(as.vector(temp))
ntemp = matrix(CDF(temp), 35, end-start+1)
summer_precip = colSums(prec)

## LOOCV comparison
cv.error.dr = c()
cv.error.kde = c()
cv.error.fr = c()
for(i in 1:nrow(ntemp)){
  print(i)
  fiti = GPfit(temp[-i,], summer_precip[-i], sigma=1, nu=Inf, sigK=1, l=15, k=1)
  cv.error.dr[i] = mean(predict(fiti, temp[i,])$E) - summer_precip[i]
  df = data.frame(precip=summer_precip[-i])
  df$temp = temp[-i,]
  fiti = pfr(precip ~ lf(temp, presmooth="bspline"), data=df)
  cv.error.fr[i] = predict(fiti, list(temp=matrix(temp[i,], 1, ncol(temp)))) - summer_precip[i]
  df = data.frame(precip=summer_precip[-i])
  df$temp = temp_dens[-i,]
  fiti = pfr(precip ~ lf(temp, presmooth="bspline", k=30), data=df)
  cv.error.kde[i] = predict(fiti, list(temp=matrix(temp_dens[i,], 1, ncol(temp_dens)))) - summer_precip[i]
}

mean(abs(cv.error.dr)^2)
mean(abs(cv.error.kde)^2)
mean(abs(cv.error.fr)^2)

# plot(fit1)
#
# predicted = c()
# for(i in 1:nrow(ntemp)){
#   val = predict(fit1, ntemp[i,])$E
#   predicted[i] = mean(val)
# }

1 - sum((summer_precip - predicted)^2) / sum((summer_precip - mean(summer_precip))^2)

df = data.frame(precip=summer_precip)
df$temp = temp
fit2 = pfr(precip ~ lf(temp, presmooth="fpca.sc"), data=df)
predict(fit2, df)
1 - sum((summer_precip - fit2$fitted.values)^2) / sum((summer_precip - mean(summer_precip))^2)

}