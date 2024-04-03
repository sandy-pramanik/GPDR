library(tidyverse)
# source("bdr.R")
library(glmnet)
library(kernlab)

## load NHANES data

## Friday has the highest correlation with Age
sat_data = Act_Analysis %>% dplyr::filter(Act_Analysis$WEEKDAY==5)
#### activity data
sat_P = as.matrix(sat_data[,paste0("MIN", 1:1440)])
sat_P[is.na(sat_P)] = 0
#### down sample to every 20 min average count data
sat_P = array(sat_P, dim=c(nrow(sat_P), 20, 1440 / 20))
sat_P = apply(sat_P, c(1, 3), mean)
### log transformation (distribution regression is invariant to it)
sat_P = log(pmax(sat_P, 0.025))
#### target variables
sat_y = c(scale(sat_data$Age))
# sat_y = sat_data$Age


## kernel embedding in landmark points
minp = min(sat_P); maxp = max(sat_P)
l = maxp - minp
landmarks = seq(minp, maxp, length.out=10)

phi = array(
  t(apply(sat_P, 1, function(s) kernelMatrix(rbfdot(.5 / l^2), s, landmarks))),
  dim = c(nrow(sat_P), ncol(sat_P), length(landmarks))
)
mu_10 = apply(phi, c(1, 3), mean)

## kernel embedding in landmark points (k = 50)
landmarks = seq(minp, maxp, length.out=50)

phi = array(
  t(apply(sat_P, 1, function(s) kernelMatrix(rbfdot(.5 / l^2), s, landmarks))),
  dim = c(nrow(sat_P), ncol(sat_P), length(landmarks))
)
mu_50 = apply(phi, c(1, 3), mean)

## 5-fold Cross Validation ----
idx = c(rep(1, 543), rep(2, 543), rep(3, 542), rep(4, 542), rep(5, 542))

cv.blr_k10 = c()
cv.blr_k50 = c()
for(it in 1:500){
  idx_it = idx[sample(nrow(sat_P))]
  cv_it_10 = c()
  cv_it_50 = c()
  for(k in 1:5){
    Xtrain10 = mu_10[idx_it != k, ]
    Xtrain50 = mu_50[idx_it != k, ]
    ytrain = sat_y[idx_it != k]
    Xtest10 = mu_10[idx_it == k, ]
    Xtest50 = mu_50[idx_it == k, ]
    ytest = sat_y[idx_it == k]
    
    # 10 points
    fit0 = cv.glmnet(Xtrain10, ytrain, alpha=0)
    pred = predict(fit0, Xtest10, s="lambda.min")
    cv_it_10[k] = mean((ytest - pred)^2) / var(ytest)
    
    # 50 points
    fit1 = cv.glmnet(Xtrain50, ytrain, alpha=0)
    pred = predict(fit1, Xtest50, s="lambda.min")
    cv_it_50[k] = mean((ytest - pred)^2) / var(ytest)
  }
  print(c(it, mean(cv_it_10), mean(cv_it_50)))
  cv.blr_k10 = rbind(cv.blr_k10, cv_it_10)
  cv.blr_k50 = rbind(cv.blr_k50, cv_it_50)
}

save(cv.blr_k10, cv.blr_k50, file="BLR_AgeExperimentResults.RData")

mean(1 - rowMeans(cv.blr_k10)); quantile(1 - rowMeans(cv.blr_k10), c(0.025, 0.975))

mean(1 - rowMeans(cv.blr_k50)); quantile(1 - rowMeans(cv.blr_k50), c(0.025, 0.975))
