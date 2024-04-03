# Install and load required packages
# install.packages(c("rstan", "ggplot2"))
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the Gaussian process regression model in Stan
gp_dist_model_code_nl_unknown <- '
data {
  int<lower=1> N;        // Number of subjects
  matrix[N,N] C;         // Input valuesGP dist cov matrix
  vector[N] y;           // Output values
  real<lower=0> s;      // data sd
  }

parameters {
  vector[N] f;       // GP latent function values
  real<lower=0,upper=s> tau; // Noise level
  real b0;
  real b2;
  real b3;
  }

model {
  // Prior for latent function values
  f ~ multi_normal(rep_vector(0, N), C);
  
  // Likelihood
  for(i in 1:N){
  y[i] ~ normal(b0+f[i]+b2*f[i]^2+b3*f[i]^3,tau);
  //y[i] ~ normal(1+f[i]+f[i]^2/2+f[i]^3/6,tau);
  }

  b2 ~ normal(0,1);
  b3 ~ normal(0,1);
}'

# Compile the Stan model
gp_dist_model_nl <- stan_model(model_code = gp_dist_model_code_nl_unknown)

# Generate synthetic data
set.seed(123)
N <- 20
M <- 25
R <- M*N
xorig <- seq(0, 1, length.out = N)
x <- pmin(1,pmax(0,rnorm(N*M,rep(xorig,each=M),sd=0.1)))
xmat=matrix(x,nrow=N,byrow=T)

## visualize the densities
df=cbind(as.data.frame(xmat),data.frame("id"=1:N))
dfnew=df %>% 
  pivot_longer(cols=!id) %>%
  as.data.frame()

dfnew %>% 
  ggplot() +
  geom_density(aes(x=value,group=as.factor(id),col=as.factor(id)))

y_true=rep(0,N)
for(i in 1:N){
  y_true[i] <- rnorm(1, mean=(1)*exp(mean(sin(2 * pi * xmat[i,]))), sd = 0.1)
  }

# Set up data list
#rho=0.005
dmat=outer(x, x, "-")

A=kronecker(diag(N),t(rep(1,M)))/M
rho=0.08
C=A%*%exp(-0.5*(dmat/rho)^2)%*%t(A)

gp_data <- list(N = N, y = y_true, C=C, s=sd(y_true))

# Fit the non-linear model to the data
gp_fit_nl <- sampling(gp_dist_model_nl, data = gp_data, chains = 4, iter = 10000)

# Extract fitted parameters
posterior_samples_nl <- as.data.frame(gp_fit_nl)
summary(posterior_samples_nl)

set.seed(1)
subindex=sample(nrow(posterior_samples_nl),100)
posterior_samples_sub_nl=posterior_samples_nl[subindex,]

#posterior_samples_sub$rho[i]

# Predict using the fitted model
x_pred <- seq(0, 1, length.out = 50)
f_pred_nl <- Reduce('rbind',lapply(1:length(subindex), function(i) {
  #rho=posterior_samples_sub$rho[i]
  t_data=exp(-0.5 * (outer(x, x, "-") / rho)^2)
  t_cross=exp(-0.5 * (outer(x_pred, x, "-") / rho)^2)
  as.vector(t_cross %*% t(A) %*% solve(A%*%t_data%*%t(A),as.numeric(posterior_samples_sub_nl[i,1:N])))
}),)

# Plot the true and fitted functions
ggplot() +
  geom_line(aes(x,  sin(2 * pi * x) , color = 'true f'), size = 1) +
  #geom_point(aes(x,  y_true ), color = "green", size = 3, alpha = 0.7) +
  geom_line(aes(x_pred, apply(f_pred_nl,2,median), color = 'fitted f'), size = 1) +
  geom_ribbon(aes(x = x_pred, ymin = apply(f_pred_nl, 2, quantile, 0.025), 
                  ymax = apply(f_pred_nl, 2, quantile, 0.975)), 
              fill = "red", alpha = 0.2) +
  labs(title = "GP Distribution Regression (non-linear): Coefficent funtion",x = "x",y = "f(x)=sim(2*pi*x)") +
  scale_color_manual(name='',
                     breaks=c('true f', 'fitted f'),
                     values=c('true f'='blue', 'fitted f'='red'))+
 theme(legend.title=element_text(size=20),
       legend.text=element_text(size=14),legend.position = "bottom")

t_pred <- seq(0, 1, length.out = 50)
phi_pred_nl <- Reduce('rbind',lapply(1:length(subindex), function(i) {
  t=t_pred
  posterior_samples_sub_nl[i,"b0"]+t+posterior_samples_sub_nl[i,"b2"]*t^2+posterior_samples_sub_nl[i,"b3"]*t^3
}),)

# Plot the true and fitted functions
ggplot() +
  geom_line(aes(t_pred,  exp(t_pred) , color = 'true phi'), size = 1) +
  #geom_point(aes(x,  y_true ), color = "green", size = 3, alpha = 0.7) +
  geom_line(aes(t_pred, apply(phi_pred_nl,2,median), color = 'fitted phi'), size = 1) +
  geom_ribbon(aes(x = x_pred, ymin = apply(phi_pred_nl, 2, quantile, 0.025), 
                  ymax = apply(phi_pred_nl, 2, quantile, 0.975)), 
              fill = "deeppink", alpha = 0.2) +
  labs(title = "GP Distribution Regression  (non-linear): Link function",x = "t",y = "phi(t)=exp(t)") +
  scale_color_manual(name='',
                     breaks=c('true phi', 'fitted phi'),
                     values=c('true phi'='darkgreen', 'fitted phi'='deeppink'))+
 theme(legend.title=element_text(size=20),
       legend.text=element_text(size=14),legend.position = "bottom")

################################################################################################
################################################################################################

#### Define the Gaussian process regression model in Stan
gp_dist_model_code <- '
data {
  int<lower=1> N;        // Number of subjects
  matrix[N,N] C;         // Input valuesGP dist cov matrix
  vector[N] y;           // Output values
  //real<lower=1e-6> rho;   // GP length scale
  }

parameters {
  vector[N] f;       // GP latent function values
  //real<lower=1e-6> sigma; // Spatial Noise level
  //real<lower=0> sigma; // Noise level
  real<lower=0> tau; // Noise level
  }

model {
  // Prior for latent function values
  f ~ multi_normal(rep_vector(0, N), C);
  
  // Likelihood
  for(i in 1:N){
  y[i] ~ normal(f[i], tau);
  }
  //sigma ~ inv_gamma(2, s);
  //tau ~ inv_gamma(2, s);
  //rho ~ uniform(0.05,0.35);
}'

gp_dist_model <- stan_model(model_code = gp_dist_model_code)

# Fit the linear model to the data
gp_fit <- sampling(gp_dist_model, data = gp_data, chains = 4, iter = 3000)

# Extract fitted parameters
posterior_samples <- as.data.frame(gp_fit)
summary(posterior_samples)

set.seed(1)
subindex=sample(nrow(posterior_samples),100)
posterior_samples_sub=posterior_samples[subindex,]

#posterior_samples_sub$rho[i]

# Predict using the fitted model
f_pred <- Reduce('rbind',lapply(1:length(subindex), function(i) {
  #rho=posterior_samples_sub$rho[i]
  t_data=exp(-0.5 * (outer(x, x, "-") / rho)^2)
  t_cross=exp(-0.5 * (outer(x_pred, x, "-") / rho)^2)
  as.vector(t_cross %*% t(A) %*% solve(A%*%t_data%*%t(A),as.numeric(posterior_samples_sub[i,1:N])))
}),)

# Plot the true and fitted functions
ggplot() +
  geom_line(aes(x,  sin(2 * pi * x) , color = 'true f'), size = 1, alpha = 0.7) +
  #geom_point(aes(x,  y_true ), color = "green", size = 3, alpha = 0.7) +
  geom_line(aes(x_pred, colMeans(f_pred), color = 'fitted f'), size = 1) +
  geom_ribbon(aes(x = x_pred, ymin = apply(f_pred, 2, quantile, 0.025), 
                  ymax = apply(f_pred, 2, quantile, 0.975)), 
              fill = "red", alpha = 0.2) +
  labs(title = "GP Distribution Regression (linear)",x = "x",y = "f(x)") +
  scale_color_manual(name='',
                     breaks=c('true f', 'fitted f'),
                     values=c('true f'='blue', 'fitted f'='red'))+
 theme(legend.title=element_text(size=20),
       legend.text=element_text(size=14),legend.position = "bottom")
