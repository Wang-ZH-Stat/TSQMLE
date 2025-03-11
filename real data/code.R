library(dplyr)
library(lubridate)
library(stabledist)
library(numDeriv)
library(Rsolnp)
library(robustbase)
library(car)
library(evd)

data_raw <- read.csv("./yield-curve-rates-1990-2024.csv")
data_raw <- data_raw[,c(1,4)]
data_raw$Date <- mdy(data_raw$Date)
data_raw <- data.frame(data_raw)
names(data_raw)[2] <- "Yield"
data_raw <- data_raw %>% arrange(desc(row_number()))

## calculate the end of month value
monthly_last <- data_raw %>%
group_by(year = year(Date), month = month(Date)) %>%
  summarise(EndOfMonthYield = last(Yield)) %>%
  ungroup() %>%
  mutate(Date = make_date(year, month, 1)) %>%
  select(Date, EndOfMonthYield)

## get data
start_year <- 1990
start_month <- 1
end_year <- 2024
end_month <- 12

data <- ts(monthly_last$EndOfMonthYield, 
           start=c(start_year, start_month), frequency=12)
y <- diff(data)

par(mfrow=c(1,2), las = 1)
plot(data, type="l", xlab = "", main = "", 
     ylab = expression(x[t]))
title("(a)", adj = 0)
plot(y, type="l", xlab = "", main = "", 
     ylab = expression(y[t]))
title("(b)", adj = 0)

f_fun <- function(x){
  x <- abs(x)
  fx <- exp(-x)/(1+exp(-x))^2
  
  return(fx)
}

######################
### DAR
######################
d <- 4
n <- length(y) - 1

sigma_fun <- function(yy, theta){
  phi0 <- theta[1]
  phi1 <- theta[2]
  alpha0 <- theta[3]
  alpha1 <- theta[4]
  sigma_t <- (alpha0 + alpha1*yy^2)^(1/2)
  
  return(sigma_t)
}

g_fun <- function(yy, theta){
  phi0 <- theta[1]
  phi1 <- theta[2]
  alpha0 <- theta[3]
  alpha1 <- theta[4]
  g_t <- phi0 + phi1*yy
  
  return(g_t)
}

q_fun <- function(theta){
  n <- length(y) - 1
  ell_t <- -log(sigma_fun(y[1:n], theta)) + 
    log(f_fun((y[2:(n+1)]-g_fun(y[1:n], theta))/sigma_fun(y[1:n], theta)))
  
  return(ell_t)
}

L_log <- function(theta){
  L <- -sum(q_fun(theta))
  
  return(L)
}
theta_DAR <- try(optim(c(0.5,0.5,0.5,0.5), fn = L_log)$par)
## parameters
round(theta_DAR, 4)

A_hat <- matrix(0, nrow = d, ncol = d)
B_hat <- matrix(0, nrow = d, ncol = d)
for(t in c(1:n)){
  qt_fun <- function(theta){
    n <- length(y) - 1
    ell <- -log(sigma_fun(y[1:n], theta)) + 
      log(f_fun((y[2:(n+1)]-g_fun(y[1:n], theta))/sigma_fun(y[1:n], theta)))
    ell_t <- ell[t]
    
    return(ell_t)
  }
  st_fun <- grad(qt_fun, theta_DAR)
  ht_fun <- -hessian(qt_fun, theta_DAR)
  A_hat <- A_hat + ht_fun
  B_hat <- B_hat + st_fun%*%t(st_fun)
}
A_hat <- A_hat / n
B_hat <- B_hat / n

## ASD
ASD_DAR <- solve(A_hat)%*%B_hat%*%solve(A_hat)/n
round(sqrt(diag(ASD_DAR)), 3)

## p_value
p_value_DAR <- 1 - pnorm(abs(theta_DAR / sqrt(diag(ASD_DAR))))
round(2*p_value_DAR, 3)

## Lyapunov exponent
eta_t_DAR <- (y[2:(n+1)]-g_fun(y[1:n], theta_DAR))/sigma_fun(y[1:n], theta_DAR)
gamma_DAR <- mean(log(abs(theta_DAR[2]+eta_t_DAR*theta_DAR[4]^(1/2))))
round(gamma_DAR, 4)

## Hill estimators
sorted_eta <- sort(eta_t_DAR)
k <- 100
top_k_values <- sorted_eta[(length(sorted_eta) - k + 1):length(sorted_eta)]
hill_DAR <- mean(log(top_k_values)) - log(sorted_eta[length(sorted_eta) - k])
round(hill_DAR, 4)

##  Log-Likelihood
LL_DAR <- -L_log(theta_DAR)
round(LL_DAR, 3)

## AIC
AIC_DAR <- 2*d + 2*L_log(theta_DAR)
round(AIC_DAR, 3)



######################
### GARCH
######################
d <- 3

sigma_fun <- function(y, theta){
  n <- length(y) - 1
  alpha0 <- theta[1]
  alpha1 <- theta[2]
  beta1 <- theta[3]
  sigma_t <- c()
  sigma_t[1] <- 0
  for(i in c(2:(n+1))){
    sigma_t[i] <- (alpha0 + alpha1*y[i-1]^2 + beta1*sigma_t[i-1]^2)^(1/2)
  }
  
  return(sigma_t[2:(n+1)])
}

q_fun <- function(theta){
  n <- length(y) - 1
  ell_t <- -log(sigma_fun(y, theta)) + 
    log(f_fun(y[2:(n+1)]/sigma_fun(y, theta)))
  
  return(ell_t)
}

L_log <- function(theta){
  L <- -sum(q_fun(theta))
  
  return(L)
}

theta_GARCH <- try(optim(rep(0.5, d), fn = L_log, method = "L-BFGS-B",
                         lower = rep(0.001, d))$par)
## parameters
round(theta_GARCH, 4)

A_hat <- matrix(0, nrow = d, ncol = d)
B_hat <- matrix(0, nrow = d, ncol = d)
for(t in c(1:n)){
  print(t)
  qt_fun <- function(theta){
    n <- length(y) - 1
    ell <- -log(sigma_fun(y, theta)) + 
      log(f_fun(y[2:(n+1)]/sigma_fun(y, theta)))
    ell_t <- ell[t]
    
    return(ell_t)
  }
  st_fun <- grad(qt_fun, theta_GARCH)
  ht_fun <- -hessian(qt_fun, theta_GARCH)
  A_hat <- A_hat + ht_fun
  B_hat <- B_hat + st_fun%*%t(st_fun)
}
A_hat <- A_hat / n
B_hat <- B_hat / n

## ASD
ASD_GARCH <- solve(A_hat)%*%B_hat%*%solve(A_hat)/n
round(sqrt(diag(ASD_GARCH)), 3)

## p_value
p_value_GARCH <- 1 - pnorm(abs(theta_GARCH / sqrt(diag(ASD_GARCH))))
round(2*p_value_GARCH, 3)

## Lyapunov exponent
eta_t_GARCH <- y[2:(n+1)]/sigma_fun(y, theta_GARCH)
gamma_GARCH <- mean(log(abs(eta_t_GARCH*theta_GARCH[2]^(1/2) + theta_GARCH[3])))
round(gamma_GARCH, 4)

## Hill estimators
sorted_eta <- sort(eta_t_GARCH)
k <- 100
top_k_values <- sorted_eta[(length(sorted_eta) - k + 1):length(sorted_eta)]
hill_GARCH <- mean(log(top_k_values)) - log(sorted_eta[length(sorted_eta) - k])
round(hill_GARCH, 4)

##  Log-Likelihood
LL_GARCH <- -L_log(theta_GARCH)
round(LL_GARCH, 3)

## AIC
AIC_GARCH <- 2*d + 2*L_log(theta_GARCH)
round(AIC_GARCH, 3)



######################
### ARMA-GARCH
######################
d <- 6

g_fun <- function(y, theta){
  n <- length(y) - 1
  phi0 <- theta[1]
  phi1 <- theta[2]
  psi1 <- theta[3]
  alpha0 <- theta[4]
  alpha1 <- theta[5]
  beta1 <- theta[6]
  varep_t <- c()
  g_t <- c()
  varep_t[1] <- 0
  g_t[1] <- 0
  for(t in c(2:(n+1))){
    varep_t[t] <- y[t] - phi0 - phi1*y[t-1] - psi1*varep_t[t-1]
    g_t[t] <- phi0 + phi1*y[t-1] + psi1*varep_t[t-1]
  }
  
  return(g_t[2:(n+1)])
}

sigma_fun <- function(y, theta){
  n <- length(y) - 1
  phi0 <- theta[1]
  phi1 <- theta[2]
  psi1 <- theta[3]
  alpha0 <- theta[4]
  alpha1 <- theta[5]
  beta1 <- theta[6]
  varep_t <- c()
  sigma_t <- c()
  varep_t[1] <- 0
  sigma_t[1] <- 0
  for(t in c(2:(n+1))){
    varep_t[t] <- y[t] - phi0 - phi1*y[t-1] - psi1*varep_t[t-1]
    sigma_t[t] <- (alpha0 + alpha1*varep_t[t-1]^2 + beta1*sigma_t[t-1]^2)^(1/2)
  }
  
  return(sigma_t[2:(n+1)])
}

q_fun <- function(theta){
  n <- length(y) - 1
  ell_t <- -log(sigma_fun(y, theta)) + 
    log(f_fun((y[2:(n+1)]-g_fun(y, theta))/sigma_fun(y, theta)))
  
  return(ell_t)
}

L_log <- function(theta){
  L <- -sum(q_fun(theta))
  
  return(L)
}
theta_AG <- try(optim(rep(0.5, d), fn = L_log, method = "L-BFGS-B",
                      lower = c(-Inf, -Inf, -Inf, 0.001, 0.001, 0.001))$par)
## parameters
round(theta_AG, 4)

A_hat <- matrix(0, nrow = d, ncol = d)
B_hat <- matrix(0, nrow = d, ncol = d)
for(t in c(1:n)){
  print(t)
  qt_fun <- function(theta){
    n <- length(y) - 1
    ell <- -log(sigma_fun(y, theta)) + 
      log(f_fun((y[2:(n+1)]-g_fun(y, theta))/sigma_fun(y, theta)))
    ell_t <- ell[t]
    
    return(ell_t)
  }
  st_fun <- grad(qt_fun, theta_AG)
  ht_fun <- -hessian(qt_fun, theta_AG)
  A_hat <- A_hat + ht_fun
  B_hat <- B_hat + st_fun%*%t(st_fun)
}
A_hat <- A_hat / n
B_hat <- B_hat / n

## ASD
ASD_AG <- solve(A_hat)%*%B_hat%*%solve(A_hat)/n
round(sqrt(diag(ASD_AG)), 3)

## p_value
p_value_AG <- 1 - pnorm(abs(theta_AG / sqrt(diag(ASD_AG))))
round(2*p_value_AG, 3)

## Lyapunov exponent
eta_t_AG <- (y[2:(n+1)]-g_fun(y, theta_AG))/sigma_fun(y, theta_AG)
h <- numeric(length(eta_t_AG))
h[1] <- var(eta_t_AG)
epsilon <- eta_t_AG - theta_AG[1]
J <- matrix(0, nrow = length(eta_t_AG), ncol = 2)
J[1,] <- c(1, 0)

for (t in 2:length(eta_t_AG)) {
  h[t] <- theta_AG[4] + theta_AG[5] * epsilon[t-1]^2 + theta_AG[6] * h[t-1]
  epsilon[t] <- eta_t_AG[t] - theta_AG[1] - theta_AG[2] * eta_t_AG[t-1] - 
    theta_AG[3] * epsilon[t-1]
  J[t,] <- c(theta_AG[5] * epsilon[t-1]^2 + theta_AG[6] * J[t-1,1], theta_AG[6] * J[t-1,2])
}
gamma_AG <- mean(log(abs(J[,1])))
round(gamma_AG, 4)

## Hill estimators
sorted_eta <- sort(eta_t_AG)
k <- 100
top_k_values <- sorted_eta[(length(sorted_eta) - k + 1):length(sorted_eta)]
hill_AG <- mean(log(top_k_values)) - log(sorted_eta[length(sorted_eta) - k])
round(hill_AG, 4)

##  Log-Likelihood
LL_AG <- -L_log(theta_AG)
round(LL_AG, 3)

## AIC
AIC_AG <- 2*d + 2*L_log(theta_AG)
round(AIC_AG, 3)


par(mfrow=c(1,3), las = 1)

hist(eta_t_DAR, breaks = 16, xlim = c(-10,10),
     ylab = "", xlab="Residuals", main = "",
     freq = FALSE)
title("(a) DAR", adj = 0)

hist(eta_t_GARCH, breaks = 16, xlim = c(-10,10),
     ylab = "", xlab="Residuals", main = "",
     freq = FALSE)
title("(b) GARCH", adj = 0)

hist(eta_t_AG, breaks = 16, xlim = c(-10,10),
     ylab = "", xlab="Residuals", main = "",
     freq = FALSE)
title("(c) ARMA-GARCH", adj = 0)


save.image("results.Rdata")
# load("results.Rdata")

