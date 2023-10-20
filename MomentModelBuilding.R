##LIBRARY INVOKATION##

library(Hmisc) # accessing describe()
library(tidyverse) # ggplot, dplyr, and other functionality
library(scales) # aesthetic scaling for ggplots
library(lubridate) # date/timestamp data handling
library(forecast) # forecasting tools
library(qlcal) # trading calendar tools
library(stats) # stats/modelling tools e.g. AIC/BIC
library(highfrequency) # estimating HAR model
library(zoo) # rollmean()
library(xts) # used for highfrequency package's HAR model estimation
library(tseries) # jarque-bera test
library(leaps) #subset selection
library(dgof) #K-S test
library(quarks) # VaR and ES backtesting tools
library(moments) #Tools for calculating higher sample moments
library(PDQutils) #Use for CFE VaR estimation.
library(optimx)
library(stats4) #use for mle()
library(bbmle) #mle2()
library(kableExtra) #nice tables


##SETTING PARAMETERS##

maxTTM <- 252



#default settings
# trainDateRange <- c(as.Date("2013-01-01"), as.Date("2021-04-30"))
# testDateRange <- c(as.Date("2021-05-01"), as.Date("2023-05-01"))

#alternate
# trainDateRange <- c(as.Date("2020-01-01"), as.Date("2023-04-30"))
# testDateRange <- c(as.Date("2017-05-01"), as.Date("2019-12-31"))

#Creating ICE London Business Calendar
setCalendar("UnitedKingdom")


modData <- function(data, RM = 'RV', K = 0, onestep = TRUE) { #adds lagged + smoothed RV estimates and fourier series
  
  mod.data <- data %>% select(c(Date, !!RM, DaystoMaturity, logRet))
  
  #log transform realised vol. & kurt.
  
  if (RM %in% c("RV", "RK", "m4")) {
    
    mod.data[[RM]] <- with(mod.data, log(get(RM)))
    
  }
  
  if (K > 0) {
    
    for (i in 1:K) {
      # create fourier series (current and smoothed) columns up to Kth harmonic. most of these just useful for forecasting later.
      
      mod.data[[paste0("S", as.character(i))]] <- with(mod.data, sin(DaystoMaturity*2*pi*i/maxTTM))
      mod.data[[paste0("S", as.character(i), "_D")]] <- with(mod.data, lag(get(paste0("S", as.character(i)))))
      mod.data[[paste0("S", as.character(i), "_W")]] <- with(mod.data, lag(rev(rollmean(rev(get(paste0("S", as.character(i)))), 5, fill = NA, align = "left"))))
      mod.data[[paste0("S", as.character(i), "_M")]] <- with(mod.data, lag(rev(rollmean(rev(get(paste0("S", as.character(i)))), 22, fill = NA, align = "left"))))
      mod.data[[paste0("C", as.character(i))]] <- with(mod.data, cos(DaystoMaturity*2*pi*i/maxTTM))
      mod.data[[paste0("C", as.character(i), "_D")]] <- with(mod.data, lag(get(paste0("C", as.character(i)))))
      mod.data[[paste0("C", as.character(i), "_W")]] <- with(mod.data, lag(rev(rollmean(rev(get(paste0("C", as.character(i)))), 5, fill = NA, align = "left"))))
      mod.data[[paste0("C", as.character(i), "_M")]] <- with(mod.data, lag(rev(rollmean(rev(get(paste0("C", as.character(i)))), 22, fill = NA, align = "left"))))
    }
    
  }
  
  mod.data[["DaystoMaturity_D"]] <- with(mod.data, lag(DaystoMaturity))
  mod.data[["DaystoMaturity_W"]] <- with(mod.data, lag(rev(rollmean(rev(DaystoMaturity), 5, fill = NA, align = "left"))))
  mod.data[["DaystoMaturity_M"]] <- with(mod.data, lag(rev(rollmean(rev(DaystoMaturity), 22, fill = NA, align = "left"))))
  
  if (onestep) {
    
    mod.data[[paste0(RM, "_D")]] <- with(mod.data, as.double(lag(get(RM)))) # for some reason need to manually stop it from changing to string
    mod.data[[paste0(RM, "_W")]] <- with(mod.data, as.double(lag(rev(rollmean(rev(get(RM)), 5, fill = NA, align = "left")))))
    mod.data[[paste0(RM, "_M")]] <- with(mod.data, as.double(lag(rev(rollmean(rev(get(RM)), 22, fill = NA, align = "left"))))) #assumes 5 trading days in a week and 22 days in a month -> are there big issues when this isn't the case?
    
  }
  
  
  return(mod.data[-c(1:22),]) #remove first 22 rows so no columns are na after monthly smoothing
}


### MODEL FITTING FUNCTIONS ###

TTM.HAR_RM.fit <- function(data, RM = 'RV', K = 0, onestep = TRUE) { #fits TTM-HAR-RM, TTM, HAR models
  #Fits TTM-HAR-RM model
  
  mod.data <- modData(data = data, RM = RM, K = K, onestep = onestep) %>% filter(between(Date, trainDateRange[1], trainDateRange[2]))
  cols <- colnames(mod.data)[-c(grep("^S\\d*_\\w", colnames(mod.data)), grep("^C\\d*_\\w", colnames(mod.data)))] #remove the smoothed trig columns
  
  if (length(cols) == 0) {cols <- colnames(mod.data)}
  
  cols <- cols[!cols %in% c(RM, 'Date', 'DaystoMaturity_D', 'DaystoMaturity_W', 'DaystoMaturity_M', 'logRet')] #remove smoothed TTM + return columns
  
  formula <- as.formula(paste(RM, paste(cols, collapse = '+'), sep = '~'))
  if (onestep) { #if onestep, then modData() has included smoothed RM variables for fitting
    
    models <- lm(formula = formula, data = mod.data, na.action = na.omit)
    
  }
  
  else { #if twostep, RM smoothed variables not included, separately fit TTM model then harmodel onto residuals
    trigmodel <- lm(formula = formula, data = mod.data, na.action = na.omit)
    # tbl <- data.frame(cbind(residuals = trigmodel$residuals, pull(mod.data, paste0(RM, '_D')), pull(mod.data, paste0(RM, '_W')), pull(mod.data, paste0(RM, '_M'))))
    harmodel <- HARmodel(xts(trigmodel$residuals, order.by = mod.data$Date, type = 'HAR', inputType = 'RM')) #check
    models <- c()
    models$trigmodel <- trigmodel
    models$harmodel <- harmodel
    
  }
  
  return(models)
  
}


##TTM-HAR-EPSILON MODEL FITTING FUNCTIONS##

HAR_eps_residuals <- function(data, par, RM, K) #Given parameters, calculates residuals
{
  #assuming par = c(betad, betaw, betam, alpha0, alpha1, alpha2, ..., alpha{2K+1}), 2K + 5 total coefs
  
  if (K > 0) {
  
    resids <- with(data=data, get(RM) - (
      par[4] + par[5] * DaystoMaturity + rowSums(t(par[6:(5+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2)), function(x){get(x)}))))
      + par[1] * (get(paste0(RM, "_D")) - (par[4] + par[5] * DaystoMaturity_D + rowSums(t(par[6:(5+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2), "_D"), function(x){get(x)}))))))
      + par[2] * (get(paste0(RM, "_W")) - (par[4] + par[5] * DaystoMaturity_W + rowSums(t(par[6:(5+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2), "_W"), function(x){get(x)}))))))
      + par[3] * (get(paste0(RM, "_M")) - (par[4] + par[5] * DaystoMaturity_M + rowSums(t(par[6:(5+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2), "_M"), function(x){get(x)}))))))
    ))
  
  } else {
    
    resids <- with(data=data, get(RM) - (
      par[4] + par[5] * DaystoMaturity
      + par[1] * (get(paste0(RM, "_D")) - (par[4] + par[5] * DaystoMaturity_D))
      + par[2] * (get(paste0(RM, "_W")) - (par[4] + par[5] * DaystoMaturity_W))
      + par[3] * (get(paste0(RM, "_M")) - (par[4] + par[5] * DaystoMaturity_M))
    ))
    
  }
  
  return(resids)
  
}

HAR_eps_LL <- function(data, par, RM, K) #loglik function. Calls HAR_eps_residuals, finds RMSE, finds LL
{
  u <- HAR_eps_residuals(data=data, par=par, RM=RM, K=K)
  N <- length(u)
  rmse <- sqrt(sum(u^2)/N)
  
  LL <- -N/2 * log(2*pi) - N*log(rmse) - sum(u^2)/(2*rmse^2)
  return(-LL)
  
  
}

HAR_eps_LL_deriv <- function(data, par, RM, K) #take score matrix to calculate score vector.
  
{
  return (colSums(HAR_eps_LL_score(data=data, par=par, RM=RM, K=K)))
  
}

HAR_eps_LL_score <- function(data, par, RM, K) #calculates a list of dl/dtheta, given the thetas in parameter space defined by K.
{
  u <- HAR_eps_residuals(data=data, par=par, RM=RM, K=K)
  N <- length(u)
  rmse <- sqrt(sum(u^2)/N)
  
  if (K > 0) {
    
    deriv <- matrix(data = c("RV1" = -rmse^(-2) * u * with(data=data, get(paste0(RM, "_D")) - (par[4] + par[5] * DaystoMaturity_D
                                                                                               + rowSums(t(par[6:(5+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2), "_D"), function(x){get(x)})))))),
                             "RV5" = -rmse^(-2) * u * with(data = data, get(paste0(RM, "_W")) - (par[4] + par[5] * DaystoMaturity_W
                                                                                                 + rowSums(t(par[6:(5+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2), "_W"), function(x){get(x)})))))),
                             "RV22" = -rmse^(-2) * u * with(data = data, get(paste0(RM, "_M")) - (par[4] + par[5] * DaystoMaturity_M
                                                                                                 + rowSums(t(par[6:(5+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2), "_M"), function(x){get(x)})))))),
                             "(Intercept)" = -rmse^(-2) * u,
                             "DaystoMaturity" = -rmse^(-2) * u * data$DaystoMaturity
    ), ncol = 5)
    
    for (i in 1:K) {
      
      deriv <- cbind(deriv, -rmse^(-2) * u * with(data=data, get(paste0('S', as.character(i)))))
      deriv <- cbind(deriv, -rmse^(-2) * u * with(data=data, get(paste0('C', as.character(i)))))
      
    }
    
  } else {
    
    deriv <- matrix(data = c("RV1" = -rmse^(-2) * u * with(data=data, get(paste0(RM, "_D")) - (par[4] + par[5] * DaystoMaturity_D)),
                             "RV5" = -rmse^(-2) * u * with(data = data, get(paste0(RM, "_W")) - (par[4] + par[5] * DaystoMaturity_W)),
                             "RV22" = -rmse^(-2) * u * with(data = data, get(paste0(RM, "_M")) - (par[4] + par[5] * DaystoMaturity_M)),
                             "(Intercept)" = -rmse^(-2) * u,
                             "DaystoMaturity" = -rmse^(-2) * u * data$DaystoMaturity
    ), ncol=5)
    
  }
  
  return(deriv)
  
  
}

LL_trig <- function(data, par, RM, K) #LL for TTM-HAR-eps with constraint of alpha_i = 0. Just use TTM.HAR_RM.fit(onestep=TRUE)$trigmodel instead, this is kept for validation of equivalence.
{
  n <- data %>% pull(RM) %>% length()
  #required data:
  # - RV_{t} (predict this)
  # - RV_{t-1} (d, m, w smoothed)
  # - TTM_{t}
  # - TTM_{t-1} (d, m, w smoothed)
  # - S_{i, t}, C_{i, t}
  # - S_{i, t-1}, C_{i, t-1} (d, w, m smoothed)
  #assuming par = c(sigma, alphad, alphaw, alpham, gamma0, gamma1, gamma2, ..., gamma{2K+1}), 2K + 6 total coefs
  #minimising negative loglikelihood function
  
  if (K == 0) {
    
    with(data, (n) * log(par[1]) + (n/2) * log(2*pi) + 1/(2*par[1]^2) * sum((get(RM) - (
      par[5] + par[6] * DaystoMaturity
    ))^2
    ))
    
  }
  
  else {
    
    with(data, (n) * log(par[1]) + (n/2) * log(2*pi) + 1/(2*par[1]^2) * sum((get(RM) - (
      par[5] + par[6] * DaystoMaturity + rowSums(t(par[7:(6+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2)), function(x){get(x)}))))
    ))^2
    ))
    
  }
  
}



TTM.HAR_eps.fit <- function(data, RM = 'RV', K = 0, method = 'L-BFGS-B', control = c(maxit = 1000, reltol=1e-12), lower = -Inf, upper = Inf) #uses predefined LL functions to fit TTM-HAR-eps model
  
{
  # print(paste0("Fitting ", RM, " model @ K = ", as.character(K)))
  startmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = FALSE)
  
  alpha <- startmodel$trigmodel$coefficients
  beta <- startmodel$harmodel$coefficients
  init <- c(beta[-1], beta[1] + alpha[1], alpha[-1]) #initial parameter estimates from previous ols.
  
  
  if (length(lower) == 1){lower <- rep(lower, length(init))}
  if (length(upper) == 1){upper <- rep(upper, length(init))}
  
  #if initial param estimates are outside constrained ranges, set params to equal the lower bound.
  for (i in 1:length(init)) {
    
    if (!between(init[i],lower[i], upper[i])){init[i] <- lower[i]}
    
  }
  
  mod.data <- modData(data = data, RM = RM, K = K, onestep = TRUE) %>% filter(between(Date, trainDateRange[1], trainDateRange[2])) #fit model on training data set
  
  if (method == 'L-BFGS-B') {
    model <- optim(par = init, fn = HAR_eps_LL, gr = HAR_eps_LL_deriv, data = mod.data, RM = RM, K = K, method = method, lower = lower, control = control, hessian=TRUE)
    
  }
  
  else {model <- optim(par = init, fn = HAR_eps_LL, gr = HAR_eps_LL_deriv, data = mod.data, RM = RM, K = K, method = method, control = control, hessian=TRUE)}
  
  model$var <- RM
  score <- HAR_eps_LL_score(data=mod.data, par=model$par, RM=RM, K=K)
  model$bread <- solve(model$hessian)
  model$hessSE <- sqrt(diag(model$bread))
  # model$OPG <- sqrt(diag(solve(t(score) %*% score)))
  model$sandwich <- sqrt(diag(model$bread %*% t(score) %*% score %*% model$bread))

  # print(paste0("Convergence Value: ", as.character(model$convergence)))
  return(model)
  
}

## MODEL CALIBRATION FUNCTIONS ##

fitmeasure <- function(Model, RM = 'RV', onestep = TRUE) {
  #Generates matrix of AIC/SIC/LL values against K for TTM-HAR-RM and constrained cases
  if (onestep) {
    k <- Model$rank
    K <- (k - 5)/2
    n <- length(Model$fitted.values)
    RSS <- sum(Model$residuals^2)
    loglik <- logLik(Model)
    RSq <- summary(Model)$r.squared
    AdjRSq <- summary(Model)$adj.r.squared
    SIC <- k * log(n) - 2*loglik
    AIC <- 2*k - 2*loglik
    
  }
  
  else {
    K <- (Model$trigmodel$rank - 2)/2
    k <- Model$trigmodel$rank + 3
    n <- length(Model$trigmodel$fitted.values)
    RSS <- sum(Model$harmodel$residuals^2)
    TSS <- sum(Model$trigmodel$residuals^2)/(1 - summary(Model$trigmodel)$r.squared) # RSq = 1 - RSS/TSS -> TSS = RSS/(1-RSq)
    RSq <- 1 - RSS/TSS    # using formula R^2 = 1 - RSS/TSS
    AdjRSq <- 1 - (1 - RSq)*(n - 1)/(n - k - 1)
    loglik <- logLik(Model$trigmodel) + logLik(Model$harmodel)
    SIC <- k*log(n) - 2*loglik
    AIC <- 2*k - 2*loglik
    
  }
  
  output <- c(K = K, k = k, n = n, RSS = RSS, loglikelihood = loglik, RSq = RSq, AdjRSq = AdjRSq , SIC = SIC, AIC = AIC)
  return(output)
}

optHarmonics <- function(data, RM = 'RV', onestep = TRUE) {
  #Extends fitmeasure() to create plots of interest and reformats output
  fitmat <- matrix(data = unlist(lapply(X = lapply(X = 0:30, FUN = TTM.HAR_RM.fit, data = data, RM = RM, onestep = onestep), FUN = fitmeasure, onestep = onestep)), byrow = TRUE, ncol = 9)
  colnames(fitmat) <- c('K', 'k', 'n', 'RSS', 'LogLik', 'RSq', 'AdjRSq', 'SIC', 'AIC')
  
  # plot(x = fitmat[,1], y = fitmat[,8], xlab = 'K', ylab = 'SIC', type = 'b', main = 'SIC vs K')
  # plot(x = fitmat[,1], y = fitmat[,9], xlab = 'K', ylab = 'AIC', type = 'b', main = 'AIC vs K')
  # plot(x = fitmat[,1], y = fitmat[,7], xlab = 'K', ylab = 'AdjRSq', type = 'b', main = 'Adj R Squared vs K')
  
  return(fitmat)
}


fitmeasure.eps <- function(Model) {
  #Same as fitmeasure() but for TTM-HAR-eps model.
  k <- length(Model$par)
  K <- (k - 5)/2
  n <- FirstDecNearby %>% modData() %>% filter(between(Date, trainDateRange[1], trainDateRange[2])) %>% pull('RV') %>% length()
  loglik <- -Model$value
  SIC <- k * log(n) - 2*loglik
  AIC <- 2*k - 2*loglik
  RSS <- n*exp((AIC - 2*k)/n) # from AIC = 2k + nln(RSS/n)
  
  output <- c(K = K, k = k, n = n, loglikelihood = loglik, SIC = SIC, AIC = AIC)
  print(output)
  return(output)
  
}

optHarmonics.eps <- function(data, RM = 'RV', KVals = 0:10, method = 'L-BFGS-B', control = c(maxit = 10000, reltol=1e-12))
{
  #Same as optHarmonics() but for TTM-HAR-eps
  fitmat <- matrix(data = unlist(lapply(X = lapply(X = KVals, FUN = TTM.HAR_eps.fit, data = data, RM = RM, method = method, control = control), FUN = fitmeasure.eps)), byrow = TRUE, ncol = 6)
  colnames(fitmat) <- c('K', 'k', 'n', 'LogLik', 'SIC', 'AIC')
  
  # plot(x = fitmat[,1], y = fitmat[,5], xlab = 'K', ylab = 'SIC', type = 'b', main = 'SIC vs K')
  # plot(x = fitmat[,1], y = fitmat[,6], xlab = 'K', ylab = 'AIC', type = 'b', main = 'AIC vs K')
  
  return(fitmat)
  
}

### MOMENT FORECASTING FUNCTIONS ###

MomForecast <- function(model, data, window, RM, spec = 'TTM-HAR.eps')
  
{
  #Forecasts realised moments given some window length and model specification
  if (spec == 'TTM-HAR.eps') {
    
    
    beta <- model$par[1:3]
    alpha <- model$par[-c(1:3)]
    K <- (length(alpha) - 2)/2


    data <- modData(data, K = K, RM = RM) %>% mutate(mu = alpha[1] + alpha[2] * DaystoMaturity
            + if (K > 0) {rowSums(t(alpha[3:(2+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2)), function(x){get(x)}))))} else {0}) #create mu_t column
    
    for (i in 1:22) { #create 22 lagged mu_{t-k} and RM_{t-k} columns
      data[[paste0('mu_', as.character(i))]] <- with(data = data, lag(mu, n = i))
      data[[paste0(RM, '_', as.character(i))]] <- with(data = data, lag(get(RM), n = i))
    }
    
    data <- data[,c('Date', 'DaystoMaturity', 'logRet', 'mu', paste0(RM, '_', c('D', 'W', 'M')), paste0('mu_', as.character(1:22)), paste0(RM, '_', as.character(22:1)), RM)] #52 initial columns
    
    #use the forecast from k days ahead (and previous lagged RM realisations) to create the k+1 day ahead forecast recursively. forecast_{i,t} = E[RM_{t|t-i}]
    beta_coefs <- c(beta[1] + beta[2]/5 + beta[3]/22, rep(beta[2]/5 + beta[3]/22, 4), rep(beta[3]/22, 17))
    
    for (i in 1:window) {
      
      col <- with(data = data, mu - ((beta[1] + beta[2]/5 + beta[3]/22) * lag(mu) + (beta[2]/5 + beta[3]/22) * rowSums(sapply(X = 2:5, FUN = function(a){lag(x = mu, n=a)})) + beta[3]/22 * rowSums(sapply(X = 6:22, FUN = function(a){lag(x = mu, n=a)})))
                  
                  + if (i == 1) {colSums(beta_coefs * t(sapply(X = c(1:22), FUN = function(a){lag(x = get(RM), n=a)})))}
                  
                  else if (i<=22) {colSums(beta_coefs[1:min(i-1, 22)] * t(sapply(X = c(1:min(i-1, 22)), FUN = function(a){lag(x = data[,52+i-a], n=a)}))) + colSums(beta_coefs[-c(1:min(i-1, 22))] * t(sapply(X = c(i:22), FUN = function(a){lag(x = get(RM), n=a)})))}
                  
                  else{colSums(beta_coefs * t(sapply(X = c(1:22), FUN = function(a){lag(x = data[,52+i-a], n=a)})))}
                  
          )
      
      data <- cbind(data, col)
      
    }
    colnames(data)[53:(52+window)] <- c(paste0('forecast_', 1:window))
    
    
  } else if (spec == 'TTM-HAR-RM') {
    
    K <- (length(model$coefficients) - 5)/2
    alpha <- model$coefficients[1:(2*K + 2)]
    beta <- model$coefficients[-c(1:(2*K + 2))]
    
    #create mu_t column
    data <- modData(data, K = K, RM = RM) %>% mutate(mu = alpha[1] + alpha[2] * DaystoMaturity
            + if (K > 0) {rowSums(t(alpha[3:(2+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2)), function(x){get(x)}))))} else {0}) 
    
    data <- data[,c('Date', 'DaystoMaturity', 'mu', 'logRet', RM)]
    for (i in 1:22) { #create 22 lagged RM_{t-k} columns
      data[[paste0(RM, '_', as.character(i))]] <- with(data = data, lag(get(RM), n = i))
    }
    
    beta_coefs <- c(beta[1] + beta[2]/5 + beta[3]/22, rep(beta[2]/5 + beta[3]/22, 4), rep(beta[3]/22, 17))
    
    for (i in 1:window) {
      
      col <- with(data = data, mu + if (i == 1) {colSums(beta_coefs * t(sapply(X = c(1:22), FUN = function(a){lag(x = get(RM), n=a)})))}
                  
                  else if (i<=22) {colSums(beta_coefs[1:(i-1)] * t(sapply(X = c(1:(i-1)), FUN = function(a){lag(x = data[,(27+i-a)], n=a)}))) + colSums(beta_coefs[-c(1:(i-1))] * t(sapply(X = c(i:22), FUN = function(a){lag(x = get(RM), n=a)})))}
                  
                  else{colSums(beta_coefs * t(sapply(X = c(1:22), FUN = function(a){lag(x = data[,(27+i-a)], n=a)})))}
      )
      
      data <- cbind(data, col)
      colnames(data)[27+i] <- paste0('forecast_', i)
      
    }
  
  } else if (spec == 'HAR') {
    
    
    alpha <- model$coefficients[1]
    beta <- model$coefficients[-1]
    
    data <- modData(data, K = 0, RM = RM)
    data <- data[,c('Date', 'DaystoMaturity', 'logRet', RM)]
    for (i in 1:22) { #create 22 lagged RM_{t-k} columns
      data[[paste0(RM, '_', as.character(i))]] <- with(data = data, lag(get(RM), n = i))
    }
    
    beta_coefs <- c(beta[1] + beta[2]/5 + beta[3]/22, rep(beta[2]/5 + beta[3]/22, 4), rep(beta[3]/22, 17))

    for (i in 1:window) {
      
      col <- with(data = data, alpha[1] + if (i == 1) {colSums(beta_coefs * t(sapply(X = c(1:22), FUN = function(a){lag(x = get(RM), n=a)})))}
                  
                  else if (i<=22) {colSums(beta_coefs[1:(i-1)] * t(sapply(X = c(1:(i-1)), FUN = function(a){lag(x = data[,(26+i-a)], n=a)}))) + colSums(beta_coefs[-c(1:(i-1))] * t(sapply(X = c(i:22), FUN = function(a){lag(x = get(RM), n=a)})))}
                  
                  else{colSums(beta_coefs * t(sapply(X = c(1:22), FUN = function(a){lag(x = data[,26+i-a], n=a)})))}
      )
      
      data <- cbind(data, col)
      colnames(data)[26+i] <- paste0('forecast_', i)
      
    }
    
  } else if (spec == 'TTM') {
    
    K <- (length(model$coefficients) - 2)/2
    alpha <- model$coefficients
    
    #create mu_t column
    if (RM %in% c("RV", "RK", "m4")) {
    data <- modData(data, K = K, RM = RM) %>% mutate(forecast_1 = exp(alpha[1] + alpha[2] * DaystoMaturity #note we only do one forecast as this function is deterministic.
                   + if (K > 0) {rowSums(t(alpha[3:(2+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2)), function(x){get(x)}))))} else{0}))
    data[[RM]] <- with(data = data, exp(get(RM)))
    }
    else {
      data <- modData(data, K = K, RM = RM) %>% mutate(forecast_1 = alpha[1] + alpha[2] * DaystoMaturity #note we only do one forecast as this function is deterministic.
                  + if (K > 0) {rowSums(t(alpha[3:(2+2*K)] * t(sapply(paste0(c('S', 'C'), rep(1:K, each = 2)), function(x){get(x)}))))} else{0})
    }
    
    data <- data[,c('Date', 'DaystoMaturity', 'logRet', RM, 'forecast_1')]
    
  }
  
  
  if (RM %in% c("RV", "RK", "m4") & spec != 'TTM') #undo the log transform to get raw RV/RK
  {
    for (i in 1:window) {
      
      data[[paste0('forecast_', i)]] <- with(data = data, exp(get(paste0('forecast_', i))))
      
    }
    
    data[[RM]] <- with(data = data, exp(get(RM)))
  }
  
  return(data[-c(1:22),])
  
}


### HYPOTHESIS (LR) TESTING FUNCTIONS ###

hypTest1 <- function(data, RM = 'RV', K=0)
  
{
  # 1: Do TTM Effects Provide Significant Benefit Over Intercept Term? LR-TEST:
  # H0: TTM coefs = 0
  # H1: TTM coefs =/= 0
  
  trigmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = FALSE)$trigmodel
  minmodel <- modData(data = data, RM = RM, K = 0, onestep = TRUE) %>% filter(between(Date, trainDateRange[1], trainDateRange[2])) %>%
    lm(formula = formula(get(RM) ~ 1))
  teststat <- as.numeric(-2 * (logLik(minmodel) - logLik(trigmodel)))
  
  df <- length(trigmodel$coefficients) - length(minmodel$coefficients)
  
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail = FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
  
}

hypTest2 <- function(data, RM = 'RV', K=0)
  
{
  # 2: Do HAR effects provide significant benefit over intercept term? LR-TEST:
  # H0: HAR effects provide little substantial benefit over intercept term
  # H1: HAR effects do provide substantial benefit
  
  minmodel <- modData(data = data, RM = RM, K = 0, onestep = TRUE) %>% lm(formula = formula(get(RM) ~ 1))
  
  harmodel <- modData(data, RM = RM, K = K, onestep = TRUE) %>%
    filter(between(Date, trainDateRange[1], trainDateRange[2])) %>%
    lm(formula = formula(get(RM) ~ get(paste0(RM, '_D')) + get(paste0(RM, '_W')) + get(paste0(RM, '_M'))))
  
  teststat <- as.numeric(-2 * (logLik(minmodel) - logLik(harmodel)))
  
  df <- length(harmodel$coefficients) - length(minmodel$coefficients)
  
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail = FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
  
}


hypTest3 <- function(data, RM = 'RV', K=0)
  
{
  
  # 3: Do the TTM factors improve the fit over just the HAR model? HAR-RM vs. TTM + HAR-RM LR TEST:
  # H0: HAR model fit reflects best i.e. trig + linear coefs = 0
  # H1: maximal model reflects best i.e. trig + linear coefs =/= 0
  
  maxmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = TRUE)
  harmodel <- modData(data, RM = RM, K = K, onestep = TRUE) %>%
    filter(between(Date, trainDateRange[1], trainDateRange[2])) %>%
    lm(formula = formula(get(RM) ~ get(paste0(RM, '_D')) + get(paste0(RM, '_W')) + get(paste0(RM, '_M'))))
  teststat <- as.numeric(-2 * (logLik(harmodel) - logLik(maxmodel)))
  df <- maxmodel$rank - harmodel$rank
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail=FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
  
}

hypTest4a <- function(data, RM = 'RV', K=0, method='L-BFGS-B')
  
{
  
  # 4a: Do the TTM factors improve the fit over just the HAR-eps model? HAR-eps vs. TTM-HAR-eps LR TEST:
  # H0: HAR-eps model fit reflects best i.e. trig + linear coefs = 0
  # H1: maximal model reflects best i.e. trig + linear coefs =/= 0
  
  maxmodel <- TTM.HAR_eps.fit(data = data, RM = RM, K = K, method=method)
  mod.data <- modData(data, RM = RM, K = K, onestep = TRUE) %>% filter(between(Date, trainDateRange[1], trainDateRange[2]))
  # harmodel <- HARmodel(xts((pull(mod.data,RM) - mean(pull(mod.data,RM))), order.by = mod.data$Date), type = 'HAR', inputType = 'RM')
  harmodel <- mod.data %>% lm(formula = formula(get(RM) ~ get(paste0(RM, '_D')) + get(paste0(RM, '_W')) + get(paste0(RM, '_M'))))
  #we should fit har-eps using mean-detrended RM. however this is pretty much the same anyway and HARmodel has different performance to lm for some reason.
  teststat <- as.numeric(-2 * (logLik(harmodel)[[1]] - (-maxmodel$value)))
  df <- length(maxmodel$par) - harmodel$rank
  
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail=FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
  
}


hypTest4b <- function(data, RM = 'RV', K=0)
  
{
  # 4b: Does the seasonal effect significantly improve the model over the linear TTM effect (excluding HAR model)? LR TEST:
  # H0: trig coefs = 0
  # H1: trig coefs =/= 0
  
  maxmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = FALSE)$trigmodel
  linmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = 0, onestep = FALSE)$trigmodel
  teststat <- as.numeric(-2 * (logLik(linmodel) - logLik(maxmodel)))
  df <- maxmodel$rank - linmodel$rank
  
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail=FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
}



hypTest5 <- function(data, RM = 'RV', K=0)
  
{
  #5: How about in the presence of HAR? Use the TTM-HAR-RM model.
  maxmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = TRUE)
  linmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = 0, onestep = TRUE)
  teststat <- as.numeric(-2 * (logLik(linmodel) - logLik(maxmodel)))
  df <- maxmodel$rank - linmodel$rank
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail=FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
  
}

hypTest6 <- function(data, RM = 'RV', K=0, method='L-BFGS-B')
  
{
  #6: Conduct the same test again but with the TTM-HAR-EPS model.
  maxmodel <- TTM.HAR_eps.fit(data = data, RM = RM, K = K, method=method)
  constrmodel <- TTM.HAR_eps.fit(data = data, RM = RM, K = 0, method=method)
  teststat <- as.numeric(2 * (constrmodel$value - maxmodel$value)) # mult. by 2, not -2 as the "value" is negative loglikelihood
  df <- 2*K #only 1
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail=FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
  
}

hypTest7 <- function(data, RM = 'RV', K=0)
  
{
  #7: Does the TTM-HAR-RM model provide significant improvement over the trig model? LR TEST:
  maxmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = TRUE)
  constrmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = FALSE)$trigmodel
  teststat <- as.numeric(-2 * (logLik(constrmodel) - (logLik(maxmodel))))
  df <- 3
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail=FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
  
}



hypTest8 <- function(data, RM = 'RV', K=0, method='L-BFGS-B') {
  
  # 8: Additionally, does the TTM-HAR-EPS model provide significant improvement over the trig model?
  # LR test: Max model vs. TTM model. Are HAR coefficients significantly non-zero? 
  maxmodel <- TTM.HAR_eps.fit(data = data, RM = RM, K = K, method=method)
  constrmodel <- TTM.HAR_RM.fit(data = data, RM = RM, K = K, onestep = FALSE)$trigmodel
  teststat <- as.numeric(-2 * (logLik(constrmodel) - (-maxmodel$value)))
  df <- 3
  return(c('teststat' = round(teststat, digits=2), 'critval99' = qchisq(0.01, df, lower.tail=FALSE), 'pval' = round(pchisq(teststat, df, lower.tail = FALSE), digits=4)))
}


hypTest9 <- function(data, RM = 'RV', K=0, method = 'L-BFGS-B') {
  
  
  # 9: Is TTM-HAR-EPS MLE significantly better/worse than TTM-HAR-RM? How about original 2step TTM-HAR-EPS? Relative AIC Ratio:
  model1 <- TTM.HAR_RM.fit(data = data, RM = RM, K = K)
  AIC1 <- 2*(5 + 2*K) - 2*logLik(model1)[1]
  model2 <- TTM.HAR_eps.fit(data = data, RM = RM, K = K, method=method)
  AIC2 <- 2*(5 + 2*K) - 2*(-model2$value)
  teststat <- exp(0.5 * -abs(AIC1 - AIC2))
  return(c('AIC1' = AIC1, 'AIC2' = AIC2, 'teststat' = teststat))
}