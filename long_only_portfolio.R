##################################################
## Project: Financial case study APG: Covariance estimation methods
## Script purpose: Estimate, and evaluate covariance estimation methods
## Methods include: DCC-GARCH, DECO-GARCH, Sample estimates, Addaptive tresholding, 1/N
## Evaluation include: Sharpe-ratios, Volatility, Turnover
## Date: 28-01-21
## Author: Karan Samlal
##################################################

# Clear variables
rm(list = ls())

#----------------------------------------------------------
#----Setting up Packages-----------------------------------
#----------------------------------------------------------

library(xtable)
library(psych)
library(rmgarch)
library(rugarch)
library(readxl)
library(fitdistrplus)
library(fGarch)
library(fitHeavyTail)
library(fMultivar)
library(cramer)
library(sn)
library(fBasics)
library(SimDesign)
library(Peacock.test)
library(latex2exp)
library(xdcclarge)
library(MSwM)
library(zoo)
library(forecast)
library(caret)
library(CovTools)
library(PortfolioAnalytics)
library(xts)
library(Quandl)
library(Rsolnp)

# ommit scientific notation
options(scipen=999)

#----------------------------------------------------------
#----Data analysis----------------------------------------
#----------------------------------------------------------

# read in data
All_data <- readRDS(file = "2021_01_comm_data_corrected.RDS")
All_data

# read sp500 data
sp500.data <- read.csv("^GSPC (5).csv")

# daily return data
ret_data <- All_data$returns
num_assets <- 27
ret <- ret_data[,1:num_assets] * 100

# obtain the in-sample length
start.date.in = as.Date("2016-01-04")
end.date.in = as.Date("2019-12-31")
in.sample_period = ret[which(rownames(data.frame(ret)) == start.date.in):
                         which(rownames(data.frame(ret)) == end.date.in),]
in_len  = length(in.sample_period$W)
in_len

# subset the data to the full sample w.r.t. the panel you have
start.date.full = as.Date("2016-01-04")
end.date.full = as.Date("2020-12-31")

# define full sample
ret_sample = ret[which(rownames(data.frame(ret)) == start.date.full):
                   which(rownames(data.frame(ret)) == end.date.full),]

# the full sample length
size_sample <- length(ret_sample$W)

# define out-of-sample sp500 returns
sp500.data.sample = sp500.data[(which(rownames(data.frame(ret)) == start.date.full) - 1):
                                 which(rownames(data.frame(ret)) == end.date.full),]
sp500.returns.sample = diff(log(sp500.data.sample$Adj.Close)) * 100

# get index returns
index_weights_sample = All_data$weights[which(rownames(data.frame(All_data$weights)) == start.date.full):
                                          which(rownames(data.frame(All_data$weights)) == end.date.full),]
index.returns <- c()
for (t in 1:length(ret_sample$W)){
  index.returns[t] = (index_weights_sample[t,]) %*% t(ret_sample[t,])
}

# read data about covariance matrix of Factor model 
cov.factor = read_excel("Cov_period1.xlsx", col_names = FALSE)
cov.factor = data.frame(cov.factor)
names(cov.factor) <- colnames(ret_sample)
cov.factor

#----------------------------------------------------------
#----DCC-GARCH----------------------------------------
#----------------------------------------------------------

# in-sample and number of out-of-sample observations
in_sample_size <- in_len
num_forc <- size_sample - in_sample_size

# define the out-of-sample return series
out_of_sample = ret_sample[(in_sample_size + 1):size_sample, ]

# out-of-sample dates
out_sample_dates = as.Date(rownames(data.frame(out_of_sample)))

# univariate normal GARCH(1,1) for each series
# use GJR-GARCH and a t-distribution
garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                          variance.model = list(garchOrder = c(1,1), 
                                                model = "sGARCH"), 
                          distribution.model = "norm")

# dcc specification - GARCH(1,1) for conditional correlations
dcc.garch11.spec = dccspec(uspec = multispec( replicate(num_assets, garch11.spec) ), 
                           dccOrder = c(1,1), 
                           distribution = "mvnorm")

dcc.fit = dccfit(dcc.garch11.spec, data = ret_sample, out.sample = num_forc, fit.control = list(stationarity = 1,
                                                                                                scale=1),
                 solver = "solnp")
dcc.fcst = dccforecast(dcc.fit, n.ahead=1, n.roll = (num_forc - 1))

# obtain the forecasted covariance and correlation matrices from the DCC model
# store them here
dcc_cov <- array(dim = c(num_assets,num_assets,num_forc)) 
dcc_cor <- array(dim = c(num_assets,num_assets,num_forc)) 
for (t in 1:num_forc){
  dcc_cov[,,t] = dcc.fcst@mforecast$H[[t]]
  dcc_cor[,,t] = dcc.fcst@mforecast$R[[t]]
}

# point where we start forecasting (where the out-of-sample starts)
T <- dcc.fit@model$modeldata$T
T

min_var <- function(w, Sigma){
  result = t(w) %*% Sigma %*%  w
  return(result)
}
# make sure weight sum up to one
equal <- function(w, Sigma) {
  z1 <- (sum(w))
  return(z1)
}
dcc_garch_gmv_multistep <- function(horizon, forc.cov.mat){
  
  # Store portfolio returns
  portfolio_returns.gmv <- c()
  
  # Initialize a vector of ones
  vOnes <- rep(1, num_assets)
  
  # keep track of the out-of-sample day
  day = 1
  
  # store the estimated weights
  dcc.weights <- matrix(nrow = (num_forc), ncol = num_assets)
  
  # the t at which we are going to rebalance the portfolio
  rebalancing_periods = seq(from = 1, to = (num_forc), by=horizon)
  
  # store absolute weight deviation when to rebalance the portfolio
  # so the adjustment  we have to make to get to the rebalanced portolio
  weights_diff = matrix(nrow = (length(rebalancing_periods)), ncol = num_assets)
  
  # track the rebalancing days
  rebalance_day = 1
  
  for (t in seq(from = 1, to = (num_forc), by=horizon)){
    Cov.mat <- forc.cov.mat[,,t]
    #pi.gmv.start <- (solve(Cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(Cov.mat) %*% vOnes)[1])
    result.longonly <-  solnp(pars =  rep(1/27, 27), Sigma = Cov.mat,
                              fun = min_var,
                              #function to optimise
                              eqfun=equal, #equality function 
                              eqB=1,   #the equality constraint
                              LB=rep(0, 27), #lower bound for parameters i.e. greater than zero
                              UB=rep(10, 27), control = list(trace = 0, tol = 1e-6)) #upper bound for parameters (I just chose 100 randomly)
    pi.gmv.start = result.longonly$pars
    print(result.longonly$convergence)
    
    if (t > 1){
      weights_diff[rebalance_day,] = t(abs(pi.gmv - pi.gmv.start ))
      rebalance_day = rebalance_day + 1
    }
    
    # so if we are at day 1, we just set up the portfolio instead of rebalancing
    else{
      weights_diff[rebalance_day,] = rep(0, 27)
      rebalance_day = rebalance_day + 1
    }
    # start with de optmized portolio 
    pi.gmv = pi.gmv.start
    
    for (i in 1:horizon){
      # check if do not are out of bound
      if ((T + t + i - 1) > size_sample){
        break
      }
      
      # store the weights
      dcc.weights[day,] = pi.gmv
      
      # calculate portfolio return using the dcc one step ahead forecasted covariance matrix
      out_sample_obs <- ret_sample[(T + t + i - 1), ]
      ret.gmv <- t(pi.gmv) %*% t(out_sample_obs)
      #print(ret.gmv)
      portfolio_returns.gmv[day] = sum(ret.gmv)
      
      # update the weights of the portfolio
      # they change troughout time
      pi.gmv = pi.gmv * t(1  + out_sample_obs/ 100)
      pi.gmv = pi.gmv / (sum(pi.gmv)) # scale to one
      
      # adjust the out-of-sample-day
      day = day + 1
      
    }
  }
  # discard NA
  #portfolio_returns.gmv <- portfolio_returns.gmv[!is.na(portfolio_returns.gmv)]
  return(list(portfolio_returns.gmv, dcc.weights, weights_diff ))
}

# obtain estimated weighs, and returns 
rebalance_horizon = 252
dcc.result = dcc_garch_gmv_multistep(rebalance_horizon, dcc_cov)
dcc.port.ret = dcc.result[[1]]
dcc.port.weights = dcc.result[[2]]
dcc.weight.diff = dcc.result[[3]]
dcc.weight.diff.sum = rowSums(dcc.weight.diff)

#----------------------------------------------------------
#----DECO-GARCH----------------------------------------
#----------------------------------------------------------

# calculate equicorrelation, obtain the mean correlation mu rho
equicorr_fun <- function(){
  equicorr <- t(dcc_cor[,,1][!upper.tri(dcc_cor[,,1], diag = TRUE)])
  for (i in 2:num_forc){
    new_equicorr = t(dcc_cor[,,i][!upper.tri(dcc_cor[,,i], diag = TRUE)])
    equicorr =  rbind(equicorr, new_equicorr)
    
  }
  return(equicorr)
}
equicorr <- equicorr_fun()
mean_rho <- apply(equicorr,1, mean) # for later

# store DECO and BLOCK-DECO covariances out-of-sample
dec_cov <- array(dim = c(num_assets,num_assets,num_forc)) # make room for the covariances of DECO
bdec_cov <- array(dim = c(num_assets,num_assets,num_forc)) # make room for the covariances of DECO

# extract the forecasted dcc covariance and correlation out-of-sample
dcccov <- rcov(dcc.fcst) 
dcccor <- rcor(dcc.fcst)

# check out 

# apply DECO method to the dcc results
for (i in 1:num_forc){
  # forecasted dcc covariance matrix
  cov.mat.dcc = dcccov[[i]]
  cov.mat.dcc = cov.mat.dcc[,,1]
  
  # diagonal matrix of dcc covariance
  dt <- sqrt(cov.mat.dcc*diag(num_assets))
  
  # store the equicorrelations
  eq_cor <- diag(num_assets)
  
  # set the equicorrelation to the mean rho at time i
  eq_cor[upper.tri(eq_cor,diag=F)] <- mean_rho[i]
  eq_cor[lower.tri(eq_cor,diag=F)] <- mean_rho[i]
  
  # estimate the DECO covariance matrix 
  dec_cov[,,i] <- dt %*% eq_cor %*% dt
  
  # estimate the BLOCK-DECO covariance matrix
  cor.mat.dcc = dcccor[[i]]
  cor.mat.dcc = cor.mat.dcc[,,1]
  block_deco_eqcorr = BLOCK.DECO.GARCH(cor.mat.dcc)
  bdec_cov[,,i] <- dt %*% block_deco_eqcorr %*% dt
}

# obtain estimated weighs, and returns for DECO-GARCH 
deco.result = dcc_garch_gmv_multistep(rebalance_horizon, dec_cov)
deco.port.ret = deco.result[[1]]
deco.port.weights = deco.result[[2]]
deco.weight.diff = deco.result[[3]]
deco.weight.diff.sum = rowSums(deco.weight.diff)

# obtain estimated weighs, and returns for DECO-GARCH 
bdeco.result = dcc_garch_gmv_multistep(rebalance_horizon, bdec_cov)
bdeco.port.ret = bdeco.result[[1]]
bdeco.port.weights = bdeco.result[[2]]
bdeco.weight.diff = bdeco.result[[3]]
bdeco.weight.diff.sum = rowSums(bdeco.weight.diff)


### In-sample performances

# DCC in-sample covariace matrices
dcc.cov.insample = dcc.fit@mfit$H
dcc.cor.insample = dcc.fit@mfit$R

# calculate equicorrelation, obtain the mean correlation mu rho
equicorr_fun_insample <- function(cor){
  equicorr <- t(cor[[1]][!upper.tri(cor[[1]], diag = TRUE)])
  for (i in 2:T){
    new_equicorr = t(cor[[i]][!upper.tri(cor[[i]], diag = TRUE)])
    equicorr =  rbind(equicorr, new_equicorr)
    
  }
  return(equicorr)
}
equicorr <- equicorr_fun_insample(dcc.cor.insample)
mean_rho <- apply(equicorr,1, mean) # for later

# DECO and BLOCK-DECO covariance matrices
dec_cov_insample <- array(dim = c(num_assets,num_assets,T)) # make room for the covariances of DECO
bdec_cov_insample <- array(dim = c(num_assets,num_assets,T)) # make room for the covariances of DECO

# apply DECO method to the dcc results
for (i in 1:T){
  # in-sample DCC covariance matrix
  cov.mat.dcc <- dcc.cov.insample[,,i]
  
  # diagonal matrix of dcc covariance
  dt <- sqrt(cov.mat.dcc*diag(num_assets))
  
  # store the equicorrelations
  eq_cor <- diag(num_assets)
  
  # set the equicorrelation to the mean rho at time i
  eq_cor[upper.tri(eq_cor,diag=F)] <- mean_rho[i]
  eq_cor[lower.tri(eq_cor,diag=F)] <- mean_rho[i]
  
  # # In-sample DECO covariance matrix
  dec_cov_insample[,,i] <- dt %*% eq_cor %*% dt
  
  # In-sample BLOCK-DECO covariance matrix
  temp = dcc.cor.insample[[i]]
  
  # make sure the correlation matrix has row, and colum names for the block-deco function (which are the commodities names)
  names(temp) <- colnames(ret_data)
  rownames(temp) = colnames(ret_data)
  
  block_deco_eqcorr = BLOCK.DECO.GARCH(temp)
  bdec_cov_insample[,,i] <- dt %*% block_deco_eqcorr %*% dt
}

#----------------------------------------------------------
#----Addaptive tresholding and Sample Estimates------------
#----------------------------------------------------------

# store Addaptive tresholding forecasted covariance
addaptive_cov <- array(dim = c(num_assets,num_assets,num_forc)) 

# store Addaptive tresholding forecasted covariance
sample_cov <- array(dim = c(num_assets,num_assets,num_forc))

cross_validation_addaptive <- function(lambda_range, nCV, est_sample){
  # split data into test and training set
  train_samples = c(1, round(T/6), 2*round(T/6), 3*round(T/6), 4*round(T/6), 5*round(T/6))
  vOnes <- rep(1, num_assets)
  lambda_cv_scores = c()
  ind_lambda = 1
  
  for (lambda in lambda_range){
    
    # store volatiltiies
    vols = c()
    
    for (grid in 1:nCV){
      # store returns
      returns = c()
      
      # define the train dataset for the first split
      train = est_sample[1:train_samples[(grid+1)],]
      
      # estimate the covariance matrix on the train dataset using the current lambda
      cov.mat = CovEst.adaptive(train, thr= lambda) 
      cov.mat = cov.mat$S
      
      # compute validation returns without re-estimating the covariance matrix
      # if we have a new train dataset, we will estimate a covariance matrix on that and evaluate on new test set
      for (t in 1:(round((T / 6)))){
        if ((train_samples[(grid+1)] + t) > T){
          break
        }
        out_sample = est_sample[(train_samples[(grid+1)] + t),]
        #gmv_weights = (solve(cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(cov.mat) %*% vOnes)[1])
        
        result.longonly <-  solnp(pars =  rep(1/27, 27), Sigma = cov.mat,
                                  fun = min_var,
                                  #function to optimise
                                  eqfun=equal, #equality function 
                                  eqB=1,   #the equality constraint
                                  LB=rep(0, 27), #lower bound for parameters i.e. greater than zero
                                  UB=rep(4, 27), control = list(trace = 0, tol = 1e-6)) #upper bound for parameters (I just chose 100 randomly)
        gmv_weights = result.longonly$pars
        returns[t] = t(gmv_weights) %*%  t(out_sample)
        
      }
      # store the volatility obtained from the validation set 
      vols[grid] = sd(returns) * sqrt(252)
    }
    # track which lambda we used
    lambda_cv_scores[ind_lambda] = mean(vols)
    ind_lambda = ind_lambda + 1
  }
  return(lambda_cv_scores)
}

# estimation window for the static sparse estimators
estimation_window_static = T

# define the grid we will search the best threshold value from
mthr <- seq(from=0.01,to=0.99,length.out=20)

# estimate the covariance matrices and get the static estimators
est_sample = ret_sample[1:T,]

# use 5-fold cross validation with a time-series split
lambda_cv_scores = cross_validation_addaptive(mthr, 5, est_sample)
treshold.star = mthr[which.min(lambda_cv_scores)]

# estimate the static (sparse) covariance matrices
est_cov.mat.addaptive = CovEst.adaptive(est_sample, thr= treshold.star)
est_cov.mat.addaptive = est_cov.mat.addaptive$S
est_cov.mat.sample = cov(est_sample)

# store them
for (i in 1:num_forc){
  # store the estimated covariance matrices without including new data
  addaptive_cov[,,i] <- est_cov.mat.addaptive
  sample_cov[,,i] <- est_cov.mat.sample
  
}

# store results for the static methods 
addaptive.result = dcc_garch_gmv_multistep(rebalance_horizon, addaptive_cov)
addaptive.port.ret = addaptive.result[[1]]
addaptive.port.weights = addaptive.result[[2]]
addaptive.weight.diff = addaptive.result[[3]]
addaptive.weight.diff.sum = rowSums(addaptive.weight.diff)

# obtain estimated weighs, and returns for sample estimates
sample.result = dcc_garch_gmv_multistep(rebalance_horizon, sample_cov)
sample.port.ret = sample.result[[1]]
sample.port.weights = sample.result[[2]]
sample.weight.diff = sample.result[[3]]
sample.weight.diff.sum = rowSums(sample.weight.diff)


#----------------------------------------------------------
#----1/N---------------------------------------------------
#----------------------------------------------------------
equally_weighted_multistep <- function(horizon){
  
  # Store portfolio returns
  portfolio_returns.gmv <- c()
  
  # Initialize a vector of ones
  vOnes <- rep(1, num_assets)
  
  # keep track of the out-of-sample day
  day = 1
  
  # store the estimated weights
  dcc.weights <- matrix(nrow = (num_forc), ncol = num_assets)
  w0 = rep((1/num_assets), num_assets)
  
  # the t at which we are going to rebalance the portfolio
  rebalancing_periods = seq(from = 1, to = (num_forc), by=horizon)
  
  # store weight differences
  weights_diff = matrix(nrow = (length(rebalancing_periods)), ncol = num_assets)
  
  # track the rebalancing days
  rebalance_day = 1
  
  for (t in seq(from = 1, to = (num_forc), by=horizon)){
    if (t > 1){
      weights_diff[rebalance_day,] = t(abs(pi.equal - w0))
      rebalance_day = rebalance_day + 1
    }
    
    # so if we are at day 1, we just set up the portfolio instead of rebalancing
    else{
      weights_diff[rebalance_day,] = rep(0, 27)
      rebalance_day = rebalance_day + 1
    }
    # we do rebalancing such that the weights are 
    # now 1/N again
    pi.equal = w0
    for (i in 1:horizon){
      # check if do not are out of bound
      if ((T + t + i - 1) > size_sample){
        break
      }
      # store the weights of each day
      dcc.weights[day,] = (pi.equal)
      
      # calculate portfolio return using the dcc one step ahead forecasted covariance matrix
      out_sample_obs <- ret_sample[(T + t + i - 1), ] # a single out-of-sample obsevation
      ret.gmv <- t(pi.equal) %*% t(out_sample_obs) # obtained portfolio returns
      portfolio_returns.gmv[day] = ret.gmv # store portfolio returns
      
      # update the weights of the portfolio
      pi.equal = pi.equal * t(1  + out_sample_obs/ 100)
      pi.equal = pi.equal / (sum(pi.equal)) # scale to one
      
      # adjust the out-of-sample-day
      day = day + 1
      
    }
  }
  # discard NA
  #portfolio_returns.gmv <- portfolio_returns.gmv[!is.na(portfolio_returns.gmv)]
  return(list(portfolio_returns.gmv, dcc.weights, weights_diff))
}
equal.result = equally_weighted_multistep(rebalance_horizon)
equal.port.ret = equal.result[[1]]
equal.port.weights = equal.result[[2]]
equal.weight.diff = equal.result[[3]]
equal.weight.diff.sum = rowSums(equal.weight.diff)

#----------------------------------------------------------
#----Forecast Combinations----------------------------
#----------------------------------------------------------

# calculate in-sample portfolio returns from dynamic methods
in_sample_portfolio_ret <- function(Cov){
  vOnes <- rep(1, num_assets)
  portfolio.ret = c()
  for (t in 1:T){
    Cov.mat <- Cov[,,t]
    #pi.gmv <- (solve(Cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(Cov.mat) %*% vOnes)[1])
    result.longonly <-  solnp(pars =  rep(1/27, 27), Sigma = Cov.mat,
                              fun = min_var,
                              #function to optimise
                              eqfun=equal, #equality function 
                              eqB=1,   #the equality constraint
                              LB=rep(0, 27), #lower bound for parameters i.e. greater than zero
                              UB=rep(4, 27), control = list(trace = 0, tol = 1e-6)) #upper bound for parameters (I just chose 100 randomly)
    pi.gmv = result.longonly$pars
    
    portfolio.ret[t] = t(pi.gmv) %*%  t(ret_sample[t,])
    
  }
  return(portfolio.ret)
}

# calculate in-sample portfolio returns from static sparse estimators
in_sample_portfolio_ret_static <- function(Cov){
  vOnes <- rep(1, num_assets)
  portfolio.ret = c()
  for (t in 1:T){
    Cov.mat <- Cov
    #pi.gmv <- (solve(Cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(Cov.mat) %*% vOnes)[1])
    result.longonly <-  solnp(pars =  rep(1/27, 27), Sigma = Cov.mat,
                              fun = min_var,
                              #function to optimise
                              eqfun=equal, #equality function 
                              eqB=1,   #the equality constraint
                              LB=rep(0, 27), #lower bound for parameters i.e. greater than zero
                              UB=rep(4, 27), control = list(trace = 0, tol = 1e-6)) #upper bound for parameters (I just chose 100 randomly)
    pi.gmv = result.longonly$pars
    portfolio.ret[t] = t(pi.gmv) %*%  t(ret_sample[t,])
    
  }
  return(portfolio.ret)
}
# compute the in-sample returns
in_sample_ret.dcc = in_sample_portfolio_ret(dcc.cov.insample)
in_sample_ret.deco = in_sample_portfolio_ret(dec_cov_insample)
in_sample_ret.bdeco = in_sample_portfolio_ret(bdec_cov_insample)
in_sample_ret.addaptive = in_sample_portfolio_ret_static(addaptive_cov[,,1])
in_sample_ret.sample = in_sample_portfolio_ret_static(sample_cov[,,1])

# create full sample returns
all_ret.dcc = c(in_sample_ret.dcc, dcc.port.ret)
all_ret.deco = c(in_sample_ret.deco, deco.port.ret)
all_ret.bdeco = c(in_sample_ret.bdeco, bdeco.port.ret)
all_ret.addaptive = c(in_sample_ret.addaptive, addaptive.port.ret)
all_ret.sample = c(in_sample_ret.sample, sample.port.ret)

lambda_forecast_combinations <- function(){
  # store the weights
  lambda.mat = matrix(ncol = 5, nrow = (num_forc))
  
  # store the forecast combination estimator
  FC.cov = array(dim = c(num_assets,num_assets,num_forc))
  for (t in 1:num_forc){
    # compute the volatility based on all observation untill time t which is known
    vol.dcc = sd(all_ret.dcc[1:(T+t-1)])
    vol.deco = sd(all_ret.deco[1:(T+t-1)])
    vol.bdeco = sd(all_ret.bdeco[1:(T+t-1)])
    vol.addaptive = sd(all_ret.addaptive[1:(T+t-1)])
    vol.sample = sd(all_ret.sample[1:(T+t-1)])
    
    # compute the weights for each method
    total_inverse_vol = 1/vol.dcc + 1/vol.deco + 1/vol.bdeco + 1/vol.addaptive + 1/vol.sample
    lambda.mat[t,1] = (1/vol.dcc) / total_inverse_vol
    lambda.mat[t,2] = (1/vol.deco) / total_inverse_vol
    lambda.mat[t,3] = (1/vol.bdeco) / total_inverse_vol
    lambda.mat[t,4] = (1/vol.addaptive) / total_inverse_vol
    lambda.mat[t,5] = (1/vol.sample) / total_inverse_vol
    
    # compute the forecast combinations covariance estimator
    FC.cov[,,t] = lambda.mat[t,1] * dcc_cov[,,t] + lambda.mat[t,2] * dec_cov[,,t] +
      lambda.mat[t,3] * bdec_cov[,,t] + lambda.mat[t,4] * addaptive_cov[,,t] + 
      lambda.mat[t,5] * sample_cov[,,t]
    
    
  }
  lambda.mat = data.frame(lambda.mat)
  names(lambda.mat) <- c("DCC", "DECO", "BDECO", "Adaptive", "sample")
  return(list(lambda.mat, FC.cov))
}
# get forecast combination results
forecast.combination.result = lambda_forecast_combinations()
lambda.FC = forecast.combination.result[[1]]
combined_cov = forecast.combination.result[[2]]

# portfolio results
combined.result = dcc_garch_gmv_multistep(rebalance_horizon, combined_cov)
combined.port.ret = combined.result[[1]]
combined.port.weights = combined.result[[2]]
combined.weight.diff = combined.result[[3]]
combined.weight.diff.sum = rowSums(combined.weight.diff)


#----------------------------------------------------------
#----Factor model-----------------------------------------
#----------------------------------------------------------

# store the factor model covariance matrices
factor_cov <- array(dim = c(num_assets,num_assets,num_forc))

for(t in 1:num_forc){
  factor_cov[,,t] = as.matrix(cov.factor)
}

# only applied if we use monthly rebalancing or 2 months ..  etc.
# portfolio results
factor.result = dcc_garch_gmv_multistep(rebalance_horizon, factor_cov)
factor.port.ret = factor.result[[1]]
factor.port.weights = factor.result[[2]]
factor.weight.diff = factor.result[[3]]
factor.weight.diff.sum = rowSums(factor.weight.diff)



#----------------------------------------------------------
#----Portfolio performance evaluation---------------------
#----------------------------------------------------------

# check correlation with the sp500
corr.dcc = cor(sp500.returns.sample[(T+1):size_sample], dcc.port.ret)
corr.deco = cor(sp500.returns.sample[(T+1):size_sample], deco.port.ret)
corr.bdeco = cor(sp500.returns.sample[(T+1):size_sample], bdeco.port.ret)
corr.addaptive = cor(sp500.returns.sample[(T+1):size_sample], addaptive.port.ret)
corr.sample = cor(sp500.returns.sample[(T+1):size_sample], sample.port.ret)
corr.equal = cor(sp500.returns.sample[(T+1):size_sample], equal.port.ret)
corr.combined = cor(sp500.returns.sample[(T+1):size_sample], combined.port.ret)
corr.factor = cor(sp500.returns.sample[(T+1):size_sample], factor.port.ret)

#C = c(corr.dcc, corr.deco, corr.bdeco, corr.addaptive, corr.sample, corr.combined, corr.equal)
C = c(corr.dcc, corr.deco, corr.bdeco, corr.addaptive, corr.sample, corr.combined, corr.equal, corr.factor)


# Annualized daily mean all methods 
mean.dcc = mean(dcc.port.ret) * 252
mean.deco = mean(deco.port.ret) * 252
mean.addaptive = mean(addaptive.port.ret) * 252
mean.sample = mean(sample.port.ret) * 252
mean.equal = mean(equal.port.ret) * 252
mean.combined = mean(combined.port.ret) * 252
mean.bdeco = mean(bdeco.port.ret) * 252
mean.factor = mean(factor.port.ret) * 252
#All_mean = c(mean.dcc, mean.deco, mean.bdeco, mean.addaptive , mean.sample, mean.combined, mean.equal)
All_mean = c(mean.dcc, mean.deco, mean.bdeco, mean.addaptive , mean.sample, mean.combined, mean.equal, mean.factor)

# Annualized daily turnover all methods 
TO.dcc = (1 / (num_forc)) * sum(dcc.weight.diff.sum) 
TO.deco =(1 / (num_forc)) * sum(deco.weight.diff.sum) 
TO.addaptive = (1 / (num_forc)) * sum(addaptive.weight.diff.sum) 
TO.sample = (1 / (num_forc)) * sum(sample.weight.diff.sum) 
TO.equal = (1 / (num_forc)) * sum(equal.weight.diff.sum) 
TO.bdeco = (1 / (num_forc)) * sum(bdeco.weight.diff.sum)
TO.combined = (1 / (num_forc)) * sum(combined.weight.diff.sum)
TO.factor = (1 / (num_forc)) * sum(factor.weight.diff.sum)
#All_TO = c(TO.dcc, TO.deco, TO.bdeco, TO.addaptive , TO.sample, TO.combined,  TO.equal) * 252
All_TO = c(TO.dcc, TO.deco, TO.bdeco, TO.addaptive , TO.sample, TO.combined,  TO.equal, TO.factor) * 252


# Annualized daily Sharpe-ratios for all methods 
SR.dcc = (mean(dcc.port.ret) / sd(dcc.port.ret)) * sqrt(252)
SR.deco = (mean(deco.port.ret) / sd(deco.port.ret)) * sqrt(252)
SR.addaptive = (mean(addaptive.port.ret) / sd(addaptive.port.ret)) * sqrt(252)
SR.sample = (mean(sample.port.ret) / sd(sample.port.ret)) * sqrt(252)
SR.equal = (mean(equal.port.ret) / sd(equal.port.ret)) * sqrt(252)
SR.bdeco = (mean(bdeco.port.ret) / sd(bdeco.port.ret)) * sqrt(252)
SR.combined = (mean(combined.port.ret) / sd(combined.port.ret)) * sqrt(252)
SR.factor = (mean(factor.port.ret) / sd(factor.port.ret)) * sqrt(252)
#All_SR = c(SR.dcc, SR.deco, SR.bdeco, SR.addaptive , SR.sample, SR.combined, SR.equal)
All_SR = c(SR.dcc, SR.deco, SR.bdeco, SR.addaptive , SR.sample, SR.combined, SR.equal, SR.factor)

# Annualized daily volatility for all methods
vol.dcc = sd(dcc.port.ret) * sqrt(252)
vol.deco = sd(deco.port.ret) * sqrt(252)
vol.addaptive = sd(addaptive.port.ret) * sqrt(252)
vol.sample = sd(sample.port.ret) * sqrt(252)
vol.equal = sd(equal.port.ret) * sqrt(252)
vol.bdeco = sd(bdeco.port.ret) * sqrt(252)
vol.combined = sd(combined.port.ret)* sqrt(252)
vol.factor = sd(factor.port.ret)* sqrt(252)
#All_vol = c(vol.dcc, vol.deco, vol.bdeco, vol.addaptive , vol.sample,vol.combined,  vol.equal )
All_vol = c(vol.dcc, vol.deco, vol.bdeco, vol.addaptive , vol.sample, vol.combined,  vol.equal, vol.factor)


# Daily Value at Risk for all methods based on the historical method
VaR.dcc = quantile(dcc.port.ret, 0.05)
VaR.deco = quantile(deco.port.ret, 0.05)
VaR.addaptive = quantile(addaptive.port.ret, 0.05)
VaR.sample = quantile(sample.port.ret, 0.05)
VaR.equal = quantile(equal.port.ret, 0.05)
VaR.bdeco = quantile(bdeco.port.ret, 0.05)
VaR.combined = quantile(combined.port.ret, 0.05)
VaR.factor = quantile(factor.port.ret, 0.05)
All_VaR = c(VaR.dcc, VaR.deco, VaR.addaptive , VaR.sample, VaR.equal, VaR.bdeco, VaR.combined)
#All_VaR = c(VaR.dcc, VaR.deco, VaR.addaptive , VaR.sample, VaR.equal, VaR.bdeco, VaR.combined, VaR.factor)

# Daily Expected shortfall for all methods based on the historical method
ES <- function(R, VaR){
  ES = R[which(R <= VaR)]
  return(mean(ES))
}
ES.dcc = ES(dcc.port.ret, VaR.dcc)
ES.deco = ES(deco.port.ret, VaR.deco)
ES.addaptive = ES(addaptive.port.ret, VaR.addaptive)
ES.sample = ES(sample.port.ret, VaR.sample)
ES.equal = ES(equal.port.ret, VaR.equal)
ES.bdeco = ES(bdeco.port.ret, VaR.bdeco)
ES.combined = ES(combined.port.ret, VaR.combined)
ES.factor = ES(factor.port.ret, VaR.factor)
#All_ES = c(ES.dcc, ES.deco,ES.bdeco,  ES.addaptive , ES.sample, ES.combined,ES.equal)
All_ES = c(ES.dcc, ES.deco, ES.addaptive , ES.sample, ES.equal, ES.bdeco, ES.combined, ES.factor)

# create dataframe with performance measures
#methods <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Sample Estimates", "Forecast Combinations", "1/N")
methods <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Sample Estimates", "Forecast Combinations", "1/N",
            "Factor model")
measures <- c("Annualized mean", "Annualized vol", "Annualized Sharpe-ratio", "Annualized turnover", "Daily ES 95%", "Correlation S&P 500")
performances = data.frame(All_mean, All_vol, All_SR, All_TO, All_ES, C)
names(performances) <- measures
rownames(performances) <- methods

# performances with a 1-day rebalancing frequency
round(performances, 2)

# create directly a latex table 
xtable(performances, 2)


# make a plot of the cummulative returns
# note that if we would do monthly or other types of rebalancing adjust the title
plot(out_sample_dates, cumsum(dcc.port.ret), type = "l", 
     main = "Cummulative out-of-sample portfolio returns daily rebalancing",
     col = "blue", xlab = "Time", ylab = "Cummulative return", ylim = c(-50,20))
lines(out_sample_dates, cumsum(deco.port.ret), col = "red")
lines(out_sample_dates, cumsum(addaptive.port.ret), col = "green")
lines(out_sample_dates, cumsum(sample.port.ret), col = "orange")
lines(out_sample_dates, cumsum(equal.port.ret), col = "yellow")
lines(out_sample_dates, cumsum(bdeco.port.ret), col = "black")
lines(out_sample_dates, cumsum(combined.port.ret), col = "pink")
lines(out_sample_dates, cumsum(index.returns[(T+1):size_sample]), col = "darkgrey")
legend("bottomleft", legend=c("Sample estimates", "DECO-GARCH", "DCC-GARCH", "Adaptive Thresholding", "1/N",
                              "BLOCK-DECO",
                              "Forecast combinations", "Index-Benchmark"),
       col=c("orange", "red", "blue", "green", "yellow", "black", "pink", "grey"), lty=1:1, cex=0.7)



#----------------------------------------------------------
#----Evaluation on all Rebalancing frequencies-------------
#----------------------------------------------------------

#### Create the plot of the sharpe=ratios and portfolio volatility given a rebalancing frequency
# frequencies = {1, .. 250}

# store sharpe-ratios  for each rebalancing horizon
dcc <- c()
deco <- c()
addaptive <- c()
sample <- c()
equal <- c()
combined <- c()
bdeco <- c()

# store portfolio volatilties for each rebalancing horizon
dcc2 <- c()
deco2 <- c()
addaptive2 <- c()
sample2 <- c()
equal2 <- c()
combined2 <- c()
bdeco2 <- c()
for (rebalance_freq in 1:250){
  # obtain results for all methods
  bdeco.res = dcc_garch_gmv_multistep(rebalance_freq, bdec_cov)
  dcc.res = dcc_garch_gmv_multistep(rebalance_freq, dcc_cov)
  deco.res = dcc_garch_gmv_multistep(rebalance_freq, dec_cov)
  addaptive.res= dcc_garch_gmv_multistep(rebalance_freq, addaptive_cov)
  sample.res = dcc_garch_gmv_multistep(rebalance_freq, sample_cov)
  equal.res = equally_weighted_multistep(rebalance_freq)
  combined.res = dcc_garch_gmv_multistep(rebalance_freq, combined_cov)
  
  # extract sharpe ratios for each method 
  dcc[rebalance_freq] = (mean(dcc.res[[1]]) / sd(dcc.res[[1]])) * sqrt(252)
  deco[rebalance_freq] = (mean(deco.res[[1]]) / sd(deco.res[[1]])) * sqrt(252)
  addaptive[rebalance_freq] = (mean(addaptive.res[[1]]) / sd(addaptive.res[[1]])) * sqrt(252)
  sample[rebalance_freq] = (mean(sample.res[[1]]) / sd(sample.res[[1]])) * sqrt(252)
  equal[rebalance_freq] = (mean(equal.res[[1]]) / sd(equal.res[[1]])) * sqrt(252)
  combined[rebalance_freq] = (mean(combined.res[[1]]) / sd(combined.res[[1]])) * sqrt(252)
  bdeco[rebalance_freq] = (mean(bdeco.res[[1]]) / sd(bdeco.res[[1]])) * sqrt(252)
  
  # calculate portfolio volatlity for each method
  dcc2[rebalance_freq] = (sd(dcc.res[[1]])) * sqrt(252)
  deco2[rebalance_freq] = (sd(deco.res[[1]])) * sqrt(252)
  addaptive2[rebalance_freq] = (sd(addaptive.res[[1]])) * sqrt(252)
  sample2[rebalance_freq] = (sd(sample.res[[1]])) * sqrt(252)
  equal2[rebalance_freq] = (sd(equal.res[[1]])) * sqrt(252)
  combined2[rebalance_freq] =  sd(combined.res[[1]]) * sqrt(252)
  bdeco2[rebalance_freq] = sd(bdeco.res[[1]]) * sqrt(252)
}
addaptive2 - dcc2
# create dataframe with sharpe-ratios for all rebalance frequencies
methods <- c("Sample estimates", "DECO-GARCH", "DCC-GARCH", "Adaptive Thresholding", "1/N",
             "BLOCK-DECO", "Forecast combinations")
sharpe_rebalance_freq = data.frame(sample, deco, dcc, addaptive, equal, bdeco, combined)
names(sharpe_rebalance_freq) <- methods
sharpe_rebalance_freq

# plot sharpe-ratios for all rebalancing frequencies for all methods
plot(c(1:250), sharpe_rebalance_freq$`DCC-GARCH`, type = "l", 
     main = "Annualized Portfolio Volatility given a rebalancing frequency",
     col = "blue", xlab = "Rebalancing frequency", ylab = "Annualized Portfolio Volatility", lwd = 2, ylim = c(-2, 2))
lines(c(1:250), sharpe_rebalance_freq$`DECO-GARCH`, col = "red", lwd = 2)
lines(c(1:250), sharpe_rebalance_freq$`Adaptive Thresholding`, col = "green", lwd = 2)
lines(c(1:250), sharpe_rebalance_freq$`Sample estimates`, col = "orange", lwd = 2)
lines(c(1:250), sharpe_rebalance_freq$`1/N`, col = "yellow", lwd = 2)
lines(c(1:250), sharpe_rebalance_freq$`BLOCK-DECO`, col = "black", lwd = 2)
lines(c(1:250), sharpe_rebalance_freq$`Forecast combinations`, col = "pink", lwd = 2)
legend("bottomleft", legend= methods,
       col=c("orange", "red", "blue", "green", "yellow", "black", "pink"), lty=1:1, cex=0.7)


# 
# #----------------------------------------------------------
# #----Factor Model----------------------------
# #----------------------------------------------------------
# library(PerformanceAnalytics)
# library(readxl)
# library(xts)
# library(ggplot2)
# options(scipen=999)
# factor.data = read_excel("factor_data.xlsx")
# 
# ### other factors
# 
# # VIX
# VIX.data <- read.csv("^VIX.csv")
# VIX.dates = as.Date(VIX.data$Date[-1])
# VIX.returns = diff(log(VIX.data$Adj.Close))
# VIX.ret.cleaned <- xts(x=VIX.returns, order.by=VIX.dates)
# Monthly.VIX = apply.monthly(VIX.ret.cleaned , (Return.cumulative))
# Monthly.VIX = Monthly.VIX * 100
# Monthly.VIX[1]
# 
# ### treasury yield
# 
# # # read in data
# # Tbill.data <- read.csv("^IRX.csv")
# # 
# # # drop null
# # Tbill.data<-Tbill.data[!(Tbill.data$Adj.Close=="null"),]
# # Tbill.data$Adj.Close = as.numeric(Tbill.data$Adj.Close)
# # 
# # # get dates
# # Tbill.dates = as.Date(Tbill.data$Date[-1])
# # 
# # # continue
# # Tbill.returns = Tbill.data$Adj.Close[-1] / 100
# # Tbill.ret.cleaned <- xts(x=Tbill.returns, order.by=Tbill.dates)
# # Monthly.Tbill = apply.monthly(Tbill.ret.cleaned , (Return.cumulative))
# # Monthly.Tbill = Monthly.Tbill
# # Monthly.Tbill[248]
# 
# ### treasury yield
# 
# # read in data
# Tbill.data <- read.csv("^GSPC (4).csv")
# Tbill.data
# # drop null
# #Tbill.data<-Tbill.data[!(Tbill.data$Adj.Close=="null"),]
# #Tbill.data$Adj.Close = as.numeric(Tbill.data$Adj.Close)
# 
# # get dates
# Tbill.dates = as.Date(Tbill.data$Date[-1])
# 
# # continue
# Tbill.returns = diff(log(Tbill.data$Adj.Close))
# Tbill.ret.cleaned <- xts(x=Tbill.returns, order.by=Tbill.dates)
# Monthly.Tbill = apply.monthly(Tbill.ret.cleaned , (Return.cumulative))
# Monthly.Tbill = Monthly.Tbill * 100
# Monthly.Tbill
# 
# 
# ### Bond-10y
# 
# # read in data
# Bond.data <- read.csv("^TNX.csv")
# Bond.data
# 
# # drop null
# Bond.data<-Bond.data[!(Bond.data$Adj.Close=="null"),]
# Bond.data$Adj.Close = as.numeric(Bond.data$Adj.Close)
# 
# # get dates
# Bond.dates = as.Date(Bond.data$Date[-1])
# 
# # continue
# Bond.returns = diff(log(Bond.data$Adj.Close))
# Bond.ret.cleaned <- xts(x=Bond.returns, order.by=Bond.dates)
# Monthly.Bond = apply.monthly(Bond.ret.cleaned , (Return.cumulative))
# Monthly.Bond = Monthly.Bond * 100
# Monthly.Bond
# 
# ###  USD
# 
# # read in data
# USD.data <- read.csv("DX-Y.NYB.csv")
# USD.data
# 
# # drop null
# USD.data<-USD.data[!(USD.data$Adj.Close=="null"),]
# USD.data$Adj.Close = as.numeric(USD.data$Adj.Close)
# 
# # get dates
# USD.dates = as.Date(USD.data$Date[-1])
# 
# # continue
# USD.returns = diff(log(USD.data$Adj.Close))
# USD.ret.cleaned <- xts(x=USD.returns, order.by=USD.dates)
# Monthly.USD = apply.monthly(USD.ret.cleaned , (Return.cumulative))
# Monthly.USD= Monthly.USD * 100
# Monthly.USD
# 
# 
# # regression
# r2 = c()
# for (i in 1:27){
#   # obtain regression inputs
#   test = ret[,i] / 100
#   MonthlyReturns <- apply.monthly(test, (Return.cumulative))
#   size.factors = length(factor.data$Time)
#   size.monthly.ret = length(MonthlyReturns)
#   MonthlyReturns = MonthlyReturns * 100
#   
#   # define regression 
#   y = MonthlyReturns[1:(size.monthly.ret-2)] 
#   x1 = factor.data$`INFLATION (%)`[-1] # inflation
#   X2 = Monthly.VIX[1:(size.monthly.ret-2)] # VIX
#   x3 = Monthly.Tbill[1:(size.monthly.ret-2)] # SP500
#   x4 = Monthly.Bond[1:(size.monthly.ret-2)] # bond returns
#   #x5 = factor.data$EXCHANGE_RATE[-1]
#   x5 = Monthly.USD[1:(size.monthly.ret-2)] # USD currency
#   
#   # fama-french 
#   SMB = 1/3 * (MonthlyReturns)
#   
#   # estimate model
#   model = lm(y ~ x1 + x3 + x4 + x5)
#   res = summary(model)
#   print(res)
#   r2[i] = res$r.squared
# }
# ### Plot result
# 
# # Create data
# data <- data.frame(
#   name=names(ret) ,  
#   value=r2
# )
# # Barplot
# p <- ggplot(data, aes(x=name, y=value), col = "blue") + 
#   geom_bar(stat = "identity")
# p + ggtitle("R-squares factor model of all commodities")+
#   xlab("Commodity") + ylab("R2")
# 
# 
# 
# #----------------------------------------------------------
# #----Dynamic Rebalance strategy----------------------------
# #----------------------------------------------------------
# 
# 
# #### Dynamic rebalancing
# 








