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
library(remotes)
library(auxPort)

# ommit scientific notation
options(scipen=999)

#----------------------------------------------------------
#----Data----------------------------------------
#----------------------------------------------------------

# read in commodity data
All_data <- readRDS(file = "2021_01_comm_data_corrected.RDS")

# read in sp500 data
sp500.data <- read.csv("^GSPC (13).csv")
ret_sp500 = c(0, diff(sp500.data$Adj.Close)/sp500.data$Adj.Close[-length(sp500.data$Adj.Close)] * 100)
sp500.data$returns = ret_sp500

# daily return datar
ret_data <- All_data$returns
num_assets <- 27
ret <- ret_data[,1:num_assets] * 100

# obtain the in-sample length
start.date.in = as.Date("2000-01-04")
end.date.in = as.Date("2014-12-31")
in.sample_period = ret[which(rownames(data.frame(ret)) == start.date.in):
                         which(rownames(data.frame(ret)) == end.date.in),]
in_len  = length(in.sample_period$W)

# subset the data to the full sample w.r.t. the panel you have
start.date.full = as.Date("2000-01-04")
end.date.full = as.Date("2020-12-31")

# define full sample
ret_sample = ret[which(rownames(data.frame(ret)) == start.date.full):
                   which(rownames(data.frame(ret)) == end.date.full),]

# the full sample length
size_sample <- length(ret_sample$W)

# get index returns
index_weights_sample = All_data$weights[which(rownames(data.frame(All_data$weights)) == start.date.full):
                                          which(rownames(data.frame(All_data$weights)) == end.date.full),]
index.returns <- c()
for (t in 1:length(ret_sample$W)){
  index.returns[t] = (index_weights_sample[t,]) %*% t(ret_sample[t,])
}


# read in data of factor covariance matrix
cov.factor = read_excel("Cov_period2.xlsx", col_names = FALSE)
cov.factor = data.frame(cov.factor)
names(cov.factor) <- colnames(ret_sample)

#----------------------------------------------------------
#----Calculate Portfolio returns----------------------------------------
#----------------------------------------------------------

# calculates portfolio returns out of sample
portfolio_returns_out_sample <- function(horizon, forc.cov.mat){
  # set up daily transaction costs of 8.6 bps
  TC = 0.01 * 8.6 
  
  # Store portfolio returns
  portfolio_returns.gmv <- rep(0, num_forc)
  
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
    # get the covariance matrix for the out-of-sample observation t
    Cov.mat <- forc.cov.mat[,,t]
    
    # calculate the weights
    pi.gmv.start <- (solve(Cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(Cov.mat) %*% vOnes)[1])
    
    # if are rebalancing then add transaction costs
    if (t > 1){
      # store the absolute weight difference
      weights_diff[rebalance_day,] = t(abs(pi.gmv - pi.gmv.start ))
      rebalance_day = rebalance_day + 1
      
      # subtract transaction costs
      portfolio_returns.gmv[(day)] = portfolio_returns.gmv[(day)] - TC * sum(t(abs(pi.gmv - pi.gmv.start )))
    }
    
    # at day 0, first allocate the portfolio and our weights are all 0
    else{
      weights_diff[(rebalance_day),] = rep(0, 27)
      rebalance_day = rebalance_day + 1
    }
    # start with de optimized portolio 
    pi.gmv = pi.gmv.start
    
    # check if we do not already reached to full out-of-sample period
    for (i in 1:horizon){
      # check if do not are out of bound
      if ((T + t + i - 1) > size_sample){
        break
      }
      
      # store the weights
      dcc.weights[day,] = pi.gmv
      
      # calculate obtained portfolio return
      out_sample_obs <- ret_sample[(T + t + i - 1), ]
      ret.gmv <- t(pi.gmv) %*% t(out_sample_obs)
      portfolio_returns.gmv[day] = portfolio_returns.gmv[day] + ret.gmv
      
      # update the weights of the portfolio
      # they change throughout time
      pi.gmv = pi.gmv * t(1  + out_sample_obs/ 100)
      pi.gmv = pi.gmv / (sum(pi.gmv)) # scale to one
      
      # adjust the out-of-sample-day
      day = day + 1
      
    }
  }
  return(list(portfolio_returns.gmv, dcc.weights, weights_diff ))
}

# calculates portfolio returns out of sample
portfolio_returns_in_sample <- function(horizon, forc.cov.mat, n_size, Return_series){
  # set up daily transaction costs of 8.6 bps
  TC = 0.01 * 8.6 
  
  # Store portfolio returns
  portfolio_returns.gmv <- rep(0, n_size)
  
  # Initialize a vector of ones
  vOnes <- rep(1, num_assets)
  
  # keep track of the out-of-sample day
  day = 1
  
  for (t in seq(from = 1, to = (n_size), by=horizon)){
    # get the covariance matrix for the out-of-sample observation t
    Cov.mat <- forc.cov.mat[,,t]
    
    # calculate the weights
    pi.gmv.start <- (solve(Cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(Cov.mat) %*% vOnes)[1])
    
    # if are rebalancing then add transaction costs
    if (t > 1){
      # subtract transaction costs
      portfolio_returns.gmv[(day)] = portfolio_returns.gmv[(day)] - TC * sum(t(abs(pi.gmv - pi.gmv.start )))
    }
    # start with de optimized portolio 
    pi.gmv = pi.gmv.start
    
    # check if we do not already reached to full out-of-sample period
    for (i in 1:horizon){
      # check if do not are out of bound
      if ((t + i) > n_size){
        break
      }
      
      # calculate obtained portfolio return
      out_sample_obs <- Return_series[(t + i), ]
      ret.gmv <- t(pi.gmv) %*% t(out_sample_obs)
      portfolio_returns.gmv[day] = portfolio_returns.gmv[day] + ret.gmv
      
      # update the weights of the portfolio
      # they change throughout time
      pi.gmv = pi.gmv * t(1  + out_sample_obs/ 100)
      pi.gmv = pi.gmv / (sum(pi.gmv)) # scale to one
      
      # adjust the out-of-sample-day
      day = day + 1
      
    }
  }
  return(portfolio_returns.gmv)
}

#----------------------------------------------------------
#----DCC-GARCH----------------------------------------
#----------------------------------------------------------

# in-sample and number of out-of-sample observations
in_sample_size <- in_len
num_forc <- size_sample - in_sample_size
in_len

# define the out-of-sample return series
out_of_sample = ret_sample[(in_sample_size + 1):size_sample, ]

# out-of-sample dates
temp_dates = ret_sample[(in_sample_size):size_sample, ]
out_sample_dates = as.Date(rownames(data.frame(temp_dates)))

# subset sp500 data
sp500.data.sample <- sp500.data[which(sp500.data$Date %in% rownames(data.frame(out_of_sample))),]

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
#dcc.fit@mfit$convergence

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


# In-sample covariances (BLOCK) DECO

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
        gmv_weights = (solve(cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(cov.mat) %*% vOnes)[1])
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

# store Addaptive tresholding forecasted covariance
addaptive_cov <- array(dim = c(num_assets,num_assets,num_forc)) 

# store Addaptive tresholding forecasted covariance
sample_cov <- array(dim = c(num_assets,num_assets,num_forc))

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

# get insample covariance matrices
addaptive_cov_insample <- array(dim = c(num_assets,num_assets,T)) 
sample_cov_insample <- array(dim = c(num_assets,num_assets, T ))
for (i in 1:T){
  # store the estimated covariance matrices without including new data
  addaptive_cov_insample[,,i] <- est_cov.mat.addaptive
  sample_cov_insample [,,i] <- est_cov.mat.sample
  
}




#----------------------------------------------------------
#----1/N---------------------------------------------------
#----------------------------------------------------------
equally_weighted_portfolio_returns <- function(horizon){
  # set up daily transaction costs of 8.6 bps
  TC = 0.01 * 8.6 
  
  # Store portfolio returns
  portfolio_returns.gmv <- rep(0, num_forc)
  
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
      
      # substract transaction costs
      portfolio_returns.gmv[(day)] = portfolio_returns.gmv[(day)] - TC * sum(t(abs(pi.equal - w0 )))
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
      portfolio_returns.gmv[day] = portfolio_returns.gmv[day] +  ret.gmv # store portfolio returns
      
      # update the weights of the portfolio
      pi.equal = pi.equal * t(1  + out_sample_obs/ 100)
      pi.equal = pi.equal / (sum(pi.equal)) # scale to one
      
      # adjust the out-of-sample-day
      day = day + 1
      
    }
  }
  return(list(portfolio_returns.gmv, dcc.weights, weights_diff))
}



#----------------------------------------------------------
#----Forecast Combinations----------------------------
#----------------------------------------------------------
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



#----------------------------------------------------------
#----Factor model-----------------------------------------
#----------------------------------------------------------

# store the factor model covariance matrices
factor_cov <- array(dim = c(num_assets,num_assets,num_forc))
for(t in 1:num_forc){
  factor_cov[,,t] = as.matrix(cov.factor)
}



#----------------------------------------------------------
#----Calculate Cummulative returns-----------------------------------------
#----------------------------------------------------------
cum_ret <- function(portfolio_returns){
  cummulative_returns = c()
  cummulative_returns[1] = 1
  for (t in 2:(length(portfolio_returns)+1)){
    cummulative_returns[t] = (1 + (portfolio_returns[(t-1)] / 100))
    
  }
  final_ret = cumprod(cummulative_returns) - 1
  final_ret = final_ret * 100
  return(final_ret)
}


#----------------------------------------------------------
#----Plotting average weights----------------------------------------
#----------------------------------------------------------

barplot_weights <- function(All.weights, title){
  data <- data.frame(
    name= colnames(All_data$returns),  
    value= colMeans(All.weights)
  )
  
  # Barplot
  ggplot(data, aes(x=name, y=value)) + 
    geom_bar(stat = "identity") + scale_x_discrete(name ="Commodity") + scale_y_continuous(name ="Average weight") +ggtitle(title)
  
}

minmax_weights <- function(){
  weights.min.max = matrix(ncol = 4, nrow = 7)
  weights.min.max[1,] = c(min(colMeans(dcc.port.weights)), colnames(All_data$returns)[which.min(colMeans(dcc.port.weights))],
                          max(colMeans(dcc.port.weights)), colnames(All_data$returns)[which.max(colMeans(dcc.port.weights))] )
  weights.min.max[2,] = c(min(colMeans(deco.port.weights)), colnames(All_data$returns)[which.min(colMeans(deco.port.weights))],
                          max(colMeans(deco.port.weights)), colnames(All_data$returns)[which.max(colMeans(deco.port.weights))] )
  weights.min.max[3,] = c(min(colMeans(bdeco.port.weights)), colnames(All_data$returns)[which.min(colMeans(bdeco.port.weights))],
                          max(colMeans(bdeco.port.weights)), colnames(All_data$returns)[which.max(colMeans(bdeco.port.weights))] )
  weights.min.max[4,] = c(min(colMeans(addaptive.port.weights)), colnames(All_data$returns)[which.min(colMeans(addaptive.port.weights))],
                          max(colMeans(addaptive.port.weights)), colnames(All_data$returns)[which.max(colMeans(addaptive.port.weights))] )
  weights.min.max[5,] = c(min(colMeans(sample.port.weights)), colnames(All_data$returns)[which.min(colMeans(sample.port.weights))],
                          max(colMeans(sample.port.weights)), colnames(All_data$returns)[which.max(colMeans(sample.port.weights))] )
  weights.min.max[6,] = c(min(colMeans(combined.port.weights)), colnames(All_data$returns)[which.min(colMeans(combined.port.weights))],
                          max(colMeans(combined.port.weights)), colnames(All_data$returns)[which.max(colMeans(combined.port.weights))] )
  
  weights.min.max[7,] = c(min(colMeans(factor.port.weights)), colnames(All_data$returns)[which.min(colMeans(factor.port.weights))],
                          max(colMeans(factor.port.weights)), colnames(All_data$returns)[which.max(colMeans(factor.port.weights))] )

  weights.min.max = data.frame(weights.min.max)
  colnames(weights.min.max) = c("Min weight", "Name of min weight", "Max weight", "Name of max weight")
  row.names(weights.min.max) = c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Sample Estimates",
                                 "Forecast Combinations", "Factor model") #"Factor model"
  print(weights.min.max)
}
max_secondmax_weights <- function(){
  weights.max = matrix(ncol = 2, nrow = 7)
  weights.max[1,] = c(sort(colMeans(dcc.port.weights), TRUE)[2],
                          max(colMeans(dcc.port.weights)) )
  weights.max[2,] = c(sort(colMeans(deco.port.weights), TRUE)[2],
                          max(colMeans(deco.port.weights)) )
  weights.max[3,] = c(sort(colMeans(bdeco.port.weights), TRUE)[2],
                          max(colMeans(bdeco.port.weights)) )
  weights.max[4,] = c(sort(colMeans(addaptive.port.weights), TRUE)[2],
                          max(colMeans(addaptive.port.weights)) )
  weights.max[5,] = c(sort(colMeans(sample.port.weights), TRUE)[2],
                          max(colMeans(sample.port.weights)) )
  weights.max[6,] = c(sort(colMeans(combined.port.weights), TRUE)[2],
                          max(colMeans(combined.port.weights)) )
  
  weights.max[7,] = c(sort(colMeans(factor.port.weights), TRUE)[2],
                          max(colMeans(factor.port.weights)) )
  
  weights.max = data.frame(weights.max)
  colnames(weights.max) = c("Second max weight", "Max weight")
  row.names(weights.max) = c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Sample Estimates",
                                 "Forecast Combinations", "Factor model")
  print(weights.max)
}

 

#----------------------------------------------------------
#----Main-----------------------------------------
#----------------------------------------------------------
rebalance_horizon = 252


### DCC result
dcc.result = portfolio_returns_out_sample(rebalance_horizon, dcc_cov)
dcc.port.ret = dcc.result[[1]]
dcc.port.weights = dcc.result[[2]]
dcc.weight.diff = dcc.result[[3]]
dcc.weight.diff.sum = rowSums(dcc.weight.diff)

### DECO result
deco.result = portfolio_returns_out_sample(rebalance_horizon, dec_cov)
deco.port.ret = deco.result[[1]]
deco.port.weights = deco.result[[2]]
deco.weight.diff = deco.result[[3]]
deco.weight.diff.sum = rowSums(deco.weight.diff)

### BLOCK-DECO result
bdeco.result = portfolio_returns_out_sample(rebalance_horizon, bdec_cov)
bdeco.port.ret = bdeco.result[[1]]
bdeco.port.weights = bdeco.result[[2]]
bdeco.weight.diff = bdeco.result[[3]]
bdeco.weight.diff.sum = rowSums(bdeco.weight.diff)

### Adaptive Tresholding result
addaptive.result = portfolio_returns_out_sample(rebalance_horizon, addaptive_cov)
addaptive.port.ret = addaptive.result[[1]]
addaptive.port.weights = addaptive.result[[2]]
addaptive.weight.diff = addaptive.result[[3]]
addaptive.weight.diff.sum = rowSums(addaptive.weight.diff)

### Sample estimates result
sample.result = portfolio_returns_out_sample(rebalance_horizon, sample_cov)
sample.port.ret = sample.result[[1]]
sample.port.weights = sample.result[[2]]
sample.weight.diff = sample.result[[3]]
sample.weight.diff.sum = rowSums(sample.weight.diff)

### Equally weighted result
equal.result = equally_weighted_portfolio_returns(rebalance_horizon)
equal.port.ret = equal.result[[1]]
equal.port.weights = equal.result[[2]]
equal.weight.diff = equal.result[[3]]
equal.weight.diff.sum = rowSums(equal.weight.diff)

### Forecast combination results
# in-sample returns given a rebalancing frequency
in_sample_ret.dcc = portfolio_returns_in_sample(rebalance_horizon,  dcc.cov.insample, T, ret_sample[1:T,])
in_sample_ret.deco = portfolio_returns_in_sample(rebalance_horizon,  dec_cov_insample, T, ret_sample[1:T,])
in_sample_ret.bdeco = portfolio_returns_in_sample(rebalance_horizon,  bdec_cov_insample, T, ret_sample[1:T,])
in_sample_ret.addaptive = portfolio_returns_in_sample(rebalance_horizon,  addaptive_cov_insample, T, ret_sample[1:T,])
in_sample_ret.sample = portfolio_returns_in_sample(rebalance_horizon,  sample_cov_insample, T, ret_sample[1:T,])
# stack returns
all_ret.dcc = c(in_sample_ret.dcc, dcc.port.ret)
all_ret.deco = c(in_sample_ret.deco, deco.port.ret)
all_ret.bdeco = c(in_sample_ret.bdeco, bdeco.port.ret)
all_ret.addaptive = c(in_sample_ret.addaptive, addaptive.port.ret)
all_ret.sample = c(in_sample_ret.sample, sample.port.ret)
# get forecast combination results
forecast.combination.result = lambda_forecast_combinations()
lambda.FC = forecast.combination.result[[1]]
combined_cov = forecast.combination.result[[2]]
# portfolio results
combined.result = portfolio_returns_out_sample(rebalance_horizon, combined_cov)
combined.port.ret = combined.result[[1]]
combined.port.weights = combined.result[[2]]
combined.weight.diff = combined.result[[3]]
combined.weight.diff.sum = rowSums(combined.weight.diff)


### Factor model
factor.result = portfolio_returns_out_sample(rebalance_horizon, factor_cov)
factor.port.ret = factor.result[[1]]
factor.port.weights = factor.result[[2]]
factor.weight.diff = factor.result[[3]]
factor.weight.diff.sum = rowSums(factor.weight.diff)


#----------------------------------------------------------
#----Final result-----------------------------------------
#----------------------------------------------------------
performances_measures_excl_factor()
performances_measures_incl_factor()



#----------------------------------------------------------
#----Bar plots of average weights-----------------------------------------
#----------------------------------------------------------
barplot_weights(dcc.port.weights, "DCC unrestricted GMV")
barplot_weights(deco.port.weights, "DECO unrestricted  GMV")
barplot_weights(bdeco.port.weights, "BLOCK-DECO unrestricted GMV")
barplot_weights(addaptive.port.weights, "Adaptive Tresholding unrestricted  GMV")
barplot_weights(sample.port.weights, "Sample covariance unrestricted  GMV")
barplot_weights(combined.port.weights, "Forecast Combinations unrestricted  GMV")
barplot_weights(factor.port.weights, "Factor model unrestricted  GMV")

#----------------------------------------------------------
#----Min/max weights----------------------------------------
#----------------------------------------------------------
minmax_weights()
max_secondmax_weights()


#----------------------------------------------------------
#----Forecast combination average weights-----------------------------------------
#----------------------------------------------------------
colMeans(lambda.FC)


#----------------------------------------------------------
#----Calculate performance given a rebalancing frequency at the end------
#----------------------------------------------------------
reb.frequency.res = rebalancing_frequencies_performance()
vol.reb.freq = reb.frequency.res[[1]]
SR.reb.freq =  reb.frequency.res[[2]]

# Plot sharpe ratios and volatilty against rebalancing frequencies
plot_vol(vol.reb.freq)
plot_sharpe(SR.reb.freq)

#----------------------------------------------------------
#----Result calculation functions-----------------------------------------
#----------------------------------------------------------

### performance evaluation excluding the factor model
performances_measures_excl_factor <- function(){
  # get correlation between portfolio returns and the SP500
  corr.dcc =  cor(dcc.port.ret, sp500.data.sample$returns)
  corr.deco = cor(deco.port.ret, sp500.data.sample$returns)
  corr.bdeco = cor(bdeco.port.ret, sp500.data.sample$returns)
  corr.addaptive = cor(addaptive.port.ret, sp500.data.sample$returns)
  corr.sample = cor(sample.port.ret, sp500.data.sample$returns)
  corr.equal =  cor(equal.port.ret, sp500.data.sample$returns)
  corr.combined = cor(combined.port.ret, sp500.data.sample$returns)
  
  # Annualized daily mean 
  mean.dcc = mean(dcc.port.ret) * 252
  mean.deco = mean(deco.port.ret) * 252
  mean.addaptive = mean(addaptive.port.ret) * 252
  mean.sample = mean(sample.port.ret) * 252
  mean.equal = mean(equal.port.ret) * 252
  mean.combined = mean(combined.port.ret) * 252
  mean.bdeco = mean(bdeco.port.ret) * 252
  
  # Annualized daily turnover 
  TO.dcc = (1 / (num_forc)) * sum(dcc.weight.diff.sum) 
  TO.deco =(1 / (num_forc)) * sum(deco.weight.diff.sum) 
  TO.addaptive = (1 / (num_forc)) * sum(addaptive.weight.diff.sum) 
  TO.sample = (1 / (num_forc)) * sum(sample.weight.diff.sum) 
  TO.equal = (1 / (num_forc)) * sum(equal.weight.diff.sum) 
  TO.bdeco = (1 / (num_forc)) * sum(bdeco.weight.diff.sum)
  TO.combined = (1 / (num_forc)) * sum(combined.weight.diff.sum)
  
  # Annualized daily Sharpe-ratios 
  SR.dcc = (mean(dcc.port.ret) / sd(dcc.port.ret)) * sqrt(252)
  SR.deco = (mean(deco.port.ret) / sd(deco.port.ret)) * sqrt(252)
  SR.addaptive = (mean(addaptive.port.ret) / sd(addaptive.port.ret)) * sqrt(252)
  SR.sample = (mean(sample.port.ret) / sd(sample.port.ret)) * sqrt(252)
  SR.equal = (mean(equal.port.ret) / sd(equal.port.ret)) * sqrt(252)
  SR.bdeco = (mean(bdeco.port.ret) / sd(bdeco.port.ret)) * sqrt(252)
  SR.combined = (mean(combined.port.ret) / sd(combined.port.ret)) * sqrt(252)
  
  # Annualized daily volatility for all methods
  vol.dcc = sd(dcc.port.ret) * sqrt(252)
  vol.deco = sd(deco.port.ret) * sqrt(252)
  vol.addaptive = sd(addaptive.port.ret) * sqrt(252)
  vol.sample = sd(sample.port.ret) * sqrt(252)
  vol.equal = sd(equal.port.ret) * sqrt(252)
  vol.bdeco = sd(bdeco.port.ret) * sqrt(252)
  vol.combined = sd(combined.port.ret)* sqrt(252)
  
  # Value at Risk 
  VaR.dcc = quantile(dcc.port.ret, 0.05)
  VaR.deco = quantile(deco.port.ret, 0.05)
  VaR.addaptive = quantile(addaptive.port.ret, 0.05)
  VaR.sample = quantile(sample.port.ret, 0.05)
  VaR.equal = quantile(equal.port.ret, 0.05)
  VaR.bdeco = quantile(bdeco.port.ret, 0.05)
  VaR.combined = quantile(combined.port.ret, 0.05)
  
  # Expected shortfall 
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
  
  # Stack result
  All_mean = c(mean.dcc, mean.deco, mean.bdeco, mean.addaptive , mean.sample, mean.combined, mean.equal)
  All_TO = c(TO.dcc, TO.deco, TO.bdeco, TO.addaptive , TO.sample, TO.combined,  TO.equal) * 252
  All_SR = c(SR.dcc, SR.deco, SR.bdeco, SR.addaptive , SR.sample, SR.combined, SR.equal)
  All_vol = c(vol.dcc, vol.deco, vol.bdeco, vol.addaptive , vol.sample, vol.combined,  vol.equal)
  C = c(corr.dcc, corr.deco, corr.bdeco, corr.addaptive, corr.sample, corr.combined, corr.equal)
  All_ES = c(ES.dcc, ES.deco,ES.bdeco,  ES.addaptive , ES.sample, ES.combined, ES.equal)
  
  # create dataframe with performance measures
  methods <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Sample Estimates", "Forecast Combinations", "1/N")
  measures <- c("Annualized mean", "Annualized vol", "Annualized Sharpe-ratio", "Annualized turnover", "Daily ES 95%", "Correlation S&P 500")
  performances = data.frame(All_mean, All_vol, All_SR, All_TO, abs(All_ES), C)
  names(performances) <- measures
  rownames(performances) <- methods
  
  # create directly a latex table 
  print(xtable(performances, 2))
  
  ### Statistical tests
  # Sharpe-ratio test
  dcc.sr.test = hac_infer(cbind(dcc.port.ret, sample.port.ret), digits = 4, type = "SR")
  deco.sr.test = hac_infer(cbind(deco.port.ret, sample.port.ret), digits = 4, type = "SR")
  bdeco.sr.test = hac_infer(cbind(bdeco.port.ret, sample.port.ret), digits = 4, type = "SR")
  addaptive.sr.test = hac_infer(cbind(addaptive.port.ret, sample.port.ret), digits = 4, type = "SR")
  combined.sr.test = hac_infer(cbind(combined.port.ret, sample.port.ret), digits = 4, type = "SR")
  equal.sr.test = hac_infer(cbind(equal.port.ret, sample.port.ret), digits = 4, type = "SR")

  # Volatility test
  dcc.vol.test = hac_infer(cbind(dcc.port.ret, sample.port.ret), digits = 4, type = "Var")
  deco.vol.test = hac_infer(cbind(deco.port.ret, sample.port.ret), digits = 4, type = "Var")
  bdeco.vol.test = hac_infer(cbind(bdeco.port.ret, sample.port.ret), digits = 4, type = "Var")
  addaptive.vol.test = hac_infer(cbind(addaptive.port.ret, sample.port.ret), digits = 4, type = "Var")
  combined.vol.test = hac_infer(cbind(combined.port.ret, sample.port.ret), digits = 4, type = "Var")
  equal.vol.test = hac_infer(cbind(equal.port.ret, sample.port.ret), digits = 4, type = "Var")
  
  SR_test = matrix(ncol = 6, nrow = 3)
  Vol_test = matrix(ncol = 6, nrow = 3)
  SR_test[1,] = c(dcc.sr.test$Difference, deco.sr.test$Difference, bdeco.sr.test$Difference, addaptive.sr.test$Difference
                  , combined.sr.test$Difference, equal.sr.test$Difference)
  SR_test[2,] = c(dcc.sr.test$Standard.Errors[1], deco.sr.test$Standard.Errors[1],
                  bdeco.sr.test$Standard.Errors[1], addaptive.sr.test$Standard.Errors[1]
    , combined.sr.test$Standard.Errors[1], equal.sr.test$Standard.Errors[1])
  SR_test[3,] = c(dcc.sr.test$p.Values[1], deco.sr.test$p.Values[1],
                  bdeco.sr.test$p.Values[1], addaptive.sr.test$p.Values[1]
                  , combined.sr.test$p.Values[1], equal.sr.test$p.Values[1])

  Vol_test[1,] = c(dcc.vol.test$Difference, deco.vol.test$Difference, bdeco.vol.test$Difference, addaptive.vol.test$Difference
                  , combined.vol.test$Difference, equal.vol.test$Difference)
  Vol_test[2,] = c(dcc.vol.test$Standard.Errors[1], deco.vol.test$Standard.Errors[1],
                  bdeco.vol.test$Standard.Errors[1], addaptive.vol.test$Standard.Errors[1]
                  , combined.vol.test$Standard.Errors[1], equal.vol.test$Standard.Errors[1])
  Vol_test[3,] = c(dcc.vol.test$p.Values[1], deco.vol.test$p.Values[1],
                  bdeco.vol.test$p.Values[1], addaptive.vol.test$p.Values[1]
                  , combined.vol.test$p.Values[1], equal.vol.test$p.Values[1])
  SR_test = data.frame(SR_test)
  row.names(SR_test) <- c("Delta.diff", "HAC se", "pvalue")
  names(SR_test) <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Forecast Combinations", "1/N")
  
  Vol_test = data.frame(Vol_test)
  row.names(Vol_test) <- c("Delta.diff", "HAC se", "pvalue")
  names(Vol_test) <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Forecast Combinations", "1/N")
  
  print(xtable(SR_test, 2))
  print(xtable(Vol_test, 2))

  # make a plot of the cummulative returns
  cum_ret.dcc = cum_ret(dcc.port.ret)
  cum_ret.deco = cum_ret(deco.port.ret)
  cum_ret.bdeco = cum_ret(bdeco.port.ret)
  cum_ret.sample = cum_ret(sample.port.ret)
  cum_ret.addaptive = cum_ret(addaptive.port.ret)
  cum_ret.combined = cum_ret(combined.port.ret)
  cum_ret.equal = cum_ret(equal.port.ret)
  cum_ret.index = cum_ret(index.returns[(T+1):size_sample])

  # note that if we would do monthly or other types of rebalancing adjust the title
  plot(out_sample_dates, cum_ret.dcc, type = "l",
       col = "blue", xlab = "Time", ylab = "Cumulative return (%)", ylim = c(-60,35), lwd = 2)
  lines(out_sample_dates, cum_ret.deco, col = "red", lwd = 2)
  lines(out_sample_dates, cum_ret.addaptive, col = "green", lwd = 2)
  lines(out_sample_dates, cum_ret.sample, col = "orange", lwd = 2)
  lines(out_sample_dates, cum_ret.equal, col = "yellow", lwd = 2)
  lines(out_sample_dates, cum_ret.bdeco, col = "black", lwd = 2)
  lines(out_sample_dates, cum_ret.combined, col = "pink", lwd = 2)
  lines(out_sample_dates, cum_ret.index, col = "darkgrey", lwd = 2)
  legend("bottomleft", legend=c("Sample estimates", "DECO-GARCH", "DCC-GARCH", "Adaptive Thresholding", "1/N",
                                "BLOCK-DECO",
                                "Forecast combinations", "Index-Benchmark"),
         col=c("orange", "red", "blue", "green", "yellow", "black", "pink", "grey"), lty=1:1, cex=0.6)


}


### performance evaluation including factor model
performances_measures_incl_factor <- function(){
  # get correlation between portfolio returns and the SP500
  corr.dcc =  cor(dcc.port.ret, sp500.data.sample$returns)
  corr.deco = cor(deco.port.ret, sp500.data.sample$returns)
  corr.bdeco = cor(bdeco.port.ret, sp500.data.sample$returns)
  corr.addaptive = cor(addaptive.port.ret, sp500.data.sample$returns)
  corr.sample = cor(sample.port.ret, sp500.data.sample$returns)
  corr.equal =  cor(equal.port.ret, sp500.data.sample$returns)
  corr.combined = cor(combined.port.ret, sp500.data.sample$returns)
  corr.factor = cor(factor.port.ret, sp500.data.sample$returns)
  
  # Annualized daily mean 
  mean.dcc = mean(dcc.port.ret) * 252
  mean.deco = mean(deco.port.ret) * 252
  mean.addaptive = mean(addaptive.port.ret) * 252
  mean.sample = mean(sample.port.ret) * 252
  mean.equal = mean(equal.port.ret) * 252
  mean.combined = mean(combined.port.ret) * 252
  mean.bdeco = mean(bdeco.port.ret) * 252
  mean.factor = mean(factor.port.ret) * 252
  
  # Annualized daily turnover 
  TO.dcc = (1 / (num_forc)) * sum(dcc.weight.diff.sum) 
  TO.deco =(1 / (num_forc)) * sum(deco.weight.diff.sum) 
  TO.addaptive = (1 / (num_forc)) * sum(addaptive.weight.diff.sum) 
  TO.sample = (1 / (num_forc)) * sum(sample.weight.diff.sum) 
  TO.equal = (1 / (num_forc)) * sum(equal.weight.diff.sum) 
  TO.bdeco = (1 / (num_forc)) * sum(bdeco.weight.diff.sum)
  TO.combined = (1 / (num_forc)) * sum(combined.weight.diff.sum)
  TO.factor= (1 / (num_forc)) * sum(factor.weight.diff.sum)
  
  # Annualized daily Sharpe-ratios 
  SR.dcc = (mean(dcc.port.ret) / sd(dcc.port.ret)) * sqrt(252)
  SR.deco = (mean(deco.port.ret) / sd(deco.port.ret)) * sqrt(252)
  SR.addaptive = (mean(addaptive.port.ret) / sd(addaptive.port.ret)) * sqrt(252)
  SR.sample = (mean(sample.port.ret) / sd(sample.port.ret)) * sqrt(252)
  SR.equal = (mean(equal.port.ret) / sd(equal.port.ret)) * sqrt(252)
  SR.bdeco = (mean(bdeco.port.ret) / sd(bdeco.port.ret)) * sqrt(252)
  SR.combined = (mean(combined.port.ret) / sd(combined.port.ret)) * sqrt(252)
  SR.factor = (mean(factor.port.ret) / sd(factor.port.ret)) * sqrt(252)
  
  # Annualized daily volatility for all methods
  vol.dcc = sd(dcc.port.ret) * sqrt(252)
  vol.deco = sd(deco.port.ret) * sqrt(252)
  vol.addaptive = sd(addaptive.port.ret) * sqrt(252)
  vol.sample = sd(sample.port.ret) * sqrt(252)
  vol.equal = sd(equal.port.ret) * sqrt(252)
  vol.bdeco = sd(bdeco.port.ret) * sqrt(252)
  vol.combined = sd(combined.port.ret)* sqrt(252)
  vol.factor = sd(factor.port.ret)* sqrt(252)
  
  # Value at Risk 
  VaR.dcc = quantile(dcc.port.ret, 0.05)
  VaR.deco = quantile(deco.port.ret, 0.05)
  VaR.addaptive = quantile(addaptive.port.ret, 0.05)
  VaR.sample = quantile(sample.port.ret, 0.05)
  VaR.equal = quantile(equal.port.ret, 0.05)
  VaR.bdeco = quantile(bdeco.port.ret, 0.05)
  VaR.combined = quantile(combined.port.ret, 0.05)
  VaR.factor = quantile(factor.port.ret, 0.05)
  
  # Expected shortfall 
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
  
  # Stack result
  All_mean = c(mean.dcc, mean.deco, mean.bdeco, mean.addaptive , mean.sample, mean.combined, mean.equal, mean.factor)
  All_TO = c(TO.dcc, TO.deco, TO.bdeco, TO.addaptive , TO.sample, TO.combined,  TO.equal, TO.factor) * 252
  All_SR = c(SR.dcc, SR.deco, SR.bdeco, SR.addaptive , SR.sample, SR.combined, SR.equal, SR.factor)
  All_vol = c(vol.dcc, vol.deco, vol.bdeco, vol.addaptive , vol.sample, vol.combined,  vol.equal, vol.factor)
  C = c(corr.dcc, corr.deco, corr.bdeco, corr.addaptive, corr.sample, corr.combined, corr.equal, corr.factor)
  All_ES = c(ES.dcc, ES.deco,ES.bdeco,  ES.addaptive , ES.sample, ES.combined, ES.equal, ES.factor)
  
  # create dataframe with performance measures
  methods <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Sample Estimates", "Forecast Combinations", "1/N" ,"Factor model")
  measures <- c("Annualized mean", "Annualized vol", "Annualized Sharpe-ratio", "Annualized turnover", "Daily ES 95%", "Correlation S&P 500")
  performances = data.frame(All_mean, All_vol, All_SR, All_TO, abs(All_ES), C)
  names(performances) <- measures
  rownames(performances) <- methods
  
  # create directly a latex table 
  print(xtable(performances, 2))
  
  ### Statistical tests
  # Sharpe-ratio test
  dcc.sr.test = hac_infer(cbind(dcc.port.ret, sample.port.ret), digits = 4, type = "SR")
  deco.sr.test = hac_infer(cbind(deco.port.ret, sample.port.ret), digits = 4, type = "SR")
  bdeco.sr.test = hac_infer(cbind(bdeco.port.ret, sample.port.ret), digits = 4, type = "SR")
  addaptive.sr.test = hac_infer(cbind(addaptive.port.ret, sample.port.ret), digits = 4, type = "SR")
  combined.sr.test = hac_infer(cbind(combined.port.ret, sample.port.ret), digits = 4, type = "SR")
  equal.sr.test = hac_infer(cbind(equal.port.ret, sample.port.ret), digits = 4, type = "SR")
  factor.sr.test = hac_infer(cbind(factor.port.ret, sample.port.ret), digits = 4, type = "SR")
  
  # Volatility test
  dcc.vol.test = hac_infer(cbind(dcc.port.ret, sample.port.ret), digits = 4, type = "Var")
  deco.vol.test = hac_infer(cbind(deco.port.ret, sample.port.ret), digits = 4, type = "Var")
  bdeco.vol.test = hac_infer(cbind(bdeco.port.ret, sample.port.ret), digits = 4, type = "Var")
  addaptive.vol.test = hac_infer(cbind(addaptive.port.ret, sample.port.ret), digits = 4, type = "Var")
  combined.vol.test = hac_infer(cbind(combined.port.ret, sample.port.ret), digits = 4, type = "Var")
  equal.vol.test = hac_infer(cbind(equal.port.ret, sample.port.ret), digits = 4, type = "Var")
  factor.vol.test = hac_infer(cbind(factor.port.ret, sample.port.ret), digits = 4, type = "Var")
  
  SR_test = matrix(ncol = 7, nrow = 3)
  Vol_test = matrix(ncol = 7, nrow = 3)
  SR_test[1,] = c(dcc.sr.test$Difference, deco.sr.test$Difference, bdeco.sr.test$Difference, addaptive.sr.test$Difference
                  , combined.sr.test$Difference, equal.sr.test$Difference, factor.sr.test$Difference)
  SR_test[2,] = c(dcc.sr.test$Standard.Errors[1], deco.sr.test$Standard.Errors[1],
                  bdeco.sr.test$Standard.Errors[1], addaptive.sr.test$Standard.Errors[1]
                  , combined.sr.test$Standard.Errors[1], equal.sr.test$Standard.Errors[1], factor.sr.test$Standard.Errors[1])
  SR_test[3,] = c(dcc.sr.test$p.Values[1], deco.sr.test$p.Values[1],
                  bdeco.sr.test$p.Values[1], addaptive.sr.test$p.Values[1]
                  , combined.sr.test$p.Values[1], equal.sr.test$p.Values[1], factor.sr.test$p.Values[1])
  
  Vol_test[1,] = c(dcc.vol.test$Difference, deco.vol.test$Difference, bdeco.vol.test$Difference, addaptive.vol.test$Difference
                   , combined.vol.test$Difference, equal.vol.test$Difference, factor.vol.test$Difference)
  Vol_test[2,] = c(dcc.vol.test$Standard.Errors[1], deco.vol.test$Standard.Errors[1],
                   bdeco.vol.test$Standard.Errors[1], addaptive.vol.test$Standard.Errors[1]
                   , combined.vol.test$Standard.Errors[1], equal.vol.test$Standard.Errors[1],
                   factor.vol.test$Standard.Errors[1])
  Vol_test[3,] = c(dcc.vol.test$p.Values[1], deco.vol.test$p.Values[1],
                   bdeco.vol.test$p.Values[1], addaptive.vol.test$p.Values[1]
                   , combined.vol.test$p.Values[1], equal.vol.test$p.Values[1],
                   factor.vol.test$p.Values[1])
  SR_test = data.frame(SR_test)
  row.names(SR_test) <- c("Delta.diff", "HAC se", "pvalue")
  names(SR_test) <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Forecast Combinations", "1/N", "Factor model")
  
  Vol_test = data.frame(Vol_test)
  row.names(Vol_test) <- c("Delta.diff", "HAC se", "pvalue")
  names(Vol_test) <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Forecast Combinations", "1/N", "Factor model")
  
  print(xtable(SR_test, 2))
  print(xtable(Vol_test, 2))
  
  
  # make a plot of the cummulative returns
  cum_ret.dcc = cum_ret(dcc.port.ret)
  cum_ret.deco = cum_ret(deco.port.ret)
  cum_ret.bdeco = cum_ret(bdeco.port.ret)
  cum_ret.sample = cum_ret(sample.port.ret)
  cum_ret.addaptive = cum_ret(addaptive.port.ret)
  cum_ret.combined = cum_ret(combined.port.ret)
  cum_ret.equal = cum_ret(equal.port.ret)
  cum_ret.index = cum_ret(index.returns[(T+1):size_sample])
  cum_ret.factor = cum_ret(factor.port.ret)
  
  # note that if we would do monthly or other types of rebalancing adjust the title
  plot(out_sample_dates, cum_ret.dcc, type = "l",
       col = "blue", xlab = "Time", ylab = "Cumulative return (%)", ylim = c(-60,35), lwd = 2)
  lines(out_sample_dates, cum_ret.deco, col = "red", lwd = 2)
  lines(out_sample_dates, cum_ret.addaptive, col = "green", lwd = 2)
  lines(out_sample_dates, cum_ret.sample, col = "orange", lwd = 2)
  lines(out_sample_dates, cum_ret.equal, col = "yellow", lwd = 2)
  lines(out_sample_dates, cum_ret.bdeco, col = "black", lwd = 2)
  lines(out_sample_dates, cum_ret.combined, col = "pink", lwd = 2)
  lines(out_sample_dates, cum_ret.index, col = "darkgrey", lwd = 2)
  lines(out_sample_dates, cum_ret.factor, col = "darkred", lwd = 2)
  legend("bottomleft", legend=c("Sample estimates", "DECO-GARCH", "DCC-GARCH", "Adaptive Thresholding", "1/N",
                                "BLOCK-DECO",
                                "Forecast combinations", "Index-Benchmark", "Factor model"),
         col=c("orange", "red", "blue", "green", "yellow", "black", "pink", "grey", "darkred"), lty=1:1, cex=0.6)

  
  
  
}



rebalancing_frequencies_performance <- function(){
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
    bdeco.res = portfolio_returns_out_sample(rebalance_freq, bdec_cov)
    dcc.res = portfolio_returns_out_sample(rebalance_freq, dcc_cov)
    deco.res = portfolio_returns_out_sample(rebalance_freq, dec_cov)
    addaptive.res= portfolio_returns_out_sample(rebalance_freq, addaptive_cov)
    sample.res = portfolio_returns_out_sample(rebalance_freq, sample_cov)
    equal.res = equally_weighted_portfolio_returns(rebalance_freq)
    
    # in-sample returns given a rebalancing frequency
    in_sample_ret.dcc = portfolio_returns_in_sample(rebalance_freq,  dcc.cov.insample, T, ret_sample[1:T,])
    in_sample_ret.deco = portfolio_returns_in_sample(rebalance_freq,  dec_cov_insample, T, ret_sample[1:T,])
    in_sample_ret.bdeco = portfolio_returns_in_sample(rebalance_freq,  bdec_cov_insample, T, ret_sample[1:T,])
    in_sample_ret.addaptive = portfolio_returns_in_sample(rebalance_freq,  addaptive_cov_insample, T, ret_sample[1:T,])
    in_sample_ret.sample = portfolio_returns_in_sample(rebalance_freq,  sample_cov_insample, T, ret_sample[1:T,])
    
    # create full sample returns
    all_ret.dcc = c(in_sample_ret.dcc, dcc.res[[1]])
    all_ret.deco = c(in_sample_ret.deco, deco.res[[1]])
    all_ret.bdeco = c(in_sample_ret.bdeco, bdeco.res[[1]])
    all_ret.addaptive = c(in_sample_ret.addaptive, addaptive.res[[1]])
    all_ret.sample = c(in_sample_ret.sample, sample.res[[1]])
    
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
    combined.res = portfolio_returns_out_sample(rebalance_freq, combined_cov)
    
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
  # create dataframe with sharpe-ratios for all rebalance frequencies
  methods <- c("DCC-GARCH", "DECO-GARCH", "BLOCK-DECO", "Adaptive Thresholding", "Sample Estimates", "Forecast Combinations", "1/N")
  volatility_reb_freq = data.frame(dcc2, deco2, bdeco2, addaptive2, sample2, combined2, equal2)
  names(volatility_reb_freq) <- methods
  
  sharpe_reb_freq = data.frame(dcc, deco, bdeco, addaptive, sample, combined, equal)
  names(sharpe_reb_freq) <- methods
  return(list(volatility_reb_freq, sharpe_reb_freq))
}


plot_vol <- function(df.vol){
  methods = c("DCC-GARCH", "DECO-GARCH","BLOCK-DECO",
              "Adaptive Thresholding","Sample estimates", "Forecast combinations", "1/N")
  # plot sharpe-ratios for all rebalancing frequencies for all methods
  plot(c(1:250), df.vol$`DCC-GARCH`, type = "l",
       col = "blue", xlab = "Rebalancing frequency", ylab = "Annualized Volatility (%)", lwd = 2, ylim = c(8, 22))
  lines(c(1:250), df.vol$`DECO-GARCH`, col = "red", lwd = 2)
  lines(c(1:250), df.vol$`BLOCK-DECO`, col = "black", lwd = 2)
  lines(c(1:250), df.vol$`Adaptive Thresholding`, col = "green", lwd = 3)
  lines(c(1:250), df.vol$`Sample Estimates`, col = "orange", lwd = 1)
  lines(c(1:250), df.vol$`Forecast Combinations`, col = "pink", lwd = 2)
  lines(c(1:250), df.vol$`1/N`, col = "yellow", lwd = 2)
  legend("topleft", legend= methods,
         col=c("blue", "red", "black", "green", "orange", "pink", "yellow"), lty=1:1, cex=0.55)
}


plot_sharpe <- function(df.sharpe){
  methods = c("DCC-GARCH", "DECO-GARCH","BLOCK-DECO",
              "Adaptive Thresholding","Sample estimates", "Forecast combinations", "1/N")
  # plot sharpe-ratios for all rebalancing frequencies for all methods
  plot(c(1:250), df.sharpe$`DCC-GARCH`, type = "l",
       col = "blue", xlab = "Rebalancing frequency", ylab = "Annualized Sharpe Ratio", lwd = 2, ylim = c(-2, 0))
  lines(c(1:250), df.sharpe$`DECO-GARCH`, col = "red", lwd = 2)
  lines(c(1:250), df.sharpe$`BLOCK-DECO`, col = "black", lwd = 2)
  lines(c(1:250), df.sharpe$`Adaptive Thresholding`, col = "green", lwd = 3)
  lines(c(1:250), df.sharpe$`Sample Estimates`, col = "orange", lwd = 1)
  lines(c(1:250), df.sharpe$`Forecast Combinations`, col = "pink", lwd = 2)
  lines(c(1:250), df.sharpe$`1/N`, col = "yellow", lwd = 2)
  legend("topright", legend= methods,
         col=c("blue", "red", "black", "green", "orange", "pink", "yellow"), lty=1:1, cex=0.6)
}



    
    
    
    
    
  



