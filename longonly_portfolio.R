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
library(Rsolnp)


#----------------------------------------------------------
#----Variance function and equality constraint function----
#----------------------------------------------------------

min_var <- function(w, Sigma){
  result = t(w) %*% Sigma %*%  w
  return(result)
}
equal.fun <- function(w, Sigma) {
  z1 <- (sum(w))
  return(z1)
}


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
    
    # calculate the weights with restrictions
    #pi.gmv.start <- (solve(Cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(Cov.mat) %*% vOnes)[1])
    result.longonly <-  solnp(pars =  rep(1/27, 27), Sigma = Cov.mat,
                              fun = min_var,
                              #function to optimise
                              eqfun=equal.fun, #equality function 
                              eqB=1,   #the equality constraint
                              LB=rep(0, 27), #lower bound for parameters i.e. greater than zero
                              UB=rep(4, 27), control = list(trace = 0, tol = 1e-6)) #upper bound for parameters (I just chose 100 randomly)
    pi.gmv.start = result.longonly$pars
    
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
    #pi.gmv.start <- (solve(Cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(Cov.mat) %*% vOnes)[1])
    result.longonly <-  solnp(pars =  rep(1/27, 27), Sigma = Cov.mat,
                              fun = min_var, #function to optimise
                              eqfun=equal.fun, #equality function 
                              eqB=1,   #the equality constraint
                              LB=rep(0, 27), #lower bound for parameters i.e. greater than zero
                              UB=rep(4, 27), control = list(trace = 0, tol = 1e-6)) #upper bound for parameters (I just chose 100 randomly)
    pi.gmv.start = result.longonly$pars
    
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
        #gmv_weights = (solve(cov.mat) %*% vOnes) / ((t(vOnes) %*% solve(cov.mat) %*% vOnes)[1])
        result.longonly <-  solnp(pars =  rep(1/27, 27), Sigma = cov.mat,
                                  fun = min_var,
                                  #function to optimise
                                  eqfun=equal.fun, #equality function 
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

# store Addaptive tresholding forecasted covariance
addaptive_cov <- array(dim = c(num_assets,num_assets,num_forc)) 

# estimation window for the static sparse estimators
estimation_window_static = T

# define the grid we will search the best threshold value from
mthr <- seq(from=0.01,to=0.99,length.out=15)

# estimate the covariance matrices and get the static estimators
est_sample = ret_sample[1:T,]

# use 5-fold cross validation with a time-series split
lambda_cv_scores = cross_validation_addaptive(mthr, 5, est_sample)
treshold.star = mthr[which.min(lambda_cv_scores)]

# estimate the static (sparse) covariance matrices
est_cov.mat.addaptive = CovEst.adaptive(est_sample, thr= treshold.star)
est_cov.mat.addaptive = est_cov.mat.addaptive$S

# store them
for (i in 1:num_forc){
  # store the estimated covariance matrices without including new data
  addaptive_cov[,,i] <- est_cov.mat.addaptive
  
}

# get insample covariance matrices
addaptive_cov_insample <- array(dim = c(num_assets,num_assets,T)) 
for (i in 1:T){
  # store the estimated covariance matrices without including new data
  addaptive_cov_insample[,,i] <- est_cov.mat.addaptive
  
}







