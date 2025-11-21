
# Library with various utility based portfolio choice models
#
# Costas Xiouros (2022)


###############################################################################
# FUNCTION:                DESCRIPTION:
#  portCRRA                 Returns portfolio weights of risky assets that 
#                           maximizes expected CRRA utility, with risk-free asset
#  portCRRA_noRf            Returns portfolio weights that maximizes expected
#                           CRRA utility, without risk-free asset
#  portCRRA_logn_p          Returns portfolio weights assuming log-normally
#                           distributed returns by providing the moments

###############################################################################
portCRRA <- function(data, Rf=1, model="normal", gamma=2, minlev=0, maxlev=1, lower=0, upper=1, n_opts=5)
{
  #ARGUMENTS:
  # data              Excess return data (T x n, n: number of risky assets)
  # Rf                Gross risk-free rate
  # model             "normal", "lognorm" or "hist"
  # gamma             Coefficient of relative risk aversion
  # minlev            Minimum leverage on the portfolio of risky assets, where
  #                   (1 - leverage) is the position on the risk-free asset
  # maxlev            Maximum leverage on the portfolio of risky assets
  # lower             Lower portfolio bounds (scalar or n x 1 vector)
  # upper             Upper portfolio bounds (scalar or n x 1 vector)
  # n_opts            Number of optimization runs with different random starts
  
  #number of assets
  n <- ncol(data) 
  
  #Normality assumption
  #====================
  if (model=="normal")
  {
    #mean vector and covariance matrix
    Mu <- colMeans(data_eret)
    Sigma <- cov(data_eret)
    
    #Find the optimal tangent portfolio
    tgW <- conTangent(Mu,Sigma,lower,upper)
    
    #compute expected excess return and standard deviation
    tgStats <- c(tgW %*% Mu, sqrt(tgW %*% Sigma %*% tgW))
    
    #univariate Gauss-Hermite quadrature
    source("variousTools.R")
    pts <- gauss.hermite(20)
    p <- pts[,2]
    Re <- tgStats[1] + tgStats[2] * pts[,1]
    
    #objective function
    certeq <- function(theta){
      util <- (Rf + theta * Re)^(1-gamma)
      (util %*% p)^(1/(1-gamma))
    }
    
    #expected utility maximization
    optPort <- optimize(certeq, c(minlev, maxlev), maximum=TRUE)
    theta <- optPort$maximum
    w <- theta * tgW
    
  } else {
    
    #risky asset lower and upper bounds
    if (length(upper)==1){
      ub <- rep(upper, n)
    } else {
      ub <- upper
    }
    if (length(lower)==1){
      lb <- rep(lower, n)
    } else {
      lb <- lower
    }
    
    #portfolio constraints:
    # ui %*% theta - ci >= 0
    #=======================
    #Unit vector and matrix
    uVn <- rep(1, times=n)
    uMn <- diag(n)
    
    #ui matrix: k x n
    uiM <- rbind(uVn, -uVn, uMn, -uMn)
    
    #ci vector: k x 1
    ciV <- c(minlev, -maxlev, lb, -ub)
    
    #objective function (theta: (n+1) x 1, Re: N x n)
    certeq <- function(theta){
      util <- (Rf + Re %*% theta)^(1-gamma)
      (p %*% util)^(1/(1-gamma))
    }
    
    #Log-normality assumption
    #========================
    if (model=="lognorm"){
      
      #log-returns: not exact since use Rf
      r_data = log(Rf + data)
      Mu <- colMeans(r_data)
      Sigma <- cov(r_data)
      
      #multivariate Gauss-Hermite quadrature
      source("variousTools.R")
      if (n==2){
        prn <- 0.1
        N <- 15
      } else if (n <= 4) {
        prn <- 0.2
        N <- 10
      } else {
        prn <- 0.4
        N <- round(10000^(1/n), digits=0)
      }
      pts <- mgauss.hermite(N, mu=Mu, sigma=Sigma, prune=prn)
      Re <- exp(pts$points) - Rf
      p <- pts$weights
    }
    
    #Non-parametric: historical simulations
    #======================================
    if (model=="hist"){
      Re <- as.matrix(data)
      T <- length(data[,1])
      p <- rep(1/T, times=T)
    }
    
    #multistart random
    max_val <- 0
    for (i in 1:n_opts){
      #random initial point for portfolio optimization
      # draw random points until constraints are satisfied
      min_c <- -1
      while (min_c < 0){
        x0 <- lb + runif(n, 0, 1)*(ub - lb)
        lev <- runif(1, minlev, maxlev)
        x0 <- lev * x0 / sum(x0)
        min_c <- min(uiM %*% x0 - ciV)
      }
      
      #optimization
      optPort <- constrOptim(x0, certeq, NULL,  ui=uiM, ci=ciV, control=list(fnscale=-1))
      
      #check if better result
      if (optPort$value > max_val){
        w <- optPort$par
        max_val <- optPort$value
      }
    }
    
  }
  
  #add column names
  names(w) <- names(data)
  
  #return value
  w
}

###############################################################################
portCRRA_noRf <- function(data, model="hist", gamma=2, lower=0, upper=1, n_opts=5)
{
  #ARGUMENTS:
  # data              Gross return data (T x n, n:number of assets)
  # model             "normal" or "hist" or "lognorm"
  # N                 Number of simulations
  # gamma             Coefficient of relative risk aversion
  # lower             Lower portfolio bounds (scalar or n x 1 vector)
  # upper             Upper portfolio bounds (scalar or n x 1 vector)
  # n_opts             Number of optimizations with different starting values
  
  #number of assets
  n <- ncol(data) 
  
  #portfolio constraints
  if (length(upper)==1){
    ub <- rep(upper, n)
  } else {
    ub <- upper
  }
  if (length(lower)==1){
    lb <- rep(lower, n)
  } else {
    lb <- lower
  }
  
  if (model != "hist"){
    #multivariate Gauss-Hermite quadrature
    source("variousTools.R")
    if (n==2){
      prn <- 0.1
      N <- 15
    } else {
      prn <- 0.2
      N <- 10
    }
    
  }
  
  #Simulations with normal distribution
  if (model=="normal")
  {
    Mu <- colMeans(data)
    Sigma <- cov(data)
    
    #multivariate Gauss-Hermite quadrature
    pts <- mgauss.hermite(N, mu=Mu, sigma=Sigma, prune=prn)
    R <- pts$points
    p <- pts$weights
  }
  
  #Non-parametric: historical simulations
  if (model=="hist"){
    R <- data
    T <- length(data[,1])
    p <- rep(1/T, times=T)
  }
  
  #Simulations with log-normal distribution
  #========================================
  if (model=="lognorm"){
    #log-returns
    r_data = log(data)
    Mu <- colMeans(r_data)
    Sigma <- cov(r_data)
    
    #multivariate Gauss-Hermite quadrature
    pts <- mgauss.hermite(N, mu=Mu, sigma=Sigma, prune=prn)
    R <- pts$points
    p <- exp(pts$weights)
  }
  
  #objective function (theta: n x 1, Re: N x n)
  certeq <- function(theta){
    util <- (R %*% theta)^(1-gamma)
    (p %*% util)^(1/(1-gamma))
  }
  
  #portfolio constraints:
  # ui %*% theta - ci >= 0
  #=======================
  #Unit vector and matrix
  uVn <- rep(1, times=n)
  uMn <- diag(n)
  
  #ui matrix: k x n
  uiM <- rbind(-uVn, uMn, -uMn)
  
  #ci vector: k x 1
  ciV <- c(-1, lb, -ub)
  
  #multistart random
  max_val <- 0
  for (i in 1:n_opts){
    #random initial point for portfolio optimization
    # draw random points until constraints are satisfied
    min_c <- -1
    while (min_c < 0){
      x0 <- lb + runif(n, 0, 1)*(ub - lb)
      x0 <- 0.95 * x0 / sum(x0)
      min_c <- min(uiM %*% x0 - ciV)
    }
    
    #optimization
    optPort <- constrOptim(x0, certeq, NULL,  ui=uiM, ci=ciV, control=list(fnscale=-1))
    
    #check if better result
    if (optPort$value > max_val){
      w <- optPort$par
      max_val <- optPort$value
    }
  }
  
  #add column names
  names(w) <- names(data)
  
  #return value
  w
}

###############################################################################
portCRRA_LT <- function(data, model="logn", T=1, rf=NA, gamma=2, minlev=0, maxlev=1, lower=0, upper=1, n_opts=5, N=1000)
{
  #ARGUMENTS:
  # data              Gross return data (L x n, n:number of assets)
  # model             Log-normal ("logn") or historical ("hist")
  # T                 Investment horizon
  # rf                Log risk-free rate. NA for no risk-free asset
  # gamma             Coefficient of relative risk aversion
  # minlev            Minimum leverage on the portfolio of risky assets, where
  #                   (1 - leverage) is the position on the risk-free asset
  # maxlev            Maximum leverage on the portfolio of risky assets
  # lower             Lower portfolio bounds (scalar or n x 1 vector)
  # upper             Upper portfolio bounds (scalar or n x 1 vector)
  # n_opts            Number of optimizations with different starting values
  # N                 Number of scenarios
  
  #number of risky assets
  n <- length(data[1,])
  
  #scenarios of long-term returns
  if (model=="hist"){
    #historical simulations
    L <- length(data[,1])
    p <- rep(1/N, times=N)
    
    #initialize matrix to store simulations
    R <- matrix(NA, nrow=N, ncol=n)
    
    #generate simulations
    for (i in seq(1,N)){
      #sample sets of returns
      I <- sample(L,T,replace=TRUE)
      
      #compute cumulative returns
      R[i,] = apply(data[I,], 2, prod)
    }
    
  } else {
    #log-normal distribution
    source("variousTools.R")
    
    #means and covariance matrix of log-returns
    data <- log(data)
    Mu <- colMeans(data) * T
    Sigma <- cov(data) * T
    
    #quadrature points
    if (n==1){
      #univariate Gauss-Hermite quadrature
      pts <- gauss.hermite(20)
      p <- pts[,2]
      R <- exp(Mu + c(sqrt(Sigma)) * pts[,1])
    } else {
      #multivariate Gauss-Hermite quadrature
      if (n==2){
        prn <- 0.1
        n_pts <- 15
      } else {
        prn <- 0.2
        n_pts <- ceiling((N / (1-prn))^(1/n))
      }
      pts <- mgauss.hermite(n_pts, mu=Mu, sigma=Sigma, prune=prn)
      R <- exp(pts$points)
      p <- pts$weights
    }
  }
  
  #risky asset weight bounds
  if (length(upper)==1){
    ub <- rep(upper, n)
  } else {
    ub <- upper
  }
  if (length(lower)==1){
    lb <- rep(lower, n)
  } else {
    lb <- lower
  }
  
  
  if (is.na(rf)){
    
    #objective function (theta: n x 1, R: N x n)
    certeq <- function(theta){
      util <- (R %*% theta)^(1-gamma)
      (p %*% util)^(1/(1-gamma))
    }
    
    #portfolio constraints:
    # ui %*% theta - ci >= 0
    #=======================
    #Unit vector and matrix
    uVn <- rep(1, times=n)
    uMn <- diag(n)
    
    #ui matrix: k x n
    uiM <- rbind(-uVn, uMn, -uMn)
    
    #ci vector: k x 1
    ciV <- c(-1, lb, -ub)
    
    #multistart random
    max_val <- 0
    for (i in 1:n_opts){
      #random initial point for portfolio optimization
      # draw random points until constraints are satisfied
      min_c <- -1
      while (min_c < 0){
        x0 <- lb + runif(n, 0, 1)*(ub - lb)
        x0 <- 0.95 * x0 / sum(x0)
        min_c <- min(uiM %*% x0 - ciV)
      }
      
      #optimization
      optPort <- constrOptim(x0, certeq, NULL,  ui=uiM, ci=ciV, control=list(fnscale=-1))
      
      #check if better result
      if (optPort$value > max_val){
        w <- optPort$par
        max_val <- optPort$value
      }
    }
    
  } else {
    
    #Gross risk-free rate
    Rf <- exp(rf*T)
    
    #excess returns
    Re <- R - Rf
    
    #optimal portfolios
    if (n==1){
      
      #objective function (theta: n x 1, Re: N x n)
      certeq <- function(theta){
        util <- (Rf + Re * theta)^(1-gamma)
        (util %*% p)^(1/(1-gamma))
      }
      
      #optimization
      optPort <- optimize(certeq, c(minlev, maxlev), maximum=TRUE)
      w <- optPort$maximum
      max_val <- optPort$objective
      
    } else {
      
      #portfolio constraints:
      # ui %*% theta - ci >= 0
      #=======================
      #Unit vector and matrix
      uVn <- rep(1, times=n)
      uMn <- diag(n)
      
      #ui matrix: k x n
      uiM <- rbind(uVn, -uVn, uMn, -uMn)
      
      #ci vector: k x 1
      ciV <- c(minlev, -maxlev, lb, -ub)
      
      #objective function (theta: (n+1) x 1, Re: N x n)
      certeq <- function(theta){
        util <- (Rf + Re %*% theta)^(1-gamma)
        (p %*% util)^(1/(1-gamma))
      }
      
      #multistart random
      max_val <- 0
      for (i in 1:n_opts){
        #random initial point for portfolio optimization
        # draw random points until constraints are satisfied
        min_c <- -1
        while (min_c < 0){
          x0 <- lb + runif(n, 0, 1)*(ub - lb)
          lev <- runif(1, minlev, maxlev)
          x0 <- lev * x0 / sum(x0)
          min_c <- min(uiM %*% x0 - ciV)
        }
        
        #optimization
        optPort <- constrOptim(x0, certeq, NULL,  ui=uiM, ci=ciV, control=list(fnscale=-1))
        
        #check if better result
        if (optPort$value > max_val){
          w <- optPort$par
          max_val <- optPort$value
        }
      }
    }
    
  }
  
  #add column names
  names(w) <- names(data)
  
  #return value
  res <- list(w, max_val)
  names(res) <- c("w","ce")
  res
  
}