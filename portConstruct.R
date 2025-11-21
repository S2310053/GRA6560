
# Library with various portfolio construction methods
#
# Costas Xiouros (2023)


###############################################################################
# FUNCTION:                DESCRIPTION:
#  portMinVar               Returns the unconstrained or constrained minimum 
#                           variance portfolio
#  uncTangent               Returns the unconstrained mean-variance tangent
#                           portfolio
#  conTangent               Returns the constrained mean-variance tangent
#                           portfolio
#  uncMeanVar_ra            Unconstrained mean-variance optimal portfolio for
#                           given risk aversion parameter "gamma"
#  uncMeanVar_mu            Unconstrained mean-variance optimal portfolio for
#                           a target mean
#  uncMeanVar_sd            Unconstrained mean-variance optimal portfolio for
#                           a target standard deviation
#  conMeanVar               Constrained mean-variance optimal portfolio for a
#                           given risk-aversion, or for target mean, or for 
#                           a target volatility
#  naiveVolRiskBudget       Returns portfolio weights based on naive volatility
#                           risk budgeting
#  naiveVarRiskBudget       Returns portfolio weights based on naive variance
#                           risk budgeting
#  trueVolRiskBudget        Returns portfolio weights based on true risk
#                           budgeting with respect to volatility contributions
#  maxDiversify             Returns portfolio weights of the unconstrained
#                           maximum diversification portfolio

library("matlib")

#=================================================
portMinVar <- function(Sigma, lower=NULL, upper=NULL)
{
  #number of assets
  n <- length(Sigma[,1]) 
  
  if (missing(lower)){
    #no constraints
    #==============
    #unit vector of size n
    uVn <- rep(1, times=n)
  
    w <- c(inv(Sigma) %*% uVn)
  } else {
    #upper and lower portfolio bounds
    #================================
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
    
    #objective function
    pvar <- function(theta){
      theta %*% Sigma %*% theta
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
    
    #optimization
    x0 <- lb
    rem <- 0.99 - sum(lb)
    x0 <- x0 + rep(rem/n, times=n)
    optPort <- constrOptim(x0, pvar, NULL,  ui=uiM, ci=ciV)
    w <- optPort$par
  }

  #return value
  names(w) <- colnames(Sigma)
  w / sum(w)
}

#=================================================
uncMeanVar_ra <- function(Mu, Sigma, gamma)
{
  #number of assets
  n <- length(Mu) 
  
  #unit vector of size n
  uVn <- rep(1, times=n)
  
  #lambda (Langrange multiplier)
  a <- c(uVn %*% inv(Sigma) %*% Mu)
  b <- c(uVn %*% inv(Sigma) %*% uVn)
  lambda <- (a - gamma)/b
  
  #Mu should be a column vector
  Mu_e <- Mu - lambda*uVn

  w <- c(inv(Sigma) %*% Mu_e) / gamma
  names(w) <- colnames(Sigma)
  
  #return value
  w / sum(w)
}

#=================================================
uncMeanVar_mu <- function(Mu, Sigma, mu)
{
  w1 <- uncMeanVar_ra(Mu, Sigma, 1)
  w2 <- uncMeanVar_ra(Mu, Sigma, 2)
  
  mu1 <- c(w1 %*% Mu)
  mu2 <- c(w2 %*% Mu)
  ww1 <- (mu - mu2)/(mu1 - mu2)
  
  w <- ww1*w1 + (1-ww1)*w2
  
  #return value
  w / sum(w)
}

#=================================================
uncMeanVar_sd <- function(Mu, Sigma, sd)
{
  #target variance
  v <- sd^2
  
  #two optimal portfolios
  w1 <- uncMeanVar_ra(Mu, Sigma, 1)
  w2 <- uncMeanVar_ra(Mu, Sigma, 2)

  #their means
  mu1 <- c(w1 %*% Mu)
  mu2 <- c(w2 %*% Mu)
  
  #their variances and covariance
  v1 <- c(w1 %*% Sigma %*% w1)
  v2 <- c(w2 %*% Sigma %*% w2)
  v12 <- c(w1 %*% Sigma %*% w2)
  
  #combinations to achieve target variance
  a <- v1 + v2 - 2*v12
  b <- 2*(v2 - v12)
  c <- v2 - v
  
  if (b^2 < 4*a*c){
    print("The target portfolio is infeasible!")
    w <- rep(0, times=length(Mu))
  } else {
    wp <- c((b + sqrt(b^2 - 4*a*c)) / (2*a))
    wm <- c((b - sqrt(b^2 - 4*a*c)) / (2*a))

    #choose the one with the highest mean
    mu_p <- c(wp*mu1 + (1-wp)*mu2)
    mu_m <- c(wm*mu1 + (1-wm)*mu2)
    if(mu_p > mu_m) {
      w <- wp*w1 + (1-wp)*w2
    } else {
      w <- wm*w1 + (1-wm)*w2
    }
  
    #return value
    names(w) <- names(Mu)
    w / sum(w)
  }
}

#=================================================
conMeanVar <- function(Mu, Sigma, target="ra", t_val, lower=0, upper=1)
{
  #number of assets
  #Mu should be a column vector
  n <- length(Mu) 
  
  dev <- function(rf){
    #hypothetical tangent portfolio for risk-free rate rf
    x <- conTangent(Mu - rf, Sigma, lower, upper)
    mu_x <- c(x %*% Mu)
    vr_x <- c(x %*% Sigma %*% x)
    
    #match leverage depending on the target
    if (target == "ra"){
      #maximize expected utility
      t_val - (mu_x - rf) / vr_x
    } else if (target == "mu") {
      #minimize volatility
      1 - t_val/mu_x
    } else {
      #maximize expected return
      1 - t_val/sqrt(vr_x)
    }
  }
  
  #minimum variance portfolio
  x_mv <- portMinVar(Sigma, lower, upper)
  mu_mv <- c(x_mv %*% Mu)
  
  #approximate max expected return to set interval for rf
  x_max <- conTangent(Mu - mu_mv, Sigma, lower, upper)
  mu_max <- c(x_max %*% Mu)
  
  #solve for hypothetical rf
  root <- uniroot(dev, interval = c(-0.5, mu_max), tol=1.0e-9)
  rf <- root$root
  w <- conTangent(Mu - rf, Sigma, lower, upper)

  #add column names
  names(w) <- names(Mu)

  #return value
  w / sum(w)
}

#=================================================
uncTangent <- function(Mu, Sigma)
{
  #Mu should be a column vector
  w <- c(inv(Sigma) %*% Mu)
  names(w) <- colnames(Sigma)

  #return value
  w / sum(w)
}

#=================================================
conTangent <- function(Mu, Sigma, lower=0, upper=1)
{
  #number of assets
  #Mu should be a column vector
  n <- length(Mu) 
  
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
  
  #objective function
  SR <- function(theta){
    vol <- sqrt(theta %*% Sigma %*% theta)
    ere <- c(theta %*% Mu)
    ere/vol;
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
  
  #minimum variance portfolio
  w_mv <- portMinVar(Sigma, lower, upper)
  mu_mv <- c(w_mv %*% Mu)
  
  #optimization
  x0 <- lb
  rem <- 0.99 - sum(lb)
  x0 <- x0 + rep(rem/n, times=n)
  optPort <- constrOptim(x0, SR, NULL,  ui=uiM, ci=ciV, control=list(fnscale=-1))
  w <- optPort$par

  #add column names
  names(w) <- names(Mu)
  
  #return value
  w / sum(w)
}

#=================================================
naiveVolRiskBudget <- function(Sigma, b=1)
{
  #extract standard deviations from the diagonal
  vols <- sqrt(diag(Sigma))
  
  #check length of b if provided to match number of assets
  n_b <- length(b) 
  if(n_b > 1){
    if(n_b != length(vols)){
      stop("Length of b vector does not match the size of Sigma")
    }
  }

  #estimate weights (not normalized) based on risk budgets b
  w <- b / vols
  
  #return value
  w / sum(w)
}

#=================================================
naiveVarRiskBudget <- function(Sigma, b=1)
{
  #extract standard deviations from the diagonal
  vars <- diag(Sigma)
  
  #check length of b if provided to match number of assets
  n_b <- length(b) 
  if(n_b > 1){
    if(n_b != length(vars)){
      stop("Length of b vector does not match the size of Sigma")
    }
  }
  
  #estimate weights (not normalized) based on risk budgets b
  w <- b / vars
  
  #return value
  w / sum(w)
}

#=================================================
trueVolRiskBudget <- function(Sigma, b=1)
{
  #read number of assets
  n <- length(Sigma[,1])

  #set equal risk budgets if b not provided
  if (length(b) == 1){
    #create column vector
    b <- as.matrix(rep(1/n, n))
  } else {
    if(length(b) != n){
      stop("Length of b vector does not match the size of Sigma")
    }
    b <- as.matrix(b)
  }
  
  #Employ Newton's method
  I <- t(as.matrix(rep(1, n)))   #row vector of ones
  W <- b;
        
  tolw <- 1 
  lambda <- 2.5
  while (tolw > 1.0e-8) {
    #Compute function
    boW <- b / W
    F <- rbind(Sigma %*% W - lambda * boW, I %*% W - 1)

    #construct the Jacobian
    J = cbind(rbind(Sigma + lambda * diag(c(b / (W^2))), I), rbind(-boW, 0))

    #Newton step
    W_n <- rbind(W, lambda) - inv(J) %*% F

    #check convergence
    tolw <- sum(abs(W_n - rbind(W, lambda))) + sum(abs(F))

    #update solution
    W <- as.matrix(W_n[1:n])
    lambda <- W_n[1+n]
  }
  
  #return value
  w <- c(W)
  names(w) <- colnames(Sigma)
  w
}

#=================================================
maxDiversify <- function(Sigma)
{
  #extract vector of standard deviations
  sigma <- sqrt(diag(Sigma))
  
  #compute portfolio weights (not normalized)
  w <- c(sigma %*% inv(Sigma))
  names(w) <- colnames(Sigma)
  
  #return value
  w / sum(w)
}
