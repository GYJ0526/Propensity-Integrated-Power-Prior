# Required libraries
library(mvtnorm)
library(invgamma)
library(MatchIt)
library(rstantools)
library(dplyr)
library(MASS)

PS_BBPPP <- function(
    rep,
    n,            # Current sample size
    p,            # Number of covariates (without intercept)
    n0,           # Historical sample size
    b,            # True coefficients (theta = b[1])
    A,            # Nominal patients to be borrowed
    var_y,        # Variance of current data
    var_y0,       # Variance of historical data
    mu,  		  # Mean vector for current data
    mu0,          # Mean vector for historical data (no trt var)
    var_x,        # Variance for current study covariates
    var_x0,       # Variance for historical study covariates
    K,            # Number of strata
    iter,         # Number of MCMC iterations
    B,            # Number of Burn-in samples
    seed=NULL) {
  
### Storage
  uw.trt.effect.rep <- numeric(rep)
  uw.CP.trt <- numeric(rep)
  total_n0_vec <- numeric(rep)
  uw.width.CI.rep <- numeric(rep)
  
  for(r in 1:rep){
    set.seed(seed + r)
    
    cov_matrix <- matrix(0, p, p)
    diag(cov_matrix) <- var_x #variance-covariance matrix of current data covariates
    X=rmvnorm(n, mean=mu, sigma=cov_matrix) #design matrix of current data (without intercept)
    colnames(X) =paste0("V", 1:p)
    X <- as.data.frame(X) %>%
      mutate(V1 = ifelse(V1 < median(V1), 0, 1)) #converting the continuous X2 covariate into binary
    e=rnorm(n,0, sd=sqrt(var_y)) #measure of uncertainty (error) of current data
    Y=as.matrix(X)%*%b+e #response vector of current data
    
    cov_matrix0 <- matrix(0, p-1, p-1)
    diag(cov_matrix0) <- var_x0 #variance-covariance matrix of historical data covariates
    X0=rmvnorm(n0, mean=mu0, sigma=cov_matrix0) #design matrix of historical data (without intercept)
    X0= cbind(rep(0, n0), X0) #treatment (only control arm) variable is included
    colnames(X0) =paste0("V", 1:p)
    e0 <- rnorm(n0, 0, sd=sqrt(var_y0))  #measure of uncertainty (error) of historical data
    Y0=as.matrix(X0)%*%b+e0 #response vector of historical data
    
    study=c(rep(1,n), rep(0,n0)) #outcome for propensity score model; 1 for current data, 0 for historical data
    fulldata=data.frame(cbind(rbind(X,X0),study)) #combined data (both current and historical)
    colnames(fulldata)=c(paste0("X", 1:p),"study")
    
    p.score <- fitted(glm(study ~ .-X1, data = fulldata, family = binomial)) 
    ps <- matchit(study ~ .-X1, data = fulldata, distance = p.score, method="subclass", subclass=K) #stratification; ensures all the strata have equal number of current subjects
    data.ps1=match.data(ps)
    t=table(data.ps1$study, data.ps1$subclass) #checking the samples within each stratum
    data.ps1$Y=rbind(Y,Y0)
    
### Trimming
    summary(p.score);summary(p.score[1:n]);summary(p.score[(n+1):(n+n0)]) #checking the summary of propensity scores
    lower_bound <- min(p.score[1:n]) #lower bound of the current study propensity score
    upper_bound <- max(p.score[1:n]) #upper bound of the current study propensity score
    data.ps <- data.ps1 %>% filter(distance >= lower_bound & distance <= upper_bound) #excluding extreme historical subjects 
    t=table(data.ps$study, data.ps$subclass) #checking the samples within each stratum
    
    total_n0 <- sum(t[1, ]) #total number of historical controls after stratification
    total_n0_vec[r] <- total_n0
    
### Computing overlapping coefficient for mean component
    mdist_list <- list()
    r_k01_list <- vector("list", K)  # r_k01 (overlapping coefficient for mean component) is a list
    
    for (k in 1:K) {
      subset_data <- data.ps[data.ps$subclass == k, 1:(p+1)]     # Subset data for kth strata
      mean_diff <- numeric(length = ncol(subset_data) - 1) # empty vector
      mean_diff_flag <- FALSE  # Track if mean_diff contains 1
      
      for (i in 1:(ncol(subset_data) - 1)) { # Exclude study column
        cur.stdy <- subset_data[subset_data$study == 1, i] # for current study
        ext.stdy <- subset_data[subset_data$study == 0, i] # for historical study
        
        if (length(ext.stdy) == 0) {
          mean_diff[i] <- 1 #if there is no subs in the historical study, the mean difference is set to 1, meaning the covariates in the current data 
          #and historical data are completely different. If you use mean difference = 0 that means there is no difference
          mean_diff_flag <- TRUE  # Mark that at least one covariate has mean_diff = 1
        } else {
          mean_diff[i] <- mean(cur.stdy) - mean(ext.stdy)
        }
      }
      
      subset_data_cov <- data.ps[data.ps$subclass == k, 1:p]
      cov_mat <- cov(subset_data_cov) # covariance matrix of the polled data for the kth strata
      
### Compute Mahalanobis distance and ensure it's stored as a scalar
      mdist <- as.numeric(sqrt(mean_diff %*% solve(cov_mat) %*% mean_diff))
      mdist_list[[k]] <- mdist  # Keep the sequence of strata intact
      
### If any mean_diff was set to 1, r_k01 = 0 for this k as there is no similarity
      if (mean_diff_flag == TRUE) {
        r_k01_list[[k]] <- 0
      }
    }
    
    total_mdist <- sum(unlist(mdist_list)) # Compute total Mahalanobis distance for normalization
    
### Normalize and compute similarity scores only for strata without mean_diff = 1
    for (k in seq_len(K)) {
      if (is.null(r_k01_list[[k]])) {  # Only compute if not pre-set to 0
        r_k01_list[[k]] <- 1 - (mdist_list[[k]] / total_mdist)
      }
    }
    
    r_k01 <- r_k01_list  # Ensure output remains a list
    
### Computing overlapping coefficient for variance component
    r_k02_list=list()
    for (k in 1:K){
      ps.cur=subset(data.ps, study==1 & subclass==k)[,"distance"] #propensity score of current data for the kth strata
      ps.ext=subset(data.ps, study==0 & subclass==k)[,"distance"] #propensity score of historical data for the kth strata
      if (length(ps.ext) < 2) { # to calculate the variance we need at least two observation in the historical data
        bk = 0
      } else {
        bk=var(ps.ext)/(var(ps.cur)+var(ps.ext))
      }
      r.val=1-2*abs(bk-0.5)
      r_k02_list=c(r_k02_list, r.val)
    }
    r_k02=r_k02_list #overlapping coefficient for variance component
    
### Calculating power parameter
    N_0k <- as.numeric(t[1, ]) #historical data within each stratum- only control
    lam_part01=A*unlist(r_k01)/sum(unlist(r_k01)) #first part of equation 9 in the main paper (for mean)
    lam_part02=A*unlist(r_k02)/sum(unlist(r_k02)) #first part of equation 9 in the main paper (for variance)
    
    ak01=list()
    ak02=list()
    for (k in 1:K) {
      lam01 = min(lam_part01[k], N_0k[k])
      power_par01 = ifelse(N_0k[k] == 0, 0, lam01 / N_0k[k]) #if lam01 = N_0k = 0, lam01/N_0k=undefined, so it's replaced by 0 (equation 9)
      ak01 = c(ak01, power_par01)
      
      lam02 = min(lam_part02[k], N_0k[k])
      power_par02 = ifelse(N_0k[k] == 0, 0, lam02 / N_0k[k])
      ak02 = c(ak02, power_par02)
    }
    
    sigma_sq=2.1 #initial vale of sigma_sq
    beta=rep(1,p) #initial value of beta
        
    MC.sample <- vector("list", length = K) #empty matrix
    
    for (k in 1:K) {
      #current data
      df = subset(data.ps, study==1 & subclass==k)
      Xk=as.matrix(df[,1:p])
      Yk=df[, "Y"]
      
### Historical data
      df0 = subset(data.ps, study==0 & subclass==k)
      Xk0=as.matrix(df0[,1:p])
      Yk0=df0[, "Y"]
      
      hatbeta=solve(t(Xk)%*%Xk)%*%t(Xk)%*%Yk #estimates using current data only
      if (nrow(Xk0) > p) { #to calculate hatbeta0, we need at least p number of obs in the historical data
        hatbeta0 = ginv(t(Xk0) %*% Xk0) %*% t(Xk0) %*% Yk0 #estimates using historical data only
      } else {
        hatbeta0=matrix(0, nrow=p) #This is just for computational benefit, will not be used later, dimension is p*1
      }
      
### Empty vectors and matrices
      posterior_sample=matrix(NA, iter, (p+1))
      
### MCMC running samples
      for(i in 1:iter){
        if (nrow(Xk0) > p) {
          beta.var=sigma_sq*solve(t(Xk)%*%Xk+ak01[[k]]*t(Xk0)%*%Xk0)
          beta.mean=solve(t(Xk)%*%Xk+ak01[[k]]*t(Xk0)%*%Xk0)%*%(t(Xk)%*%Xk%*%hatbeta+ak01[[k]]*t(Xk0)%*%Xk0%*%hatbeta0)
          beta=as.vector(rmvnorm(1, mean=beta.mean, sigma=beta.var))
          
          sigma_sq_shape=(nrow(df)+p*ak01[[k]]+ak02[[k]]*(nrow(df0)-p))/2
          sigma_sq_rate=(t(Yk-Xk%*%hatbeta)%*%(Yk-Xk%*%hatbeta)+ak02[[k]]*t(Yk0-Xk0%*%hatbeta0)%*%(Yk0-Xk0%*%hatbeta0))/2
          sigma_sq=rinvgamma(1, shape=sigma_sq_shape, rate=sigma_sq_rate)
        } else {
          beta.var=sigma_sq*solve(t(Xk)%*%Xk)
          beta.mean=solve(t(Xk)%*%Xk)%*%(t(Xk)%*%Xk%*%hatbeta)
          beta=as.vector(rmvnorm(1, mean=beta.mean, sigma=beta.var))
          
          sigma_sq_shape=(nrow(df))/2 
          sigma_sq_rate=(t(Yk-Xk%*%hatbeta)%*%(Yk-Xk%*%hatbeta))/2
          sigma_sq=rinvgamma(1, shape=sigma_sq_shape, rate=sigma_sq_rate)
        }
        posterior_sample[i,] = c(beta, sigma_sq)
      } 
      MC.sample[[k]] <- posterior_sample
    }
    
### Unweighted overall estimate of treatment effect
    true.trt=b[1] #based on the current study only
    
    uw.MC.sample <- matrix(NA, iter, K)
    for(k in 1:K){
      uw.new.MC.sample=MC.sample[[k]][,1]
      uw.MC.sample[,k]=uw.new.MC.sample
    }
    uw.trt.effect.vec=apply(uw.MC.sample[-(1:B),], 1, mean) #ignoring first B as burn in
    uw.trt.effect=mean(uw.trt.effect.vec) #overall treatment effect
    uw.trt.effect.rep[r] <- uw.trt.effect
    uw.quantiles.trt=quantile(uw.trt.effect.vec, probs=c(0.025, 0.975)) #95% credible intervals
    uw.width.CI.rep[r] <- diff(uw.quantiles.trt) #width of 95% CI
    #uw.width.CI.rep=c(uw.width.CI.rep, uw.width.CI)
    
    if (true.trt >= uw.quantiles.trt[1] && true.trt <= uw.quantiles.trt[2]) {
      uw.CP.trt[r] <- 1
    } else {
      uw.CP.trt[r] <- 0
    }
    cat("Running iteration",r,"\n")
  }
  
### Output lists
  list(
    mean_TE = mean(uw.trt.effect.rep),
    bias = mean(true.trt - uw.trt.effect.rep),
    abs_bias = mean(abs(true.trt - uw.trt.effect.rep)),
    RMSE = sqrt(mean((true.trt - uw.trt.effect.rep)^2)),
    SE = sd(uw.trt.effect.rep),
    CP = mean(uw.CP.trt),
    mean_CI_width = mean(uw.width.CI.rep),
    used_hist_sub= mean(total_n0_vec)
  )
}

### Running the developed function to run
results <- PS_BBPPP(rep = 100,n = 100, p = 4, n0 = 1000, b = c(2, 1, 1.5, -1.3), A = 30,var_y = 1,var_y0 = 1,         
                    mu = c(1,1.2,1.5,1.6),mu0 = c(1,1,1),var_x = 1,var_x0 = 1,K = 5,iter = 5000,B = 2000, seed = 1234)
