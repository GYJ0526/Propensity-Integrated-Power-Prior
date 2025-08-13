# README Document for PS_BBPPP

## The R script for PS_BBPPP (Propensity Score–Based Borrowing-by-Parts Power Prior) can be stored in a folder named PS_BBPPP. The method can be run by sourcing the script and setting the required input parameters in the R environment:
PS_BBPPP(rep = , n = , p = , n0 = , b = , A = , var_y = , var_y0 = , mu = , mu0 = , var_x =, var_x0 =, K =, iter = , B = , seed =)

# Inputs for PS_BBPPP
	rep: Number of simulation replications to run. Required.
	n: Sample size of the current study. Required.
	p: Number of covariates (excluding the intercept). Required.
	n0: Sample size of the historical dataset. Required.
	b: Vector of true parameter values for data generation. The first element corresponds to the true treatment effect (θ). Required.
	A: Nominal number of historical subjects to be borrowed (for mean and variance components). Required.
	var_y, var_y0: Variances of the outcome in the current and historical datasets, respectively. Required.
	mu: Mean vector of covariates for the current study (length = p). Required.
	mu0: Mean vector of covariates for the historical study (length = p − 1, excluding treatment variable). Required.
	var_x: Variance for current study covariates (dimension = p*p). Required.
	var_x0: Variance for historical study covariates (dimension = p-1*p-1). Required.
	K: Number of propensity score strata. Required.
	iter: Number of MCMC iterations per stratum. Required.
	B: Number of burn-in samples to discard from each MCMC chain. Required.
	seed: Indicates the initial seed to generate sets of seed numbers in generating data and running MCMC. Optional.

# Outputs for PS_BBPPP
This function automatically produces the following outputs: 
	mean_TE: Mean treatment effect estimate across r repetition.
	bias: Mean bias of the treatment effect estimates.
	abs_bias: Mean absolute bias of the treatment effect estimates.
	RMSE: Mean root mean squared error (RMSE) of the treatment effect estimates.
	SE: Mean simulation error (SE) of the treatment effect estimates.
	CP: Mean coverage probability (CP) of the 95% credible interval.
	mean_CI_width: Mean width of the 95% credible interval.
	used_hist_sub: Mean number of historical subjects retained after trimming.

# Running Time for PS_BBPPP
It takes approximately 14 seconds to complete a single repetition (rep=1) with 5,000 MCMC iterations when n=100, n_0=1000, p=4, and A=30.
