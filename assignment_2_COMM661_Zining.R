#####################################################################################################################################################
# Assignment 2 of COMM661
# Zining Wang
# Nov 5, 2019
#####################################################################################################################################################

#install.packages("MCMCpack")
#install.packages("matrixcalc")
require(data.table)
require(ggplot2)
require(MCMCpack)
require(matrixcalc) #  vec operator

# install.packages("monomvn")
# write down your functions for 
# multivariate normal
# gamma/inverse gamma
# wishart/inverse wishart
# (check the density fct)
# generate mean and variance to confirm
# y_i = x_i' \beta + epsilon_i
# \pi(\beta | y_i, \sigma^2)
# \pi(\sigma^2 | \beta, y_i)


# input function: 
# rnorm(n, 0, 1)
# runif(n, 0, 1)
#####

#####################################################################################################################################################
# multivariate distributon
#####################################################################################################################################################
mvrnorm_fct <-
  function(n, mu, Sigma, tol=1e-8)
  {
    d <- length(mu)
    if(!all(dim(Sigma) == c(d,d))) stop("incompatible arguments")
    em <- eigen(Sigma, symmetric = TRUE)
    ev <- em$values
    
    if(!all(ev >= tol)) stop("covariance matrix is not positive definite") # use a tolarance value slightly larger than 0 for computational purposes
    
    X <- matrix(rnorm(d * n), n) # gen vector/matrix of raondom numbers from the normal distirbution
    
    X <- mu + em$vectors %*% diag(sqrt(pmax(ev, 0)), d) %*% t(X) # another choice: Cholesky factorization 
    #X<- mu + chol(Sigma) %*% X
    # https://stats.stackexchange.com/questions/83850/confused-about-cholesky-and-eigen-decomposition
    if(n == 1) X else t(X)
  }
#####################################################################################################################################################
# data check
sim_data_mvrnorm = mvrnorm_fct(10000, c(0.3, 0.4), Sigma = diag(2))
sim_data_mvrnorm = as.data.frame(sim_data_mvrnorm)
colnames(sim_data_mvrnorm) = c("var_1: 0.3", "var_2: 0.4")

p_mvrnorm_1 <-ggplot(sim_data_mvrnorm, aes(x=`var_1: 0.3`)) + geom_histogram() + ggtitle("Multivariate: var_1")
p_mvrnorm_1

p_mvrnorm_2 <-ggplot(sim_data_mvrnorm, aes(x=`var_2: 0.4`)) + geom_histogram()+ ggtitle("Multivariate: var_2")
p_mvrnorm_2
#####################################################################################################################################################
# gamma and inverse gamma distribuion
#####################################################################################################################################################
# the inverse of CDF failed--
# try the acceptance rejection method
# Generate U, V and W as iid uniform (0, 1] variates.
# Gamma(n + \delta, 1)

set.seed =(1234)

gamma_fun = function(alpha, beta){
  
  k = alpha
  theta = 1/beta
  n=floor(k)
  delta = k - n
  
  r = runif(n,  min = 0, max = 1)
  z = -sum(log(r)) # z~Gamma(n, 1)
  
  eta = 1
  threshold = 0
  #ind = 0
  #  Ahrens-Dieter acceptance–rejection method: source: wkipedia
  while (eta > threshold){
    U = runif(1, min = 0, max = 1)
    V = runif(1, min = 0, max = 1)
    W  =runif(1, min = 0, max = 1)
    
    if (U <= exp(1)/ (exp(1) + delta)){
      xi = V^{1/delta}
      eta = W*xi^(delta-1)
      
    }else{
      xi = 1 - log(V)
      eta = W*exp(1)^(-xi)
      
    }
    threshold = xi^(delta-1)*exp(-xi)
  }
  #xi_vec[i] = xi
  x = theta*(xi+z)
  x_inv =1/x
  c(x, x_inv) # return x and x_inv
  
}

#####################################################################################################################################################
# data check

sim_data_gamma = matrix(0, nrow = 10000, ncol = 2)

for (i in 1:10000){
  sim_data_gamma[i,] = gamma_fun(5,1)
}

hist(sim_data_gamma[,1])
hist(sim_data_gamma[,2])
sim_data_gamma = as.data.frame(sim_data_gamma)
colnames(sim_data_gamma) = c("Gamma", "Inverse Gamma")
p_gamma = ggplot(sim_data_gamma, aes(x=Gamma)) + geom_histogram() + ggtitle("Histagram: Gamma Distribution with alpha=5, beta=1")
p_gamma
p_invgamma = ggplot(sim_data_gamma, aes(x=`Inverse Gamma`)) + geom_histogram() + ggtitle("Histagram: Inverse Gamma Distribution with alpha=5, beta=1")
p_invgamma
#####################################################################################################################################################
# wishart and inverse wishart distribuion
#####################################################################################################################################################
# Wishart distribution is a generalization to multiple dimensions of the gamma distribution
#  G is a p × n matrix, each column of which is independently drawn from a p-variate normal distribution
#  In Bayesian statistics, the Wishart distribution is the conjugate prior of the inverse covariance-matrix of a multivariate-normal random-vector.
# install.packages("matlib")
# require(matlib)
wishart_fun <- function(V, p, n){
  G = mvrnorm_fct(n, mu = rep(0, p), Sigma = V)
  S = G%*%t(G)
  # S_inv = G_inv%*%t(G_inv)
}

invwishart_fun <- function(V, p, n){
  G = mvrnorm_fct(n, mu = rep(0, p), Sigma = V)
  G_inv = mvrnorm_fct(n, mu = rep(0, p), Sigma = solve(V)) # solve is the base R fct to get matrix inverse
  S_inv = G_inv%*%t(G_inv)
}

#####################################################################################################################################################
# data check

sim_data_wishart = wishart_fun(V=diag(3), p=3, n=1)
sim_data_invwishart = invwishart_fun(V=diag(3), p=3, n=1)


#####################################################################################################################################################
# MCMC
#####################################################################################################################################################
# code the Gibbs updating steps, using conjugate priors

# 1 univariate: only one y
# Data Simulation

nsim = 1000 # of simulation
N = 1000 # of obs
betas_true = c(0.4, 0.5)

v_0 = 5 # hyper prior of the inverse gamma \sigma^2 ~ IG(V_0/2, \delta_0/2)
s_0 = 1 # hyper prior of the inverse gamma \
intcpt = 1
X = mvrnorm_fct(N, c(0, 0), Sigma = diag(2))
Y = intcpt + X %*% betas_true+ rnorm(n=N, mean=0, sd=sigma)

mean=c(1,1,1) # beta_0
sigma=diag(3) # B_0
alhpa=1 # inverse gamma
beta=1  # inverse gamma

univariate_reg <- function(X, Y, N, mean, sigma, alhpa, beta, nsim){ # N: number of iterartion
  
  nparams=ncol(X)+2 # IV + intcpt + sigma
  
  X=cbind(rep(1, N),X)
  
  sims=matrix(0,nrow=nsim,ncol=nparams)
  
  sims[1,]=1 # to start
  
  #Gibbs Sampler
  for(i in 2:nsim){
    sims[i,1]=gamma_fun(v_0/2, v_0*s_0^2/2)[2] # draw sigma: sigma^2|y, to draw from the IG, use uninformative priors
    
    mean_p = solve(t(X) %*% X /sims[i,1]+solve(sigma)) %*% (t(X)%*%Y/sims[i,1]+solve(sigma) %*% mean)
    sigma_p = sims[i,1] * solve(t(X) %*% X + solve(sigma/sims[i,1]) ) # A = sigma/sims[i,1]
    
    sims[i,2:nparams]=mvrnorm_fct(1, mean_p, sigma_p) # draw beta: from multivariate normal: Pi(beta|y, sigma^2) <- N_k (beta_hat, B_hat)
  }
  return(sims)
}


output = as.data.frame(univariate_reg(X, Y, N, mean, sigma, alhpa, beta,  nsim))
colnames(output) = c("sigma", "intcpt: 1", "beta_1: 0.4", "beta_2: 0.5")
output = as.data.table(output)
output[,draws:=seq_len(nrow(output))]
output_melt = melt(output, id.vars = "draws", 
                   measure.vars = c("sigma", "intcpt: 1", "beta_1: 0.4", "beta_2: 0.5"))
output_melt
#####################################################################################################################################################
# plot the posteriors

p_posterior <- ggplot(output_melt[variable!="sigma",], aes( x =draws , y = value ) ) +geom_point()  + facet_wrap(vars(variable), ncol=3) + ggtitle("Posteriors: Univariate") 
p_posterior 

p_posterior_sigma <- ggplot(output_melt[variable=="sigma",], aes( x =draws , y = value ) ) +geom_point() + ggtitle("Posteriors of sigma: Univariate") 
p_posterior_sigma

#####################################################################################################################################################

# 2.1 multivariate: multiple y's
# Y = XB+E (n*m, n*k, k*m, n*m)
# Sigma ~ IW(v_0, V_0) # Sigma is a m*m matrix
# beta|Sigma ~ N(beta_bar, Sigma \otimes inv(A))
# posteriors:
  # Sigma|Y, X ~ IW(v_0 + n, V_0 + S)
  # beta|Y, X, Sigma ~ N(beta_tilda, Sigma \otimes inv(X'X + A)), in which
  # beta_tilda = vec(B_tilda)
  # B_tilda = inv(X'X + A)(X'XB_hat + AB_bar)
  # S = (Y - XB_tilda)'(Y-XB_tilda) + (B_tilda - B_bar)'A(B_tilda - B_bar)

# data simulation

nsim = 1000
m = 2
N = 1000
betas_true = matrix(c(0.3, 0.4, 0.2, 0.6), nrow = 2, byrow = F) 
E = matrix(rnorm(n=N*m, mean=0, sd=1), nrow = N); E
intcpt = matrix(rep(1, N*m), nrow=N); intcpt
X = mvrnorm_fct(N, c(0, 0), Sigma = diag(2))
Y = intcpt + X %*% betas_true+ E ;Y
v_0 = 1
V_0 = diag(2)
A = matrix(rep(1, (m+1)^2), nrow = m+1)

multivariate_reg <- function(X, Y, N, v_0, V_0, A, betas_true, nsim){ # N: number of iterartion
  
  nparams = m*(ncol(X)+m+1) # m*(IV + intcpt + epsilon)
  sims = matrix(0, nrow=nsim, ncol=nparams)
  X=cbind(rep(1, N), X)
  
  betas = rbind(c(1, 1), betas_true) # add intercept to the beta vector
  B_tilda = solve(t(X) %*% X + A) %*% ( t(X) %*% Y + A %*% betas) # replace XB_hat with Y
  S = t(Y-X %*% B_tilda) %*% (Y-X %*% B_tilda) + t(B_tilda - betas) %*% A %*% (B_tilda - betas)
  
  #Gibbs Sampler
  
  for(i in 1:nsim){
    # wishart_fun
    Sigma_p = riwish( v_0 + N, V_0 + S) # update Sigma
    sims[i,1:m^2 ] = as.vector(Sigma_p) # flat by col and store the updated values

    beta_tilda = vec(B_tilda)
    Var = kronecker(Sigma_p, solve(t(X) %*% X + A))

    betas_p = mvrnorm_fct(1, beta_tilda, Var)
    sims[i, (m^2+1):nparams] = betas_p # intcpt1, beta_1.1, beta2.1, intcpt2, beta_1.2, beta2.2
  }
  return(sims)
}

output = as.data.frame(multivariate_reg(X, Y, N, v_0, V_0, A, betas_true, nsim))
colnames(output) = c("sigma_1: 1", "sigma_2: 0", "sigma_3: 0", "sigma_4: 1", "intcpt_1.1: 1", "beta_1.1: 0.3", "beta_2.1: 0.4","intcpt_1.2: 1", "beta_1.2: 0.2", "beta_2.2: 0.6")

output = as.data.table(output)
output[,draws:=seq_len(nrow(output))]
output_melt = melt(output, id.vars = "draws", 
                   measure.vars = c("sigma_1: 1", "sigma_2: 0", "sigma_3: 0", "sigma_4: 1", "intcpt_1.1: 1", "beta_1.1: 0.3", "beta_2.1: 0.4","intcpt_1.2: 1", "beta_1.2: 0.2", "beta_2.2: 0.6"))
output_melt
#####################################################################################################################################################
# plot the posteriors

p_posterior <- ggplot(output_melt, aes( x =draws , y = value ) ) +geom_point()  + facet_wrap(vars(variable), ncol=3) + ggtitle("Posteriors: Multivariate") 
p_posterior 



