# Recreates Figures 8 and 10 in Cotsakis et al. (2022)

# PARAMETERS
sigma = c(2,0.5) # The singular values of the anisotropy matrix.
us = c(0, 0.25, 0.5) # Threshold level (Standard Gaussian quantile)
sc = 1/30 # Domain has height and width 1/sc
nu = 2.5 # Smoothness parameter > 0 of the Matern covariance function. Set to Inf for a Gaussian covariance function.
N_sim = 500 # Number of simulations per value of n. 500 in the original paper
N = 1024 # the height and width of the image measured in pixels
m = 8
theta = pi/4

# Generate the random fields and compute the excursion perimeters
source("excursions.R")
library(RandomFields)
U = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol=2) # unitary matrix
A = diag(sigma) %*% t(U)
model <- RPgauss(RMmatern(nu = nu, scale = sc*N, var = 1, Aniso=A))

estimates = c()
for (i in 1:N_sim){
  f <- RFsimulate(model, x=1:N, y=1:N)$variable1
  f = matrix(f,nrow=N)
  new_row = c()
  for (u in us){
    bf = f>u
    new_row = c(new_row, perimBF(bf, box_size = 8)/(N-1))
  }
  estimates = rbind(estimates, new_row)
  
  if(i %% 10 == 0){
    print(paste(i, "/", N_sim, "replications completed for n =", m))
    imageBF(bf)
    save(estimates, file = "figure8.RData")
  }
}
estimates = data.frame(estimates)