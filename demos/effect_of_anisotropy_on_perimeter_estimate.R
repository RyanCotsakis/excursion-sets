# Recreates Figure 6 in Cotsakis et al. (2022)

# Generate the random fields and compute the excursion perimeters
source("excursions.R")
library(RandomFields)

N = 256 # the height and width of the image measured in pixels
u = 0.5 # Threshold level (Standard Gaussian quantile)
sc = 0.2 # Domain has height and width 1/sc
sigma = c(2,0.5) # The singular values of the anisotropy matrix.
nu = 2.5 # Smoothness parameter > 0 of the Matern covariance function. Set to Inf for a Gaussian covariance function.
N_sim = 500

estimates = c()
for (theta in seq(0,pi/2+0.001,pi/32)){
  U = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol=2) # unitary matrix
  A = diag(sigma) %*% t(U)
  model <- RPgauss(RMmatern(nu = nu, scale = sc*N, var = 1, Aniso=A))
  
  for (i in 1:N_sim){
    f <- RFsimulate(model, x=1:N, y=1:N)$variable1
    f = matrix(f,nrow=N)
    bf = f>u
    estimates = rbind(estimates, c("Phat1"=perimBF(bf, p=1)/(N-1),
                                   "Phat2"=perimBF(bf, box_size = 8)/(N-1),
                                   "truth"=perimCF(f,u)/(N-1),
                                   "theta"=theta))
    if(i %% 10 == 0){
      print(paste(i, "/", N_sim, "replications completed for theta =", theta))
      imageBF(bf)
    }
  }
  save(estimates, file="figure6.RData")
}
estimates = data.frame(estimates)

source("demos/utils.R")
expect = c1star(diag(sigma), u, nu)*2/sc^2



#' REFERENCES
#'
#' Ryan Cotsakis, Elena Di Bernardino, Thomas Opitz. On the perimeter estimation
#' of pixelated excursion sets of 2D anisotropic random fields. 2022. hal-03582844v2