# Recreates Figure 6 in Cotsakis et al. (2022)

# setwd("path/to/excursion-sets")

source("excursions.R")
source("demos/utils.R")
library(RandomFields)

N = 256 # the height and width of the image measured in pixels
u = 0.5 # Threshold level (Standard Gaussian quantile)
sc = 0.2 # Domain has height and width 1/sc
sigma = c(2,0.5) # The singular values of the anisotropy matrix.
nu = 2.5 # Smoothness parameter > 0 of the Matern covariance function. Set to Inf for a Gaussian covariance function.
N_sim = 200
thetas = seq(0,pi/2+0.001,pi/32)

estimates = c()
for (theta in thetas){
  U = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol=2) # unitary matrix
  A = diag(sigma) %*% t(U)
  model <- RPgauss(RMmatern(nu = nu, scale = sc*(N-1), var = 1, Aniso=A))
  
  for (i in 1:N_sim){
    f <- RFsimulate(model, x=1:N, y=1:N)$variable1
    f = matrix(f,nrow=N)
    bf = f>u
    estimates = rbind(estimates, c("Phat1"=perimBF(bf, p=1)/(N-1)/sc,
                                   "Phat2"=perimBF(bf, box_size = 11)/(N-1)/sc,
                                   "truth"=perimCF(f,u)/(N-1)/sc,
                                   "theta"=theta))
    if(i %% 10 == 0){
      print(paste(i, "/", N_sim, "replications completed for theta =", theta))
      imageBF(bf)
    }
  }
}
estimates = data.frame(estimates)

# Compute mean estimates
df = data.frame()
for(theta in thetas){
  inds = which(estimates$theta == theta)
  if(length(inds) == 0) break
  new_row = colMeans(data.matrix(estimates[inds,]))
  df = rbind(df, new_row)
}
names(df) = names(estimates)

# plotting
library(ggplot2)
expect = c1star(diag(sigma), u, nu)*2/sc^2
labels = c("0",
           expression(frac(pi,8)),
           expression(frac(pi,4)),
           expression(frac(3*pi,8)),
           expression(frac(pi,2)))
p = ggplot(df, aes(x=thetas)) +
  scale_x_continuous(breaks = thetas[seq(1,17,4)], labels = labels) +
  geom_abline(slope=0,intercept = 0,color="red",alpha=0.9) +
  geom_line(aes(y=Phat1*pi/4 - truth), size=1, linetype=3, color=3) +
  geom_point(aes(y=Phat1*pi/4 - truth), size=2, shape=15, color=3) +
  geom_line(aes(y=Phat2 - truth), size=1, linetype=3, color=4) +
  geom_point(aes(y=Phat2 - truth), size=2, shape=16, color=4) +
  theme(legend.position = "none",
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major= element_line(colour = "lightgrey"),
        panel.grid.minor=element_line(colour = "lightgrey")) +
  ylab("Average signed error") + xlab(expression("Angle of anisotropy, "* theta))
p


#' REFERENCES
#'
#' Ryan Cotsakis, Elena Di Bernardino, Thomas Opitz. On the perimeter estimation
#' of pixelated excursion sets of 2D anisotropic random fields. 2022. hal-03582844v2