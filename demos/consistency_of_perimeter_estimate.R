# Recreates Figures 7 and 9 in Cotsakis et al. (2022)

# setwd("path/to/excursion-sets")

source("excursions.R")
source("demos/utils.R")
library(RandomFields)

# PARAMETERS
anisotropy = 2 # The amount of anisotropy. Set to 1 for an isotropic model
u = 0.5 # Threshold level (Standard Gaussian quantile)
sc = 0.2 # Domain has height and width 1/sc
nu = 2.5 # Smoothness parameter > 0 of the Matern covariance function. Set to Inf for a Gaussian covariance function.
N_sim = 50 # Number of simulations per value of n. 500 in the original paper

A = diag(c(1/anisotropy, anisotropy))
ms = 2:16
estimates = c()
for (m in ms){
  N = floor(10*m^(3/2))
  model <- RPgauss(RMmatern(nu = nu, scale = sc*(N-1), var = 1, Aniso=A))
  
  for (i in 1:N_sim){
    f <- RFsimulate(model, x=1:N, y=1:N)$variable1
    f = matrix(f,nrow=N)
    bf = f>u
    estimates = rbind(estimates, c("Phat1"=perimBF(bf, p=1)/(N-1)/sc,
                                   "Phat2"=perimBF(bf, box_size = m)/(N-1)/sc,
                                   "truth"=perimCF(f,u)/(N-1)/sc,
                                   "m"=m))
    if(i %% 10 == 0){
      print(paste(i, "/", N_sim, "replications completed for n =", m))
      imageBF(bf)
    }
  }
}
estimates = data.frame(estimates)

# Compute mean absolute errors
df = data.frame()
for(m in ms){
  inds = which(estimates$m == m)
  if(length(inds) == 0) break
  
  errors1 = abs(pi/4*estimates$Phat1[inds] - estimates$truth[inds])
  df = rbind(df, c(mean(errors1), m, 1))
  errors2 = abs(estimates$Phat2[inds] - estimates$truth[inds])
  df = rbind(df, c(mean(errors2), m, 2))
}
names(df) = c("mae", "m", "p")

# Plot mean absolute errors
library(ggplot2)
p = ggplot(df, aes(x=m, color=as.factor(p))) +
  scale_x_continuous(breaks = ms[which(ms %% 2== 0)]) +
  geom_line(aes(y=mae), size=1,linetype=3) +
  geom_point(aes(y=mae, shape=as.factor(p)),size=3) +
  scale_shape_manual(values=c(15, 16))+
  scale_color_manual("",values = c(3, 4)) +
  theme(legend.position = "none",
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major= element_line(colour = "lightgrey"),
        panel.grid.minor=element_line(colour = "lightgrey")) +
  ylab("MAE") + xlab("n") + 
  geom_abline(slope=0,intercept = 0,color="red",alpha=0.4)
p


#' REFERENCES
#'
#' Ryan Cotsakis, Elena Di Bernardino, Thomas Opitz. On the perimeter estimation
#' of pixelated excursion sets of 2D anisotropic random fields. 2022. hal-03582844v2
