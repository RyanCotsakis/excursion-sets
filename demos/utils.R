# Install necessary packages
if (!is.element("RandomFieldsUtils", installed.packages()[, "Package"])){
  if (tolower(readline(
    "'RandomFieldsUtils' and 'RandomFields' are not installed. Would you like to install them? [y/n]: "
    )) == "y"){
    install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_0.5.3.tar.gz", repos = NULL, type = "source")
    install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.8.tar.gz", repos = NULL, type = "source")
  }
} else if (!is.element("RandomFields", installed.packages()[, "Package"])){
  if (tolower(readline(
    "'RandomFields' is not installed. Would you like to install it? [y/n]: "
    )) == "y"){
    install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.8.tar.gz", repos = NULL, type = "source")
  }
}

# compute the perimeter of an ellipse with major and minor axis (radii) a and b.
ellipse = function(a, b){
  r = ((a-b)/(a+b))^2
  C = pi*(a+b)*(3*r/(sqrt(-3*r+4)+10)+1)
  return(C)
}

# the theoretical density of half perimeter length of excursion set at u
# for the affine random field X(s) = Y(As), where Y ia a stationary, isotropic,
# standard Gaussian random field with Matern covariance function with smoothness
# parameter nu.
c1star = function(A, u, nu){
  if (is.matrix(A)) sigma = sqrt(eigen(t(A)%*%A)$values)
  else if(length(A) == 1) sigma = c(A,A)
  else sigma = A[1:2]
  
  lambda2 = nu/(nu-1)
  return(sqrt(lambda2)/(8*pi)*exp(-u^2/2)*ellipse(sigma[1],sigma[2]))
}
