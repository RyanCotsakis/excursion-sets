library(Rcpp)
sourceCpp("binary_fields.cpp")


#' Plot a binary field as a black and white image
#'
#' @param bf: a binary field represented by a matrix of 0's and 1's
#' @return A black and white image of the binary field.
imageBF = function(bf){
  image(t(bf), col=gray.colors(2,start=1,end=0.5))
}


#' Computes the area of a binary field.
#'
#' @param bf: The binary field as a matrix of 1s and 0s.
#' @return The number of ones in the field.
areaBF <- function(bf){
  bf <- as.matrix(bf)*1
  return(sum(bf))
}


#' Gives a summary of the connected components of a binary field
#' 
#' @param bf A matrix of 1s and 0s representing the binary field.
#' @param connectivity An integer indicating the connectivity of the components.
#'   - If 8, components connected by corners are considered connected.
#'   - If 6, these situations are considered connected with probability 1/2 (stochastic output).
#'   - If 4 or anything else, components are not connected by corners.
#' @param incl_perim A logical value indicating whether to compute the perimeter of each component.
#'   The perimeter will not include the boundary in any case.
#' @param incl_area A logical value indicating whether to compute the area of each component.
#' @return A list containing the following variables:
#'   - ncc: The number of connected components.
#'   - nholes: The number of holes.
#'   - euler: The Euler characteristic (ncc - nholes).
#'   - labeled: A matrix like bf that has its ones changed to unique labels
#'     corresponding to connected components.
#'   - holes_labeled: A matrix like bf that has its ones changed to unique labels
#'     corresponding to holes.
#'   - hole_labels: A numeric vector of the labels used in holes_labeled.
#'   - table: A data.frame with one column corresponding to the labels in labeled,
#'     and the other columns corresponding to the area if incl_area and the perimeter
#'     of that component if incl_perim.
ccBF <- function(bf, connectivity = 6, incl_perim = TRUE, incl_area = TRUE){
  bf = as.matrix(bf)*1
  cc = cpp_connected_components(bf, connectivity)
  labels = setdiff(cc$cc, 0)
  ncc = length(labels)
  hole_labels = setdiff(cc$ch, 0)
  nholes = length(hole_labels)
  euler = ncc - nholes
  m = ncol(bf)
  n = nrow(bf)
  
  tab = data.frame("label" = labels)
  if (incl_area){
    area_from_label = function(label) return(areaBF(cc$cc == label))
    areas = sapply(labels, area_from_label)
    tab = cbind(tab, "area" = areas)
  }
  if (incl_perim){
    perimeter_from_label = function(label) return(perimBF(cc$cc == label, p = 1))
    perims = sapply(labels, perimeter_from_label)
    tab = cbind(tab, "perim" = perims)
  }
  on_border = labels %in% c(cc$cc[1,], cc$cc[,1], cc$cc[nrow(cc$cc),], cc$cc[,ncol(cc$cc)])
  tab = cbind(tab, "border" = on_border)
  summary = list("ncc" = ncc,
                 "nholes" = nholes,
                 "euler" = euler,
                 "labeled" = cc$cc,
                 "holes_labeled" = cc$ch,
                 "table" = tab)
  return(summary)
}


#' Computes the full perimeter of a binary field.
#' 
#' @param bf: The binary field as a matrix of 1s and 0s.
#' @param p: (Optional) The p-norm used to estimate perimeter. Default is 2.
#' @param incl_boundary: (Optional) If TRUE, include the perimeter of the excursions that is
#' on the boundary of the matrix.
#' @param box_size: (Optional) Use "auto" or a positive integer to choose the side length of
#' the boxes used in our estimator measured in pixels. Default is "auto".
#' 
#' @return The full perimeter estimate.
#'
perimBF <- function(bf, p = 2, incl_boundary = FALSE, box_size = "auto"){
  bf = as.matrix(bf)*1
  if (box_size == 0 | p == 1){
    per = cpp_perimeter_bierme(bf)
  } else if (box_size == "auto"){
    cc = ccBF(bf, connectivity = 8, incl_perim = FALSE, incl_area = FALSE)
    if (cc$ncc == 0) return(0)
    bs = max(floor((length(bf)/(cc$nholes + cc$ncc))^(1/3)/3), 1)
    per = cpp_perimeter_cotsakis(bf, bs)
  } else {
    per = cpp_perimeter_cotsakis(bf, box_size)
  }
  if (!incl_boundary){
    return(per)
  }
  border_count = sum(c(bf[1,], bf[,1], bf[nrow(bf),], bf[,ncol(bf)]))
  return(per+border_count)
}


#' Computes the full perimeter of the excursion set of a continuous random field
#' at a given threshold.
#' 
#' @param f: The random field as a matrix.
#' @param u: The threshold value.
#' @param incl_boundary: Logical indicating whether or not to include the perimeter
#' on the boundary of the matrix.
#' 
#' @return The full perimeter estimate.
perimCF <- function(f, u, incl_boundary = FALSE){
  f = as.matrix(f)
  per = cpp_perimeter_unbiased(f,u)
  if (!incl_boundary){
    return(per)
  }
  bf = as.matrix(f > u)
  border_count = sum(c(bf[1,], bf[,1], bf[nrow(bf),], bf[,ncol(bf)]))
  return(per+border_count)
}


#' Produce a dilated version of a binary image. The boundary of the input is
#' assumed to be padded with 0's.
#'
#' @param M: The binary field as a matrix of 1s and 0s.
#' @param r: The dilation distance (inclusive). Negative dilations correspond
#' to erosion.
#' @return: The dilated image as a binary field.
dilateBF <- function(M, r){
  M = as.matrix(M)*1
  # get translations
  translations = c()
  for (i in -ceiling(abs(r)):ceiling(abs(r))){
    for (j in -ceiling(abs(r)):ceiling(abs(r))){
      if(norm(c(i,j),type = "2") <= abs(r) & (i | j)){
        translations = rbind(translations, c(i,j))
      }
    }
  }
  # decide if it's a dilation or erosion
  if (r == 0) return(M)
  if (r > 0) return(cpp_dilation(M,translations,1))
  return(1-cpp_dilation(1-M,translations,1))
}


#' Computes the estimator for the largest radius such that the closing operation of a binary
#' set is equivalent to the identity operation (hat{rnv}).
#' The estimate provides an upper bound for the true value.
#'
#' @param M: the binary field. A matrix of 1s and 0s
#' @param r0: the initial guess for the estimator. The default is 10.
#' @param tol: the tolerance for the estimator. The default is 0.01.
#' @return the estimator for the largest radius.
rconvBF <- function(M, r0=10, tol=0.01){
  if(r0 < 5){
    r0 = 5
  }
  max_guess = 8*r0
  epsilon = sqrt(2)/2
  r = r0
  n = dim(M)[1]
  m = dim(M)[2]
  prev_min = r
  prev_max = r
  while(TRUE){
    # Update r
    padding = ceiling(r)+1
    n_padded = n + 2*padding
    m_padded = m + 2*padding
    M_padded = matrix(rep(0,n_padded*m_padded),nrow = n_padded)
    M_padded[(1:n)+padding,(1:m)+padding] = M
    M_closed = dilateBF(dilateBF(M_padded,r-epsilon),-(r+epsilon))
    too_big = sum(M_closed > M_padded) > 0
    if(too_big){
      prev_max = r
      if(r <= prev_min){ # still dividing by 2
        r = r/2
        prev_min = r
      } else {
        r = (r + prev_min)/2
      }
    } else {
      prev_min = r
      if(r >= prev_max){ # still multiplying by 2
        r = r*2
        prev_max = r
      } else {
        r = (r + prev_max)/2
      }
    }
    if(prev_max - prev_min < tol){
      return(prev_max)
    }
    if(r > max_guess){
      return(max_guess)
    }
  }
}


#' Computes the reach (hat{rch}) estimator for a square grid of samples of a set A
#' 
#' @param M: the binary field
#' @param max_guess: the estimate will not exceed this value
#' @return the reach estimate.
reachBF <- function(M, max_guess) {
  epsilon <- sqrt(1.25)
  Mup <- cpp_dilation(M, matrix(c(0,1), nrow=1), 0)
  Mdown <- cpp_dilation(M, matrix(c(0,-1), nrow=1), 0)
  Mleft <- cpp_dilation(M, matrix(c(-1,0), nrow=1), 0)
  Mright <- cpp_dilation(M, matrix(c(1,0), nrow=1), 0)
  M_boundary <- M * (Mup + Mdown + Mleft + Mright > 0)
  boundary_indices <- which(M_boundary == 1)
  guess <- Inf
  g <- function(x, alpha) {
    return(alpha^2/(8*x) + x/2)
  }
  for(i in 1:(length(boundary_indices)-1)) {
    coord_i <- get_coord(boundary_indices[i], nrow(M))
    my_neighbours <- c()
    for(j in (i+1):length(boundary_indices)) {
      coord_j <- get_coord(boundary_indices[j], nrow(M))
      if(norm(coord_i-coord_j, type="2") <= 2*max_guess) {
        coord_mid <- (coord_i + coord_j)/2
        flr_mid <- floor(coord_mid)
        if(M[flr_mid[1], flr_mid[2]]) {
          next
        }
        radius <- max(abs(coord_mid-coord_i))
        r_min <- max(floor(coord_mid[1] - radius), 1)
        r_max <- min(ceiling(coord_mid[1] + radius), nrow(M))
        c_min <- max(floor(coord_mid[2] - radius), 1)
        c_max <- min(ceiling(coord_mid[2] + radius), ncol(M))
        min_distance <- Inf
        for (k in r_min:r_max) {
          for (l in c_min:c_max) {
            if(M[k,l] && norm(coord_mid - c(k,l), type="2") < min_distance) {
              min_distance <- norm(coord_mid - c(k,l), type="2")
            }
          }
        }
        if (min_distance > epsilon) {
          guess <- min(guess, g(min_distance-epsilon, norm(coord_i-coord_j, type="2")))
        }
      }
    }
  }
  return(min(guess, max_guess))
}
