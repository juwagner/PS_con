

#####-------------------------------------------------
#####  B-spline Basis Matrix (equidistant knots)
bspline_matrix = function(x, m, q, Omega){
  
  trunc_pol = function(x,t,q) (x-t)^q * (x>t)
  h <- (Omega[2]-Omega[1]) / (m+1)
  knots <- seq(Omega[1]-q*h, Omega[2]+q*h, by = h)
  H <- outer(x, knots, trunc_pol, q)
  D <- diff(diag(dim(H)[2]), diff = q+1) / (gamma(q+1)*h^q)
  Phi <- (-1)^(q+1) * tcrossprod(H,D)
  Phi[Phi<1e-10] <- 0
  return(Phi)
  
}

#####-------------------------------------------------
#####  Gramian Matrix of B-splines with respect to L2
L2_norm_matrix = function(m, q, l, Omega){
  
  # Parameter
  K <- m+q+1
  h <- (Omega[2]-Omega[1]) / (m+1)
  
  # Gauss-Quadratur
  n_gauss <- q+1
  gauss_legendre <- legendre.quadrature.rules(n_gauss)[[n_gauss]]
  w_ref <- gauss_legendre$w[order(gauss_legendre$x)]
  points_ref <- (h/2)*sort(gauss_legendre$x) + (Omega[1]+(h/2))
  points <- as.vector( sapply( 1:(m+1), function(j) points_ref+(j-1)*h ) )
  w <- rep(w_ref,(m+1))
  Phi <- bspline_matrix(points, m, (q-l), Omega)
  G <- (h/2)*crossprod(Phi, Phi*w)
  
  # Gramian Matrix
  if(l==0){
    Lambda <- G
  } else{
    Delta <- diff(diag(K),diff=l)
    Lambda <- (1 / (h^(2*l))) * crossprod(Delta,G%*%Delta)
  }
  return(Lambda)
  
}


#####-------------------------------------------------
##### Shape-Constraints Matrix
shape_constraints_matrix = function(x,m,q,r,Omega){
  
  # Parameter
  K <- m+q+1
  h <- (Omega[2]-Omega[1]) / (m+1)
  
  # Matrix
  Phi <- bspline_matrix(x, m, (q-r), Omega)
  if(r==0){
    Gamma <- Phi
  } else{
    C <- (1/(h^r))*diff(diag(K),diff=r)
    Gamma <- Phi%*%C
  }
  return(Gamma)  
  
}