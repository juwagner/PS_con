#####-----------------------------------------------------------
##### preamble

rm(list=ls())
library(gaussquad)
library(quadprog)
library(Matrix)
source("./functions.R")

#####-----------------------------------------------------------
##### simulated data

set.seed(123)
n <- 100                          
x <- sort( runif(n) )
fx <- sin(pi*x)                                                  # arbitrary test function
D <- 10                                                          # number of random intercepts
area <- sample(1:10, n, replace=T)                               # intercept indicator
y <- fx + rnorm(n, sd=0.1) + rnorm(length(area), sd=0.5)[area]   # test function + random error + area specific intercept

plot(x,y)
lines(x, fx, col="red", lwd=2)

#####-----------------------------------------------------------
##### spline matrices

Omega <- c(0,1)                           # or c(min(x), max(x))
m <- 40                                   # number of spline knots
q <- 3                                    # spline degree
K <- m+q+1                                # number of basis functions
Phi <- bspline_matrix(x, m, q, Omega )    # B-spline basis matrix at the covariates
l <- 2                                    # degree of penalty
#Delta <- diff( diag(K), diff=l)          # difference penalty
Delta <- L2_norm_matrix(m, q, l, Omega)   # curvature penalty

x_grid <- seq(Omega[1],Omega[2],length=200)       # grid values for visualization and shape constraints
Phi_grid <- bspline_matrix(x_grid, m, q, Omega )  # B-spline basis matrix at the gridded values
lambda <- 0.000001                                 # weight of the penalty term (manually choosen)

#####-----------------------------------------------------------
##### common P-spline

A <- crossprod(Phi) + lambda*crossprod(Delta)    # coefficient matrix
b <- crossprod(Phi,y)                            # right-hand site
alpha <- solve(A,b)                              # coefficient vector as solution of the linear system

lines(x_grid, Phi_grid%*%alpha, col="green", lwd=2)

#####-----------------------------------------------------------
##### P-spline with shape constraints

C1 <- shape_constraints_matrix(x_grid,m,q,r=0,Omega)      # nonnegativity constraint
C2 <- -shape_constraints_matrix(x_grid,m,q,r=2,Omega)     # concavity constraint
C <- rbind(C1,C2)                                         

A <- crossprod(Phi) + lambda*crossprod(Delta)             # matrix for the quadratic program 
b <- crossprod(Phi,y)                                     # vector for the quadratic program
sol <- solve.QP(A, b, t(C))                               # solve the quadratic program
alpha <- sol$solution                                 

lines(x_grid, Phi_grid%*%alpha, col="blue", lwd=2)


#####-----------------------------------------------------------
##### P-spline with shape constraints and random intercept

W <- diag(D)[area,]                                     # intercept indicator matrix
B <- cbind(Phi, W)                                      # extension of the basis matrix
P <- bdiag( crossprod(Delta) , diag(D))                 # extension of the penaly matrix
C_ext <- cbind(C, matrix(0, nrow=dim(C)[1], ncol=D))    # extension of the shape constraint matrix

A <- crossprod(B) + lambda*P            # matrix for the quadratic program 
b <- crossprod(y,B)                     # matrix for the quadratic program 
sol <- solve.QP(A, b, t(C_ext))         # solve the quadratic program
alpha <- sol$solution[1:K]
u <- sol$solution[(K+1):(K+D)]

lines(x_grid,Phi_grid%*%alpha, col="purple", lwd=2)


# legend("topright", 
#        legend=c("true function", "P-spline", "shape constrained P-spline", "shape constrained P-spline with intercept"),
#        col=c("red", "green", "blue", "purple"),
#        lwd=rep(2,4)
# )



