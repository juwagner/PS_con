# PS_con: Shape-Constrained Penalized Splines
Penalized splines is a popular method for function estimation under the assumption of smoothness.
The provided code extends the well known P-spline method to constraints on its shape, such as monotonicity or convexity, and allows for the consideration of additional random effects within the model.
Further detail can be found in [Wagner et al., 2017](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssa.12295).

## Manual
For selected parameters, the B-spline basis matrix and the matrix of the penalty term are assembled, where the difference penalty or the curvature penalty (default) can be used.
```{r}
Phi <- bspline_matrix(x, m, q, Omega)     # B-spline basis matrix at the covariates
#Delta <- diff(diag(K), diff=l)           # difference penalty
Delta <- L2_norm_matrix(m, q, l, Omega)   # curvature penalty (default)
```


```{r}
A <- crossprod(Phi) + lambda*crossprod(Delta)    # coefficient matrix
b <- crossprod(Phi,y)                            # right-hand site
alpha <- solve(A,b)                              # coefficient vector as solution of the linear system
```

```{r}
C1 <- shape_constraints_matrix(x_grid,m,q,r=0,Omega)      # nonnegativity constraint
C2 <- -shape_constraints_matrix(x_grid,m,q,r=2,Omega)     # concavity constraint
C <- rbind(C1,C2)                                         
```

```{r}
A <- crossprod(Phi) + lambda*crossprod(Delta)             # matrix for the quadratic program 
b <- crossprod(Phi,y)                                     # vector for the quadratic program
sol <- solve.QP(A, b, t(C))                               # solve the quadratic program
alpha <- sol$solution   
```

```{r}
W <- diag(D)[area,]                                     # intercept indicator matrix
B <- cbind(Phi, W)                                      # extension of the basis matrix
P <- bdiag( crossprod(Delta) , diag(D))                 # extension of the penaly matrix
C_ext <- cbind(C, matrix(0, nrow=dim(C)[1], ncol=D))    # extension of the shape constraint matrix
```

```{r}
A <- crossprod(B) + lambda*P            # matrix for the quadratic program 
b <- crossprod(y,B)                     # matrix for the quadratic program 
sol <- solve.QP(A, b, t(C_ext))         # solve the quadratic program
alpha <- sol$solution[1:K]
```

## Output



