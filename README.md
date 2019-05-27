# PS_con: Shape-Constrained Penalized Splines
Penalized splines is a popular method for function estimation under the assumption of smoothness.
The provided code extends the well known P-spline method to constraints on its shape, such as monotonicity or convexity, and allows for the consideration of additional random effects within the model.
Further detail can be found in [Wagner et al., 2017](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssa.12295).

## Manual
For selected parameters, the B-spline basis matrix and the matrix of the penalty term are assembled, where the difference penalty or the curvature penalty (default) can be used.
```{r}
Phi <- bspline_matrix(x, m, q, Omega )    # B-spline basis matrix at the covariates
#Delta <- diff( diag(K), diff=l)          # difference penalty
Delta <- L2_norm_matrix(m, q, l, Omega)   # curvature penalty (default)
```


## Output
