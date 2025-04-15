data = data.frame(
  'gen' = factor(rep(paste0("G",1:5),2)),
  'rept' = factor(rep(paste0("R",1:2),each= 5)),
  'y' = c(18.36, 8.23, 16, 18.25, 9.95, 21.54, 7.25,10, 20, 10.01)
)


X = model.matrix(y ~-1 + rept, data = data)
Z = model.matrix(y ~-1 + gen, data= data)


Xly = crossprod(X, data$y)

Zly = crossprod(Z, data$y)


sigma2g = 26.42
sigma2e = 6.17

G = sigma2g*diag(1,5)
R = sigma2e*diag(1,10)

H = tcrossprod(Z %*% G, Z) + R
Hinv = solve(H)

BLUE = solve(crossprod(X, Hinv) %*% X) %*% (crossprod(X, Hinv)%*%na.omit(data$y))


BLUP = (tcrossprod(G, Z) %*% Hinv) %*% (data$y - X %*% BLUE)
