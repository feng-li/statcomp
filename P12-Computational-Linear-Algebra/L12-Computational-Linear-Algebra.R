###----------------------------------------------------------------------------
### Determinant
###----------------------------------------------------------------------------

A  =   matrix(rnorm(100), 10, 10)
det(A)

B = matrix(c(1:4), 2, 2)
det(B)

determinant(B, logarithm = TRUE)


det(solve(A))

1/det(A)

###----------------------------------------------------------------------------
### Eigenvalues
###----------------------------------------------------------------------------

A

eigen(B)


sum(eigen(B)$values)  == sum(diag(B))

prod(eigen(B)$values)  == det(B)

all.equal(prod(eigen(B)$values), det(B) )

abs(prod(eigen(B)$values)- det(B))< 1e-8

## Trace of hat matrix
beta <- matrix(c(1, 2, 3, 4))
n = 1000
X = matrix(rnorm(n*4), n)
epsilon <- rnorm(n)
Y = X%*%beta + epsilon

beta_hat <- solve(t(X)%*%X)%*%t(X)%*%Y
H = X%*%solve(t(X)%*%X)%*%t(X)
sum(eigen(H)$values)
sum(diag(H))


F <- matrix(c(1, 2, 1, 2), 2)
det(F)
solve(F)
eigen(F)

###----------------------------------------------------------------------------
### Upper triangular matrix
###----------------------------------------------------------------------------


n <- 25
nrow <- sqrt(n)

whichcols <- ceiling((1:n)/nrow)
whichrows <- (1:n)%%nrow
whichrows[whichrows == 0] <- nrow


T <- matrix(rnorm(25), 5, 5)

T[whichrows<whichcols]
T
###----------------------------------------------------------------------------
###
###----------------------------- Outer product-----------------------------------------------

x1 <- 1:5
outer(x1, x1, "/")

y1 <- 1:5
outer(x1, y1, "+")

###----------------------------------------------------------------------------
###
###----------------------------- General Inverse-----------------------------------------------

library("MASS")
Y = matrix(c(1, 2, 1, 2), 2)

YP <- ginv(Y)


a <- matrix(rnorm(100), 10)


###----------------------------------------------------------------------------
### Linear regression with QR decomposition
###----------------------------------------------------------------------------

X <- cbind(1, 1:5, matrix(rnorm(30), 5))
y <- matrix(c(7.97, 10.2, 14.2, 16.0, 21.2))


## QR decomposition of X matrix
QR <- qr(X)

## Extract Q and R matrices
R <- qr.R(QR, complete = TRUE)
Q <- qr.Q(QR, complete = TRUE)

betaHat <- solve(Q)%*%t(Q)%*%y


betaHat <- qr.solve(X, y)


X1 <- cbind(X, 3)
QR <- qr(X1)

## Extract Q and R matrices
R <- qr.R(QR)
Q <- qr.Q(QR)
betaHat <- solve(R)%*%t(Q)%*%y


###----------------------------------------------------------------------------
### SVD
###----------------------------------------------------------------------------

a <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1), 9, 4)

svd(a)

###----------------------------------------------------------------------------
###
###----------------------------------------------------------------------------

library(png)
x <- readPNG("flower.png")
dim(x)

plot.new()
rasterImage(x, 100, 50, 20, 70)
