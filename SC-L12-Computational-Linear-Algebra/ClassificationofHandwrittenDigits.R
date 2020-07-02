###############################################################################
### 1
###############################################################################

## The image matrix for training sample. 256x1707
azip <- read.table("azip.dat")

## The true digits given in the training sample. lenght = 1707
dzip <- as.numeric(read.table("dzip.dat"))

## The testing image matrix. 256x2007
testzip <- read.table("testzip.dat")

## The true digits for the testing sample. length = 2007
dtest <- read.table("dtest.dat")

## Display the image
i <- 120
image(matrix(azip[ , i], ncol = 16)[ , 16:1], col = gray(255:0/255))

###----------------------------------------------------------------------------
### The naive method
###----------------------------------------------------------------------------

## The mean of training sample of a single digit
digits <- 0:9 # The possible digits in the US postal code
img.mean <- matrix(0, 256, length(digits))

for(i in digits)
  {
    idx <- (i == dzip) # the location indicator for the ith digit
    imgi <- azip[, idx, drop = FALSE]
    imgi.mean <- rowMeans(imgi)
    img.mean[,  i+1] <- imgi.mean
  }


## Plot the mean image
par(mfrow = c(2, 5))
for(i in 1:10)
  {
    image(matrix(img.mean[, i], 16)[, 16:1],
          col = gray(255:0/255))
  }


## Sketch a distance function to compute the Euclidean distance between two
## matrices in row wise.
rdist <- function(X, Y)
  {
    dim.X <- dim(X)
    dim.Y <- dim(Y)
    sum.X <- matrix(rowSums(X^2), dim.X[1], dim.Y[1])
    sum.Y <- matrix(rowSums(Y^2), dim.X[1], dim.Y[1], byrow = TRUE)
    dist0 <- sum.X + sum.Y - 2*tcrossprod(X, Y)
    out <- sqrt(dist0)
    return(out)
  }

## For an unknown testing digit image, compare the distance to the means
test.sample <- 1:5

## Let's first plot those testing image
par(mfcol = c(ceiling(length(test.sample)/5), 5)) # five columns

for(i in test.sample)
  {
    image(matrix(testzip[, i], 16)[, 16:1], col = gray(255:0/255))
  }

## Calculate the distance from testing sample to the mean in the training
## sample.
img.dist <- rdist(t(testzip[, test.sample]), t(img.mean))


## The classification results by the naive method
apply(img.dist, 1, which.min)-1

###----------------------------------------------------------------------------
### The SVD method
###----------------------------------------------------------------------------

## Compute the singular matrix of a single digit in the training sample
digit <- 9

## Subtract the matrix for that digit
img.mat <- azip[, digit == dzip, drop = FALSE]
img.matSVD <- svd(img.mat)


X = img.matSVD$u
y = testzip[, 1]


img.mat <- azip[, dzip == 7, drop = FALSE]
img.matSVD <- svd(img.mat)
X = img.matSVD$u


## Plot the singular matrix under different basis.
par(mfrow = c(2, 5))
for(i in 1:10)
  {
    image(matrix(img.matSVD$u[, i], 16)[, 16:1], col = gray(255:0/255), main =
          paste("singular image ", i, sep = ""))
  }

## Do the least square method with different basis and find the minimal
## residuals.

## The testing digit matrix
test.idx <- 3
image(matrix(testzip[, test.idx], 16)[, 16:1], col = gray(255:0/255), main =
      paste("Testing digit"))

resid.norm <- matrix(NA, 10, 1, dimnames = list(0:9,"resid"))
for(i in 0:9)
  {
    img.mat <- azip[, i == dzip, drop = FALSE]
    img.matSVD <- svd(img.mat)
    ## basis.max <- ncol(img.matSVD$u)
    basis.max <- 10

    resid.norm[i+1, ] <- norm(matrix(lm(testzip[, test.idx]~
                                        0+img.matSVD$u[, 1:basis.max])$resid), "F")
  }
resid.norm

## Find the optimum basis function.

test.idx <- 1:100
basis.max <- 1:20
nCorrectOut <- numeric(length(test.idx))
resid.norm <- matrix(0, 10, 1, dimnames = list(0:9,"resid"))
for(b in basis.max)
  {
    nCorrect <- 0
    for(k in test.idx)
      {
        for(i in 0:9)
          {
            img.mat <- azip[, i == dzip, drop = FALSE]
            img.matSVD <- svd(img.mat)
            resid.norm[i+1, ] <- norm(matrix(lm(testzip[, test.idx[k]]~
                                                0+img.matSVD$u[, 1:basis.max[b]])$resid), "F")
          }
        if((which.min(resid.norm)-1)  == dtest[k])
          {
            nCorrect <- nCorrect+1
          }
      }
    nCorrectOut[b] <- nCorrect
    print(b)
  }

###----------------------------------------------------------------------------
### Regression with QR,  SVD
###----------------------------------------------------------------------------

y <- matrix(1:5)
X <- cbind(1, rnorm(5), runif(5), 1+runif(5)/1000)

## Traditional method
betaHat1 <- solve(t(X)%*%X)%*%t(X)%*%y

## QR
XQR <- qr(X)
Q <- qr.Q(XQR)
R <- qr.R(XQR)

XPlus2 <- solve(R)%*%t(Q)

XPlus2%*%X

betaHat2 <- XPlus2%*%y


## SVD

XSVD <- svd(X)
U <- XSVD$u
V <- XSVD$v
Sigma <- diag(XSVD$d)
Sigma1 <- diag(1/XSVD$d)

XPlus3 <- V%*%Sigma1%*%t(U)

betaHat3 <- XPlus3%*%y


cbind(betaHat1, betaHat2, betaHat3)
