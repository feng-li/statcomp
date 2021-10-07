###----------------------------------------------------------------------------
### Newton-Raphson Algorithm to Find the Root
###----------------------------------------------------------------------------
## The target function
f <- function(x, b)
    {
        x^2 -b
    }
## The gradient wrt x
df <- function(x, b)
    {
        2*x
    }
## Loop to find the best x
n <- 50
x <- rep(NA, n)
b <- 3
x[1] <- -100 # The initial value
for(i in 1:n)
{
    x[i+1] <- x[i]- (df(x[i], b))^(-1)*f(x[i], b)
}
## Check the convergence
par(mfrow = c(3, 1))
plot(x, type = "l", col = "red")
plot(f(x, b), type = "l", col = "red")
plot(df(x, b), type = "l", col = "red")


###----------------------------------------------------------------------------
### Newton-Raphson Algorithm to Maximize a Function
###----------------------------------------------------------------------------

## The target function
f <- function(x, a, b, c)
    {
        a*x^2 +b*x +c
    }

## The gradient wrt x
df <- function(x, a, b, c)
    {
        2*a*x +b
    }
## The Hessian wrt x
d2f <- function(x, a, b, c)
    {
        2*a
    }

## Loop to find the best x
n <- 50
x <- rep(NA, n)
a <- -1
b <- 3
c <- 10
x[1] <- 300 # The initial value
for(i in 1:n)
{
    x[i+1] <- x[i]- (d2f(x[i], a, b, c))^(-1)*df(x[i], a, b, c)
}
## Check the convergence
par(mfrow = c(4, 1))
plot(x, type = "l", col = "red")
plot(f(x, a, b, c), type = "l", col = "red")
plot(df(x, a, b, c), type = "l", col = "red")
plot(rep(d2f(x, a, b, c), n), type = "l", col = "red")
