###----------------------------------------------------------------------------
### Inverse Method
###----------------------------------------------------------------------------

u2x <- function(u, lambda)
    {
        x <- - log(1-u)/lambda
        return(x)
    }

rmyexp <- function(n, lambda)
    {

        ## n <- 1000
        ## lambda <- 0.5
        u <- runif(n, 0, 1)
        x <- u2x(u, lambda)
        return(x)
    }

n = 1000
lambda = 2

par(mfrow = c(1, 2))
ylim <- c(0, 2)
xlim <- c(0, 5)
hist(rmyexp(n, lambda),
     breaks = 20, freq = FALSE,
     xlim = xlim, ylim = ylim)

xax <- seq(0, max(x), 0.01)
xDens <- lambda*exp(-lambda*xax)

lines(xax, xDens, type = "l",
     col = "blue", lwd = 3)


hist(rexp(n, rate = lambda),
     breaks = 20, freq = FALSE,
     xlim = xlim, ylim = ylim)
lines(xax, xDens, type = "l",
     col = "blue", lwd = 3)

###----------------------------------------------------------------------------
### Reject Method
###---------------------- ------------------------------------------------------

i <- 1
y <- NULL
while(i <1000)
    {
        u <- runif(1, 0, 1)
        v <- runif(1, -5, 5)
        M <- 1/(sqrt(2*pi)*0.1)

        test <- dnorm(v, 0, 1)/(0.1*M)
        if(u<test)
            {
                y[i] <- v
                i <- i+1
            }
}

###----------------------------------------------------------------------------
### Metropolis Algorithm
###----------------------------------------------------------------------------
