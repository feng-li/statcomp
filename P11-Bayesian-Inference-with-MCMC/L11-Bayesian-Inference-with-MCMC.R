Rprof()

## Bayesian Regression with MCMC

###----------------------------------------------
### The data
###----------------------------------------------
## DGP: Data generating process
n <- 1000
p <- 3

epsilon <- rnorm(n, 0, 0.1)
beta <- matrix(c(-3, 1, 3, 5))
X <- cbind(1, matrix(runif(n*p), n))
y <- X%*%beta + epsilon

###----------------------------------------------
### Write down the log posterior function
###----------------------------------------------

library(mvtnorm)
LogPosterior <- function(beta, sigma2, y, X)
    {
        p <- length(beta)

        ## The log likelihood function
        LogLikelihood <- sum(dnorm(x = y,
                                   mean = X%*%beta,
                                   sd = sqrt(sigma2),
                                   log = TRUE))
        ## browser(text = "I am in log Posterior")
        ## The priors (use your expertize)
        LogPrior4beta <- dmvnorm(
                x = t(beta),
                mean = matrix(0, 1, p),
                sigma = diag(p), log = TRUE)

        LogPrior4sigma2 <- dchisq(x = sigma2, df = 10,
                                  log = TRUE)

        LogPrior <- LogPrior4beta + LogPrior4sigma2

        ## The log posterior
        LogPosterior <- LogLikelihood + LogPrior

        ## Return
        return(LogPosterior)
    }

###----------------------------------------------
### The Metropolis algorithm with Gibbs
###----------------------------------------------

nIter <- 10000
p <- ncol(X)

## Reserve space
beta <- matrix(NA, nIter, p)
sigma2 <- matrix(NA, nIter, 1)
acc <- matrix(NA, nIter, 2)

## Initial values
beta[1, ] <- c(-2, 3, 4, 1)
sigma2[1, ] <- 0.5


for(i in 2:nIter)
{
    beta.current <- beta[i-1, ]
    sigma2.current <- sigma2[i-1,]

        ## browser(text = paste("I am in a big loop:", i))

    ## The Metropolis Algorithm For beta
    beta.prop <- rmvnorm(
            n = 1,
            mean = matrix(beta.current, 1),
            sigma = diag(p)) # FV

    logPosterior.beta.prop <- LogPosterior(
        beta = t(beta.prop),
        sigma2 = sigma2.current,
        y = y, X = X)

    logPosterior.beta.current <- LogPosterior(
        beta = beta.current,
        sigma2 = sigma2.current,
        y = y, X = X)

    logratio <- (logPosterior.beta.prop -
                 logPosterior.beta.current)

    acc.prob.beta <- min(exp(logratio), 1)

    if(!is.na(acc.prob.beta) &
       runif(1) < acc.prob.beta) # accept the proposal
        {
            beta.current <- beta.prop
        }
    acc[i, 1] <- acc.prob.beta
    beta[i, ] <- beta.current


    ## The Metropolis Algorithm For sigma2, Check: is this right?
    sigma2.prop <- rgamma(1, 1) # FV
    logPosterior.sigma2.prop <- LogPosterior(
        beta = matrix(beta.current),
        sigma2 = sigma2.prop,
        y = y, X = X)

    logPosterior.sigma2.current <- LogPosterior(
        beta = matrix(beta.current),
        sigma2 = sigma2.current,
        y = y, X = X)

    logratio <- logPosterior.sigma2.prop -
        logPosterior.sigma2.current

    acc.prob.sigma2 <- min(exp(logratio), 1)

    if(!is.na(acc.prob.sigma2)&
       runif(1, 0, 1) < acc.prob.sigma2) # accept the proposal
        {
            sigma2.current <- sigma2.prop
        }

    ## Update the matrix
    acc[i, 2] <- acc.prob.sigma2
    sigma2[i, ] <- sigma2.current

    if ((i %% 100) == 0) cat(i, " of ", nIter, "iterations done!\n")
}

## Summarize beta and sigma2
apply(beta, 2, mean)
apply(beta, 2, sd)

mean(sigma2)
sd(sigma2)

Rprof(NULL)
print(summaryRprof())
