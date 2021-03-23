###----------------------------------------------------------------------------
### Debugging
###----------------------------------------------------------------------------

func <- function(x, y)
    {
        x^3-3*x+y^3-3*y
    }

integLoop <- function(func, xint, yint, n)
{
    local_sum <- 0
    xincr <- ( xint[2]-xint[1] ) / n
    yincr <- ( yint[2]-yint[1] ) / n
    browser()
    for(xi in seq(xint[1], xint[2],length.out = n)){
        for(yi in seq(yint[1], yint[2],length.out = n)){

            box <- func(xi, yi) * xincr * yincr
            local_sum <- local_sum + box

        }
    }
    return(local_sum)
}

for(i in 1:1000)
    {
        A <- matrix(runif(100), 10, 10)
        A.try <- try(solve(A), silent = TRUE)
        if(is(A.try, "try-error"))
            {
                A.inv <- NA
            }
        else
            {
                A.inv <- A.try
            }
    }


###----------------------------------------------------------------------------
### Profiling
###----------------------------------------------------------------------------

library(MASS)
library(boot)

storm.fm <- nls(Time ~ b*Viscosity/(Wt - c), stormer,
                start = c(b=30.401, c=2.2183))
st <- cbind(stormer, fit=fitted(storm.fm))
storm.bf <- function(rs, i)
{
    st$Time <- st$fit + rs[i]
    tmp <- nls(Time ~ (b * Viscosity)/(Wt - c), st,
               start = coef(storm.fm))
    tmp$m$getAllPars()
}
rs <- scale(resid(storm.fm), scale = FALSE) # remove the mean


## Speed profiling
Rprof("boot.out")
storm.boot <- boot(rs, storm.bf, R = 4999) # slow enough to profile
Rprof(NULL)
summaryRprof("boot.out")

## Memory profiling
Rprof("boot.out", memory.profiling = TRUE)
storm.boot <- boot(rs, storm.bf, R = 4999) # slow enough to profile
Rprof(NULL)
summaryRprof("boot.out")
