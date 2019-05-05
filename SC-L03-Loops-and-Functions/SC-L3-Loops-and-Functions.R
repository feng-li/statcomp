maxdiff  =  function(a, b = a)
    {
        if(!is.numeric(a) ||  !is.numeric(b))
        {
          stop("a and be must be numeric.")
        }


        out <- max(a) -min(b)

        return(out)
    }


dat  = matrix(rnorm(1000), 200, 5)
nRow = dim(dat)[1]
out <- rep(NA, nRow)

for(i in 1:nRow)
  {
    a = dat[i, ]
    outi = maxdiff(a)
    out[i] <- outi
  }

out <- rep(NA, nRow)
i <- 1
while(i <= nRow)
  {
    a = dat[i, ]
    outi = maxdiff(a)
    out[i] <- outi

    i  = i + 1

    if(i == 10) break

  }


#######################
x <- c(1, 2, 3)
y <- c(2, 4, 5)

if(x<y)
    {
        print("Good!")
        print("very good")
    }

n = 10000
A = matrix(1:n, n, 1)
B =  matrix(1:n, n, 1)

a = Sys.time()
C = matrix(NA, n, 1)
for(i in 1:n)
    {
       C[i]  =  A[i] + B[i]
    }
b = Sys.time()-a


a = Sys.time()
## C = matrix(NA, n, 1)
C  = NULL
for(i in 1:n)
    {
       C[i]  =  A[i] + B[i]
    }
b = Sys.time()-a


i = 1
while(TRUE)
    {
        print(i)
        C[i, ]  =  A[i, ] + B[i, ]
        i  =  i + 1
    }



#####

n = 100000
B = matrix(n,n,1)
C = matrix(n,n,1)

tic = proc.time()

## A = NULL
A  = matrix(NA, n, 1)
for(i in 1:n)
{
A[i] = B[i] + C[i]
}

toc = proc.time()
toc -tic


myFun2 =  function(y = 5)
    {
     if(is.numeric(y) == FALSE | y<0)
         {
             stop("Negative y is not allowed!")
         }


      out = log(y)
      return(out)
    }

myFun = function(x)
    {
        logx = myFun2(x)
        logx.length = length(logx)
        return(logx.length)
    }


matrixInv <- function(x)

    {
        if(is.numeric(x) & nrow(x)  == ncol(x))
            {

              xTry =  try(solve(x), silent = TRUE)
              if(is(xTry, "try-error"))
                  {
                      return(NA)
                  }


              #xInv =  solve(x)
              #xLen = length(x)
              return(xTry)
            }
        else
            {
                stop("X must be a numeric squared matrix")
            }
    }






quaroot <- function(a, b, c)
  {
     if(b^2-4*a*c>=0)
      {
          x1 <- (-b+sqrt(b^2-4*a*c))/(2*a)
          x2 <- (-b-sqrt(b^2-4*a*c))/(2*a)
      }
     else{
         stop("Something wrong")
     }
    out <- c(x1, x2)
    return(out)
  }

n <- 100
out <- matrix(0, n, 2)
for(i in 1:n)
    {
        a = rnorm(1)
        b = rnorm(1)
        c = rnorm(1)
        out[i,] <- quaroot(a, b, c)

        ## if(is(out.try, "try-error"))
        ##     {
        ##        ##  out[i, ] <- NA
        ##         browser()
        ##     }
        ## else
        ##     {
        ##         out[i, ] <- out.try
        ##     }
    }

#####################################################################
## Lab: Quaroot
#####################################################################

quaroot(1, 4, -1)

quaroot(5, 2, 1) # warning with NaN
quaroot <- function(a, b, c)
  {
    d <- b^2-4*a*c

    if(d<0)
      {
        x1 <- NA
        x2 <- NA
      }
    else
      {
        x1 <- (-b+sqrt(d))/(2*a)
        x2 <- (-b-sqrt(d))/(2*a)
      }

    out <- c(x1, x2)
    return(out)
  }

quaroot(5, 2, 1) # NA without warnings




###----------------------------------------------------------------------------
### Lab: Leapday
###----------------------------------------------------------------------------

# 1.1
isLeapYear <- function(year)
  {
   mod4 <- year%%4
   mod100 <- year%%100
   mod400 <- year%%400

   LeapdayIndex = ((mod4  == 0 & mod100 != 0) |  mod400  == 0)


   ## if((mod4  == 0 & mod100 != 0) |  mod400  == 0)
   ##   {
   ##     out <- TRUE
   ##   }
   ## else
   ##   {
   ##     out <- FALSE
   ##   }
   ## return(out)

   out =  list(mod4 = mod4, LeapdayIndex = LeapdayIndex)
   return(out)
   ## return(LeapdayIndex)
  }

# 1.2
isLeapday(1998)
isLeapday(2000)
isLeapday(2002)
isLeapday(2100)



whichLeapday <- function(year)
  {
    if(!is.numeric(year))
      {
        stop("You must specify a numerical input.")
      }

   mod4 <- year%%4
   mod100 <- year%%100
   mod400 <- year%%400

   leapdayIndex <- ((mod4  == 0 & mod100 != 0) | mod400  == 0)

    out <- year[leapdayIndex]

    return(out)
  }

whichLeapday(c(1998:2100))

fact <- function(n, log)
{

 logsum =  sum(log(1:n))


  if(log == TRUE)
    {
      return(logsum)
    }
  else
    {
      return(exp(logsum))
    }
}

mygamma <- function(n)
{
    if(n > 1)
    {
        out =  mygamma(n-1)*(n-1)
    }
    else if (n  == 1)
    {
        out <- 1
    }

    return(out)
}



logmygamma <- function(n)
{
    if(n > 1)
    {
        # out = log(mygamma(n-1)*(n-1))
        out = logmygamma(n-1) + log(n-1)

    }
    else if (n  == 1)
    {
        out <- log(1)
    }

    return(out)
}


mygamma <- function(n, log)
{
    if(n > 1)
    {
        out.log = mygamma(n-1, log=TRUE) + log(n-1)
    }
    else if (n  == 1)
    {
        out.log <- log(1)
    }

    if(log == TRUE)
    {
        out = out.log
    }    else {
        out = exp(out.log)
    }
    return(out)
}
