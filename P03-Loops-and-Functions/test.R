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

quaroot(1, 9, 3)
quaroot(1, 2, 3)
quaroot(1, 4, 3)

myBatch <- function()
    {
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
    }

myBatch()
