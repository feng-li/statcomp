## Create a sequence

## Sequence with special patterns

## 1, 2, 3, ..., n

rep()


## seq(-5, 5, length.out = 8)
## seq(-5, 5, length = 8)
## seq(-5, 5, length.o = 8)

rep(1:5, each = 3)

cbind(rep(1:5, 5), rep(1:5, each = 5))

cbind(1:3, 2:4, 4:6) # combine vectors

mycbind <- function(...)
    {
        args <- list(...)

        print(args)

        nargs <- length(args)
        out <- NULL

        for(i in 1:nargs)
            {
                out <- cbind(out, args[[i]])
            }

        return(out)
    }


A <- matrix(1:24, 4, 6)


B <- array(1:24, c(2, 3, 4))



## Matrices

A  = matrix(1:24, 4, 12)



## Array

## Data Frame

NM = c("Zhang3", "Li4", "Wang5")
ID = c(201701, 201702, 201703)
Stat = c(90, 59, 89)


## List

Class1 = data.frame(ID=ID, Stat=Stat)
Class2 = matrix(1:24, 4, 6)
Class3 = 1:10

mylist = list(Class1, Class2, Class3)
