###----------------------------------------------------------------------------
### Read Txt format
###----------------------------------------------------------------------------


mytable <- read.table(file = "table1.txt", sep = "\t", header = TRUE)

as.matrix(mytable[, 2:5])


write.table(mytable, file = "mysavedtable.txt", sep = ";")

###----------------------------------------------------------------------------
### Read Excel
###----------------------------------------------------------------------------

install.packages("XLConnect")

## library("XLConnect")
require("XLConnect")

NameList <- loadWorkbook("NameList.xls", create = TRUE)

out <- readWorksheet(NameList, sheet = "12-1",
                     startRow = 0, endRow = 100,
                     startCol = 0, endCol = 5)

library("RODBC")

###----------------------------------------------------------------------------
### Read and write Matlab format
###----------------------------------------------------------------------------
install.packages("R.matlab")
require("R.matlab")

a <- readMat("RajanData.mat")


###----------------------------------------------------------------------------
### Read Database
###----------------------------------------------------------------------------

# RODBC Example
# import 2 tables (Crime and Punishment) from a DBMS
# into R data frames (and call them crimedat and pundat)

library(RODBC)
myconn <-odbcConnect("mydsn", uid="Rob", pwd="aardvark")
crimedat <- sqlFetch(myconn, "Crime")
pundat <- sqlQuery(myconn, "select * from Punishment")
close(myconn)



channel <- odbcConnect("test")
sqlSave(channel, USArrests, rownames = "State", verbose = TRUE)
# options(dec=".") # optional, if DBMS is not locale-aware or set to ASCII
## note case of State, Murder, Rape are DBMS-dependent,
## and some drivers need column and table names double-quoted.
sqlQuery(channel, paste("select State, Murder from USArrests",
                        "where Rape > 30 order by Murder"))
close(channel)
## End(Not run)

###----------------------------------------------------------------------------
### Load an image
###----------------------------------------------------------------------------

install.packages("jpeg")

library("jpeg")
myprofile <- readJPEG("FengLi.jpg")
dim(myprofile)
image(1:1539, 1:1154, myprofile[, , 1])
image(1:1539, 1:1154, myprofile[, , 2])
image(1:1539, 1:1154, myprofile[, , 3])

plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "")
rasterImage(myprofile, 100, 300, 150, 350, interpolate = FALSE)
rasterImage(myprofile, 100, 400, 150, 450)
rasterImage(myprofile, 200, 300, 200 + xinch(.5), 300 + yinch(.3), interpolate = FALSE)
rasterImage(myprofile, 200, 400, 250, 450, angle = 15, interpolate = FALSE)
