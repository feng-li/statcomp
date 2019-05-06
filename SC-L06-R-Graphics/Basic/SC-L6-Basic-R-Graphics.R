par(mfrow = c(1, 2))

x <- seq(-5, 5, 1)
plot(x = x, col = "purple", ylab = c(-1, 1),
     xlab = "hello world!",
     main = "My best plot美图秀秀", type = "p",
     lwd = 1)
points(x, type = "l")

n <- 100
a  = rnorm(100)
plot(a, xlab = "Day", ylab = "Return", type = "l",
     col = "red")

b  = rnorm(100, sd = 0.3)
points(b, col = "blue", type = "l", lty = "dotted",
       lwd = 10)

# points(x = sort(x), sin(x[order(x)]), type = "l")


plot(x = rnorm(100), pch = 14)



par(mfcol = c(3, 2))


#######################################
# 2.1.1 the style of axis labels      #
#######################################
par(mfrow=c(2,2))

par(las=0)
plot(0,main="par(las=0) default")


par(las=1)
plot(0,main="par(las=1)")
par(las=2)
plot(0,main="par(las=2)")
par(las=3)
plot(0,main="par(las=3)")

#######################################
# 2.1.2 the position of axis labels   #
#######################################
par(mfrow=c(2,2))
par(adj=0)
plot(0,main="par(adj=0)",ylab="ylab",xlab="xlab")
par(adj=.5)
plot(0,main="par(adj=.5)default",ylab="ylab",xlab="xlab")
par(adj=.7)
plot(1,main="par(adj=.7)",ylab="ylab",xlab="xlab")
par(adj=1)
plot(1,main="par(adj=1)",ylab="ylab",xlab="xlab")

#######################################
# 2.1.3 the type of the box          #
#######################################

par(mfrow=c(2,4))
par(bty="o")
plot(0,main="par(bty=\"o\")")
par(bty="l")
plot(0,main="par(bty=\"l\")")
par(bty="7")
plot(0,main="par(bty=\"7\")")
par(bty="n")
plot(0,main="par(bty=\"n\")")
par(bty="c")
plot(0,main="par(bty=\"c\")")
par(bty="u")
plot(0,main="par(bty=\"u\")")
par(bty="]")
plot(0,main="par(bty=\"]\")")

#######################################
# 2.1.4 figure in figure              #
#######################################

x = seq(-5, 5, 0.1)
fn = dnorm(x)
ft = dt(x, df = 4)
plot(x, fn, xlab = "X", ylab = "Density", main = "Density plot", col = "Blue", type = "l", lwd = 3)
points(x, ft, col = "Red", type = "l", lwd = 3)



x <- rt(100, df=4)                                      # give 100 random t with df 12
hist( x, col = "light blue")                             # make a histgram
op <- par(fig = c(.02,.42,.53,.99),new = TRUE)           # prepar to add a new small graph
qqnorm(x, xlab = "", ylab = "", main = "", axes = FALSE) # show the qqnorm in the small graph
qqline(x, col = "red", lwd = 2)                          # add a qq line in the small graph
box(lwd=3)                                               # tell the line width of the frame box for the small graph
par(op)                                                  # apply


#######################################
# 2.1.5 Double axises                 #
#######################################
set.seed(1)      # set the seed
par(mar=c(5,5,5,5))  # set the margin
x1<-sort(rnorm(100,mean=20,sd=5)) # asumption this Interest rate
x2<-x1^3 # this is GDP
plot(x2,axes=FALSE,type="l",col="blue",xlab="",ylab="") # plot GDP
axis(1) # add the axis on bottom
axis(2,col="blue") # add the asis on left
par(new=TRUE)  # add new plot on the current graph
plot(x1,,axes=FALSE,type="l",col="red",xlab="",ylab="") # plot Interest rate
axis(4,col="red") # add the axis on right
mtext(c("Time","GDP(mil. $)","Interest rate(%)"),side=c(1,2,4),line=3) # give the three axises labels
title("GDP with interest rate") # give the main title

#######################################
# 2.2   figures layout                #
#######################################
split.screen(c(2,1))                      # split display into two screens(1,2)
split.screen(c(1,2),2)                    # split bottom half in two(3,4)
screen(1)                                 # active screen (1), the top one
plot((1:10)^2,type="l",main="Screen(1)")  # plot something on screen 1
screen(3)                                 # active screen (3), the bottom left one
plot((1:10)^.5,type="l",main="Screen(3)") # plot something on screen 3
screen(4)                                 # active screen (4), the bottom right one
plot((1:10),type="l",main="Screen(4)")    # plot something on screen 4
close.screen(all = TRUE)                  # exit split-screen mode

## Example
split.screen(c(2,2))
split.screen(c(2,1),4)

screen(1)
plot((1:10)^2, type = "l", main = "screen1")

screen(2)
plot((1:10)^2, type = "l", main = "2")


screen(3)
plot((1:10)^2, type = "l", main = "3")

screen(5)
plot(1, type = "l", main = "screen5")

screen(6)
plot(1, type = "l", main = "screen6")

#######################################
# 2.2.1 grid in figures               #
#######################################
par(mfrow=c(2,2))
plot(1:3)
grid(nx=NA, ny=5,col="red") # grid only in y-direction
plot(1:3)
grid(nx=5,ny=NA, col="red") # only x-direction
plot(1:3)
grid(nx = 5, ny = 5,lty="dashed",col="blue") # both x and Y
plot(1:3)
grid(nx = 3, ny = 5) # x and y but different numbers

#######################################
# 2.3 Math in Graphics               #
######################################
# make a hist but use the density
hist(rnorm(100),xlim=c(-3,3),freq=FALSE)
# add a normal curve here
# 'expr' must be a function or an expression containing 'x'
x<-seq(-5,5,.01)
curve(dnorm(x),add=TRUE)
# add a density expression
text(x=2, y=.35, expression(paste(frac(1,sqrt(2*pi)*sigma),
    exp(-frac(1,2*sigma^2)*((x-mu)^2)))))
# check all symbols in R
demo(plotmath)

#######################################
# 2.4 long title line                 #
#######################################
plot(0,main=paste(strwrap("This is a really long title that I can not
     type it proprtly but we can use strwap function to wrap it", width=55)))

plot(0,main="This is a really long title that I can not type it proprtly but we can use strwap function to wrap it")

#######################################
# 2.5 Legend in graphics              #
#######################################
x<-seq(0,1,.01)
par(lwd=2)
plot(x,dbeta(x,5,1),type="l", ylim=c(0,3), col="green", lty=1,xlab="",ylab="")
lines(x,dbeta(x,.5,.5), col="red",lty=2)
lines(x,dbeta(x,1,3),col="blue",lty=3)
lines(x,dbeta(x,2,2), col="pink",lty=4)
lines(x,dbeta(x,2,5), col="black",lty=5)
lines(x,dbeta(x,1,1),col="orange",lty=6)

legend("top",
c("Beta(5,1)","Beta(.5,.5)","Beta(1,3)","Beta(2,2)","Beta(2,5)","Beta(1,1)"),
lty=c(1,2,3,4,5,6),
col=c("green","red","blue","pink","black","orange"),
ncol=3)

legend(x = 0.5, y =  1.5,
   c("Beta(5,1)","Beta(.5,.5)","Beta(1,3)","Beta(2,2)","Beta(2,5)","Beta(1,1)"),
       lty=c(1,2,3,4,5,6),
       col="green",
       lwd = c(1, 1, 5, 1, 1, 1),
       ncol=3)

#######################################
# 2.6.1  Stars                       #
#######################################
palette(rainbow(12)) ## # 12 color rainbow
stars(mtcars[, 1:7], len = 0.8, key.loc = c(12, 1.5),
      main = "Motor Trend Cars", draw.segments = TRUE)

#######################################
# 2.6.2  Faces                       #
#######################################
library("aplpack")
faces(matrix(sample(1:1000,128,),16,8),main="random faces")


#######################################
# 2.7 Time series plot                #
#######################################
# prepare the monthly time series data (8 series), from the Jan, 1961
z <- ts(matrix(rt(200 * 8, df = 3), 200, 8), start = c(1961, 1), frequency = 12)
# we can plot some of them for a period (1961.1-1969.12) in multiple frames
z1 <- window(z[,1:3], end = c(1969,12))
plot(z1, type = "b")
# or we can plot them in a single frame
plot(z, plot.type="single", lty=1:3, col=4:2)
# the leg opeator plot
# we use the monthly mean relative sunspot numbers from 1749 to 1983 in "sunspots"
# x_{t-1} aginst x_t
plot(lag(sunspots, 1), sunspots, pch = ".")
# the acf and pacf
# use the data "Luteinizing Hormone in Blood Samples"
acf(lh)
pacf(lh)


#######################################
# 2.8 grey band for standard deviation#
#######################################
set.seed(125) # set the seed
x<-1:200
y<-cumsum(rnorm(200)) # simulate an AR(1) model
int1<-y-1.96 # calculate the 95% level confident interval
int2<-y+1.96
plot(y,xlab="t",ylab=expression(y[t])) # plot the data
polygon(c(x, rev(x)), c(int1, rev(int2)),col="grey",border=NA) # show the band
lines(y,type="l",col="blue",lwd=2) # plot the line again


#######################################
# 2.9 3D plot (non-interactive)        #
#######################################
window()
x <- seq(-10, 10, .3)
y <- x
rotsinc <- function(x,y)
{
    sinc <- function(x) { y <- sin(x)/x ; y[is.na(y)] <- 1; y }
    sinc( sqrt(x^2+y^2) )
}

z <- outer(x, y, rotsinc)

## x <- seq(-10, 10, 0.1)
## y <- x
## xy <- expand.grid(x, y)
## z <- rotsinc(x = xy[, 1], y = xy[, 2])

persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightgreen",
      ltheta = 120, shade = 0.75, ticktype = "detailed",
      xlab = "X", ylab = "Y", zlab = "Z")
title(main = expression(z == sinc(sqrt(x^2 + y^2))))

#######################################
# 2.10  Project 3D onto 2D with colors#
#######################################
window()
x<-seq(-10,10,.3)
y<-x
z<-outer(x, y, rotsinc) # the rotsinc is defined in section 2.6
image(x, y, z)
filled.contour(x,y,z) # show the color legend

#######################################
# 2.11  the contour plot               #
#######################################
# first we get some data to use...
N <- 50
x <-y<- seq(-1, 1, length=N)
xx <- matrix(x, nrow=N, ncol=N)
yy<-t(xx)
z <- 1 / (1 + xx^2 + (yy + .5 * sin(5*yy))^2)
# here is the contour plot
contour(x, y, z, main = "Contour plot")
image(x, y, z)
#######################################
# 2.12  3D  interactive               #
#######################################
install.packages("rgl")
library("rgl")

demo("rgl")


## Non-interactive
X11()
png(filename = "my2x2png.png")

## plotting
par(mfrow = c(1, 2))
plot(rnorm(10), col = "blue")
plot(rt(10, df = 3), col = "red")

## Close the device
dev.off()
