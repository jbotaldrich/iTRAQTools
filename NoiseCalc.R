#Variance function estimation

# xx <- runif(10000, 1, 100000)
# yy <- rnorm(10000, 0, sd = exp(-0.9 * log(xx,10)))
# plot(xx,yy, pch = ".")
# xx <- yy^10


NoiseCalc <- function(inData, inbase = 2)
{
  inData <- log(inData,inbase)
  repratio <- inData[,1]-inData[,2]
  repaverage <- 0.5 * (inData[,1] + inData[,2])
  
  
  myout <- nlminb(c(0,10,1),minfunc, lower = rep(1e-10, 3), 
                  x = repratio , y = repaverage)
  return(myout$par)
}


varfunc <- function(mypar, x)
{
  return(mypar[1] + mypar[2]*exp(-mypar[3]* x))
}

minfunc <- function(mypar, x,y){
  return(-sum(dnorm(x, mean = 0, sd = varfunc(mypar,y),log = TRUE)))
}


PlotNoise <- function(inData, mypar, inbase = 2){
  inData <- log(inData,inbase)
  repratio <- inData[,1]-inData[,2]
  repaverage <- 0.5 * (inData[,1] + inData[,2])
  repend <- max(repaverage) 
  repstart <- min(repaverage)
  myp <- seq(repstart, repend, by =0.05)
  
  plot(repaverage, repratio, pch = ".")
  lines(myp, varfunc(mypar,myp), col = "red", pch = ".")
  lines(myp, -varfunc(mypar,myp), col = "red", pch = ".")
}

