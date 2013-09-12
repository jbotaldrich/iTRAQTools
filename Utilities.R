
#List apply with a progressbar courtesy of 
#http://ryouready.wordpress.com/2010/01/11/progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/
lapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}

PNGPlots <- function(inPlot, filename)
{
  #To do: put functionality for quick plot save
}


#Writes out data.frame tables as they are displayed
IOWriter <- function(TableToWrite, separator, Filename)
{
  
  out <- file(Filename, "w")
  header <- paste(c("rowname", names(TableToWrite)), collapse = separator)
  cat(header, "\n", file = out)
  write.table(TableToWrite, out, sep = separator, col.names=FALSE,row.names = TRUE, quote = FALSE)
  close(out)
} 

#Writes out result tables
WriteTables <- function(resulttoWrite, dir)
{
  original <- getwd()
  dir.create(file.path(original, dir), showWarnings = FALSE)
  setwd(file.path(original, dir))
  for(i in names(resulttoWrite))
  {
    if(i != "AnovaModel")
    {
      IOWriter(resulttoWrite[[i]], sep = '\t', paste0(i, ".txt"))
    }
  }
  setwd(original)
}