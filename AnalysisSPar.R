library(gridExtra)
library(reshape2)
library(scales)
library(plyr)

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

#function to normalize based on reference
#assumes log transformed
#requires a dataset and dataframe associating each reference with a sample group
PoolRefNorm <- function(indata, alias, ref = "ref")
{
  outdata <- indata
  for(i in unique(alias$Sample))
  {
    refvar <- alias[alias$Sample == i,]$Alias
    refs <- grep(ref, refvar)
    refmean <- apply(outdata[,refvar[refs]], 1, mean)
    for( j in refvar)
    {
      outdata[,j] <- outdata[,j] - refmean
    }
  }
  return(outdata)
}

#Writes out data.frame tables as they are displayed
IOWriter <- function(TableToWrite, separator, Filename)
{
  
  out <- file(Filename, "w")
  header <- paste(c("rowname", names(TableToWrite)), collapse = separator)
  cat(header, "\n", file = out)
  write.table(TableToWrite, out, sep = separator, col.names=FALSE,row.names = TRUE)
  close(out)
} 

#Runs analysis, TODO: Needs to be further parameterized for reuse
# RunAnalysis <- function(T_RowC, T_PepQuant,T_Alias, outdir)
# {
#   PairPlot(T_PepQuant,outdir, "Scatter")
#   T_RowQ <- subset(T_RowC, T_RowC$Peptide %in% rownames(T_PepQuant))
#   T_PepNorm <- T_PepQuant
#   NormPlot(log(T_PepQuant,2), outdir, filename = "Original")
#   T_PepNorm <- MedNorm(T_PepNorm)
#   T_PepNorm[,1:ncol(T_PepNorm)] <- log(T_PepNorm[,1:ncol(T_PepNorm)],2)
#   NormPlot(T_PepNorm, outdir, "Unbiased")
#   T_Melt <- PepProtMelter(T_PepNorm, T_RowQ)
#   T_Melt <- merge(T_Alias,T_Melt, by.x = "Alias", by.y="variable", all = FALSE)
#   T_Model <- GetProteinEffects(as.formula("value ~ Treatment"), as.formula("~ 1 | Peptide"), data = T_Melt)
#   ptable <- getPValueTable(T_Model)
#   ptable$padj <- p.adjust(ptable$'p-value')
#   plist <- sqldf("Select Protein from ptable WHERE padj < 0.05")
#   VolPlot(ptable$Value, ptable$padj, outdir, filename = "Volcano")
#   #test <- printTukeyTest(T_Model, printGraph = TRUE, mcpvar = "Treatment", output = outdir)
#   
#   outlist <- list("NormData" = T_PepNorm, "AnovaModel" = T_Model, "PTable" = ptable, "Tukey" = NULL, "descTable" = resultdesc(ptable))
#   return(outlist)
# }



#Median central tendency normalization
MedNorm <- function(mydf)
{
  # mydf <- mydf[apply(mydf, 1, function(x){all(x != 0) && all(!is.na(x))}),]
  mydf[,1:ncol(mydf)] <- apply(mydf, 2, function(x){x / median(x, na.rm = TRUE)})
  return(mydf)
}

#Gets the results of applying LME to a set of data


#Melts the dataset and connects the row metadata for later 
#processing
PepProtMelter <- function(T_Data, T_Row)
{
  require(reshape2)
  T_Data <- cbind(T_Data, Peptide = rownames(T_Data))
  tempData <- merge(T_Data, T_Row, by.x = "Peptide", by.y = "Peptide")
  melt_T_Data <- melt(tempData, id = c("Protein", "Peptide","MultiProtein"))
  return(melt_T_Data)
}



#Applys the LMEobject method to generate the protein effect info
GetProteinEffects <- function(model, rand, datasub)
{ 
  require(sqldf)
  require(snowfall)
  
  proteinList <- sqldf("SELECT Protein FROM datasub GROUP BY Protein")
  sfInit(parallel = TRUE, cpus = 3, type = "SOCK")
  
  GetLMEobj <- function(sprot, rand, model)
  {
    datasub <- subset(datasub, datasub$Protein %in% sprot)
    datasub$Peptide <- factor(datasub$Peptide)
    require(nlme)
    fit1 <- NULL
    rs <- NULL
    if(length(unique(datasub$Peptide))>1)
    {
      rs <- tryCatch(fit1 <- lme(model, random = rand, data = datasub, 
                                 na.action = "na.omit"), error = function(e) NULL)
    }
    if(is.null(rs))
    {
      return(NULL)
    }
    return(fit1)
  }
  sfLibrary(nlme)
  result <- sfLapply_pb(proteinList$Protein, GetLMEobj, rand , model)
  sfStop()
  names(result) <- proteinList$Protein
  result <- result[!sapply(result, is.null)]
  return(result)

}


sfLapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal + 1, envir = env)
    setTxtProgressBar(get("pb", envir = env), curVal + 1)
    FUN(...)
  }
  res <- sfLapply(X, wrapper, ...)
  close(pb)
}


#Generates a pvalue table using results of GetProteinEffects
getPValueTable <- function(data)
{
  mydf <- data.frame(matrix(ncol =4, nrow = length(data)))
  colnames(mydf) <- c("Protein", "Value", "t-value", "p-value")
  count = 1
  for(i in names(data))
  {
    result <- summary(data[[i]])$tTable
    mydf[count,1] <- i
    mydf[count,2] <- result[,"Value"][2]
    mydf[count,3] <- result[,'t-value'][2]
    mydf[count,4] <- result[,'p-value'][2]
    count = count+1  
  }
  return(mydf)
}





#Applies Tukey test and outputs a set of results, generates plots if necessary for every test
printTukeyTest <- function(alist, printGraph = FALSE, proteins = NULL, mcpvar = "variable", output = "img")
{
  if(require(multcomp))
  {
    detach("package:multcomp", unload = TRUE)
  }
  library(multcomp)
  myargs <- list("Tukey")
  names(myargs) <- mcpvar
  cmp <- do.call(mcp, myargs)
  pb <- txtProgressBar(min = 0, max = length(alist), style = 3)
  tukeyobj <- summary(glht(alist[[1]], linfct = cmp))
  tukeyNames <- names(tukeyobj$test$coefficients)
  tmelt <- data.frame(Protein = NA, Comparison = NA, mean = NA, upper95Conf = NA, lower95Conf = NA, pval = NA)
  pcount = 0
  options(warn = -1)
  mydir <- getwd()
  if(printGraph)
  {
    dir.create(file.path(mydir, output), showWarnings = FALSE)
    setwd(file.path(file.path(mydir, output)))
  }
  if(is.null(proteins))
  {
    proteins = names(alist)
  }
  
  for(i in proteins)
  {
    pcount = pcount + 1
    setTxtProgressBar(pb, pcount)
    tt <- summary(glht(alist[[i]], linfct = cmp))
    pvalCount = 1
    for(j in tukeyNames)
    {
      
      tmelt <- rbind(tmelt, c(i, j, tt$test$coef[j][[1]], tt$test$coef[j][[1]]-tt$test$sigma[j][[1]],
                              tt$test$coef[j][[1]]+tt$test$sigma[j][[1]],  tt$test$pvalues[[pvalCount]]))
      pvalCount =pvalCount + 1
    }
    if(printGraph)
    {
      png(file = paste0(getwd(), "/", sub("|", "", i,fixed = TRUE),"_Tukey.png"), width = 640)
      par(mar =c(4, 20, 4, 4))
      plot(tt)
      dev.off()
    }
  }
  if(printGraph)
  {
    setwd(mydir)
  }
  return(tmelt)
}

#Generates a plot to show the normalization effect, breaks glht for some reason
NormPlot <- function (T_Data, output, filename = "default") {
  originald <- getwd()
  dir.create(file.path(originald, output), showWarnings = FALSE)
  setwd(file.path(file.path(originald, output)))
  
  png(file = paste0(getwd(), "/", filename, ".png"), height = 1000,  width = 640)
  T_Data.melt <- melt(T_Data)
  T_Data.melt[T_Data.melt == -Inf] <- NA
  require(lattice)
  p1 <- histogram( ~ value | variable, T_Data.melt, ylab = "Abundance", xlab = "Sample")
  
  p2 <- bwplot(value ~ variable, T_Data.melt, scales=list(x=list(rot=45)), ylab = "Abundance", xlab = "Sample")
  
  grid.arrange(p1, p2, ncol = 1, nrow = 2)
  dev.off()
  setwd(originald)
}

#prints out teh results of applying an anova test
PrintAnovaResults <- function(aovObjects, modelParams)
{
  result_Table <- matrix(ncol = length(modelParams) * 2, nrow = length(aovObjects))
  rownames(result_Table) <-  names(aovObjects) 
  result_cols <- c()
  for(i in modelParams)
  {
    result_cols <- c(result_cols, paste0(i, "_F value"), paste0(i, "_Pr(>F)"))
  }
  colnames(result_Table) <- result_cols    
  for(i in names(aovObjects))
  {  
    pepstats <- summary(aovObjects[[i]])[[1]]
    for(k in modelParams)
    {
      result_cols2 <- c(paste0(k, "_F value"), paste0(k, "_Pr(>F)"))
      if(c("F value", "Pr(>F)") %in% colnames(pepstats))
      {
        result_Table[i,result_cols2] <- c(pepstats[k,]$'F value', pepstats[k,]$'Pr(>F)')
      }
    }
    
  }
  return(result_Table)
}

#Generates a neat little heat map
QuickHeatMap <- function(mydata)
{
  require(ggplot2)
  data.m <- ddply(mydata, .(Peptide), transform, rescale = rescale(value))
  p <- ggplot(data.m, aes(Alias, Peptide)) + geom_tile(aes(fill = rescale),colour = "white")
  p + scale_fill_gradient(low = "white", high = "steelblue")
}

#Generates a nicely formatted volcano plot
VolPlot <- function(inDatax, inDatay, output, maxy = NULL, filename = "default")
{
  originald <- getwd()
  dir.create(file.path(originald, output), showWarnings = FALSE)
  setwd(file.path(file.path(originald, output)))
  
  png(file = paste0(getwd(), "/", filename, ".png"), width = 640)
  if(is.null(maxy))
  {
    maxy = rev(range(inDatay))[2]  
  }
  maxx <- ceiling(max(abs(inDatax)))
  plot(inDatax, inDatay,  log = "y", ylim = c(1, maxy), xlim = c(-maxx, maxx))
  grid(lwd = 3)
  abline(h = 0.05, col = "red")
  dev.off()
  setwd(originald)
}


PairPlot <- function(inData, output, filename = "default")
{
  originald <- getwd()
  dir.create(file.path(originald, output), showWarnings = FALSE)
  setwd(file.path(file.path(originald, output)))
  inData[inData == 0] <- NA
  inData <- log(inData,2)
  png(file = paste0(getwd(), "/", filename, ".png"), width = 640)
  pairs(inData, upper.panel = panel.cor)
  dev.off()
  setwd(originald)
}

#For use with the pairs function
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


resultdesc <- function(myPTable)
{
  library(sqldf)
  result <- data.frame(Identifications = c("Total", "Significant", "Up", "Down"), Amount = NA)
  result[1,2] = sqldf("SELECT COUNT(Protein) FROM myPTable")
  result[2,2] = sqldf("SELECT COUNT(Protein) FROM myPTable WHERE padj < 0.05")
  result[3,2] = sqldf("SELECT COUNT(Protein) FROM myPTable WHERE padj < 0.05 AND Value > 0")
  result[4,2] = sqldf("SELECT COUNT(Protein) FROM myPTable WHERE padj < 0.05 AND Value < 0")
  return(result)
  
}


PNGPlots <- function(inPlot, filename)
{
  #To do: put functionality for quick plot save
}