
library(gridExtra)
library(reshape2)
library(scales)
library(plyr)
library(fastmatch)

#######################   Normalize    ##############################

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
    refmean <- apply(as.matrix(outdata[,refvar[refs]]), 1, mean)
    for( j in refvar)
    {
      outdata[,j] <- outdata[,j] / refmean
    }
  }
  return(outdata)
}

CorrectZeros <- function(T_Data, T_Alias)
{
  kk <- as.character(unique(T_Alias$Sample))
  for(i in kk)
  { 
    xx <- c(as.character(T_Alias[T_Alias$Sample == i, "Alias"]))
    T_Data[which(apply(T_Data[,xx], 1, function(x) {sum(x, na.rm = TRUE)==0})),xx] <- NA
  }
  return(T_Data)
}

#No longer used.  Normalize the reporter ions to the pooled reference
PoolNormalize <- function (T_In,poolColumn, colsToNormalize) {
  for(i in 1:nrow(T_PepNorm))
  {
    for(j in colsToNormalize)
    {
      T_In[i,j] = T_In[i, j] / T_In[i,poolColumn]
    }
  }
  return(T_In[,colsToNormalize])
}

#Median central tendency normalization
MedNorm <- function(mydf)
{
 # mydf <- mydf[apply(mydf, 1, function(x){all(x != 0) && all(!is.na(x))}),]
  mydf[,1:ncol(mydf)] <- apply(mydf, 2, function(x){x / median(x, na.rm = TRUE)})
  return(mydf)
}

ProteinQuant <- function (Tdata, Trow) {
  T_PepTemp <- Tdata
  T_PepTemp[sapply(T_PepTemp, is.infinite)] <- NA
  T_PepTemp$Peptide <- rownames(T_PepTemp)
  T_PepTemp <- merge(T_PepTemp, Trow, by.x = "Peptide", by.y = "Peptide")
  T_ProtTab <- aggregate(2^T_PepTemp[,2:(ncol(T_PepTemp)-2)], list(Protein = T_PepTemp$Protein), sum, na.rm = TRUE)
  T_ProtTab <- data.frame(log(T_ProtTab[,2:(ncol(T_PepTemp)-2)],2), row.names = T_ProtTab$Protein)
  T_ProtTab[sapply(T_ProtTab, is.infinite)] <- NA
  return(T_ProtTab)
}


#Melts the dataset and connects the row metadata for later 
#processing
PepProtMelter <- function(T_Data, T_Row)
{
  require(reshape2)
  T_Data <- cbind(T_Data, Peptide = rownames(T_Data))
  tempData <- merge(T_Data, T_Row, by.x = "Peptide", by.y = "Peptide")
  melt_T_Data <- melt(tempData, id = c("Protein", "Peptide","MultiProtein"))
  melt_T_Data[is.infinite(melt_T_Data$value), "value"] <-NA
  return(melt_T_Data)
}

####################   LME  #################################### 

#Applys the LMEobject method to generate the protein effect info
GetProteinEffectsLME <- function(model, rand, data)
{ 
  require(sqldf)
  proteinList <- sqldf("SELECT Protein FROM data GROUP BY Protein")
  resultList <- list()
  pb <- txtProgressBar(min = 0, max = nrow(proteinList), style = 3)
  pcount <- 0
  for(i in proteinList$Protein)
  {  
    setTxtProgressBar(pb, pcount) 
    pcount = pcount + 1
    tempData <- subset(data, data$Protein %in% i)
    tempData$Peptide <- factor(tempData$Peptide)
    fit <- NULL
    if(length(unique(tempData$Peptide))>1)
    {
      fit <- GetLMEobj(tempData, rand, model)
    }
    if(!is.null(fit))
    {
      resultList[[i]] <- fit
    }
  }
  return(resultList)
}


#Gets the results of applying LME to a set of data
GetLMEobj <- function(datasub, rand, model)
{
  require(nlme)
  rs <- tryCatch(fit1 <- lme(model, random = rand, data = datasub, 
                             na.action = "na.omit"), error = function(e) NULL)
  if(is.null(rs))
  {
    return(NULL)
  }
  return(fit1)
}

############################   LMER   ##################################

#Gets the results of applying LMER to a set of data
GetLMERobj <- function(datasub, model)
{
  require(lme4)
  rs <- tryCatch(fit1 <- lmer(model, data = datasub, 
                             na.action = "na.omit"), error = function(e) NULL)
  if(is.null(rs))
  {
    return(NULL)
  }
  return(fit1)
}

PrintLMERAnovaResult <- function(lmmodlist)
{
  require(car)
  outlist <- lapply(lmmodlist, Anova)
  return(outlist)
}

GetLMERAnovaTable <- function(lmanovlist, infactors)
{
  outframe <- data.frame(Protein = rep(names(lmanovlist), each = length(infactors)), 
                         Type = rep(infactors, each = length(lmanovlist)), "Pr(>Chisq)" = NA, check.names = FALSE)
  for(i in names(lmanovlist))
  {
    for(j in infactors)
    {
      outframe[outframe$Protein == i & outframe$Type == j, 'Pr(>Chisq)'] <- lmanovlist[[i]][j, 'Pr(>Chisq)']
    }
  }
  return(outframe)
}



GetProteinEffectsLMER <- function(model, data)
{ 
  require(sqldf)
  proteinList <- sqldf("SELECT Protein FROM data GROUP BY Protein")
  resultList <- list()
  pb <- txtProgressBar(min = 0, max = nrow(proteinList), style = 3)
  pcount <- 0
  for(i in proteinList$Protein)
  {  
    setTxtProgressBar(pb, pcount) 
    pcount = pcount + 1
    tempData <- subset(data, data$Protein %in% i)
    tempData$Peptide <- factor(tempData$Peptide)
    fit <- NULL
    if(length(unique(tempData$Peptide))>1)
    {
      fit <- GetLMERobj(tempData, model)
    }
    if(!is.null(fit))
    {
      resultList[[i]] <- fit
    }
  }
  return(resultList)
}

#########################  PVALUEs!!!!    ######################

#Generates a pvalue table using results of GetProteinEffects
getPValueTable <- function(data)
{
  mydf <- data.frame(matrix(ncol =5, nrow = length(data)))
  colnames(mydf) <- c("Protein", "Variable","Value", "t-value", "p-value")
  count = 1
  for(i in names(data))
  {
    result <- summary(data[[i]])$tTable
    for(j in rownames(result)){
      mydf[count,1] <- i
      mydf[count,2] <- j
      mydf[count,3] <- result[j,"Value"]
      mydf[count,4] <- result[j,'t-value']
      mydf[count,5] <- result[j,'p-value']
      count = count+1  
    }
  }
  return(mydf)
}

#Benjamini Hochberg p-value correction for multiple testing error
AdjustPvalues <- function(pvalueframe, infactors)
{
  pvalueframe$padj <- NA
  for(i in infactors)
  {
    pvalueframe[pvalueframe$Type == i,"padj"] <- 
      p.adjust(pvalueframe[pvalueframe$Type == i,'Pr(>Chisq)'])
  }
  return(pvalueframe)
}

###########################  TUKEY   ####################################

#performs a list application of glht to a single comparison variable. 
# outputs list
ApplyTukeyTest <- function(inlist, proteins = NULL, mcpvar = "variable")
{
  require(multcomp)
  myargs <- list("Tukey")
  names(myargs) <- mcpvar
  cmp <- do.call(mcp, myargs)
  outlist <- lapply_pb(inlist[proteins], glht, linfct = cmp)
  return(outlist)
  
}









