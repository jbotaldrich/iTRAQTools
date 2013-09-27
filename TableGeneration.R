

#Generates the table output from the tukey test
GenTukeyTable <- function(tukeylist, CL = FALSE)
{
  if(CL)
  {
    tresult <- data.frame(Protein = NA, Comparison = NA, mean = NA, StdError = NA, lwr = NA,  upr = NA, pval = NA)
  }
  else
  {
    tresult <- data.frame(Protein = NA, Comparison = NA, mean = NA, StdError = NA, zvalue = NA,   pval = NA)
  }
  temp <- lapply_pb(tukeylist, summary)
  tukeyNames <- names(temp[[1]]$test$coefficients)
  tester = 0
  for(i in names(temp))
  {
    tukeyNames <- names(temp[[i]]$test$coefficients)
    pcount = 1
    for(j in tukeyNames)
    {
      if(CL){
        tt <- confint(temp[[i]])
        tresult <- rbind(tresult, c(i,j, temp[[i]]$test$coef[j][[1]],
          temp[[i]]$test$sigma[j][[1]],
          tt$confint[pcount, 2],
          tt$confint[pcount, 3],
          temp[[i]]$test$pvalues[[pcount]]))
      }
      else{
        tresult <- rbind(tresult, c(i, j, temp[[i]]$test$coef[j][[1]], 
                                    temp[[i]]$test$sigma[j][[1]],
                                    temp[[i]]$test$tstat[j][[1]],
                                    temp[[i]]$test$pvalues[[pcount]]))
      }
      pcount= pcount+1
    }
    
  }
  tresult <- na.omit(tresult)
  return(tresult)
}



#prints out the results of applying an anova test
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
    pepstats <- anova(aovObjects[[i]])
    for(k in modelParams)
    {
      result_cols2 <- c(paste0(k, "_F value"), paste0(k, "_Pr(>F)"))
      if(result_cols2 %in% colnames(result_Table))
      {
        result_Table[i,result_cols2] <- c(pepstats[k,]$'F-value', 
                                          pepstats[k,]$'p-value')
      }
    }
    
  }
  return(result_Table)
}


#summary of peptide identifications
peptideresults_desc<-function(peptable)
{
  require(sqldf)
  result <- data.frame(Identifications = c("Peptide Spectrum Matches", 
                                           "Unique Peptides"), 
                       Amount=NA)
  result[1,2] <- nrow(peptable)
  result[2,2] <- length(as.character(unique(peptable$Peptide)))
  return(result)
}

#summary of significant expression
resultdesc <- function(myPTable)
{
  require(sqldf)
  result <- data.frame(Identifications = 
                         c("Total", "Significant", "Up", "Down"), 
                       Amount = NA)
  result[1,2] = sqldf("SELECT COUNT(Protein) FROM myPTable")
  result[2,2] = sqldf("SELECT COUNT(Protein) FROM myPTable WHERE padj < 0.05")
  result[3,2] = sqldf("SELECT COUNT(Protein) FROM myPTable WHERE padj < 0.05 AND Value > 0")
  result[4,2] = sqldf("SELECT COUNT(Protein) FROM myPTable WHERE padj < 0.05 AND Value < 0")
  return(result)
  
}


#summary of identified and significant peptides.
resultdescTuk <- function(myPTable,myATable, trowmet)
{
  require(sqldf)
  result <- data.frame(Identifications = c("Total", "Significant Anova","Significant TukeyHSD"), Amount = NA)
  protcheck <- sqldf("SELECT Protein,  Count(Peptide) FROM trowmet GROUP BY PROTEIN")
  result[1,2] = length(which(protcheck$'Count(Peptide)' > 1))
  result[2,2] = sqldf("SELECT COUNT(Protein) FROM myATable WHERE padj < 0.05")
  result[3,2] = sqldf("SELECT COUNT(Protein) FROM myPTable WHERE pval < 0.05")
  return(result)
  
}
