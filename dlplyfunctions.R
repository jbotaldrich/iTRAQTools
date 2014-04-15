get_coeffs <- function(mylev)
{
  return(full_coeffs <- function (fit) {
    myfixed <- data.frame(matrix(NA, 1, length(mylev)))
    colnames(myfixed) <- mylev
    myfixed[,names(fixef(fit))] <- fixef(fit) 
    return(myfixed)
  })
  
}

getconfint <- function(myconfs)
{
  return(full_conf <- function(fit) {
    myfixed <- data.frame(matrix(NA, 1, length(myconfs)),check.names = FALSE)
    colnames(myfixed) <- myconfs
    temp <- confint(fit, method = "Wald")
    temp <- melt(temp)
    temp2 <- temp$value
    names(temp2) <- paste(temp[,1], temp[,2], sep = "_")
    myfixed[,names(temp2)] <- temp2
    return(myfixed)
  })
}

getAnova <- function(df)
{
  xx <- as.data.frame(Anova(df))
  xx$Variable <- rownames(xx)
  return(xx)
}

Heater <- function(amelt, aprot){
  xx <- dcast(amelt[amelt$Protein == aprot,], Peptide ~ Alias)
  rownames(xx) <- xx[,1]
  xx <- as.matrix(xx[,-1])
  heatmap.2(xx, trace = "none", na.rm = TRUE, distfun = na.dist, scale = "row", dendrogram = "row", Colv = FALSE)
  
}