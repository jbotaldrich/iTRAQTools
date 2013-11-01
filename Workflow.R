setwd("C:/Users/aldr699/Documents/2013/Dworkin")

#Load necessary functions
source("R/AnalysisSource.R", echo = FALSE)
source("R/Utilities.R", echo = FALSE)
source("R/PlotUtils.R", echo = FALSE)
source("R/TableGeneration.R", echo = FALSE)

#Read in data
T_filtered <- read.delim("Data/T_Filtered_Results.txt")
T_Data <- read.delim("Data/T_Data.txt", check.names=FALSE)
T_Data <- data.frame(T_Data[,2:ncol(T_Data)], row.names = T_Data$Peptide, check.names = FALSE)
T_Row_Metadata <- read.delim("Data/T_Row_Metadata.txt")
t_alias <- read.delim("Data/t_alias.txt")

isobar = 4

for(i in 1:(ncol(T_Data)/isobar))
{
  
  startp <- (i-1)*isobar+1
  stopp <- startp + (isobar-1)
  T_Data[apply(T_Data[,startp:stopp], 1, sum, na.rm = TRUE)==0,startp:stopp] <- NA
}



T_RowC <- T_Row_Metadata
T_PepQuant <- T_Data
T_Alias <- t_alias
outdir_ext = "All"

figpath = paste("Fig", outdir_ext, sep="/")

#Correlations plots of all and each samples group separately with correlation 
#values.
PairPlot(T_PepQuant,figpath, "PepScatter")
PairPlot(T_PepQuant[,1:4],figpath, "PepScatterS1")
PairPlot(T_PepQuant[,5:8],figpath, "PepScatterS2")


#Normalize and plot the normalization steps of the data
T_RowQ <- subset(T_RowC, T_RowC$Peptide %in% rownames(T_PepQuant))
T_PepNorm <- T_PepQuant
NormPlot(log(T_PepQuant,2), figpath, filename = "Original")
T_PepNorm <- MedNorm(T_PepNorm)#normalizae by median central tendency
NormPlot(log(T_PepNorm[,1:ncol(T_PepNorm)],2), figpath, "Unbiased")
T_PepNorm[,1:ncol(T_PepNorm)] <- log(T_PepNorm[,1:ncol(T_PepNorm)],2)


#Get rolled up protein values and plot
ProteinTable <- ProteinQuant(T_PepNorm, T_RowQ)
PairPlot(2^ProteinTable, figpath, "ProtScatter")

#Unpivot our table and merge to factor table for linear model
T_Melt <- PepProtMelter(T_PepNorm, T_RowQ)
T_Melt <- merge(T_Alias,T_Melt, by.x = "Alias", by.y="variable", all = FALSE)

#Use lmer to fit peptides to the model below.  Fixed Effect for treatment random
#effect for Peptide
T_Model <- GetProteinEffectsLMER(
  as.formula("value ~ Treatment + ( 1 | Peptide)"), data = T_Melt)

#Apply anova to the model and determine significance of treatment effect
#spit out some tables and adjust for multiple comparisons
T_Model.Anova <- PrintLMERAnovaResult(T_Model)
result.anovatable <- GetLMERAnovaTable(T_Model.Anova, "Treatment")
result.anovatable <- AdjustPvalues(result.anovatable, "Treatment")

#collect proteins with significant treatment effects
result.sigprotein <- result.anovatable[which(result.anovatable$padj < 0.05),]

#Apply tukeys honestly signficant difference test export plot results
test <- ApplyTukeyTest(T_Model, result.sigprotein$Protein, mcpvar = "Treatment")
PlotTukeyResult(test, output = figpath)
test.table <- GenTukeyTable(test)

#Volcano plots!
VolPlot(test.table$mean,test.table$pval, figpath, filename = "Volcano")

#Prepare tables for export
crtab <- na.omit(merge(T_RowC, T_PepNorm, by.x = "Peptide", by.y = "row.names"))
test.table$pval <- as.numeric(test.table$pval)
sigtable <- test.table[as.numeric(test.table$pval) < 0.05,]
sigtable <- sigtable[order(sigtable$pval),]

#Throw all the tables into a single data structure and send it to export files
outlist <- list("NormPeptideData" = crtab, "AnovaModel" = T_Model, 
                "Tukey" = test.table, 
                "SignificanceCountTable" = resultdescTuk(test.table, 
                                                         result.anovatable,T_RowQ), 
                "ProteinQuantTable" = ProteinTable,
                "PeptideIdentTable" = peptideresults_desc(T_filtered), 
                "ProteinAnovaTable" = result.anovatable,
                "SignificantComparisonTable" = sigtable)
WriteTables(outlist, paste("output", outdir_ext, sep="/"))
  
