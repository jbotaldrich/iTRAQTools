

projtemplate <- function(directory)
{
  setwd(directory)
  
  dir.create(file.path(getwd(), "Fig"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "Data"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "Output"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "R"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "Doc"), showWarnings = FALSE)
  download.file("https://raw.github.com/jbotaldrich/iTRAQTools/master/PlotUtils.R", file.path(directory, "R/Plotutils.R"))
  download.file("https://raw.github.com/jbotaldrich/iTRAQTools/master/AnalysisSource.R", file.path(getwd(), "R/AnalysisSource.R"))
  download.file("https://raw.github.com/jbotaldrich/iTRAQTools/master/TableGeneration.R", file.path(directory, "R/TableGeneration.R"))
  download.file("https://raw.github.com/jbotaldrich/iTRAQTools/master/Utilities.R", file.path(directory, "R/Utilities.R"))
  
}




