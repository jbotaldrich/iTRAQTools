

projtemplate <- function(directory)
{
  setwd(directory)
  
  dir.create(file.path(getwd(), "Fig"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "Data"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "Output"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "R"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "Doc"), showWarnings = FALSE)
}
