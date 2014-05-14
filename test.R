mymelter <- function(T_Meltees, T_RowQ, T_Alias)
{

  T_Meltees <- llply(T_Meltees, failwith(NULL, function(x){
    x <- PepProtMelter(x, T_RowQ)
    x <- merge(T_Alias, x, by.x = "Alias", by.y="variable", all = FALSE)
    x <- na.omit(x)
   return(x)
  }, quiet = FALSE), .progress = "tk")
  return(T_Meltees)
}


