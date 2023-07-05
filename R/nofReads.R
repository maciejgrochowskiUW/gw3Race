#' Plot of number of reads/genes in sample
#' @description Plot number of reads in data frames stored on list with gw3RACE::nofReads()
#' @param listofcsv A list of data frames, output of gw3RACE::data_fix() or gw3RACE::tail_summary()
#'
#' @return Barplot
#' @export
#'
#' @examples nofReads(dataList)
nofReads <- function(listofcsv){
  numberofreads <- numeric()
  for (n in c(1:length(listofcsv))) {
    numberofreads <- c(numberofreads,length(listofcsv[[n]][,1]))
  }
  ylabel <- numberofreads
  xlabel <- names(listofcsv)
  leftmar <- max(strwidth(ylabel, "inch") + 0.3, na.rm = TRUE)
  bottommar <- max(strwidth(xlabel, "inch") + 0.3, na.rm = TRUE)
  par(mai=c(bottommar, leftmar, 0.1, 0.1))
  barplot(numberofreads, names.arg = xlabel, las = 2)
}
