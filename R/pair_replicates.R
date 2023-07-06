#' Teach R to recognize replicates of independent samples in Genome-wide 3' RACE
#' @description pair_replicates() creates list of two character vectors, first contains names of used strains, second contains expressions that can be used to grep the said strains from the list of data frames
#' @param data A list that is an output of gw3RACE::data_fix()
#' @param grep1 an argument passed to gsub(), regular expression that can be used to erase the indication of replication of experiment from name of data table stored on the list. By default it is "_clone\\\\d"
#'
#' @return List of vectors
#' @export
#'
#' @examples strainslist <- pair_replicates(dataList)
pair_replicates <- function(data, grep1 = "_clone\\d"){
  strnames <- c()
  for (a in 1:length(data)) {
    temp <- gsub(grep1, "", names(dataList)[a])
    strnames <- c(strnames, temp)
  }
  strnames <- unique(strnames)
  namestogrep <- c()
  for (a in strnames) {
    namestogrep <- append(namestogrep, gsub("(.*)","^\\1_",a))
  }
  mylist <- list(strnames, namestogrep)
  return(mylist)
}
