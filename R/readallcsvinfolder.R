#' Upload all .csv files in one folder to R list
#' @description readallcsvinfolder() will read all .csv files in given folder and create list of them
#' @param folder A path to the folder where you store your .csv files, by default it is your current working directory
#' @param columnstokeep A vector of column names that is passed to subset function, by default it's = c("gene", "coord_R2", "RNA_type", "tail_from", "tail_len", "tail_type", "stop_R2", "distance_to_TES", "rel_distance_to_TES")
#' @param RNAtype A character or a vector of characters describing RNA types, that you want to keep in your data, must be identical to names used in RNA_type column in original file, by default it's "mRNA"
#'
#' @return List of data frames
#' @export use_pipe(export = TRUE)
#'
#' @examples dataList <- readallcsvinfolder()
readallcsvinfolder <- function(folder = getwd(), columnstokeep = c("gene", "start_R1", "stop_R1", "RNA_type", "tail_from", "tail_len", "tail_type", "stop_R2", "distance_to_TES", "rel_distance_to_TES", "UMI"), RNAtype = "mRNA"){
  setwd(folder)
  c0 <- list.files(getwd())
  c0 <- grep(".csv", c0, value = T)
  c0 <- grep(".*_clone.*.csv", c0, value = T)
  c1 <- gsub("(_clone\\d).*.csv", "\\1.csv", c0)
  c2 <- gsub("delta", "", c1)
  c3 <- gsub(".csv", "", c2)
  dflist <- list()
  options(datatable.fread.datatable=FALSE)
  for (a in 1:length(c0)) {
    temp <- data.table::fread(c0[a])
    temp <- subset(temp, select = columnstokeep)
    temp <- temp %>% dplyr::filter (RNA_type == RNAtype)
    gc()
    dflist[c3[a]] <- list(temp)
  }
  return(dflist)
}
