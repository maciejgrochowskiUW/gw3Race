#' Cleaning of the Genome-wide 3' RACE data
#' @description data_fix() clears the data frame of reads that couldn't possibly exist and are some artifacts of initial analysis
#' @param data A data frame stored on the list that is the output of gw3RACE::readallcsvinfolder() function, generally it is recommended to use this function with lapply()
#' @param excludedgenes A vector with names of genes that are supposed to be excluded from the analysis
#'
#' @return Data frame
#' @export
#'
#' @examples dataList <- lapply(dataList, data_fix)
data_fix <- function(data, excludedgenes = c()) {

  data <- data[!(data$tail_type == "jakis_bias"),]

  data <- data[!(data$tail_from == "grep" & data$tail_len < 10),]

  data <- data[!(data$gene %in% excludedgenes),]

  data[data$coord_R2 == 0,]$stop_R2 <- NA

  data[data$coord_R2 == 0,]$distance_to_TES <- NA

  data[data$coord_R2 == 0,]$rel_distance_to_TES <- NA


  return(data)

}
