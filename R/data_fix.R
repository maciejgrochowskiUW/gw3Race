#' Cleaning of the Genome-wide 3' RACE data
#' @description data_fix() clears the data frame of reads that couldn't possibly exist and are some artifacts of initial analysis
#' @param data A data frame stored on the list that is the output of gw3RACE::readallcsvinfolder() function, generally it is recommended to use this function with lapply()
#' @param excludedgenes A vector with names of genes that are supposed to be excluded from the analysis
#' @param mtgenes A vector of grep queries that will reassign genes that originally are not assigned to mitochondrial class but should be (ex. mt tRNA and mt rRNA)
#'
#' @return Data frame
#' @export
#'
#' @examples dataList <- lapply(dataList, data_fix)
data_fix <- function(data, excludedgenes = c(), mtgenesgrep = c()) {

  data <- data %>% distinct(start_R1, stop_R1, UMI, .keep_all = TRUE)

  data <- data[!(data$tail_from == "grep" & data$tail_len < 10),]

  data <- data[!(data$gene %in% excludedgenes),]


  if("coord_R2" %in% colnames(data)) {
    data[data$start_R1 == 0,]$stop_R2 <- NA

    data[data$start_R1 == 0,]$distance_to_TES <- NA

    data[data$start_R1 == 0,]$rel_distance_to_TES <- NA
  }

  if(length(data[data$RNA_type == "other_type",]$distance_to_TES) > 0) {
    data[data$RNA_type == "other_type",]$distance_to_TES <- NA
  }

  for (query in mtgenesgrep) {

    if(length(data[grep(query, data$gene),]$RNA_type) > 0) {
      data[grep(query, data$gene),]$RNA_type <- "mito_RNA"}

  }



  return(data)

}
