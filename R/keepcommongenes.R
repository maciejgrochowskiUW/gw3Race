#' Keep the genes that surpassed the cutoff for every single sample in your experiment
#' @description keepcommongenes() filters data for genes that are present in all samples
#' @param data A data frame stored on the list that is an output of gw3RACE::data_fix() or gw3RACE::tail_summary(), generally it is recommended to use this function with lapply()
#' @param querry Output of gw3RACE::tail_summary(). By default = dataList2.
#'
#' @return List of data frames
#' @export
#'
#' @examples dataList2 <- lapply (dataList2, keepcommongenes)
#' @examples dataList <- lapply (dataList, keepcommongenes)
keepcommongenes <- function(data, querry = dataList2){
  commongenes <- lapply(querry, "[",  "gene")
  commongenes <- Reduce(dplyr::intersect, commongenes)
  data <- data %>% dplyr::filter(gene %in% commongenes$gene)
  return(data)
}
