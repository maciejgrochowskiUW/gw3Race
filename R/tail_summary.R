#' Get summary of per gene analysis
#' @description tail_summary() summarises data that is stored in the initial data frames so now it holds info on (per gene) mean tail length and frequency of each tail kind
#' @param data A data frame stored on the list that is an output of gw3RACE::data_fix(), generally it is recommended to use this function with lapply()
#' @param threshold A numeric value, minimal number of reads per gene in a sample that qualifies that gene for future analysis, by default it's equal to 20
#' @param gen_len A data frame that contains lengths of all transcripts of studied organism (TRANSCRIPTS NOT ORFS), it is used to normalize notail reads to gene length. By default = gen_len.
#'
#' @return List of data frames
#' @export
#'
#' @examples dataList2 <- lapply(dataList, tail_summary, gen_len = gen_len)
tail_summary <- function(data, threshold = 20, gen_len) {
  tail_length <- data %>% dplyr::filter(tail_len > 0) %>%
    dplyr::select(gene, tail_len) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(mean = mean(tail_len), reads = dplyr::n())
  type <- data %>% dplyr::select(gene, RNA_type) %>% dplyr::distinct()
  dat1 <- data %>% dplyr::select(gene, tail_type)
  dat1 <- dat1 %>% dplyr::group_by(gene) %>% dplyr::count(tail_type)
  dat1 <- as.data.frame(reshape2::acast(dat1, gene~tail_type))
  dat1$counts <- rowSums(dat1, na.rm = TRUE)
  dat1 <- merge(dat1, gen_len, by.x =0, by.y = 1)
  dat1$no_tail <- dat1$no_tail/dat1$len
  dat1$all_reads <- rowSums(dat1[,c("mixed_tail", "no_tail", "polyA", "polyAU", "oligoU")], na.rm = TRUE)
  dat1 <- merge(dat1, tail_length, by = 1)
  dat1 <- dat1[dat1$counts >= threshold,]
  dat1 <- merge(dat1,type, by = 1)
  names(dat1)[names(dat1) == 'Row.names'] <- "gene"
  return(dat1)
}
