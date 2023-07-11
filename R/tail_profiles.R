#' Create histograms of distribution of length of tails of different classes.
#' @description tail_profiles() creates histograms of the frequency of occurrence of individual classes of tails, taking into account their length. New plots are generated within newly created directory ../profiles.
#' @param listofnames The output of gw3RACE::pair_replicates(). By default = strainslist.
#' @param data The output of gw3RACE::data_fix(). By default = dataList.
#' @param plotcolours A vector of 3 colours. By default = c("#6C91E8", "#64F069", "#C04EDC")
#' @param same_scale Whether all generated plots should have the same limits of Y axis. By default = T.
#' @param plotwidth Numeric value of width of plots that is being passed to pdf function.
#' @param plotheight Numeric value of height of plots that is being passed to pdf function.
#' @param returnplots If TRUE function will create histograms (separate for each strain) within new folder ../profiles/, if you don't like the plots and want to create your own, type FALSE and function will return a list of data frames with data ready to be plotted. By default = T.
#'
#' @return Either a list of data frames or .pdf files with histograms
#' @export
#'
#' @examples tail_profiles(plotwidth = 7, plotheight = 5)
tail_profiles <- function(listofnames = strainslist, data = dataList, plotcolours = c("#6C91E8", "#64F069", "#C04EDC"), same_scale = T, plotwidth, plotheight, returnplots = T){
  namestogrep <- listofnames[[2]]
  listofdataframes <- list()
  for (a in 1:length(namestogrep)) {
    temp <- grep(namestogrep[a], names(data), value = T)
    getmeanoftail <- function(x = temp, type_of_tail){
      sublist <- list()
      for (b in 1:length(x)) {
        temp1 <- as.data.frame(data[[x[b]]])
        temp1 <- temp1 %>% dplyr::select(c("tail_len","tail_type"))
        tempx <- temp1 %>% dplyr::filter(tail_type == c("polyA", "polyAU", "oligoU"))
        numberofreads <- length(tempx[,1])
        len_tail <- as.data.frame(table(temp1 %>% dplyr::filter(tail_type == type_of_tail) %>% dplyr::select(tail_len)))
        colnames(len_tail) <- c("tail_len", type_of_tail)
        len_tail[,2] <- len_tail[,2]/numberofreads*100
        sublist <- append(sublist, list(len_tail))
      }
      temp2 <- sublist %>% purrr::reduce(dplyr::full_join, by="tail_len")
      temp2[is.na(temp2)] <- 0
      temp2[,type_of_tail] <- rowMeans(temp2[,-1])
      temp2 <- temp2[,c("tail_len", type_of_tail)]
      return(temp2)
    }
    len_A <- getmeanoftail(type_of_tail = "polyA")
    len_AU <- getmeanoftail(type_of_tail = "polyAU")
    len_U <- getmeanoftail(type_of_tail = "oligoU")
    len_all <- merge(len_A, len_AU, by = "tail_len", all = T)
    len_all <- merge(len_all, len_U, by = "tail_len", all = T)
    len_all[is.na(len_all)] <- 0
    len_all <- len_all[1:89,]
    listofdataframes <- append(listofdataframes, list(len_all))
  }
  names(listofdataframes) <- listofnames[[1]]
  get_max_of_y <- function(x){
    localmax <- x[,2] + x[,3] + x[,2]
    localmax <- max(localmax)
    return(localmax)
  }
  max_of_y <- lapply(listofdataframes, get_max_of_y)
  max_of_y <- unlist(max_of_y)
  max_of_y <- max(max_of_y)+1
  melttable <- function(x){
    len_all <- reshape2::melt(x, id = "tail_len")
    colnames(len_all) <- c("tail_len", "tail_type", "frequency")
    len_all$tail_len <- as.numeric(len_all$tail_len)
    return(len_all)
  }
  listofdataframes <- lapply(listofdataframes, melttable)
  group.colors <- c("polyA" = plotcolours[1], "polyAU" = plotcolours[2], "oligoU" = plotcolours[3])
  tailsprofileplot <- function(data, plottitle, filename) {
    p1 <- ggplot2::ggplot(data, ggplot2::aes(fill = tail_type, x= tail_len, y=frequency)) +
      ggplot2::geom_histogram(position = ggplot2::position_stack(), stat = "identity", binwidth = 1) +
      ggplot2::scale_fill_manual(values= group.colors) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold", size=12), axis.text.y = ggplot2::element_text(face="bold", size=12)) +
      ggplot2::theme(axis.title = ggplot2::element_text(size=14,face="bold")) +
      ggplot2::xlab("Length of RNA tail") +
      ggplot2::ylab("Frequency [%]") +
      ggplot2::guides(fill=ggplot2::guide_legend(title="Tail type")) +
      ggplot2::ggtitle(plottitle) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20))

    if (same_scale == T) {
      p1 <- p1 + ggplot2::ylim(0, max_of_y)
    }

    pdf(filename, plotwidth, plotheight)
    print(p1)
    dev.off()
  }
  if (returnplots == T) {
    folder <- "profiles"
    if (file.exists(folder)) {
    } else {
      dir.create(folder)
    }
    for (a in 1:length(listofdataframes)) {
      temp <- listofdataframes[[a]]
      myname <- paste("profiles/", listofnames[[1]][a], ".pdf", sep = "")
      tailsprofileplot(data = temp, plottitle = listofnames[[1]][a], filename = myname)
    }
  } else {
    return(listofdataframes)
  }
}
