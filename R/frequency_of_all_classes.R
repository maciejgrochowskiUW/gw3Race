#' Create bar plots of frequency of tail types for tested strains
#' @description frequency_of_all_classes() creates bar plots of mean frequencies of all classes of RNA tails for every strain within newly created directory ../types.
#' @param listofnames The output of gw3RACE::pair_replicates(). By default = strainslist.
#' @param data The output of gw3RACE::tail_summary(). By default = dataList2.
#' @param keepMIX Whether mixed_tailed class should be plotted or not. By default = F.
#' @param plotcolours A vector of 5 colours. By default = c("#6C91E8", "#64F069", "#C04EDC", "#D5B54B", "#F23C28").
#' @param same_scale Whether all generated plots should have the same limits of Y axis. By default = T.
#' @param plotwidth Numeric value of width of plots that is being passed to pdf function.
#' @param plotheight Numeric value of height of plots that is being passed to pdf function.
#' @param returnplots If TRUE function will create bar plots (separate for each strain) within new folder ../types/, if you don't like the plots and want to create your own, type FALSE and function will return a list of data frames with data ready to be plotted. By default = T.
#'
#' @return Either a list of data frames or .pdf files with barplots
#' @export
#'
#' @examples frequency_of_all_classes(plotwidth = 5, plotheight = 3)
frequency_of_all_classes <- function(listofnames = strainslist, data = dataList2, keepMIX = F, plotcolours = c("#6C91E8", "#64F069", "#C04EDC", "#D5B54B", "#F23C28"), same_scale = T, plotwidth, plotheight, returnplots = T){
  namestogrep <- listofnames[[2]]
  listofdataframes <- list()
  forselect <- c("oligoU", "polyA", "polyAU", "no_tail", "mixed_tail", "all_reads")
  forrownames <- c("oligoU", "polyA", "polyAU", "no_tail", "mixed_tail")
  for (a in 1:length(namestogrep)) {
    temp <- grep(namestogrep[a], names(data), value = T)
    sublist <- list()
    for (b in 1:length(temp)) {
      temp1 <- as.data.frame(data[[temp[b]]])
      temp1 <- temp1 %>% dplyr::select(all_of(forselect))
      all_frequencies <- c()
      for (cols in 1:5) {
        class_frequency <- sum(temp1[cols], na.rm = T)/sum(temp1[6], na.rm = T)*100
        all_frequencies <- c(all_frequencies, class_frequency)
      }
      all_frequencies <- matrix(data = all_frequencies, ncol = 1)
      all_frequencies <- as.data.frame(all_frequencies)
      rownames(all_frequencies) <- forrownames
      sublist <- append(sublist, list(all_frequencies))
    }
    listofdataframes <- append(listofdataframes, list(sublist))
  }
  listofdataframes <- lapply(listofdataframes, dplyr::bind_cols)
  list_of_matrixes <- list()
  mean_and_sd <- function(x){
    x <- as.matrix(x)
    mymean <- rowMeans(x)
    mysd <- matrixStats::rowSds(x)
    x <- rbind(mymean,mysd)
    x <- data.frame(t(x))
    x <- tibble::rownames_to_column(x, "tail_class")
    return(x)
  }
  for (a in 1:length(listofdataframes)) {
    temp <- listofdataframes[[a]]
    temp <- mean_and_sd(temp)
    temp$tail_class <- factor(temp$tail_class, levels = c("mixed_tail", "no_tail", "oligoU", "polyAU", "polyA"))
    group.colors <- c("polyA" = plotcolours[1], "polyAU" = plotcolours[2], "oligoU" = plotcolours[3], "no_tail" = plotcolours[4], "mixed_tail" = plotcolours[5])
    if (keepMIX == F) {
      temp <- temp %>% dplyr::filter(tail_class != "mixed_tail")
    }
    list_of_matrixes <- append(list_of_matrixes, list(temp))
  }
  names(list_of_matrixes) <- listofnames[[1]]
  get_max_of_y <- function(x){
    localmax <- sum(max(x$mymean), max(x$mysd))
  }
  max_of_y <- lapply(list_of_matrixes, get_max_of_y)
  max_of_y <- unlist(max_of_y)
  max_of_y <- max(max_of_y)
  tailclassesplot <- function(data, plottitle, filename) {
    p1 <- ggplot2::ggplot(data, ggplot2::aes(fill = tail_class, y = mymean, x = tail_class)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values= group.colors) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mymean - mysd, ymax = mymean + mysd), width = .2) +
      ggplot2::coord_flip() +
      ggplot2::xlab("Tail type") +
      ggplot2::ylab("Fraction of mapped reads [%]") +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::ggtitle(plottitle) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
            axis.text.x = ggplot2::element_text(size = 12, face = "bold"), legend.position = "none",
            axis.text.y = ggplot2::element_text(face = "bold", size = 11),
            axis.title.x = ggplot2::element_text(face = "bold", size = 14), axis.title.y = ggplot2::element_text(face = "bold", size = 14))

    if (same_scale == T) {
      p1 <- p1 + ggplot2::ylim(0, max_of_y)
    }

    pdf(filename, plotwidth, plotheight)
    print(p1)
    dev.off()
  }
  if (returnplots == T) {
    folder <- "types"
    if (file.exists(folder)) {
    } else {
      dir.create(folder)
    }
    for (a in 1:length(list_of_matrixes)) {
      temp <- list_of_matrixes[[a]]
      myname <- paste("types/", listofnames[[1]][a], ".pdf", sep = "")
      tailclassesplot(data = temp, plottitle = listofnames[[1]][a], filename = myname)
    }
  } else {
    return(list_of_matrixes)
  }
}
