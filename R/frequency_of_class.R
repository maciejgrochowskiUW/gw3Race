#' Create bar plot of frequency of given tail type among tested strains
#' @description frequency_of_class() creates a bar plot of mean frequency of a class of RNA tails for every strain
#' @param listofnames The output of gw3RACE::pair_replicates(). By default = strainslist.
#' @param data The output of gw3RACE::tail_summary(). By default = dataList2.
#' @param tailclass Class of a RNA tails that is to be plotted (choose one these: "mixed_tail" "no_tail" "oligoU" "polyA" "polyAU")
#' @param nameofcontrol Name of strain that served as a control in the experiment, must be exactly the same as in the output of gw3RACE::pair_replicates().
#' @param returnplot If TRUE function will return a bar plot, if you don't like the plot and want to create your own, type FALSE and function will return a data frame with data ready to be plotted
#'
#' @return
#' @export
#'
#' @examples frequency_of_class(nameofcontrol = "Wild_type", tailclass = "polyAU")
frequency_of_class <- function(listofnames = strainslist, data = dataList2, tailclass, nameofcontrol, returnplot = T){
  namestogrep <- listofnames[[2]]
  listofdataframes <- list()
  forselect <- c(tailclass, "all_reads")
  for (a in 1:length(namestogrep)) {
    temp <- grep(namestogrep[a], names(data), value = T)
    sublist <- list()
    for (b in 1:length(temp)) {
      temp1 <- as.data.frame(data[[temp[b]]])
      temp1 <- temp1 %>% dplyr::select(all_of(forselect))
      temp1 <- sum(temp1[1], na.rm = T)/sum(temp1[2], na.rm = T)*100
      sublist <- append(sublist, list(temp1))
    }
    listofdataframes <- append(listofdataframes, list(sublist))
  }
  listofdataframes <- lapply(listofdataframes, unlist)
  mean_and_sd <- function(x){
    result <- c(mean = mean(x), sd = sd(x))
    return(result)
  }
  listofdataframes <- lapply(listofdataframes, mean_and_sd)
  toplot <- dplyr::bind_rows(listofdataframes)
  toplot$strain <- listofnames[[1]]
  p <- ggplot2::ggplot(toplot, ggplot2::aes(x=reorder(strain, -mean), y=mean, fill=factor(ifelse(strain== nameofcontrol, "Highlighted", "Normal")))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(name = "strain", values=c("red","grey50")) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
    ggplot2::ylab(paste(tailclass, "frequency", "[%]", sep = " ")) +
    ggplot2::xlab("Strain") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold.italic", size = 11, angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold", size = 12)) +
    ggplot2::theme(axis.title = ggplot2::element_text(size=14,face="bold")) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
  if (returnplot == T) {
    return(p)
  } else {
    return(toplot)
  }
}
