#' Create violin plot of RNA tails length
#' @description plot_of_tail_length() creates violin plot of mean tail lengths per gene measured in Genome-wide 3' RACE for all strains
#' @param listofnames The output of gw3RACE::pair_replicates(). By default = strainslist.
#' @param data The output of gw3RACE::tail_summary(). By default = dataList2.
#' @param nameofcontrol Name of strain that served as a control in the experiment, must be exactly the same as in the output of gw3RACE::pair_replicates()
#' @param returnplot if TRUE function will return a violin plot, if you don't like the plot and want to create your own, type FALSE and function will return a data frame with data ready to be plotted. By defult = T.
#'
#' @return Either a data frame or a violin plot based on it
#' @export
#'
#' @examples plot_of_tail_length(nameofcontrol = "Wild_type")
plot_of_tail_length <- function(listofnames = strainslist, data = dataList2, nameofcontrol, returnplot = T){
  namestogrep <- listofnames[[2]]
  listofdataframes <- list()
  for (a in 1:length(namestogrep)) {
    temp <- grep(namestogrep[a], names(data), value = T)
    sublist <- list()
    for (b in 1:length(temp)) {
      temp1 <- as.data.frame(data[[temp[b]]])
      temp1 <- temp1 %>% dplyr::select(c("gene","mean"))
      sublist <- append(sublist, list(temp1))
    }
    temp3 <- sublist %>% purrr::reduce(dplyr::full_join, by="gene")
    temp3$mean <- rowMeans(temp3[,-1])
    temp3 <- temp3[,c("gene","mean")]
    temp3$strain <- listofnames[[1]][a]
    listofdataframes <- append(listofdataframes, list(temp3))
  }
  names(listofdataframes) <- listofnames[[1]]
  toplot <- dplyr::bind_rows(listofdataframes)

  p <- ggplot2::ggplot(toplot, ggplot2::aes(x=reorder(strain, -mean), y=mean, fill=factor(ifelse(strain== nameofcontrol, "Highlighted", "Normal")))) +
    ggplot2::geom_violin(trim=FALSE)+
    ggplot2::scale_fill_manual(name = "strain", values=c("red","grey50")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Strain") +
    ggplot2::ylab("Tail lenght") +
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
