#' Creates histograms of the distribution of reads relative to the TES site.
#' @description distance_to_TES() creates histograms of distribution of mapping of reads assigned to given class of tails relative to the TES site. New plots are generated within newly created directory ../distance.
#' @param listofnames The output of gw3RACE::pair_replicates(). By default = strainslist.
#' @param data The output of gw3RACE::data_fix(). By default = dataList.
#' @param plotcolours A vector of 5 colours. By default = c("#6C91E8", "#64F069", "#C04EDC", "#D5B54B", "#F23C28").
#' @param same_scale Whether all generated plots (within 1 type of tail) should have the same limits of Y axis. By default = T.
#' @param plotwidth Numeric value of width of plots that is being passed to pdf function.
#' @param plotheight Numeric value of height of plots that is being passed to pdf function.
#' @param returnplots If TRUE function will create histograms (separate for each strain) within new folder ../distance/, if you don't like the plots and want to create your own, type FALSE and function will return a list of data frames with data ready to be plotted.
#'
#' @return Either a list of data frames or .pdf files with histograms
#' @export
#'
#' @examples distance_to_TES(plotwidth = 5, plotheight = 4)
distance_to_TES <- function(listofnames = strainslist, data = dataList, plotcolours = c("#6C91E8", "#64F069", "#C04EDC", "#D5B54B", "#F23C28"), same_scale = T, plotwidth, plotheight, returnplots = T){
  listofdataframes <- list()
  namestogrep <- listofnames[[2]]
  for (a in 1:length(namestogrep)) {
    temp <- grep(namestogrep[a], names(data), value = T)
    getmeanofdistance <- function(x = temp, type_of_tail){
      sublist <- list()
      for (b in 1:length(x)) {
        temp1 <- as.data.frame(data[[x[b]]])
        temp1 <- temp1 %>% dplyr::select(c("tail_type", "rel_distance_to_TES"))
        temp1[,2] <- round(temp1[,2], digits = 3)
        temp1 <- temp1 %>% dplyr::filter(rel_distance_to_TES <= 0.05 & rel_distance_to_TES >= -1)
        allreadsnumber <- length(temp1[,1])
        temp2 <- temp1 %>% dplyr::filter(tail_type == type_of_tail)
        temp2 <- as.data.frame(table(temp2$rel_distance_to_TES))
        temp2[,2] <- temp2[,2]/allreadsnumber*100
        colnames(temp2) <- c("distance", type_of_tail)
        sublist <- append(sublist, list(temp2))
      }
      temp3 <- sublist %>% purrr::reduce(dplyr::full_join, by="distance")
      temp3[is.na(temp3)] <- 0
      temp3[,type_of_tail] <- rowMeans(temp3[,-1])
      temp3 <- temp3[,c("distance", type_of_tail)]
      return(temp3)
    }
    distance_A <- getmeanofdistance(type_of_tail = "polyA")
    distance_AU <- getmeanofdistance(type_of_tail = "polyAU")
    distance_U <- getmeanofdistance(type_of_tail = "oligoU")
    distance_notail <- getmeanofdistance(type_of_tail = "no_tail")
    distance_mix <- getmeanofdistance(type_of_tail = "mixed_tail")
    distance_all <- merge(distance_A, distance_AU, by = "distance", all = T)
    distance_all <- merge(distance_all, distance_U, by = "distance", all = T)
    distance_all <- merge(distance_all, distance_notail, by = "distance", all = T)
    distance_all <- merge(distance_all, distance_mix, by = "distance", all = T)
    distance_all[is.na(distance_all)] <- 0
    listofdataframes <- append(listofdataframes, list(distance_all))
  }
  melttable <- function(x){
    distance_all <- reshape2::melt(x, id = "distance")
    colnames(distance_all) <- c("distance", "tail_type", "frequency")
    distance_all$distance <- as.numeric(levels(distance_all$distance))
    return(distance_all)
  }
  listofdataframes <- lapply(listofdataframes, melttable)
  get_max_of_y <- function(x,y){
    localmaxes <- c()
    for (a in y) {
      temp <- x
      temp <- temp %>% dplyr::filter(tail_type == a)
      localmax <- max(temp$frequency)
      localmaxes <- c(localmaxes, localmax)
    }
    return(localmaxes)
  }
  maxofy <- lapply(listofdataframes, get_max_of_y, y = c("polyA", "polyAU", "oligoU", "no_tail", "mixed_tail"))
  get_max_of_y <- function(x, tailtype){
    maxy <- c()
    for (a in 1:length(maxofy)) {
      temp <- maxofy[[a]]
      maxy <- c(maxy, temp[tailtype])
    }
    maxy <- max(maxy)
    return(maxy)
  }
  maxyA <- get_max_of_y(maxofy, 1)
  maxyAU <- get_max_of_y(maxofy, 2)
  maxyU <- get_max_of_y(maxofy, 3)
  maxyN <- get_max_of_y(maxofy, 4)
  maxyM <- get_max_of_y(maxofy, 5)
  maxy <- c("polyA" = maxyA, "polyAU" = maxyAU, "oligoU" = maxyU, "no_tail" = maxyN, "mixed_tail" = maxyM)
  group.colors <- c("polyA" = plotcolours[1], "polyAU" = plotcolours[2], "oligoU" = plotcolours[3], "no_tail" = plotcolours[4], "mixed_tail" = plotcolours[5])
  distanceplot <- function (data, plottitle,filename, color, max_of_y){
    p1 <- ggplot2::ggplot(data, ggplot2::aes(x=distance, y=frequency)) +
      ggplot2::geom_histogram(stat = "identity", fill= color, binwidth = 0.001) +
      ggplot2::xlab("Distance from TES") +
      ggplot2::ylab("Frequency [%]") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold",size=12),axis.text.y = ggplot2::element_text(face="bold", size=12)) +
      ggplot2::theme(axis.title = ggplot2::element_text(size=14,face="bold")) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))+
      ggplot2::ggtitle(plottitle)

    if (same_scale == T) {
      p1 <- p1 + ggplot2::ylim(0, max_of_y)
    }

    pdf(filename, plotwidth, plotheight)
    print(p1)
    dev.off()
  }
  if (returnplots == T) {
    folder <- "distance"
    if (file.exists(folder)) {
    } else {
      dir.create(folder)
    }
    setwd("distance")
    if (file.exists("polyA")) {
    } else {
      dir.create("polyA")
    }
    if (file.exists("polyAU")) {
    } else {
      dir.create("polyAU")
    }
    if (file.exists("oligoU")) {
    } else {
      dir.create("oligoU")
    }
    if (file.exists("no_tail")) {
    } else {
      dir.create("no_tail")
    }
    if (file.exists("mixed_tail")) {
    } else {
      dir.create("mixed_tail")
    }
    setwd("../")
    for (a in 1:length(listofdataframes)) {
      temp <- listofdataframes[[a]]
      for (b in c("polyA", "polyAU", "oligoU", "no_tail", "mixed_tail")) {
        temp1 <- temp %>% dplyr::filter(tail_type == b)
        myname <- paste("distance/", b, "/", listofnames[[1]][a], ".pdf", sep = "")
        distanceplot(data = temp1, plottitle = listofnames[[1]][a], filename = myname, color = group.colors[b], max_of_y = maxy[b])
      }
    }
  } else {
    return(listofdataframes)
  }
}
