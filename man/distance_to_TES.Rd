% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance_to_TES.R
\name{distance_to_TES}
\alias{distance_to_TES}
\title{Creates histograms of the distribution of reads relative to the TES site.}
\usage{
distance_to_TES(
  listofnames = strainslist,
  data = dataList,
  plotcolours = c("#6C91E8", "#64F069", "#C04EDC", "#D5B54B", "#F23C28"),
  same_scale = T,
  plotwidth,
  plotheight,
  returnplots = T
)
}
\arguments{
\item{listofnames}{The output of gw3RACE::pair_replicates(). By default = strainslist.}

\item{data}{The output of gw3RACE::data_fix(). By default = dataList.}

\item{plotcolours}{A vector of 5 colours. By default = c("#6C91E8", "#64F069", "#C04EDC", "#D5B54B", "#F23C28").}

\item{same_scale}{Whether all generated plots (within 1 type of tail) should have the same limits of Y axis. By default = T.}

\item{plotwidth}{Numeric value of width of plots that is being passed to pdf function.}

\item{plotheight}{Numeric value of height of plots that is being passed to pdf function.}

\item{returnplots}{If TRUE function will create histograms (separate for each strain) within new folder ../distance/, if you don't like the plots and want to create your own, type FALSE and function will return a list of data frames with data ready to be plotted.}
}
\value{
Either a list of data frames or .pdf files with histograms
}
\description{
distance_to_TES() creates histograms of distribution of mapping of reads assigned to given class of tails relative to the TES site. New plots are generated within newly created directory ../distance.
}
\examples{
distance_to_TES(plotwidth = 5, plotheight = 4)
}
