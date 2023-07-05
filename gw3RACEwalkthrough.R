install_github("maciejgrochowskiUW/gw3RACE")
library(gw3RACE)
#every function is documented, you can check additional flags in help files
#for some of the visualizations you can pick the colours of your liking. I suggest getting colour codes from https://www.rapidtables.com/web/color/RGB_Color.html



#set the directory where you hold your .csv files as your working directory
setwd("C:/Users/TEAM_2018/Desktop/doktorat/Tail-seq")

#read files and load them onto list. Exemplary .csv files are on my github.
dataList <- readallcsvinfolder()

#get rid of artifacts (this step is unnecessary for outputs of newer versions of Lidia Lipinska-Zubrycka software)
#you can use this function to eliminate bias-creating genes from analysis (sadly you probably have to go through entire script and do some data mining to figure out what are those genes)
dataList <- lapply(dataList, data_fix, excludedgenes = c("SPAC6F6.03c","SPCC1235.08c","SPCC663.10","SPAC4D7.10c", "SPCC576.08c", "SPBC18E5.06", "SPBC119.02", "SPAC22H12.04c", "SPAC4F10.14c", "SPCC1259.01c", "SPBC19G7.03c", "SPAC3G6.13c", "SPAC589.10c"))

#check number of reads of files with gw3RACE::nofReads() 
nofReads(dataList)

#In the next step a data frame with lengths of genes is necessary. Exemplary gen_len.csv for S. pombe is on my github.
gen_len <- read.csv("gen_len.csv", row.names = 1)

#Transform data so now ir presents mean values for individual genes rather then individual reads
dataList2 <- lapply(dataList, tail_summary, gen_len = gen_len)

#this time nofreads() will return barplot with number of genes that passed the given threshold
nofReads(dataList2)

#To obtain list of genes that are present in all samples after setting threshold of minimal amount of reads in gw3RACE::tail_summary() run following lines:
dataList2 <- lapply (dataList2, keepcommongenes)
dataList <- lapply (dataList, keepcommongenes)

#plot number of genes again, should be even for all samples
nofReads(dataList2)

#get list of strains/cell lines used in experiment. This is a potential pitfall as you have to know how to use REGEX. 
strainslist <- pair_replicates(dataList)

#Create violin plot of RNA tails length
plot_of_tail_length(nameofcontrol = "Wild_type")

#Create bar plot of frequency of given tail type among tested strains
frequency_of_class(nameofcontrol = "Wild_type", tailclass = "polyAU")

#Create bar plots of frequency of tail types separately for each tested strain. New folder "types" is created, this is where you can find plots.
frequency_of_all_classes(plotwidth = 5, plotheight = 3)

#Create histograms of distribution of length of tails of different classes separately for each tested strain. New folder "profiles" is created, this is where you can find plots.
tail_profiles(plotwidth = 7, plotheight = 5)

#Create histograms of the distribution of reads relative to the TES site separately for each tested strain. New folder "distance" is created, this is where you can find plots.
distance_to_TES(plotwidth = 5, plotheight = 4)