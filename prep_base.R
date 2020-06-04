args     <- commandArgs(trailingOnly=TRUE)

xcms_file <- args[1]   #xcms file
outfile  <- args[2]   #output

#######################################
#Collecting and prepare data from XCMS#
#######################################

library(readr)
library(dplyr)
library(readxl)

###########################
#Data for info for samples#
###########################

#header, mzmed, rtmed, ID - and making them tidy
#When all data should be run:

ms_info <- read_csv("/faststorage/project/forensics/02solveig/data/ms_info.csv") %>% select(name, day, batch)

##############
#Read data in#
##############

temp <- read_csv(paste("/faststorage/project/forensics/02solveig/results/base/", xcms_file, ".csv", sep=""))
#Only peaks with at least 5 samples and a retention time of at least 30:
temp <- temp %>% filter(npeaks>=5) %>% filter(rtmed>30) %>% select(-npeaks, -rtmed, -mzmed)
temp1 <- data.frame(t(temp[-1]))
colnames(temp1) <- temp$feat
data_cor <- temp1 %>% mutate("name"=rownames(temp1))
data_cor$name <- sub(".mzXML","",data_cor$name)
data_cor$name <- sub(".mzML","",data_cor$name)

# Add the information:
temp <- merge(ms_info, data_cor, by="name")

###############
#Changing NAs#
###############

#If more than 5% of a predictor is NA -> it is removed
temp1 <- temp[colSums(!is.na(temp)) >= nrow(temp)*0.95]
#temp1$day <- as.numeric(temp1$day)

#Rest of NAs as average of that day
temp <- temp1 %>% group_by(day) %>% 
  mutate_all(funs(ifelse(is.na(.), mean(., na.rm = TRUE),.))) %>% ungroup()


###################
#Row normalization#
###################

#Row normalization: number/sum(row)
row_sum <- rowSums(temp %>% select(-name, -day, -batch))
data_row <- temp %>% mutate_at(vars(-name, -day, -batch), ~ ./row_sum)
rm(row_sum)

#####################
#Batch normalization#
#####################

######Normalization - mean centering######
data_norm <- data_row
for (i in 4:length(data_norm)) {
  col_stand <- select(data_row,batch,i) %>% 
    group_by(batch) %>% 
    mutate_at(2, ~. - (mean(.) - (colSums(data_row[i])/nrow(data_row))))
  data_norm[i] <- col_stand[2]
}

rm(col_stand,i)

##############
#Collect data#
##############

full_data <- data_norm %>% 
  filter(!grepl(pattern="QC", x= day)) %>% 
  filter(batch!="2018sub") %>% 
  mutate(response=as.numeric(day)) %>% 
  select(-day,-batch,-name)

write_csv(full_data, path=paste(outfile,".csv", sep=""))

