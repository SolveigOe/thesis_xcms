args     <- commandArgs(trailingOnly=TRUE)

xcms_file <- args[1]   #xcms file
outfile  <- args[2]   #output

library(readr)
library(dplyr)

data <- read_csv(paste("/faststorage/project/forensics/02solveig/results/QC_base/", xcms_file, ".csv", sep=""))

data$rtmed <- as.numeric(data$rtmed)
data <- data %>% filter(rtmed>30) %>% filter(rtmed<1000)
temp1 <- data[rowSums(!is.na(data)) == ncol(data),]

temp1 <- temp1 %>% select(-feat, -mzmed, -rtmed, -npeaks)

pairs <- merge(colnames(temp1), colnames(temp1)) %>% 
  filter(x != y) %>% filter(case_when(grepl("b2018", y) ~ "2018",
                                      grepl("b2019", y) ~ "2019",
                                      grepl("b2017", y) ~ "2017") == 
                              case_when(grepl("b2018", x) ~ "2018",
                                        grepl("b2019", x) ~ "2019",
                                        grepl("b2017", x) ~ "2017"))

pairs <- pairs[!duplicated(t(apply(pairs, 1, sort))),] #sort row-wise and choose non-duplicated
pairs$corr <- NA
pairs$ang <- NA

pairs

for (i in 1:nrow(pairs)) {
  print(i)
  v1 <- pull(temp1[as.character(pairs$x[i])][,1])
  v2 <- pull(temp1[as.character(pairs$y[i])][,1])
  pairs$corr[i] <- cor(v1, v2, method="pearson")
  pairs$ang[i] <- acos(sum(v1*v2)/(sqrt(sum(v1*v1))*sqrt(sum(v2*v2))))
  
}

pairs <- pairs %>% mutate(batch = case_when(grepl("b2018", y) ~ "2018",
                                            grepl("b2019", y) ~ "2019",
                                            grepl("b2017", y) ~ "2017"))



batch_means <- pairs %>% group_by(batch) %>% summarise(mean_cor = mean(corr),
                                                       mean_ang = mean(ang))

means <- pairs %>% summarise(mean_cor = mean(corr), 
                             mean_ang=mean(ang)) %>% 
  bind_cols(predictors=nrow(temp1))

write_csv(batch_means, path=paste(outfile,"_batch.csv", sep=""))
write_csv(means, path=paste(outfile,".csv", sep=""))
