library(xcms)
library(mzR)
library(magrittr)
library(msdata)
library(readr)
register(bpstart(MulticoreParam(15))) #To make sure it can run

args     <- commandArgs(trailingOnly=TRUE)

parameter <- args[1]  #value of the parameter that will change
outfile  <- args[2]   #outfilebase = e.g. "outfiles/age_regression.outfile.compounds.1.ranger.1"

#parameter     <- as.integer(parameter)

debug=T

###################
#### Read file ####
###################
cat("\n\nStart:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(as.character(Sys.time()), file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)

file <- dir("/faststorage/project/forensics/02solveig/data/all_data/", all.files = T,
            full.name = TRUE,
            pattern = "ML")

s_names <- sub(basename(file), pattern = ".ML",
               replacement = "", fixed = TRUE)

s_group <- c()

for (i in s_names){
  if (grepl("QC", i, fixed=T)){
    s_group <- c(s_group, "QC")
  } else {
    s_group <- c(s_group, "SAMP")
  }
}

pd <- data.frame(sample_name = s_names,
                 sample_group = s_group,
                 stringsAsFactors = FALSE)

xcms <- readMSData(file, pdata = new("NAnnotatedDataFrame", pd), mode="onDisk")

cat("\n\nreadMSData:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(as.character(Sys.time()), file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
if(debug==T){
cat("\ndebug:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(xcms@phenoData@data[["sample_name"]], file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
}

########################################
#### Chromatographic Peak Detection ####
########################################

#Peak detection from the extracted chromatogram: Different stuffs can be changed here
#Peakwidth, snthresh, ppm, mzdiff are all specific for the machine used.
cwp <- CentWaveParam(peakwidth = c(4,30), snthresh = 6, ppm = 12, prefilter = c(3,200), mzdiff = 0.01) 
xdata <- findChromPeaks(xcms, param=cwp)

cat("\n\nChromPeaks:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(as.character(Sys.time()), file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
if(debug==T){
cat("\ndebug:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(xdata@phenoData@data[["sample_name"]], file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
}


#Merging of overlapping peaks - can be added if it has interest
#mpp <- MergeNeighboringPeaksParam(expandRt = 4)
#data_cent_pp <- refineChromPeaks(xdata, param = mpp)

#####################
##### Alignment #####
#####################

#Retention time correction
#Obiwarp method used to align the samples
#binSize can be changes
#Evt muligt at bruge subset-based alignment
# Evt kan der bruges PeakGroupsParam i stedet for obiwarp

xdata_adj <- adjustRtime(xdata, param=ObiwarpParam(binSize = 1))

cat("\n\adjustRtime:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(as.character(Sys.time()), file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
if(debug==T){
cat("\ndebug:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(xdata_adj@phenoData@data[["sample_name"]], file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
}

##########################
##### Correspondence #####
##########################

#Possible to change parameters here
#minfraction, bw, evt kan andet tilfojes
pdp <- PeakDensityParam(sampleGroups = xdata_adj$sample_group, minFraction = 0.4, bw=2.5)
xdata_cor <- groupChromPeaks(xdata_adj, param=pdp)

cat("\n\ngroupChromPeaks:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(as.character(Sys.time()), file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
if(debug==T){
cat("\ndebug:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(xdata_cor@phenoData@data[["sample_name"]], file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
}


#Fill out missing peaks, some parameters can possibly be changed here - eg fixedrt
xdata_cor1 <- fillChromPeaks(xdata_cor)
fv <- featureValues(xdata_cor1, value="into")
fd <- featureDefinitions(xdata_cor1)

names <- c("feat", "mzmed", "rtmed", "npeaks", xdata_cor@phenoData@data[["sample_name"]])
data <- cbind(feat = rownames(fd), mzmed = fd$mzmed, rtmed = fd$rtmed, npeaks = fd$npeaks ,fv)
data <- as.data.frame(data)
colnames(data) <- names

cat("\n\ndataframe:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(as.character(Sys.time()), file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
if(debug==T){
cat("\ndebug:\n", file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
cat(as.character(colnames(data)), file="/faststorage/project/forensics/02solveig/scripts/xcms_base/xcms_run.txt", append=T)
}

write_csv(data, path=paste(outfile,".csv", sep=""))