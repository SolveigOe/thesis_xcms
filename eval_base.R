args     <- commandArgs(trailingOnly=TRUE)

prep_file <- args[1]  #Which of the data files to use
run      <- args[2]   #run = 1:10
outfile   <- args[3]  #outfilebase = e.g. "outfiles/age_regression.outfile.compounds.1.ranger.1"


#########################
##EVALUATION OF RESULTS##
#########################

library(caret)
library(readr)
library(tidyr)
library(dplyr)

ms <- read_csv(paste("/faststorage/project/forensics/02solveig/results/base/", prep_file, ".csv", sep=""))

run <- as.integer(run)
set.seed(run)

# Data for model:

input_x <- as.matrix(select(ms, -response))
input_y <- ms$response

# Using CV:
tune_control <- caret::trainControl(
  method = "repeatedcv", number = 10, repeats = 10,
  verboseIter = TRUE,
  allowParallel = F,
  classProbs = ifelse(is.factor(input_y), TRUE, FALSE),
  savePredictions = "final")

#Training model:
model_cv <- caret::train(
  x = input_x,
  y = input_y,
  method = "ranger",
  trControl = tune_control,
  tuneLength  = 10,             #number of default values to try - choose e.g. 10
  metric = ifelse(is.factor(input_y), "Kappa", "RMSE"))


rmse_cv <- min(model_cv$results$RMSE)

value = sub("prep_xcms_", "", prep_file)

result <- tibble(run=run, value=value, rmse_cv)

write_csv(result, path=paste(outfile,".csv", sep=""))
