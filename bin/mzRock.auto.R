#load R commands
source('mzRock.beta.R')

#Load Samples
loadSamples();

#load existing model ( change here if you have your own model)
load('DEFAULT.MODEL7.Rdata');

#write report
writeReport(model);
