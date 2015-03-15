# Introduction #

mzRock is a machine learning toolkit for automation of mass spectroscopy data analysis. The program was designed to work with chromatograms obtained from triple-quad instruments, but it could be easily adapted to work with other type of instruments.

The main goal of the software is to mimic human calls in peak selection, grouping, and annotation of presence/absence of metabolite in a sample. From user point of view, all that needs to be done is to drop files into designated folder, and from this point on every step of the analysis is automated. The data is extracted into simple csv file format[mzCSV](mzCSV.md) and software carries out the following steps:

  * detect peaks
  * group peaks
  * calculate various features used in classification
  * use machine trained model to classify peaks and assign goodness probabilities
  * match retention times to known metabolites
  * write a report and various quality control spreadsheets

# Installation #

Before you start, you will need to install Perl and R. Perl is only necessary if you are planning to run monitorFolder.pl script, which will automatically convert mass spec files into [mzCSV format](mzCSV.md) and execute R scripts.

Unfortunately, I have not found a way to convert Thermo .RAW files without Xcalibar installation, so if you are using Thermo QQQ instruments, you will need to run on Windows with Xcalibar installed.
There are rumors that Xcalibar can be run on Linux with Wine.

### Detailed Instructions for Windows User ###

  1. Start R, and install "randomForest" and "waveslim" libararies. This can be done from command line with "install.packages command. Or if you prefer using menu options in R GUI.
  1. Download mzrock.release.zip file and unzip into C:\mzRock folder. This will create two subfolders. C:\mzRock\bin and C:\mzRock\DataSets
  1. Edit configuration options in "mzrock.pl" file to specify location of R executable.
  1. Start monitorFolder command in the bin deirectory. This monitorFolder service waits for new files to appear in the DataSets folder, converts files to mzCSV and runs R scripts.

# Training Models #

The program comes with a default model for making calls on peak goodness. This model is by no means universal, and you will probably need to generate your own model for your instrument.

First, start R, and change the work directory to a folder with mzCSV files.

Training new model is very simple.  First load sample, ask user for input about quality of peaks, generate features that will be used in classification. The sequence of commands is shown below. The two most critical commands are getQuality() and validateSamples().  getQuality() computes features that will be used in training and validateSamples() will ask you to judge if peaks are good or bad. These calls are recorded and used in training the model.

When you run validateSamples() command, you will be presented with a plot of a peak group and asked to type in one of the following responses.

  * 'g' Peak is good  press
  * 'b' Peak is bad
  * 'G' All peaks in group are good
  * 'B' All peaks in group bad bad
  * 'i' Ignore this peak
  * 'I' Ignore all peaks
  * 'u' Undo
  * 'Q' Stop traing

The sequence of commands below also includes generation of confusion matrix which can be used to asses the quality of the model.

```
#load R commands
source('../../bin/mzRock.beta.R')

#Load Samples
loadSamples();

#features that will be used in classification
Q = getQuality();

# COMMANDS TO TRAIN A MODEL
#ask user for peak labels. Press 'g' for good peaks, 'b' for bad peaks, 'i' to i
gnore, 'Q' to stop
validateSamples()

#pull out subset of peaks that were validated
V = which(!is.na(mzrock$peaks$valid))
V= data.frame( peakid=V, valid=mzrock$peaks$valid[V]);
V = V[ V$valid == "b" | V$valid == "g", ]  #keep only peaks marked as 'g' or 'b'

V$valid = factor(V$valid)

#create a training and testing subsets 
A=sample(1:nrow(V),nrow(V)*3/4);

#split 3/4 for training,  1/4 for testing
Train=V[A,]; Test=V[-A,];
#classifiation 
M=classifySamples(Q,Train,Test); 

#confusion table
table(predict(M$model,Q[Test$peakid,]),Test$valid)

#save model to a file
save(M$model, "mymodel.Rdata");

```

# Running R code manually #

Normally, the code below will be executed automatically. But sometimes it is useful
to run it within R manually for debugging purposes. Once again, your R work directory
should be set to a folder with mzCSV files.

```

#load R commands
source('../../bin/mzRock.beta.R')

#Load Samples
loadSamples();

#load existing model ( change here if you have your own model)
load('../../bin/DEFAULT.MODEL7.Rdata');

#write report
writeReport(model);


```

# Matching Metabolite Retention Times #

Matching of metabolite retention times is done in the last step of the report generation. The file containing retention times is called SRM2.csv. It is located in the top level of DataSets folder. In our lab, this folder and the SRM2 file are placed on a network shared drive writable by everyone in a group. Each user has a subfolder with DataSets subdirectory for their files, and we append new metabolites to the end of the SRM file whenever new metabolites are being monitored.

The names of the fields in SRM2 are adapted from Thermo speak.

  * polarity
  * compound = name of the compoumd
  * precursorIntensity = Parent m/z
  * collisionEnergy
  * productMz = Fragment m/z
  * expectedRt = Retention time in minutes

```
  polarity,compound,category,precursorIntensity,collisionEnergy,productMz,expectedRt
-,2-oxobutanoate,Carboxylic Acids,101,11,57.2,21.8
-,acetoacetate,Carboxylic Acids,101.05,12,57.2,14.1
-,glycerate,Carbohydrate Precursors/ Derivatives,105,15,75,1

```

# Reports #

If all goes well, mzRock will generate a report folder with number of spreadsheets and pdf files.

  * qualityControl.csv file which contains most complete information with one peak per line output.
  * highestPeaks.csv file which contains table of groups vs samples. Here the rows are groups, samples are columns, and values in cells is AreaTop measure.  AreaTop is an average of top 3 intensities in the peak for a given sample.

  * predictionReport.pdf file contains superimposed EICs from all samples
  * stackedPlots.pdf file which contains separated EICs for each group.

## note about naming of groups ##

Groups are numbered based on quality, one being the highest. When matching the retention time, the group that has the highest quality within 2 minutes of expected retention time is renamed to group zero.













