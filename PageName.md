# Introduction #

# Installation #

# Training #

Training new model is very simple.  First load sample, ask user for input about quality of peaks, generate features that will be used in classification, train model.  The sequence of commands is shown below.

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

# Details #

Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages