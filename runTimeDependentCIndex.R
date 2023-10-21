########################################################
##                                                    ##
## EXAMPLE OF RUNNING THE TimeDependentCIndex PROGRAM ##
##                                                    ##
########################################################

### Install and load the required packages
# install.packages('dplyr')
# install.packages('tidyr')
# install.packages('survival')
# install.packages('bench')
# install.packages('nloptr')
# install.packages('ggplot2')
rm(list = ls())
library(dplyr)
library(tidyr)
library(survival)
library(bench)
library(nloptr)
library(ggplot2)
### Load the program
source('TimeDependentCIndex.R')

### Generate the test data sets
source('DataReformat.R')

### Read the data
train.data<-read.csv('TrainData.csv')
validate.data<-read.csv('ValidateData.csv')

### Define the parameters
n.var<-3  ## number of predictors
B<-200  ## number of perturbations

### Calculate alpha values and perturbation alpha values
results<-TimeDependentCIndex(train.data,validate.data,n.var,B)
alpha<-results[[1]] ## alpha estimation
alpha.pb<-results[[2]] ## alpha perturbation results




