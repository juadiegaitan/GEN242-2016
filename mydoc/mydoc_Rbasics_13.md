---
title: Analysis Routine
keywords: 
last_updated: Thu Apr  7 17:52:19 2016
---

## Overview

The following exercise introduces a variety of useful data analysis utilities in R. 

## Analysis Routine: Data Import

- __Step 1__: To get started with this exercise, direct your R session to a dedicated workshop directory and download into this directory the following sample tables. Then import the files into Excel and save them as tab delimited text files.

    - [MolecularWeight_tair7.xls](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls)
    - [TargetP_analysis_tair7.xls](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls)

__Import the tables into R__

Import molecular weight table


{% highlight r %}
my_mw <- read.delim(file="MolecularWeight_tair7.xls", header=T, sep="\t") 
my_mw[1:2,]
{% endhighlight %}

Import subcelluar targeting table

{% highlight r %}
my_target <- read.delim(file="TargetP_analysis_tair7.xls", header=T, sep="\t") 
my_target[1:2,]
{% endhighlight %}

Online import of molecular weight table

{% highlight r %}
my_mw <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls", header=T, sep="\t") 
my_mw[1:2,]
{% endhighlight %}

{% highlight txt %}
##   Sequence.id Molecular.Weight.Da. Residues
## 1 AT1G08520.1                83285      760
## 2 AT1G08530.1                27015      257
{% endhighlight %}

Online import of subcelluar targeting table

{% highlight r %}
my_target <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls", header=T, sep="\t") 
my_target[1:2,]
{% endhighlight %}

{% highlight txt %}
##      GeneName Loc   cTP   mTP    SP other
## 1 AT1G08520.1   C 0.822 0.137 0.029 0.039
## 2 AT1G08530.1   C 0.817 0.058 0.010 0.100
{% endhighlight %}

## Merging Data Frames

- __Step 2__: Assign uniform gene ID column titles


{% highlight r %}
colnames(my_target)[1] <- "ID"
colnames(my_mw)[1] <- "ID" 
{% endhighlight %}

- __Step 3__: Merge the two tables based on common ID field


{% highlight r %}
my_mw_target <- merge(my_mw, my_target, by.x="ID", by.y="ID", all.x=T)
{% endhighlight %}

- __Step 4__: Shorten one table before the merge and then remove the non-matching rows (NAs) in the merged file


{% highlight r %}
my_mw_target2a <- merge(my_mw, my_target[1:40,], by.x="ID", by.y="ID", all.x=T)  # To remove non-matching rows, use the argument setting 'all=F'.
my_mw_target2 <- na.omit(my_mw_target2a) # Removes rows containing "NAs" (non-matching rows).
{% endhighlight %}

- __Homework 3D__: How can the merge function in the previous step be executed so that only the common rows among the two data frames are returned? Prove that both methods - the two step version with `na.omit` and your method - return identical results. 
- __Homework 3E__: Replace all `NAs` in the data frame `my_mw_target2a` with zeros.



## Filtering Data

- __Step 5__: Retrieve all records with a value of greater than 100,000 in 'MW' column and 'C' value in 'Loc' column (targeted to chloroplast).


{% highlight r %}
query <- my_mw_target[my_mw_target[, 2] > 100000 & my_mw_target[, 4] == "C", ] 
query[1:4, ]
{% endhighlight %}

{% highlight txt %}
##              ID Molecular.Weight.Da. Residues Loc   cTP   mTP    SP other
## 219 AT1G02730.1               132588     1181   C 0.972 0.038 0.008 0.045
## 243 AT1G02890.1               136825     1252   C 0.748 0.529 0.011 0.013
## 281 AT1G03160.1               100732      912   C 0.871 0.235 0.011 0.007
## 547 AT1G05380.1               126360     1138   C 0.740 0.099 0.016 0.358
{% endhighlight %}

{% highlight r %}
dim(query)
{% endhighlight %}

{% highlight txt %}
## [1] 170   8
{% endhighlight %}

- __Homework 3F__: How many protein entries in the `my`_mw`_target` data frame have a MW of greater then 4,000 and less then 5,000. Subset the data frame accordingly and sort it by MW to check that your result is correct.


{% highlight r %}
query2 <- my_mw_target[my_mw_target[, 2] > 4000 & my_mw_target[, 2] < 5000, ] 
dim(query2)
{% endhighlight %}

{% highlight txt %}
## [1] 38  8
{% endhighlight %}

{% highlight r %}
query2[order(query2[,2]),] 
{% endhighlight %}

{% highlight txt %}
##                ID Molecular.Weight.Da. Residues Loc   cTP   mTP    SP other
## 13211 AT2G47660.1                 4043       35   _ 0.153 0.099 0.140 0.909
## 13808 AT3G05080.1                 4068       35   * 0.091 0.200 0.231 0.704
## 17007 AT3G42090.1                 4113       37   * 0.089 0.071 0.345 0.825
## 31743 ATCG00550.1                 4117       40   S 0.003 0.065 0.953 0.366
## 24539 AT5G02502.1                 4133       35   S 0.006 0.013 0.995 0.106
## 31739 ATCG00510.1                 4134       37   S 0.005 0.025 0.992 0.066
## 31713 ATCG00080.1                 4168       36   S 0.012 0.030 0.991 0.044
## 31748 ATCG00600.1                 4204       37   S 0.010 0.054 0.981 0.058
## 27714 AT5G35480.1                 4298       38   * 0.049 0.157 0.398 0.672
## 192   AT1G02490.1                 4337       38   * 0.094 0.085 0.438 0.695
## 30410 AT5G57567.1                 4356       41   S 0.066 0.044 0.843 0.331
## 14816 AT3G13240.1                 4366       40   * 0.216 0.075 0.261 0.784
## 24276 AT4G39403.1                 4390       36   * 0.009 0.713 0.363 0.101
## 5668  AT1G59535.1                 4404       39   _ 0.122 0.132 0.135 0.895
## 31745 ATCG00570.1                 4424       39   S 0.005 0.309 0.718 0.213
## 31762 ATCG00760.1                 4461       37   * 0.110 0.542 0.074 0.334
## 31744 ATCG00560.1                 4470       38   S 0.016 0.148 0.737 0.321
## 823   AT1G07610.1                 4495       45   * 0.656 0.044 0.248 0.589
## 288   AT1G03200.1                 4573       41   * 0.089 0.358 0.119 0.647
## 822   AT1G07600.1                 4580       45   * 0.545 0.046 0.234 0.713
## 16351 AT3G25717.1                 4620       40   M 0.022 0.890 0.039 0.208
## 292   AT1G03240.1                 4644       42   * 0.114 0.317 0.132 0.597
## 25923 AT5G14560.1                 4681       40   * 0.091 0.268 0.122 0.779
## 1568  AT1G13245.1                 4719       41   * 0.024 0.746 0.067 0.306
## 31756 ATCG00700.1                 4722       43   S 0.011 0.048 0.972 0.099
## 6552  AT1G67148.1                 4734       45   _ 0.203 0.096 0.126 0.897
## 6566  AT1G67265.1                 4738       40   * 0.059 0.283 0.098 0.738
## 26710 AT5G21020.2                 4744       42   S 0.032 0.069 0.951 0.079
## 30432 AT5G57730.1                 4763       42   * 0.103 0.308 0.126 0.511
## 14906 AT3G13857.1                 4778       42   * 0.244 0.476 0.303 0.018
## 27114 AT5G24980.1                 4839       42   S 0.004 0.030 0.991 0.076
## 13050 AT2G46390.1                 4884       46   S 0.003 0.155 0.957 0.075
## 16657 AT3G28120.1                 4897       41   S 0.006 0.075 0.970 0.145
## 3797  AT1G32670.1                 4926       44   * 0.233 0.607 0.129 0.103
## 13470 AT3G02390.1                 4933       43   _ 0.135 0.118 0.147 0.877
## 9128  AT2G14460.1                 4950       41   * 0.143 0.165 0.139 0.784
## 1389  AT1G11785.1                 4983       46   _ 0.142 0.110 0.146 0.901
## 23211 AT4G31030.1                 4983       44   * 0.041 0.133 0.352 0.795
{% endhighlight %}

## String Substitutions

- __Step 6__: Use a regular expression in a substitute function to generate a separate ID column that lacks the gene model extensions.
<<label=Exercise 4.7, eval=TRUE, echo=TRUE, keep.source=TRUE>>=


{% highlight r %}
my_mw_target3 <- data.frame(loci=gsub("\\..*", "", as.character(my_mw_target[,1]), perl = TRUE), my_mw_target)
my_mw_target3[1:3,1:8]
{% endhighlight %}

{% highlight txt %}
##        loci          ID Molecular.Weight.Da. Residues Loc  cTP   mTP    SP
## 1 AT1G01010 AT1G01010.1                49426      429   _ 0.10 0.090 0.075
## 2 AT1G01020 AT1G01020.1                28092      245   * 0.01 0.636 0.158
## 3 AT1G01020 AT1G01020.2                21711      191   * 0.01 0.636 0.158
{% endhighlight %}

- __Homework 3G__: Retrieve those rows in `my_mw_target3` where the second column contains the following identifiers: `c("AT5G52930.1", "AT4G18950.1", "AT1G15385.1", "AT4G36500.1", "AT1G67530.1")`. Use the `%in%` function for this query. As an alternative approach, assign the second column to the row index of the data frame and then perform the same query again using the row index. Explain the difference of the two methods.



## Calculations on Data Frames

- __Step 7__: Count the number of duplicates in the loci column with the `table` function and append the result to the data frame with the `cbind` function.


{% highlight r %}
mycounts <- table(my_mw_target3[,1])[my_mw_target3[,1]]
my_mw_target4 <- cbind(my_mw_target3, Freq=mycounts[as.character(my_mw_target3[,1])]) 
{% endhighlight %}

- __Step 8__: Perform a vectorized devision of columns 3 and 4 (average AA weight per protein)


{% highlight r %}
data.frame(my_mw_target4, avg_AA_WT=(my_mw_target4[,3] / my_mw_target4[,4]))[1:2,5:11] 
{% endhighlight %}

{% highlight txt %}
##   Loc  cTP   mTP    SP other Freq avg_AA_WT
## 1   _ 0.10 0.090 0.075 0.925    1  115.2121
## 2   * 0.01 0.636 0.158 0.448    2  114.6612
{% endhighlight %}

- __Step 9__: Calculate for each row the mean and standard deviation across several columns


{% highlight r %}
mymean <- apply(my_mw_target4[,6:9], 1, mean)
mystdev <- apply(my_mw_target4[,6:9], 1, sd, na.rm=TRUE)
data.frame(my_mw_target4, mean=mymean, stdev=mystdev)[1:2,5:12] 
{% endhighlight %}

{% highlight txt %}
##   Loc  cTP   mTP    SP other Freq   mean     stdev
## 1   _ 0.10 0.090 0.075 0.925    1 0.2975 0.4184595
## 2   * 0.01 0.636 0.158 0.448    2 0.3130 0.2818912
{% endhighlight %}

## Plotting Example

- __Step 10__: Generate scatter plot columns: 'MW' and 'Residues' 


{% highlight r %}
plot(my_mw_target4[1:500,3:4], col="red")
{% endhighlight %}

![](Rbasics_files/plot_example-1.png)

## Export Results and Run Entire Exercise as Script

- __Step 11__: Write the data frame `my_mw_target4` into a tab-delimited text file and inspect it in Excel.


{% highlight r %}
write.table(my_mw_target4, file="my_file.xls", quote=F, sep="\t", col.names = NA) 
{% endhighlight %}

- __Homework 3H__: Write all commands from this exercise into an R script named `exerciseRbasics.R`, or download it from [here](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/exerciseRbasics.R). Then execute the script with the `source` function like this: `source("exerciseRbasics.R")`. This will run all commands of this exercise and generate the corresponding output files in the current working directory.


{% highlight r %}
source("exerciseRbasics.R")
{% endhighlight %}


