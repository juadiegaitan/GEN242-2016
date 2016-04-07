---
title: "Introduction_to_R" 
author: "Author: Thomas Girke"
date: "Last update: 07 April, 2016" 
output:
  BiocStyle::html_document:
    toc: true
    toc_depth: 3
    fig_caption: yes

fontsize: 14pt
bibliography: bibtex.bib
---
<!--
%% \VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{Overview Vignette}
%% \VignetteDepends{methods}
%% \VignetteKeywords{compute cluster, pipeline, reports}
%% \VignettePackage{longevityTools}
-->

<!---
- Compile from command-line
echo "rmarkdown::render('Rbasics.Rmd', clean=F)" | R -slave; R CMD Stangle Rbasics.Rmd; Rscript ../md2jekyll.R Rbasics.knit.md 8

- Commit to github
git commit -am "some edits"; git push -u origin master

- To customize font size and other style features, add this line to output section in preamble:  
    css: style.css
-->

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>



# Overview

## What is R?

[R](http://cran.at.r-project.org) is a powerful statistical environment and
programming language for the analysis and visualization of data.  The
associated [Bioconductor](http://bioconductor.org/) and CRAN package
repositories provide many additional R packages for statistical data analysis
for a wide array of research areas. The R software is free and runs on all
common operating systems. 

## Why Using R?
* Complete statistical environment and programming language
* Efficient functions and data structures for data analysis
* Powerful graphics
* Access to fast growing number of analysis packages
* Most widely used language in bioinformatics
* Is standard for data mining and biostatistical analysis
* Technical advantages: free, open-source, available for all OSs

## Books and Documentation
* simpleR - Using R for Introductory Statistics (John Verzani, 2004) - [URL](http://cran.r-project.org/doc/contrib/Verzani-SimpleR.pdf)
* Bioinformatics and Computational Biology Solutions Using R and Bioconductor (Gentleman et al., 2005) - [URL](http://www.bioconductor.org/help/publications/books/bioinformatics-and-computational-biology-solutions/)
* More on this see "Finding Help" section in UCR Manual - [URL](http://manuals.bioinformatics.ucr.edu/home/R\_BioCondManual\#TOC-Finding-Help)

## R Working Environments

## Working environments (IDEs) for R
<center><img title="R_Interfaces" src="images/rinterface.png"/></center>
<center> R Projects and Interfaces</center>

Some R working environments with support for syntax highlighting and utilities to send code 
to the R console: 

* [RStudio](https://www.rstudio.com/products/rstudio/features): excellent choice for beginners ([Cheat Sheet](http://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf)) 
* Basic R code editors provided by Rguis 
* [gedit](https://wiki.gnome.org/Apps/Gedit), [Rgedit](http://rgedit.sourceforge.net/), [RKWard](https://rkward.kde.org/), [Eclipse](http://www.walware.de/goto/statet), [Tinn-R](http://www.sciviews.org/Tinn-R/), [Notepad++](https://notepad-plus-plus.org/), [NppToR](http://sourceforge.net/projects/npptor/)
* [Vim-R-Tmux](http://manuals.bioinformatics.ucr.edu/home/programming-in-r/vim-r): R working environment based on vim and tmux 
* [Emacs](http://www.xemacs.org/Download/index.html) ([ESS add-on package](http://ess.r-project.org/))
	
### Example: RStudio 

New integrated development environment (IDE) for [R](http://www.rstudio.com/ide/download/). Highly functional for both beginners and 
advanced.

<center><img title="RStudio" src="images/rstudio.png"/></center>
<center> RStudio IDE</center>

Some userful shortcuts: `Ctrl+Enter` (send code), `Ctrl+Shift+C` (comment/uncomment), `Ctrl+1/2` (switch window focus)

### Example: Vim-R-Tmux

Terminal-based Working Environment for R: [Vim-R-Tmux](http://manuals.bioinformatics.ucr.edu/home/programming-in-r/vim-r)

<center><img title="Vim-R-Tmux" src="images/screenshot.png" ></center>
<center>Vim-R-Tmux IDE for R</center>

# R Package Repositories

* CRAN (>8,000 packages) general data analysis - [URL](http://cran.at.r-project.org/)
* Bioconductor (>1,100 packages) bioscience data analysis - [URL](http://www.bioconductor.org/)
* Omegahat (>90 packages) programming interfaces - [URL](https://github.com/omegahat?tab=repositories)

# Installation of R and Add-on Packages

(1.) Install R for your operating system from [CRAN](http://cran.at.r-project.org/).

(2.) Install RStudio from [RStudio](http://www.rstudio.com/ide/download).

(3.) Install CRAN Packages from R console like this:


```r
install.packages(c("pkg1", "pkg2")) 
install.packages("pkg.zip", repos=NULL)
```

(4.) Install Bioconductor packages as follows:


```r
source("http://www.bioconductor.org/biocLite.R")
library(BiocInstaller)
BiocVersion()
biocLite()
biocLite(c("pkg1", "pkg2"))
```

(5.) For more details consult the [Bioc Install page](http://www.bioconductor.org/install/)
and [BiocInstaller](http://www.bioconductor.org/packages/release/bioc/html/BiocInstaller.html) package.

# Getting Around

## Startup and Closing Behavior

* __Starting R__:
    The R GUI versions, including RStudio, under Windows and Mac OS X can be
    opened by double-clicking their icons. Alternatively, one can start it by
    typing `R` in a terminal (default under Linux). 

* __Startup/Closing Behavior__:
    The R environment is controlled by hidden files in the startup directory:
    `.RData`, `.Rhistory` and `.Rprofile` (optional). 
	
    
* __Closing R__:


```r
q()  
```
```
Save workspace image? [y/n/c]:
```
        
* __Note__:
    When responding with `y`, then the entire R workspace will be written to
    the `.RData` file which can become very large. Often it is sufficient to just
    save an analysis protocol in an R source file. This way one can quickly
    regenerate all data sets and objects. 


## Navigating directories

Create an object with the assignment operator `<-` or `=`

```r
object <- ...
```

List objects in current R session

```r
ls()
```

Return content of current working directory

```r
dir()
```

Return path of current working directory

```r
getwd()
```

Change current working directory

```r
setwd("/home/user")
```

# Basic Syntax

General R command syntax


```r
object <- function_name(arguments) 
object <- object[arguments] 
```

Finding help


```r
?function_name
```

Load a library/package


```r
library("my_library") 
```

List functions defined by a library


```r
library(help="my_library")
```

Load library manual (PDF or HTML file)


```r
vignette("my_library") 
```

Execute an R script from within R


```r
source("my_script.R")
```

Execute an R script from command-line (the first of the three options is preferred)


```sh
$ Rscript my_script.R
$ R CMD BATCH my_script.R 
$ R --slave < my_script.R 
```

# Data Types 

## Numeric data

Example: `1, 2, 3, ...`


```r
x <- c(1, 2, 3)
x
```

```
## [1] 1 2 3
```

```r
is.numeric(x)
```

```
## [1] TRUE
```

```r
as.character(x)
```

```
## [1] "1" "2" "3"
```

## Character data

Example: `"a", "b", "c", ...`


```r
x <- c("1", "2", "3")
x
```

```
## [1] "1" "2" "3"
```

```r
is.character(x)
```

```
## [1] TRUE
```

```r
as.numeric(x)
```

```
## [1] 1 2 3
```

## Complex data

Example: mix of both


```r
c(1, "b", 3)
```

```
## [1] "1" "b" "3"
```

## Logical data

Example: `TRUE` of `FALSE`


```r
x <- 1:10 < 5
x  
```

```
##  [1]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
```

```r
!x
```

```
##  [1] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
```

```r
which(x) # Returns index for the 'TRUE' values in logical vector
```

```
## [1] 1 2 3 4
```

# Data objects

## Object types

## Vectors (1D)

Example: `numeric` or `character`


```r
myVec <- 1:10; names(myVec) <- letters[1:10]  
myVec[1:5]
```

```
## a b c d e 
## 1 2 3 4 5
```

```r
myVec[c(2,4,6,8)]
```

```
## b d f h 
## 2 4 6 8
```

```r
myVec[c("b", "d", "f")]
```

```
## b d f 
## 2 4 6
```

## Factors (1D)

Example: vectors with grouping information


```r
factor(c("dog", "cat", "mouse", "dog", "dog", "cat"))
```

```
## [1] dog   cat   mouse dog   dog   cat  
## Levels: cat dog mouse
```

## Matrices (2D)

Example: two dimensional structures with data of same type


```r
myMA <- matrix(1:30, 3, 10, byrow = TRUE) 
class(myMA)
```

```
## [1] "matrix"
```

```r
myMA[1:2,]
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## [1,]    1    2    3    4    5    6    7    8    9    10
## [2,]   11   12   13   14   15   16   17   18   19    20
```

```r
myMA[1, , drop=FALSE]
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## [1,]    1    2    3    4    5    6    7    8    9    10
```

## Data Frames (2D)

Example: two dimensional objects with data of variable types


```r
myDF <- data.frame(Col1=1:10, Col2=10:1) 
myDF[1:2, ]
```

```
##   Col1 Col2
## 1    1   10
## 2    2    9
```

## Arrays

Example: data structure with one, two or more dimensions


## Lists

Example: containers for any object type


```r
myL <- list(name="Fred", wife="Mary", no.children=3, child.ages=c(4,7,9)) 
myL
```

```
## $name
## [1] "Fred"
## 
## $wife
## [1] "Mary"
## 
## $no.children
## [1] 3
## 
## $child.ages
## [1] 4 7 9
```

```r
myL[[4]][1:2] 
```

```
## [1] 4 7
```

## Functions

Example: piece of code


```r
myfct <- function(arg1, arg2, ...) { 
	function_body 
}
```

## Subsetting of data objects

__(1.) Subsetting by positive or negative index/position numbers__


```r
myVec <- 1:26; names(myVec) <- LETTERS 
myVec[1:4]
```

```
## A B C D 
## 1 2 3 4
```

__(2.) Subsetting by same length logical vectors__


```r
myLog <- myVec > 10
myVec[myLog] 
```

```
##  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z 
## 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
```

__(3.) Subsetting by field names__


```r
myVec[c("B", "K", "M")]
```

```
##  B  K  M 
##  2 11 13
```

__(4.) Subset with `$` sign__: references a single column or list component by its name 


```r
iris$Species[1:8]
```

```
## [1] setosa setosa setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica
```

# Important Utilities
	
## Combining Objects

The `c` function combines vectors and lists


```r
c(1, 2, 3)
```

```
## [1] 1 2 3
```

```r
x <- 1:3; y <- 101:103
c(x, y)
```

```
## [1]   1   2   3 101 102 103
```

```r
iris$Species[1:8]
```

```
## [1] setosa setosa setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica
```

The `cbind` and `rbind` functions can be used to append columns and rows, respecively.

```r
ma <- cbind(x, y)
ma
```

```
##      x   y
## [1,] 1 101
## [2,] 2 102
## [3,] 3 103
```

```r
rbind(ma, ma)
```

```
##      x   y
## [1,] 1 101
## [2,] 2 102
## [3,] 3 103
## [4,] 1 101
## [5,] 2 102
## [6,] 3 103
```

## Accessing Dimensions of Objects

Length and dimension information of objects


```r
length(iris$Species)
```

```
## [1] 150
```

```r
dim(iris)
```

```
## [1] 150   5
```

## Accessing Name Slots of Objects

Accessing row and column names of 2D objects

```r
rownames(iris)[1:8]
```

```
## [1] "1" "2" "3" "4" "5" "6" "7" "8"
```

```r
colnames(iris)
```

```
## [1] "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"
```

Return name field of vectors and lists

```r
names(myVec)
```

```
##  [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X"
## [25] "Y" "Z"
```

```r
names(myL)
```

```
## [1] "name"        "wife"        "no.children" "child.ages"
```




# Graphics example

Plotting example

```r
barplot(1:10, col="green")
```

![](Rbasics_files/figure-html/plot_example-1.png)\

# Session Info


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.3 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
##  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  utils     datasets  grDevices methods   base     
## 
## other attached packages:
## [1] ggplot2_2.0.0   limma_3.26.3    BiocStyle_1.8.0
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.3      codetools_0.2-14 digest_0.6.9     plyr_1.8.3       grid_3.2.3      
##  [6] gtable_0.1.2     formatR_1.2.1    magrittr_1.5     evaluate_0.8     scales_0.3.0    
## [11] stringi_1.0-1    rmarkdown_0.9.2  tools_3.2.3      stringr_1.0.0    munsell_0.4.2   
## [16] yaml_2.1.13      colorspace_1.2-6 htmltools_0.3    knitr_1.12
```

# References

