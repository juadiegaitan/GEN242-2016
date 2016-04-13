---
title: "Programming in R" 
author: "Author: Thomas Girke"
date: "Last update: 13 April, 2016" 
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
echo "rmarkdown::render('Programming_in_R.Rmd', clean=F)" | R -slave; R CMD Stangle Programming_in_R.Rmd; Rscript ../md2jekyll.R Programming_in_R.knit.md sidebartitle='Programming in R' 9

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

One of the main attractions of using the R
([http://cran.at.r-project.org](http://cran.at.r-project.org)) environment is
the ease with which users can write their own programs and custom functions.
The R programming syntax is extremely easy to learn, even for users with no
previous programming experience. Once the basic R programming control
structures are understood, users can use the R language as a powerful
environment to perform complex custom analyses of almost any type of data.


## Why Programmin in R?

* Powerful statistical environment and programming language
* Facilitates reproducible research
* Efficient data structures make programming very easy
* Ease of implementing custom functions
* Powerful graphics
* Access to fast growing number of analysis packages
* Most widely used language in bioinformatics
* Is standard for data mining and biostatistical analysis
* Technical advantages: free, open-source, available for all OSs


## R Basics 

The previous [Rbasics](http://girke.bioinformatics.ucr.edu/GEN242/mydoc/mydoc_Rbasics_01.html) provides a general introduction to the usage of the R environment and its basic command syntax.
More details can be found R & BioConductor manual [here](http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual).

## Code Editors for R

Several excellent code editors are available that provide functionalities like R syntax highlighting, auto code indenting and utilities to send code/functions to the R console.

* [RStudio](https://www.rstudio.com/products/rstudio/features/): GUI-based IDE for R
* [Vim-R-Tmux](http://manuals.bioinformatics.ucr.edu/home/programming-in-r/vim-r): R working environment based on vim and tmux
* [Emacs](http://www.xemacs.org/Download/index.html) ([ESS add-on package](http://ess.r-project.org/))
* [gedit](https://wiki.gnome.org/Apps/Gedit) and [Rgedit](https://wiki.gnome.org/Apps/Gedit)
* [RKWard](https://rkward.kde.org/)
* [Eclipse](http://www.walware.de/goto/statet)
* [Tinn-R](http://jekyll.math.byuh.edu/other/howto/tinnr/install.shtml)
* [Notepad++ (NppToR)](https://sourceforge.net/projects/npptor/)

<center> Programming in R using Vim or Emacs</center>
<center><img title="vim-r" src="images/vimR.png"/></center>

<center>Programming in R using RStudio</center>
<center><img title="R_Interfaces" src="images/rstudio.png"/></center>

## Finding Help

Reference list on R programming (selection)

* [Advanced R](http://adv-r.had.co.nz/), by Hadley Wickham
* [R Programming for Bioinformatics](http://master.bioconductor.org/help/publications/books/r-programming-for-bioinformatics/), by Robert Gentleman
* [S Programming](http://www.stats.ox.ac.uk/pub/MASS3/Sprog/), by W. N. Venables and B. D. Ripley
* [Programming with Data](http://www.amazon.com/Programming-Data-Language-Lecture-Economics/dp/0387985034), by John M. Chambers
* [R Help](http://www1.maths.lth.se/help/R/) & [R Coding Conventions](http://www1.maths.lth.se/help/R/RCC/), Henrik Bengtsson, Lund University
* [Programming in R](http://zoonek2.free.fr/UNIX/48_R/02.html) (Vincent Zoonekynd)
* [Peter's R Programming Pages](http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r), University of Warwick
* [Rtips](http://pj.freefaculty.org/R/statsRus.html), Paul Johnsson, University of Kansas
* [R for Programmers](http://heather.cs.ucdavis.edu/~matloff/r.html), Norm Matloff, UC Davis
* [High-Performance R](http://www.statistik.uni-dortmund.de/useR-2008/tutorials/useR2008introhighperfR.pdf), Dirk Eddelbuettel tutorial presented at useR-2008
* [C/C++ level programming for R](http://www.stat.harvard.edu/ccr2005/index.html), Gopi Goswami



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

# Control Structures

## Startup and Closing Behavior


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


# Session Info


```r
sessionInfo()
```

```
## R version 3.2.4 Revised (2016-03-16 r70336)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.4 LTS
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
##  [1] Rcpp_0.12.3      codetools_0.2-14 digest_0.6.9     plyr_1.8.3       grid_3.2.4      
##  [6] gtable_0.1.2     formatR_1.2.1    magrittr_1.5     evaluate_0.8     scales_0.3.0    
## [11] stringi_1.0-1    rmarkdown_0.9.2  tools_3.2.4      stringr_1.0.0    munsell_0.4.2   
## [16] yaml_2.1.13      colorspace_1.2-6 htmltools_0.3    knitr_1.12
```

# References

