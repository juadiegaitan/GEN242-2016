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
echo "rmarkdown::render('Programming_in_R.Rmd', clean=F)" | R -slave; R CMD Stangle Programming_in_R.Rmd; Rscript ../md2jekyll.R Programming_in_R.knit.md 9

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
environment to perform complex custom analyses of almost any type of data [@Gentleman2008-xo].


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


# Control Structures

## Important Operators

### Comparison operators

* `==` (equal)
* `!=` (not equal)
* `>` (greater than)
* `>=` (greater than or equal)
* `<` (less than)
* `<=` (less than or equal)

### Logical operators
		
* `&` (and)
* `|` (or) 
* `!` (not)

## Conditional Executions: `if` Statements

An `if` statement operates on length-one logical vectors.

__Syntax__

```r
if(TRUE) { 
	statements_1 
} else { 
	statements_2 
}
```

__Example__

```r
if(1==0) { 
	print(1) 
} else { 
	print(2) 
}
```

```
## [1] 2
```

## Conditional Executions: `ifelse` Statements`

The `ifelse` statement operates on vectors.

__Syntax__

```r
ifelse(test, true_value, false_value)
```
__Example__

```r
x <- 1:10 
ifelse(x<5, x, 0)
```

```
##  [1] 1 2 3 4 0 0 0 0 0 0
```

# Loops

## `for` loop

`for` loops iterate over elements of a looping vector.

__Syntax__

```r
for(variable in sequence) { 
	statements 
}
```
__Example__

```r
mydf <- iris
myve <- NULL
for(i in seq(along=mydf[,1])) {
	myve <- c(myve, mean(as.numeric(mydf[i,1:3])))
}
myve[1:8]
```

```
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

__Note:__ Inject into objecs is much faster than append approach with `c`, `cbind`, etc.

__Example__

```r
myve <- numeric(length(mydf[,1]))
for(i in seq(along=myve)) {
	myve[i] <- mean(as.numeric(mydf[i,1:3]))
}
myve[1:8]
```

```
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

### Conditional Stop of Loops

The `stop` function can be used to break out of a loop (or a function) when a condition becomes `TRUE`. In addition, an error message will be printed.

__Example__

```r
x <- 1:10
z <- NULL
for(i in seq(along=x)) { 
	if(x[i] < 5) { 
		z <- c(z, x[i]-1)  
	} else { 
		stop("values need to be < 5") 
	}
}
```

## `while` loop

Iterates as long as a condition is true.

__Syntax__

```r
while(condition) {
	statements
}
```

__Example__

```r
z <- 0
while(z<5) { 
	z <- z + 2
	print(z)  
}
```

```
## [1] 2
## [1] 4
## [1] 6
```

## The `apply` Function Family

### `apply`

__Syntax__

```r
apply(X, MARGIN, FUN, ARGs)
```

__Arguments__

* `X`: `array`, `matrix` or `data.frame`
* `MARGIN`: `1` for rows, `2` for columns
* `FUN`: one or more functions
* `ARGs`: possible arguments for functions

__Example__

```r
apply(iris[1:8,1:3], 1, mean)
```

```
##        1        2        3        4        5        6        7        8 
## 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

### `tapply`

Applies a function to vector components that are defined by a factor.

__Syntax__

```r
tapply(vector, factor, FUN)
```

__Example__

```r
iris[1:2,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
```

```r
tapply(iris$Sepal.Length, iris$Species, mean)
```

```
##     setosa versicolor  virginica 
##      5.006      5.936      6.588
```

### `sapply` and `lapply`

Both apply a function to vector or list objects. The `lapply` function always returns a list object, while `sapply` returns `vector` or `matrix` objects when it is possible. 

__Examples__

```r
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
lapply(x, mean)
```

```
## $a
## [1] 5.5
## 
## $beta
## [1] 4.535125
## 
## $logic
## [1] 0.5
```

```r
sapply(x, mean)
```

```
##        a     beta    logic 
## 5.500000 4.535125 0.500000
```

Often used in combination with a function definition:

```r
lapply(names(x), function(x) mean(x))
```

```
## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA
```

```
## [[1]]
## [1] NA
## 
## [[2]]
## [1] NA
## 
## [[3]]
## [1] NA
```

```r
sapply(names(x), function(x) mean(x))
```

```
## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA

## Warning in mean.default(x): argument is not numeric or logical: returning NA
```

```
##     a  beta logic 
##    NA    NA    NA
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

