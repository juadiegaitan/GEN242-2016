## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(limma) 
    library(ggplot2) }) 

## ----if_statement, eval=FALSE--------------------------------------------
## if(TRUE) {
## 	statements_1
## } else {
## 	statements_2
## }

## ----if_statement_example, eval=TRUE-------------------------------------
if(1==0) { 
	print(1) 
} else { 
	print(2) 
}

## ----ifelse_statement, eval=FALSE----------------------------------------
## ifelse(test, true_value, false_value)

## ----ifelse_statement_example, eval=TRUE---------------------------------
x <- 1:10 
ifelse(x<5, x, 0)

## ----for_loop, eval=FALSE------------------------------------------------
## for(variable in sequence) {
## 	statements
## }

## ----for_loop_example, eval=TRUE-----------------------------------------
mydf <- iris
myve <- NULL
for(i in seq(along=mydf[,1])) {
	myve <- c(myve, mean(as.numeric(mydf[i,1:3])))
}
myve[1:8]

## ----for_loop_inject_example, eval=TRUE----------------------------------
myve <- numeric(length(mydf[,1]))
for(i in seq(along=myve)) {
	myve[i] <- mean(as.numeric(mydf[i,1:3]))
}
myve[1:8]

## ----for_loop_stop_example, eval=FALSE-----------------------------------
## x <- 1:10
## z <- NULL
## for(i in seq(along=x)) {
## 	if(x[i] < 5) {
## 		z <- c(z, x[i]-1)
## 	} else {
## 		stop("values need to be < 5")
## 	}
## }

## ----while_loop, eval=FALSE----------------------------------------------
## while(condition) {
## 	statements
## }

## ----while_loop_example, eval=TRUE---------------------------------------
z <- 0
while(z<5) { 
	z <- z + 2
	print(z)  
}

## ----apply_loop, eval=FALSE----------------------------------------------
## apply(X, MARGIN, FUN, ARGs)

## ----apply_loop_example, eval=TRUE---------------------------------------
apply(iris[1:8,1:3], 1, mean)

## ----tapply_loop, eval=FALSE---------------------------------------------
## tapply(vector, factor, FUN)

## ----tapply_loop_example, eval=TRUE--------------------------------------
iris[1:2,]
tapply(iris$Sepal.Length, iris$Species, mean)

## ----lapply_loop_example, eval=TRUE--------------------------------------
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
lapply(x, mean)
sapply(x, mean)

## ----lapply_loop_fct_example, eval=FALSE---------------------------------
## lapply(names(x), function(x) mean(x))
## sapply(names(x), function(x) mean(x))

## ----sessionInfo---------------------------------------------------------
sessionInfo()

