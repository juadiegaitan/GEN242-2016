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

## ----install_cran, eval=FALSE--------------------------------------------
## install.packages(c("pkg1", "pkg2"))
## install.packages("pkg.zip", repos=NULL)

## ----install_bioc, eval=FALSE--------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## library(BiocInstaller)
## BiocVersion()
## biocLite()
## biocLite(c("pkg1", "pkg2"))

## ----closing_r, eval=FALSE-----------------------------------------------
## q()

## ----r_assignment, eval=FALSE--------------------------------------------
## object <- ...

## ----r_ls, eval=FALSE----------------------------------------------------
## ls()

## ----r_dirshow, eval=FALSE-----------------------------------------------
## dir()

## ----r_dirpath, eval=FALSE-----------------------------------------------
## getwd()

## ----r_setwd, eval=FALSE-------------------------------------------------
## setwd("/home/user")

## ----r_syntax, eval=FALSE------------------------------------------------
## object <- function_name(arguments)
## object <- object[arguments]

## ----r_find_help, eval=FALSE---------------------------------------------
## ?function_name

## ----r_package_load, eval=FALSE------------------------------------------
## library("my_library")

## ----r_package_functions, eval=FALSE-------------------------------------
## library(help="my_library")

## ----r_load_vignette, eval=FALSE-----------------------------------------
## vignette("my_library")

## ----r_execute_script, eval=FALSE----------------------------------------
## source("my_script.R")

## ----sh_execute_script, eval=FALSE, engine="sh"--------------------------
## $ Rscript my_script.R
## $ R CMD BATCH my_script.R
## $ R --slave < my_script.R

## ----r_numeric_data, eval=TRUE-------------------------------------------

x <- c(1, 2, 3)
x
is.numeric(x)
as.character(x)

## ----r_character_data, eval=TRUE-----------------------------------------
x <- c("1", "2", "3")
x
is.character(x)
as.numeric(x)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

