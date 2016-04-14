---
title: Building R Packages
keywords: 
last_updated: Wed Apr 13 16:38:00 2016
---

## Short Overview of Package Building Process

Automatic package building with the package.skeleton function. The given example will create a directory named `mypackage` containing the skeleton of the package for all functions, methods and classes defined in the R script(s) passed on to the `code_files` argument. The basic structure of the package directory is described [here](http://manuals.bioinformatics.ucr.edu/home/programming-in-r#Progr_pack). The package directory will also contain a file named `Read-and-delete-me` with instructions for completing the package:


{% highlight r %}
package.skeleton(name="mypackage", code_files=c("script1.R", "script2.R"))
{% endhighlight %}

Once a package skeleton is available one can build the package from the command-line (Linux/OS X). This will create a tarball of the package with its version number encoded in the file name. Subequently, the package tarball needs to be checked for errors with:


{% highlight r %}
R CMD build mypackage
R CMD check mypackage_1.0.tar.gz
{% endhighlight %}

Install package from source

{% highlight r %}
install.packages("mypackage_1.0.tar.gz", repos=NULL) 
{% endhighlight %}

For more details see [here](http://manuals.bioinformatics.ucr.edu/home/programming-in-r#TOC-Building-R-Packages)

## Exercise - build your own R package

__Task 1__: Save one or more of your functions to a file called `script.R` and build the package with the `package.skeleton` function.


{% highlight r %}
package.skeleton(name="mypackage", code_files=c("script1.R"), namespace=TRUE)
{% endhighlight %}

__Task 2__: Build tarball of the package


{% highlight r %}
system("R CMD build mypackage")
{% endhighlight %}

__Task 3__: Install and use package


{% highlight r %}
install.packages("mypackage_1.0.tar.gz", repos=NULL, type="source")
library(mypackage)
?myMAcomp # Opens help for function defined by mypackage
{% endhighlight %}

