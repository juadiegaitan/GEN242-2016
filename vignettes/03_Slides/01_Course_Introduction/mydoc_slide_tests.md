---
title: Slide Tests
last_updated: 08-Mar-16
---

## Test: slide embedding 1 for Google Slides

<iframe src="https://docs.google.com/presentation/d/19odvhFkVXjkEFQdJyH4w2gAp6Jq8yAYyogTp3jgx9nc/embed?start=false&loop=true&delayms=60000" frameborder="0" width="960" height="569" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>


## Test: slide embedding 2 for ioslides (HTML)

<iframe src="https://htmlpreview.github.io/?https://github.com/tgirke/systemPipeR/master/inst/extdata/slides/systemPipeRslides.html" frameborder="0" width="960" height="750" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

## Test: slide embedding 3 for PDF slides (e.g. Beamer or anything)

<iframe src="http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_5_8_2014/Rbasics/Introduction_into_R.pdf" frameborder="0" width="600" height="400" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>


## Code could go here

General R command syntax


{% highlight r %}
object <- function_name(arguments) 
object <- object[arguments] 
{% endhighlight %}

Finding help


{% highlight r %}
?function_name
{% endhighlight %}

Load a library/package


{% highlight r %}
library("my_library") 
{% endhighlight %}

List functions defined by a library


{% highlight r %}
library(help="my_library")
{% endhighlight %}

Load library manual (PDF or HTML file)


{% highlight r %}
vignette("my_library") 
{% endhighlight %}

Execute an R script from within R


{% highlight r %}
source("my_script.R")
{% endhighlight %}

Execute an R script from command-line (the first of the three options is preferred)


{% highlight sh %}
$ Rscript my_script.R
$ R CMD BATCH my_script.R 
$ R --slave < my_script.R 
{% endhighlight %}

