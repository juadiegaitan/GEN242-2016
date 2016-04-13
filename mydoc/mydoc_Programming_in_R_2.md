---
title: Control Structures
keywords: 
last_updated: Wed Apr 13 12:24:22 2016
---

## Startup and Closing Behavior


{% highlight r %}
q()  
{% endhighlight %}
{% highlight txt %}
Save workspace image? [y/n/c]:
{% endhighlight %}
        
* __Note__:
    When responding with `y`, then the entire R workspace will be written to
    the `.RData` file which can become very large. Often it is sufficient to just
    save an analysis protocol in an R source file. This way one can quickly
    regenerate all data sets and objects. 


## Navigating directories

Create an object with the assignment operator `<-` or `=`

{% highlight r %}
object <- ...
{% endhighlight %}

List objects in current R session

{% highlight r %}
ls()
{% endhighlight %}

Return content of current working directory

{% highlight r %}
dir()
{% endhighlight %}

Return path of current working directory

{% highlight r %}
getwd()
{% endhighlight %}

Change current working directory

{% highlight r %}
setwd("/home/user")
{% endhighlight %}

