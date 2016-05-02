---
title: Combine annotation results among samples
keywords: 
last_updated: Sun May  1 20:44:33 2016
---

To simplify comparisons among samples, the \Rfunction{combineVarReports}
function combines all variant annotation reports referenced in a
`SYSargs` instance (here `args`). At the same time the function
allows to consider only certain feature types of interest. For instance, the
below setting `filtercol=c(Consequence="nonsynonymous")` will include
only nonsysynonymous variances listed in the `Consequence` column of
the annotation reports. To omit filtering, one can use the setting
`filtercol="All"`.

## Combine results from `GATK`  


{% highlight r %}
args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_gatk_filtered.txt")
combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_gatk.xls", quote=FALSE, row.names=FALSE, sep="\t")
{% endhighlight %}

## Combine results from `BCFtools`  


{% highlight r %}
args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_sambcf.xls", quote=FALSE, row.names=FALSE, sep="\t")
{% endhighlight %}

## Combine results from `VariantTools`


{% highlight r %}
args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_vartools_filtered.txt")
combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_vartools.xls", quote=FALSE, row.names=FALSE, sep="\t")
{% endhighlight %}

