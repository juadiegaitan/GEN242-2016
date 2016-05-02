---
title: Summary statistics of variants
keywords: 
last_updated: Sun May  1 20:44:33 2016
---

The `varSummary` function counts the number of variants for each feature type
included in the anntation reports.

## Summary for `GATK`


{% highlight r %}
args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_gatk_filtered.txt")
write.table(varSummary(args), "./results/variantStats_gatk.xls", quote=FALSE, col.names = NA, sep="\t")
{% endhighlight %}

## Summary for `BCFtools`


{% highlight r %}
args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
write.table(varSummary(args), "./results/variantStats_sambcf.xls", quote=FALSE, col.names = NA, sep="\t")
{% endhighlight %}

## Summary for `VariantTools`  


{% highlight r %}
args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_vartools_filtered.txt")
write.table(varSummary(args), "./results/variantStats_vartools.xls", quote=FALSE, col.names = NA, sep="\t")
{% endhighlight %}

