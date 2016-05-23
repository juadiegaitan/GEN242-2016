---
title: R Markdown
keywords: 
last_updated: Mon May 23 14:05:23 2016
---

## Overview

R Markdown combines markdown (an easy to write plain text format)
with embedded R code chunks. When compiling R Markdown reports, the 
code components can be evaluated so that both the code and its output 
can be included in the final document. This way reports are highly 
reproducible by allowing to automatically regenerate them when the 
underlying R code or data changes. R Markdown documents can be rendered to
various formats including HTML and PDF. For rendering documents, the environment
uses `knitr` and `pandoc`. Historically, R Markdown is an extension of the 
older `Sweave/Latex` environment.


