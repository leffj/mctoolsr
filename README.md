# mctoolsr
Microbial community analysis tools in R


---
title: "mctoolsr examples"
author: "Jon Leff"
date: "July 7, 2015"
output: pdf_document
---

This document serves as a brief introduction to using **mctoolsr**. This document will go through getting the pro-package working and a few examples using the most popular functions.

### Getting and using **mctoolsr**

**mctoolsr** is available on Github at: https://github.com/leffj/mctoolsr

To use:

1. [Install git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) if you haven't already
2. Open a terminal window and clone **mctoolsr** in a directory where you keep software on your machine
```
git clone https://github.com/leffj/mctoolsr.git
```
3. When using **mctoolsr** in an R script, source the "routine_analysis_functions.R" file in the mctoolsr/R directory. For example, I include the following line at the top of all my R scripts using **mctoolsr**:
