#! /usr/bin/env Rscript


# https://github.com/IRkernel/IRkernel
install.packages('IRkernel', lib = .libPaths()[1])
IRkernel::installspec()  # to register the kernel in the current R installation
# jupyter labextension install @techrah/text-shortcuts  # for RStudio’s shortcuts
