# TODO: add in functions using it, change to ::
#' @import BiocManager
#' @import dplyr
#' @import magrittr
#' @import Biostrings
#' @import genbankr
#' @import coRdon
#' @import stringr
#' @import grDevices
#' @import reutils
#' @import pkgload
#' @import shape
#' @import diagram
#' @import gplots


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(dplyr, quietly = TRUE)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(magrittr, quietly = TRUE)){
  install.packages("magrittr")
  library(magrittr)
}
if(!require(Biostrings, quietly = TRUE)){
  BiocManager::install("Biostrings")
  library(Biostrings)
}
if(!require(genbankr, quietly = TRUE)){
  BiocManager::install("genbankr")
  library(genbankr)
}
if(!require(coRdon, quietly = TRUE)){
  BiocManager::install("coRdon")
  library(coRdon)
}
if(!require(stringr, quietly = TRUE)){
  install.packages("stringr")
  library(stringr)
}
if(!require(grDevices, quietly = TRUE)){
  install.packages("grDevices")
  library(grDevices)
}
if(!require(reutils, quietly = TRUE)){
  install.packages("reutils")
  library(reutils)
}
if(!require(pkgload, quietly = TRUE)){
  install.packages("pkgload")
  library(pkgload)
}
if(!require(shape, quietly = TRUE)){
  install.packages("shape")
  library(shape)
}
if(!require(diagram, quietly = TRUE)){
  install.packages("diagram")
  library(diagram)
}
if(!require(gplots, quietly = TRUE)){
  install.packages("gplots")
  library(gplots)
}
