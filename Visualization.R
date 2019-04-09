library(ggplot2)
library(devtools)
library(gdxtools)

setwd("C:/Users/pietr/Desktop/Oil&gas/Progetto/PeakFFmodel")

res = gdx('resultsshort.gdx')
re = res$parameters
re$PD_res
