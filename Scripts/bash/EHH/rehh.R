#!/usr/bin/env Rscript
arguments=commandArgs(trailingOnly=TRUE)
#arguments<-c(10937139,1)
##########
#
# Extended haplotype homozygosity
# 2023-06-07
#
##########

setwd("$PATH-TO-SERVER/A_thaliana/fastphase/")

#packages and custom function
library(rehh)
library(gplots)
source("$PATH-TO-SERVER/private/Athaliana_flowersize/rehh/my_calc_ehh.R")

#files
thefiles<-grep(arguments[1] ,list.files(),value = T)
haplohh<-data2haplohh(hap_file=grep("switch",thefiles,value = T),map_file = grep(paste0(arguments[1],".bimlike"),thefiles,value = T))
map<-read.table(grep(paste0(arguments[1],".bimlike"),thefiles,value = T))

#position of interes
mrk<-which(map$V3==arguments[1])

# computation and plot
pdf(file = paste0("EHH_plot.",arguments[2],".",arguments[1],".pdf"))
res.ehh<-my_calc_ehh(haplohh,mrk=mrk)
dev.off()

