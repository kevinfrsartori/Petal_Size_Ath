#!/usr/bin/env Rscript
arguments=commandArgs(trailingOnly=TRUE)
#arguments<-c("4","ab")
##########
#
# Extended haplotype homozygosity
# scan 
# 2023-06-07
#
##########

setwd("$PATH-TO-SERVER/A_thaliana/fastphase/")

#packages and custom function
library(rehh)
library(gplots)
#source("$PATH-TO-SERVER/private/Athaliana_flowersize/rehh/my_calc_ehh.R")


##### once everything computed

#files
myhaplohh<-data2haplohh(hap_file = paste0("Chrm.",arguments[1],".",arguments[2],"_hapguess_switch.out"),map_file = paste0("SNP_1001g_filtered.",arguments[1],".",arguments[2],".bimlike") )
#scan
scanhh<-scan_hh(haplohh = myhaplohh,limhaplo = 2,limehh = 0.05,limehhs = 0.05)
# save data
write.table(x = scanhh,file = paste0("scanEHH.chrm-",arguments[1],".",arguments[2],".txt"),quote = F,row.names = F,col.names = T)


