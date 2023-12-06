####################
#
# Define habitat suitability classes
# KS - 2023-11-15
#
####################

filt_acc<-read.table("Genetics/1001genomes_snp_2all_mac1_rel.95.fam")
all_acc<-read.table("Niche_Modelling/accessions_1001g_habitatsuitability.csv",h=T,sep=",")
hist(all_acc$HS)
filt_acc<-all_acc[which(all_acc$accession_name %in% filt_acc$V1),]
hist(all_acc$HS)

# IMPORTANT - check previous script and remove few weird accessions
# for now I just remove HS < .2
filt_acc<-filt_acc[-which(filt_acc$HS < .2),]

quantile(filt_acc$HS,c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
#       0%       10%       20%       30%       40%       50%       60%       70%       80%      90%      100% 
#0.2005882 0.4928905 0.5408645 0.5856892 0.6070545 0.6172905 0.6363386 0.6585770 0.6785327 0.7004914 0.8035722 

g01<-filt_acc$accession_name[which(filt_acc$HS <= quantile(filt_acc$HS,.1))]
g02<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.1) & filt_acc$HS <= quantile(filt_acc$HS,.2))]
g03<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.2) & filt_acc$HS <= quantile(filt_acc$HS,.3))]
g04<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.3) & filt_acc$HS <= quantile(filt_acc$HS,.4))]
g05<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.4) & filt_acc$HS <= quantile(filt_acc$HS,.5))]
g06<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.5) & filt_acc$HS <= quantile(filt_acc$HS,.6))]
g07<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.6) & filt_acc$HS <= quantile(filt_acc$HS,.7))]
g08<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.7) & filt_acc$HS <= quantile(filt_acc$HS,.8))]
g09<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.8) & filt_acc$HS <= quantile(filt_acc$HS,.9))]
g10<-filt_acc$accession_name[which(filt_acc$HS > quantile(filt_acc$HS,.9) & filt_acc$HS <= quantile(filt_acc$HS,1))]

write.table(x = data.frame(V1=g01,v2=g01),file = "Genetics/accessions_HS01.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g02,v2=g02),file = "Genetics/accessions_HS02.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g03,v2=g03),file = "Genetics/accessions_HS03.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g04,v2=g04),file = "Genetics/accessions_HS04.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g05,v2=g05),file = "Genetics/accessions_HS05.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g06,v2=g06),file = "Genetics/accessions_HS06.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g07,v2=g07),file = "Genetics/accessions_HS07.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g08,v2=g08),file = "Genetics/accessions_HS08.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g09,v2=g09),file = "Genetics/accessions_HS09.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g10,v2=g10),file = "Genetics/accessions_HS10.txt",quote = F,col.names = F,row.names = F)
