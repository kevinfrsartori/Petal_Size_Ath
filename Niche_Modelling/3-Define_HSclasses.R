####################
#
# Define habitat suitability classes
# KS - 2023-11-15
#
####################

filt_acc<-read.table("Genetics/1001genomes_snp_2all_mac1_rel.95.fam")
all_acc<-read.table("Genetics/accessions.txt",h=F,sep=",")
HS<-stack("Habitat_suitability_Ath_2023-04-27.grd")
