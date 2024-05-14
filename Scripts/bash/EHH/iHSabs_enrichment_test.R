#
# Perform enrichment test from Arabidopsis-male-gamete paper
# Converted to two sided test
# KS - 2024-02-06
#

ihs<-read.table("ancestry447/Genome_wide_iHS.txt",h=T)
ihs$maf<-apply(X = ihs[,9:10],MARGIN = 1,FUN = min)
ihs<-ihs[which(ihs$maf>=0.05),]
# make two sided test
ihs$iHS_abs<-abs(ihs$iHS)
# keep one value per snp
ihs<-ihs[order(ihs$iHS_abs,decreasing = T),]
ihs<-ihs[-which(duplicated(ihs$snpID)),]

#traits
filename<-grep("SNP_1001g_filtered_",list.files("$PATH-TO-SERVER/private/Athaliana_flowersize/gwas/"),value = T)
traits<-unlist(lapply(strsplit(substr(filename,start = 20,stop = 50),split = ".assoc"),"[[",1))

#matrix to fill
enrich_matrix_flower_neg<-matrix(data = NA,nrow = length(traits),ncol = 4,
                                 dimnames = list(traits,c("0.1","0.05","0.025","0.01")))

for (t in 1:length(traits)) {
  cat(t)
gwas_flower<-read.table(paste0("$PATH-TO-SERVER/private/Athaliana_flowersize/gwas/SNP_1001g_filtered_",traits[t],".assoc.txt"),h=T)

enrich_ratio_flower_0.1_neg <- numeric(0)
enrich_ratio_flower_0.05_neg <- numeric(0)
enrich_ratio_flower_0.025_neg <- numeric(0)
enrich_ratio_flower_0.01_neg <- numeric(0)

sig_snp <- gwas_flower[(-log10(gwas_flower$p_lrt) > 4)&(gwas_flower[,7]>0.05),2]

top_win <- which(ihs$snpID %in% sig_snp)

# make a two sided test
enrich_ratio_flower_0.05_neg <-  length( ihs$iHS_abs[top_win][ ihs$iHS_abs[top_win] > sort(ihs$iHS_abs,decreasing = T)[ceiling(length(ihs$iHS_abs)/20)] ])/length(ihs$iHS_abs[top_win])
enrich_ratio_flower_0.1_neg <-   length( ihs$iHS_abs[top_win][ ihs$iHS_abs[top_win] > sort(ihs$iHS_abs,decreasing = T)[ceiling(length(ihs$iHS_abs)/10)] ])/length(ihs$iHS_abs[top_win])
enrich_ratio_flower_0.025_neg <- length( ihs$iHS_abs[top_win][ ihs$iHS_abs[top_win] > sort(ihs$iHS_abs,decreasing = T)[ceiling(length(ihs$iHS_abs)/40)] ])/length(ihs$iHS_abs[top_win])
enrich_ratio_flower_0.01_neg <-  length( ihs$iHS_abs[top_win][ ihs$iHS_abs[top_win] > sort(ihs$iHS_abs,decreasing = T)[ceiling(length(ihs$iHS_abs)/100)]])/length(ihs$iHS_abs[top_win])

enrich_matrix_flower_neg[t,] <- rbind(
  enrich_ratio_flower_0.1_neg/0.1,
  enrich_ratio_flower_0.05_neg/0.05,
  enrich_ratio_flower_0.025_neg/0.025,
  enrich_ratio_flower_0.01_neg/0.01
)
}

#enrich_matrix_flower_neg[(enrich_matrix_flower_neg==0)]  <- NA


#### generate the enrichment plot 
par(mfrow=c(1,1))
cutoff <- c(1,2,3,4)
plot(cutoff,enrich_matrix_flower_neg[1:4,1],type="n",
     ylim=c(0,13),ylab="Fold enrichment",xlab="Top X % of selection scan",main="",cex=1,xaxt="n")
for(i in 1:13){
  points(cutoff,enrich_matrix_flower_neg[i,1:4],type="l",col="grey")
}
points(cutoff,enrich_matrix_flower_neg[2,1:4],type="l",col="black",lwd=2)
points(cutoff,enrich_matrix_flower_neg[3,1:4],type="l",col="black",lwd=2)
points(cutoff,enrich_matrix_flower_neg[4,1:4],type="l",col="black",lwd=2)

abline(h=1,lty=2,lwd=2)
abline(v=1,lty=2,lwd=1)
abline(v=2,lty=2,lwd=1)
abline(v=3,lty=2,lwd=1)
abline(v=4,lty=2,lwd=1)

write.table(enrich_matrix_flower_neg,file = "enrich_alltraits.txt",quote = F,row.names = T,col.names = T)
