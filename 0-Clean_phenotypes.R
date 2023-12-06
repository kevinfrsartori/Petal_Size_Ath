###################
#
# Computing genetic ogan size values
# 2023-05-04 Kevin
#
###################

organs<-read.table("U_Shaped_Data.csv",h=T,sep = ",")
# fixing some misreading of labels
# 4749 is 4779 
# 6421 is 6424 
# 7527 is 7525
# 6947 is 6945
organs[which(organs$Genotype == 4779),11:13]<-organs[which(organs$Genotype == 4749),11:13]
organs<-organs[-which(organs$Genotype == 4749),]
organs[which(organs$Genotype == 6424),11:13]<-organs[which(organs$Genotype == 6421),11:13]
organs<-organs[-which(organs$Genotype == 6421),]
organs[which(organs$Genotype == 7525),11:13]<-organs[which(organs$Genotype == 7527),11:13]
organs<-organs[-which(organs$Genotype == 7527),]
organs$Genotype[which(organs$Genotype==6947)]<-6945

metadt<-read.table("data_mani_fleur2022_cleaned_KSedited_2023-05-04.csv",h=T,sep=";",na.strings = c("NA",""),dec = ",")

# make dates
metadt$date.semis<-as.Date(metadt$date.semis,format = "%d/%m/%Y")
metadt$date.transplantation<-as.Date(metadt$date.transplantation,format = "%d/%m/%Y")
metadt$date.fleur1<-as.Date(metadt$date.fleur1,format = "%d/%m/%Y")
metadt$date.recolte<-as.Date(metadt$date.recolte,format = "%d/%m/%Y")
metadt$date.death<-as.Date(metadt$date.death,format = "%d/%m/%Y")
# make julian dates
library(timeDate)
metadt$j.semis<-julian(metadt$date.semis,origin = as.Date("01/09/2020",format = "%d/%m/%Y"))
metadt$j.transplantation<-julian(metadt$date.transplantation,origin = as.Date("01/09/2020",format = "%d/%m/%Y"))
metadt$j.fleur1<-julian(metadt$date.fleur1,origin = as.Date("01/09/2020",format = "%d/%m/%Y"))
metadt$j.recolte<-julian(metadt$date.recolte,origin = as.Date("01/09/2020",format = "%d/%m/%Y"))
metadt$j.death<-julian(metadt$date.death,origin = as.Date("01/09/2020",format = "%d/%m/%Y"))
# make durations 
metadt$flowering_time<-metadt$j.fleur1-metadt$j.semis
metadt$life_span<-metadt$j.death-metadt$j.semis
hist(metadt$flowering_time)

# merge datasets
organs<-merge(organs,metadt,by="Genotype",all.x=T)
organs$stem[which(is.na(organs$stem))]<-"unknown"
organs$mean_rank[which(is.na(organs$mean_rank))]<-mean(organs$mean_rank,na.rm = T)
# it creates duplicates
organs[ which(organs$Genotype %in% organs$Genotype[which(duplicated(organs$Genotype))]),]
# based on the comments and values, best is to keep first row per duplicate 
organs<-organs[-which(duplicated(organs$Genotype)),]

# test some confounding effects
plot(organs$Seed~organs$flowering_time)
plot(organs$Seed~organs$mean_rank)
plot(organs$Petal_Area~organs$flowering_time)
plot(organs$Petal_Area~organs$flowering_time,col=as.factor(organs$note),pch=16)
plot(organs$Petal_Area~organs$mean_rank,col=as.factor(organs$stem),pch=16)
summary(lm(Petal_Area ~ stem + mean_rank, data = organs))

# Simplify organs
organs<-organs[,c(1,14,22,23,2:13,27:33)]

# Computing genotype's values without confounding effects (essentially the petal rank) 
organs$stem<-as.factor(organs$stem)
organs$table<-as.factor(organs$table)
organs_corrected<-organs[,c(1:4)]
#  
i<-5
for (i in 5:16) {
  donnees<-na.omit(organs[,c(1:4,i)])
  rownames(donnees)<-donnees$Genotype
  colnames(donnees)[5]<-"trait"
  mod<-lm(trait ~ stem * mean_rank + table, data = donnees)
  newtrait<-data.frame(Genotype=names(mod$residuals),trait=(residuals(mod)+mean(organs[,i],na.rm=T)))
  organs_corrected<-merge(organs_corrected,newtrait,by="Genotype",all.x = T)
  names(organs_corrected)[i]<-colnames(organs)[i]
}

summary(organs_corrected)
summary(organs)

organs_corrected[,17:23]<-organs[,17:23]

plot(organs_corrected$Petal_Area ~ organs_corrected$flowering_time)
cor.test(organs_corrected$Petal_Area,organs_corrected$flowering_time)

plot(organs_corrected$Leaf_Area ~ organs_corrected$flowering_time)
cor.test(organs_corrected$Leaf_Area,organs_corrected$flowering_time)

plot(organs_corrected$Seed ~ organs_corrected$flowering_time)
cor.test(organs_corrected$Seed,organs_corrected$flowering_time)

plot(organs_corrected$Seed ~ organs_corrected$Petal_Area)
cor.test(organs_corrected$Seed,organs_corrected$Petal_Area)

write.table(x = organs_corrected,file = "U_Shaped_Data_corrected_2023-05-05.csv",quote = F,row.names = F,col.names = T)
