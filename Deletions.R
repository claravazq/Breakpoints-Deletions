library(tidyverse)
#Clearing the Global enviroment:
remove(list = ls())

#Getting the BED files:
setwd("C:/Users/cvazque/Desktop/demos")
  ##CSV file##
df<- read.csv("BED files DS split complete.csv", sep=",")%>%
  rename(donor=V,
         chr=V2,
         chrS=V3,
         chrE=V4,
         read=V5,
         score=V6,
         strand=V7,
         insS=V8,
         insE=V9,
         RGB=V10,
         blocks=V11,
         blocklength=V12,
         blockS=V13)
  ##BED file##
df<- read.table("19006_MG_act_results_sorted.bed",header = F, sep="\t",stringsAsFactors=F)%>%
  rename(chr=V1,
         chrS=V2,
         chrE=V3,
         read=V4,
         score=V5,
         strand=V6,
         insS=V7,
         insE=V8,
         RGB=V9,
         blocks=V10,
         blocklength=V11,
         blockS=V12)
#Important information for plotting the antibody locus#
  #Switch coordinates##hg19:
Sm<-c(106323092,106326646)
Sg3<-c(106239278,106240979)
Sg1<-c(106210937,106213275)
Sg2<-c(106111949,106113759)
Sg4<-c(106093938,106094977)
Sa1<-c(106175604,106178608)
Sa2<-c(106055538,106057612)
  #Set of coordinates indicating regions##hg19
IgH<-c(106303092:106097612)
swSg1<- c(106210937:106213275)
swSm<-c(106320092:106329646)
swSg2<- c(106101949:106123759)
swSg3<- c(106229278:106260979)
swSg4<- c(106091938:106096977)
swSa1<- c(106173604:106180608)
swSa2<- c(106053538:106059612)

##Random number of reads##
  #First, filter your df by sample - example:
dfAT1<-filter(df, grepl("35.3", donor))
dfAT2<-filter(df, grepl("35.4", donor))
dfL<-filter(df, grepl("36", donor))
dfN<-filter(df, grepl("37", donor))
  #Second, make R to randomly take, in this case, 3000 rows corresponding to 3000 reads
  #and add a number to track the whole read to every row:
dfAT1<-dfAT1[(sample(1:nrow(dfAT1), 3000)), ]
dfAT1<-mutate(dfAT1,"numb"=c(1:(nrow(dfAT1))), "ymin"=5*(numb-1)+2)
dfAT2<-dfAT2[(sample(1:nrow(dfAT2), 3000)), ]
dfAT2<-mutate(dfAT2,"numb"=c(1:(nrow(dfAT2))), "ymin"=5*(numb-1)+2)
dfL<-dfL[(sample(1:nrow(dfL), 3000)), ]
dfL<-mutate(dfL, "numb"=c(1:(nrow(dfL))), "ymin"=5*(numb-1)+2)
dfN<-dfN[(sample(1:nrow(dfN), 3000)), ]
dfN<-mutate(dfN,"numb"=c(1:(nrow(dfN))), "ymin"=5*(numb-1)+2)
  #third, bind in a single df the single ones from before:
df500<- rbind(dfAT1, dfAT2, dfL, dfN)

##Depict the block characteristics into real coordinates:
tab<-separate(df500, blocklength, into=c("L1", "L2","L3", "L4", "L5", "L6"), sep="[\\,]")%>%
  separate(blockS, into=c("S1", "S2","S3", "S4", "S5", "S6"), sep = "[\\,]")%>%
  arrange(desc(S6))
tab<- mutate(tab, "S1"=chrS,
             "S2"=ifelse(str_detect(S2,"\\d"),as.numeric(S2)+as.numeric(chrS),"NO"),
             "S3"=ifelse(str_detect(S3,"\\d"),as.numeric(S3)+as.numeric(chrS),"NO"),
             "S4"=ifelse(str_detect(S4,"\\d"),as.numeric(S4)+as.numeric(chrS),"NO"),
             "S5"=ifelse(str_detect(S5,"\\d"),as.numeric(S5)+as.numeric(chrS),"NO"),
             "S6"=ifelse(str_detect(S6,"\\d"),as.numeric(S6)+as.numeric(chrS),"NO"))
tab<-mutate(tab, "L1"=as.numeric(L1)+as.numeric(S1),
            "L2"=ifelse(str_detect(L2,"\\d"),as.numeric(L2)+as.numeric(S2),"NO"),
            "L3"=ifelse(str_detect(L3,"\\d"),as.numeric(L3)+as.numeric(S3),"NO"),
            "L4"=ifelse(str_detect(L4,"\\d"),as.numeric(L4)+as.numeric(S4),"NO"),
            "L5"=ifelse(str_detect(L5,"\\d"),as.numeric(L5)+as.numeric(S5),"NO"),
            "L6"=ifelse(str_detect(L6,"\\d"),as.numeric(L6)+as.numeric(S6),"NO"))
tab<-mutate(tab, "B1"= paste(tab$S1, tab$L1, sep="-"),
            "B2" = paste (tab$S2, tab$L2, sep="-"),
            "B3" = paste (tab$S3, tab$L3, sep="-"),
            "B4" = paste (tab$S4, tab$L4, sep="-"),
            "B5" = paste (tab$S5, tab$L5, sep="-"),
            "B6" = paste (tab$S6, tab$L6, sep="-"))

  #separating tables:
tab<-separate(df, blocklength, into=c("L1", "L2","L3", "L4", "L5", "L6"), sep="[\\,]")%>%
  separate(blockS, into=c("S1", "S2","S3", "S4", "S5", "S6"), sep = "[\\,]")%>%
  arrange(desc(S6))
tab<- mutate(tab, "S1"=chrS,
             "S2"=ifelse(str_detect(S2,"\\d"),as.numeric(S2)+as.numeric(chrS),"NO"),
             "S3"=ifelse(str_detect(S3,"\\d"),as.numeric(S3)+as.numeric(chrS),"NO"),
             "S4"=ifelse(str_detect(S4,"\\d"),as.numeric(S4)+as.numeric(chrS),"NO"),
             "S5"=ifelse(str_detect(S5,"\\d"),as.numeric(S5)+as.numeric(chrS),"NO"),
             "S6"=ifelse(str_detect(S6,"\\d"),as.numeric(S6)+as.numeric(chrS),"NO"))
tab<-mutate(tab, "L1"=as.numeric(L1)+as.numeric(S1),
            "L2"=ifelse(str_detect(L2,"\\d"),as.numeric(L2)+as.numeric(S2),"NO"),
            "L3"=ifelse(str_detect(L3,"\\d"),as.numeric(L3)+as.numeric(S3),"NO"),
            "L4"=ifelse(str_detect(L4,"\\d"),as.numeric(L4)+as.numeric(S4),"NO"),
            "L5"=ifelse(str_detect(L5,"\\d"),as.numeric(L5)+as.numeric(S5),"NO"),
            "L6"=ifelse(str_detect(L6,"\\d"),as.numeric(L6)+as.numeric(S6),"NO"))%>%
            mutate("numb"=c(1:(nrow(tab))), "ymin"=4*(numb-1)+2)
tab<-mutate(tab, "B1"= paste(tab$S1, tab$L1, sep="-"),
            "B2" = paste (tab$S2, tab$L2, sep="-"),
            "B3" = paste (tab$S3, tab$L3, sep="-"),
            "B4" = paste (tab$S4, tab$L4, sep="-"),
            "B5" = paste (tab$S5, tab$L5, sep="-"),
            "B6" = paste (tab$S6, tab$L6, sep="-"))
  #add donor if necessary#
taba<-select(tab, chr, read, B1, B2, B3, B4, B5, B6, numb, ymin, donor)
demo<-gather(taba,"breaki", "coor", B1, B2, B3, B4, B5, B6)%>%arrange(desc(numb))%>%mutate("g"=ifelse(grepl("NO|NA", coor), "NO", "YES"))%>%
  filter(g=="YES")%>%select(numb, breaki, coor,chr, ymin, donor)
  #Separating the coords in two different columns:
tabla<-separate(demo, coor, into=c("coor1", "coor2"), sep="[\\-]")
  #Addition of two columns for ymin and ymax for plotting purposes:
tabla<-mutate(tabla, "ymax"=ymin+2)%>%arrange(donor,numb)
  #Calcuating the deletions:
bin<-mutate(tabla, "deletion"=ifelse(tabla$numb == lead(tabla$numb, 1), lead(tabla$coor1),"NO"), "size"=as.numeric(deletion)-as.numeric(coor2))
binn<-filter(bin, chr=="chr14", coor1 %in% IgH)

##Plotting the reads with a positive deletion size with facet_wrap:
filter(binn, size>0)%>%
  ggplot()+
  #geom_rect(aes(xmin=abs(as.numeric(Sm[1])-as.numeric(PSm)), xmax=abs(as.numeric(Sm[2])-as.numeric(PSm)), ymin=0, ymax=10), color="black", alpha=0.3)+
  #geom_rect(aes(xmin=abs(as.numeric(Sg3[1])-Sm[2]), xmax=abs(as.numeric(Sg3[2])-Sm[2]), ymin=0, ymax=10), color="black", alpha=0.3)+
  #geom_rect(aes(xmin=abs(as.numeric(Sg1[1])-PSm), xmax=abs(as.numeric(Sg1[2])-PSm), ymin=0, ymax=10), color="black", alpha=0.3)+
  #geom_rect(aes(xmin=abs(as.numeric(Sg2[1])-PSm), xmax=abs(as.numeric(Sg2[2])-PSm), ymin=0, ymax=10), color="black", alpha=0.3)+
  #geom_rect(aes(xmin=abs(as.numeric(Sg4[1])-PSm), xmax=abs(as.numeric(Sg4[2])-PSm), ymin=0, ymax=10), color="black", alpha=0.3)+
  #geom_rect(aes(xmin=abs(as.numeric(Sa1[1])-PSm), xmax=abs(as.numeric(Sa1[2])-PSm), ymin=0, ymax=10), color="black", alpha=0.3)+
  #geom_rect(aes(xmin=abs(as.numeric(Sa2[1])-Sm[2]), xmax=abs(as.numeric(Sa2[2])-Sm[2]), ymin=0, ymax=10), color="black", alpha=0.3)+
  geom_histogram(aes(x=size), fill="blue", alpha=0.5, bins=500)+
  geom_density(aes(x=size), alpha=0.5)+
  scale_y_log10()+
  #xlim(c(0,300000))+
  theme(axis.text.y = element_text(size=7), legend.position = "none")
  facet_grid(donor~., scales="free_y")
