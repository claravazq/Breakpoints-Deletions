#Libraries
library(tidyverse)
#Find your BED files
setwd("C:/Users/cvazque/Desktop/demos/")

#1.- BREAKPOINTS# Choose one of the ways to import your BED file - A or B
  #>>A#If your BED file has one column with this information: chrX:zzzzzzzz-zzzzzzz
tab<-read.csv("QiangP.csv", sep=",")%>%select(sample, chr)

tab<-separate(tab, chr, into=c("chro", "begin"), sep="[\\:]")%>%
  separate(begin, into=c("end", "begin"), sep = "[\\-]")%>%
  separate(begin, into=c("begin", "X"), sep = "[\\s]")%>%
  select(sample, begin, end)
  #Choose the different samples found in your file and separate them to gather them#Example below: 
tabB<-filter(tab, sample=="BRCA1")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("BRCA1"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabL<-filter(tab, sample=="Lig4")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("Lig4"))%>%group_by(sample, key, value)%>%summarise(n=n())
tadA<-filter(tab, sample=="AID")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("AID"))%>%group_by(sample, key, value)%>%summarise(n=n())
tab6<-filter(tab, sample=="ATM106")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("ATM106"))%>%group_by(sample, key, value)%>%summarise(n=n())
tab7<-filter(tab, sample=="ATM107")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("ATM107"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabn<-filter(tab, sample=="NIPBL")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("NIPBL"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabc1<-filter(tab, sample=="Control1")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("Control1"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabc2<-filter(tab, sample=="Control2")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("Control2"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabc3<-filter(tab, sample=="Control3")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("Control3"))%>%group_by(sample, key, value)%>%summarise(n=n())
X<-rbind(tabB, tadA, tab6, tab7, tabn, tabc1, tabc2, tabc3, tabL)
X.<-mutate(X, "type"= ifelse(grepl("Control", sample), "Control", "Case"))
  
  #>>B#If your BED files uses blocks as information:
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
          ##only chr14 in IgH locus##To catch only the antibody locus:
tabl<-separate(df, blocklength, into=c("L1", "L2","L3", "L4", "L5", "L6"), sep="[\\,]")%>%
  separate(blockS, into=c("S1", "S2","S3", "S4", "S5", "S6"), sep = "[\\,]")%>%
  arrange(desc(S6))%>%filter(chr=="chr14")
tabl<- mutate(tabl, "S1"=chrS,
              "S2"=ifelse(str_detect(S2,"\\d"),as.numeric(S2)+as.numeric(chrS),"NO"),
              "S3"=ifelse(str_detect(S3,"\\d"),as.numeric(S3)+as.numeric(chrS),"NO"),
              "S4"=ifelse(str_detect(S4,"\\d"),as.numeric(S4)+as.numeric(chrS),"NO"),
              "S5"=ifelse(str_detect(S5,"\\d"),as.numeric(S5)+as.numeric(chrS),"NO"),
              "S6"=ifelse(str_detect(S6,"\\d"),as.numeric(S6)+as.numeric(chrS),"NO"))
tabl<-mutate(tabl, "L1"=as.numeric(L1)+as.numeric(S1),
             "L2"=ifelse(str_detect(L2,"\\d"),as.numeric(L2)+as.numeric(S2),"NO"),
             "L3"=ifelse(str_detect(L3,"\\d"),as.numeric(L3)+as.numeric(S3),"NO"),
             "L4"=ifelse(str_detect(L4,"\\d"),as.numeric(L4)+as.numeric(S4),"NO"),
             "L5"=ifelse(str_detect(L5,"\\d"),as.numeric(L5)+as.numeric(S5),"NO"),
             "L6"=ifelse(str_detect(L6,"\\d"),as.numeric(L6)+as.numeric(S6),"NO"))%>%
  mutate("numb"=c(1:(nrow(tabl))), "ymin"=4*(numb-1)+2)
tabl<-mutate(tabl, "B1"= paste(tabl$S1, tabl$L1, sep="-"),
             "B2" = paste (tabl$S2, tabl$L2, sep="-"),
             "B3" = paste (tabl$S3, tabl$L3, sep="-"),
             "B4" = paste (tabl$S4, tabl$L4, sep="-"),
             "B5" = paste (tabl$S5, tabl$L5, sep="-"),
             "B6" = paste (tabl$S6, tabl$L6, sep="-"))

demol<-select(tabl, chr, read, B1, B2, B3, B4, B5, B6, numb, ymin, donor)%>% gather("breaki", "coor", B1, B2, B3, B4, B5, B6)%>%arrange(desc(numb))%>%mutate("g"=ifelse(grepl("NO|NA", coor), "NO", "YES"))%>%
  filter(g=="YES")%>%select(numb, breaki, coor,chr, ymin, donor)
tablal<-separate(demol, coor, into=c("coor1", "coor2"), sep="[\\-]")%>%mutate("ymax"=ymin+2)%>%arrange(numb)

##Plots##
 #Information for feautures of the plot# 
  #Switch coordinates##hg38
Sm<-c(105856987,105860436)
Sg3<-c(105772941,105774642)
Sg1<-c(105744600,105746938)
Sg2<-c(105645612,105647422)
Sg4<-c(105627601,105628640)
Sa1<-c(105709267,105712271)
Sa2<-c(105589201,105591275)
Se<- c(105602375,105603285)
  #Switch coordinates entre puntos##hg38
SSm<-c(105856987:105860436)
SSg3<-c(105772941:105774642)
SSg1<-c(105744600:105746938)
SSg2<-c(105645612:105647422)
SSg4<-c(105627601:105628640)
SSa1<-c(105709267:105712271)
SSa2<-c(105589201:105591275)
SSe<- c(105602375:105603285)
  #Plot of all the IgH locus#
ggplot(X.)+
  geom_rect(aes(xmin=Sm[1], xmax=Sm[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg3[1], xmax=Sg3[2], ymin=0, ymax=74827),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg1[1], xmax=Sg1[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg2[1], xmax=Sg2[2], ymin=0, ymax=74827),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg4[1], xmax=Sg4[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sa1[1], xmax=Sa1[2], ymin=0, ymax=74827),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sa2[1], xmax=Sa2[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "AID"="darkgreen", 
                              "ATM106"="black",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=8, vjust=0.5))
  
#If you want to zoom in one switch region#
  #Filter first the coordinates included in one of the switches#
filter(X., value %in% SSm)%>%
  ggplot()+
  geom_rect(aes(xmin=Sm[1], xmax=Sm[2], ymin=0, ymax=74827),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=105860975, xmax=105861001, ymin=0, ymax=74827),fill="yellow", alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "Lig4"="black",
                              "AID"="darkgreen", 
                              "ATM106"="turquoise4",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints in Sm")+
  geom_hline(yintercept = -140, color="grey75")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample_f~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=7.5,face="bold", vjust=0.5))
