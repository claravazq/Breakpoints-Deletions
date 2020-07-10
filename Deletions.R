library(tidyverse)
#Clearing the Global enviroment:
remove(list = ls())
#Switch coordinates##hg19:
Sm<-c(106323092,106326646)
Sg3<-c(106239278,106240979)
Sg1<-c(106210937,106213275)
Sg2<-c(106111949,106113759)
Sg4<-c(106093938,106094977)
Sa1<-c(106175604,106178608)
Sa2<-c(106055538,106057612)
IgH<-c(106303092:106097612)
swSg1<- c(106210937:106213275)
swSm<-c(106320092:106329646)
swSg2<- c(106101949:106123759)
swSg3<- c(106229278:106260979)
swSg4<- c(106091938:106096977)
swSa1<- c(106173604:106180608)
swSa2<- c(106053538:106059612)
  #Primer binding:
PSm<-106327185

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
##Random number of reads #Limiting factor 19037-52518 reads#
dfAT1<-filter(df, grepl("35.3", donor))
dfAT2<-filter(df, grepl("35.4", donor))
dfL<-filter(df, grepl("36", donor))
dfN<-filter(df, grepl("37", donor))



dfAT1<-dfAT1[(sample(1:nrow(dfAT1), 3000)), ]
dfAT1<-mutate(dfAT1,"numb"=c(1:(nrow(dfAT1))), "ymin"=5*(numb-1)+2)
dfAT2<-dfAT2[(sample(1:nrow(dfAT2), 3000)), ]
dfAT2<-mutate(dfAT2,"numb"=c(1:(nrow(dfAT2))), "ymin"=5*(numb-1)+2)
dfL<-dfL[(sample(1:nrow(dfL), 3000)), ]
dfL<-mutate(dfL, "numb"=c(1:(nrow(dfL))), "ymin"=5*(numb-1)+2)
dfN<-dfN[(sample(1:nrow(dfN), 3000)), ]
dfN<-mutate(dfN,"numb"=c(1:(nrow(dfN))), "ymin"=5*(numb-1)+2)

df500<- rbind(dfAT1, dfAT2, dfL, dfN)

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

##separating tables:
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
#Addition of two columns for ymin and ymax:
tabla<-mutate(tabla, "ymax"=ymin+2)%>%arrange(donor,numb)
#Addition of deletions:
bin<-mutate(tabla, "deletion"=ifelse(tabla$numb == lead(tabla$numb, 1), lead(tabla$coor1),"NO"), "size"=as.numeric(deletion)-as.numeric(coor2))
binn<-filter(bin, chr=="chr14", coor1 %in% IgH)
#Plotting the size of the deletion:
filter(binn, size>160000)%>%
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

##Plotting only 3000 reads intra switch deletions:
  filter(bin, deletion != "NO", size > 0) %>%
    filter(coor1 %in% swSg4)%>%
    filter(deletion %in% swSg4)%>%
    ggplot()+
    #geom_rect(aes(xmin=as.numeric(coor2), xmax=as.numeric(deletion), ymin=ymin, ymax=ymax, color=donor),fill="black", alpha=0.5)+
    #geom_rect(aes(xmin=Sg3[1], xmax=Sg3[2], ymin=0, ymax=15000),fill="grey95", alpha=0.3)+
    #geom_rect(aes(xmin=Sg1[1], xmax=Sg1[2], ymin=0, ymax=15000), fill="grey95",alpha=0.3)+
    #geom_rect(aes(xmin=Sg2[1], xmax=Sg2[2], ymin=0, ymax=15000),fill="grey95", alpha=0.3)+
    geom_rect(aes(xmin=Sg4[1], xmax=Sg4[2], ymin=0, ymax=15000), fill="grey95",alpha=0.3)+
    #geom_rect(aes(xmin=Sa1[1], xmax=Sa1[2], ymin=0, ymax=15000),fill="grey95", alpha=0.3)+
    #geom_rect(aes(xmin=Sa2[1], xmax=Sa2[2], ymin=0, ymax=15000),fill="grey95",alpha=0.3)+
    #geom_rect(aes(xmin=Sm[1], xmax=Sm[2], ymin=0, ymax=15000), fill="grey95",alpha=0.3)+
    #geom_rect(aes(xmin=106097612, xmax=106303092, ymin=0, ymax=50000), fill="red",alpha=0.3)+
    geom_rect(aes(xmin=as.numeric(coor2), xmax=as.numeric(deletion), ymin=ymin, ymax=ymax, color=donor),fill="black", alpha=0.5)+
    #geom_rect(aes(xmin=as.numeric(coor1), xmax=as.numeric(coor2), ymin=ymin, ymax=ymax), color="grey65")+
    scale_color_manual(values=c("D1935.3"="coral", 
                                "D1936"="black",
                                "D1935.4"="turquoise4",
                                "D1937"="blue3"))+
    ggtitle("Deletions from 3000 samples Sg4")+
    geom_hline(yintercept = 0, color="black")+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
          axis.title.y=element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
          axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
          strip.text.y = element_text(size=7.5,face="bold", vjust=0.5), legend.position="none")+
    facet_grid(donor~.)
  
#Most common deletion #ANOVA#
mode<-filter(bin, coor2 != "NO", size > 0) %>%
    filter(coor1 %in% swSm)%>%
    filter(deletion %in% swSm)%>%
  select(donor, coor2)

group_by(mode, donor)%>%
  summarise(n=mean(deletion))

a<-filter(mode,grepl('37', donor))%>%rename('DS1937'='coor2')
a<-as.vector(a$DS1937)
b<-filter(mode,grepl('36', donor))%>%rename('DS1936'='coor2')
b<-as.vector(b$DS1936)
c<-filter(mode,grepl('35.3', donor))%>%rename('DS1935.3'='coor2')
c<-as.vector(c$DS1935.3)
d<-filter(mode,grepl('35.4', donor))%>%rename('DS1935.4'='coor2')
d<-as.vector(d$DS1935.4)

mode4<-as.data.frame(cbind(a,b,d,c))
mode4<-rename(mode4,D1937=a, D1936=b,D1935.4=d,D1935.3=c)
summary(mode4)
ao<-aov(coor2~donor, data=mode)  
plot(TukeyHSD(ao))


#Looking for the peaks in the histogram:
      #Density calculation#
ar<-filter(binn, grepl("37", donor))  
dn<-as.numeric(paste(ar[,"size"], sep=","))
den<-density(na.omit(dn))  
  ##Two ways to get highest Peak##
den$x[den$y==max(den$y)]  #Gives you all highest Peaks
den$x[which.max(den$y)] #Gives you the first highest Peak
  ##3 ways to get all Peaks
den$x[c(F, diff(diff(den$y)>=0)<0)] #This detects also a plateau
den$x[which(diff(sign(diff(den$y)))<0)+1]

ds1938<-den$x[which(diff(sign(diff(den$y)))<0)+1]
peaks<-rbind(ds1935.3, ds1935.4,ds1936,ds1938)
peaks<-t(peaks)%>%as.data.frame()

#as.data.frame(peaks)%>%
ggplot(peaks)+
  geom_boxplot(aes(x=ds1938), fill="red", alpha=0.3)
  #geom_histogram(aes(x=ds1936), fill="green", alpha=0.3)+
  #geom_histogram(aes(x=ds1935.3), fill="blue", alpha=0.3)+
  #geom_histogram(aes(x=ds1935.4), fill="coral", alpha=0.3)

#Distances between every peak:
peaks<-mutate(peaks, d353=lead(ds1935.3), s354=lead(ds1935.4), s38=lead(ds1938), s36=lead(ds1936))
peaks<-mutate(peaks, s1935.3=d353-ds1935.3, s1935.4=s354-ds1935.4, s1938=s38-ds1938, s1936=s36-ds1936)
peak<-select(peaks, -d353, -s354, -s38,-s36)
