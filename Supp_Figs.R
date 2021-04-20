library(ggplot2)
library(ape)
library(repr)
library("readxl")
library('gridExtra')
library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library(ggsci)
library(wesanderson)
library(RColorBrewer)
library(ggalt)


library(sf)
library(raster)
library(dplyr)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
library(ggnewscale)

data2<-read_excel('new_metadata.xlsx')
data2<-subset(data2,region=='Africa')

df_africa<-data2


df_africa$date<-as.Date(df_africa$date)

df_africa$days<-as.Date(cut(df_africa$date,breaks = "day",start.on.monday = FALSE))
df_africa$date2<-as.Date(cut(df_africa$date,breaks = "2 weeks",start.on.monday = FALSE))
df_africa$date3<-as.Date(cut(df_africa$date,breaks = "1 month",start.on.monday = FALSE))
df_africa$date_submitted<-as.Date(cut(df_africa$date_submitted,breaks = "day",start.on.monday = FALSE))

df_africa$submission_lag<-df_africa$date_submitted-df_africa$days

p_lag<-ggplot(df_africa,aes(x=reorder(country, submission_lag,median),y=submission_lag, fill=country))+
  #scale_fill_viridis_d(option='magma',alpha=0.7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
  geom_boxplot()+
  #stat_summary(fun=median, colour="blue", geom="point", 
  #             shape=18, size=3)+
  #geom_boxplot(aes(x=region,y=submission_lag))+
  scale_fill_viridis_d(alpha=0.7)+
  theme(legend.position = 'none')+
  ylab('Days from sample collection to sequence submission')+
  xlab('African countries with sequencing data')+
  scale_y_continuous(breaks=seq(0, 400, 50))
p_lag



pangolin_count<-as.data.frame(table(df_africa$pango_lineage))
pangolin_count<-pangolin_count[order(pangolin_count$Freq),]
pangolin_count_top20<-tail(pangolin_count,n=20)
pangolin_count_top20

pangolin_count_top20<-pangolin_count_top20 %>% 
  rename(
    pango_lineage = Var1,
  )

pangolin_count_top20$pangolin_africa<-pangolin_count_top20$pango_lineage



df_africa1 <- merge(df_africa,pangolin_count_top20,by="pango_lineage",all.x = TRUE)


df_africa1_DRC<-subset(df_africa1,country=='Democratic Republic of the Congo')
pangolin_count_DRC<-as.data.frame(table(df_africa1_DRC$pango_lineage))
pangolin_count_DRC<-pangolin_count_DRC[order(pangolin_count_DRC$Freq),]
pangolin_count_DRC_top5<-tail(pangolin_count_DRC,n=5)
pangolin_count_DRC_top5

df_africa1_mayotte<-subset(df_africa1,country=='Mayotte')
pangolin_count_mayotte<-as.data.frame(table(df_africa1_mayotte$pango_lineage))
pangolin_count_mayotte<-pangolin_count_mayotte[order(pangolin_count_mayotte$Freq),]
pangolin_count_mayotte_top5<-tail(pangolin_count_mayotte,n=5)
pangolin_count_mayotte_top5

df_africa1_Zim<-subset(df_africa1,country=='Zimbabwe')
pangolin_count_Zim<-as.data.frame(table(df_africa1_Zim$pango_lineage))
pangolin_count_Zim<-pangolin_count_Zim[order(pangolin_count_Zim$Freq),]
pangolin_count_Zim_top5<-tail(pangolin_count_Zim,n=5)
pangolin_count_Zim_top5

#df_africa2<-subset(df_africa1,is.na(Freq))

custom2<-c("tan4",'peachpuff3','antiquewhite2','darkseagreen2','darkseagreen4','darkolivegreen',
           #custom2<-c("darkolivegreen",'darkseagreen3','darkseagreen2','antiquewhite2','peachpuff3','tan4',
           'red3','deeppink3','hotpink2','plum2','grey',
           'thistle3','mediumpurple','mediumorchid3','purple4',
           'blue3','dodgerblue1','cadetblue3','skyblue1','slategray1', 'white')



custom1 <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
             "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352", "black","slategrey",
             "burlywood4","thistle4","lightsalmon2","coral2","lightpink3","indianred3","orange3","yellow3","yellowgreen",
             "khaki4","springgreen4","aquamarine3","lightblue2","lightblue4")


#pangolin_Moz$pangolin_africa<-pangolin_count_top20$pango_lineage

df_africa1 <- merge(df_africa,pangolin_count_top20,by="pango_lineage",all.x = TRUE)

pango_lineage<-c('C.1','C.1.1','B.1.351')
Moz_lin<-c('C.1','C.1.1','B.1.351')
Moz_lin_df<-data.frame(pango_lineage,Moz_lin)
df_africa1<-merge(df_africa1,Moz_lin_df,by='pango_lineage',all.x=TRUE)

pango_lineage<-c('B.1.525','B.1.1.7')
Nigeria_lin<-c('B.1.525','B.1.1.7')
Nigeria_lin_df<-data.frame(pango_lineage,Nigeria_lin)
df_africa1<-merge(df_africa1,Nigeria_lin_df,by='pango_lineage',all.x=TRUE)


pango_lineage<-c('B.1.1.54','B.1.1.56','C.1','B.1.351')
SA_lin<-c('B.1.1.54','B.1.1.56','C.1','B.1.351')
SA_lin_df<-data.frame(pango_lineage,SA_lin)
df_africa1<-merge(df_africa1,SA_lin_df,by='pango_lineage',all.x=TRUE)


pango_lineage<-c('A.23','A.23.1')
UG_lin<-c('A.23','A.23.1')
UG_lin_df<-data.frame(pango_lineage,UG_lin)
df_africa1<-merge(df_africa1,UG_lin_df,by='pango_lineage',all.x=TRUE)


pango_lineage<-c('B.1.1.7','B.1.525')
Ghana_lin<-c('B.1.1.7','B.1.525')
Ghana_lin_df<-data.frame(pango_lineage,Ghana_lin)
df_africa1<-merge(df_africa1,Ghana_lin_df,by='pango_lineage',all.x=TRUE)


pango_lineage<-c('B.1.351','B.1.160')
Mayotte_lin<-c('B.1.351','B.1.160')
Mayotte_lin_df<-data.frame(pango_lineage,Mayotte_lin)
df_africa1<-merge(df_africa1,Mayotte_lin_df,by='pango_lineage',all.x=TRUE)



pango_lineage<-c('B.1.351','C.2','B.1.1.29','B.1.1.111')
Zim_lin<-c('B.1.351','C.2','B.1.1.29','B.1.1.111')
Zim_lin_df<-data.frame(pango_lineage,Zim_lin)
df_africa1<-merge(df_africa1,Zim_lin_df,by='pango_lineage',all.x=TRUE)


P_lin_all<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(df_africa1, !is.na(strain)), mapping = aes(x = date2, fill=pangolin_africa))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month", limits = as.Date(c('2020-02-01','2021-03-01')))+
  theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=15))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=15))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = custom2, name='Lineage')+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=14))+
  theme(legend.position = "top")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')
#ggtitle('Africa - Top 20 circulating lineages')

P_lin_all


P_Mozambique_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(df_africa1, country=='Mozambique'), mapping = aes(x = date2, fill=Moz_lin))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('mediumorchid3','cadetblue4','cadetblue2'), name='Lineage', labels=c('B.1.351','C.1','C.1.1','Others'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Mozambique')

P_Mozambique_lin


P_Nigeria_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(df_africa1, country=='Nigeria'), mapping = aes(x = date2, fill=Nigeria_lin))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('grey','skyblue3'), name='Lineage', labels=c('B.1.1.7','B.1.525'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Nigeria')

P_Nigeria_lin


P_SA_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(df_africa1, country=='South Africa'), mapping = aes(x = date2, fill=SA_lin))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('darkseagreen3','darkorange3','mediumorchid3','cadetblue4'), name='Lineage', labels=c('B.1.1.54','B.1.1.56','B.1.351','C.1','Others'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('South Africa')

P_SA_lin


P_Mayotte_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(df_africa1, country=='Mayotte'), mapping = aes(x = date2, fill=Mayotte_lin))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black',size=0.1)+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%d-%b\n%Y",date_breaks = "2 weeks")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('darkseagreen4','mediumorchid3','cadetblue4'), name='Lineage', label=c('B.1.160','B.1.351','Others'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Mayotte')

P_Mayotte_lin


P_Botswana_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(subset(df_africa1, country=='Botswana'), !is.na(pangolin_africa)), mapping = aes(x = date3, fill=pangolin_africa=='B.1.351'))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=20,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 months")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('white','mediumorchid3'), name='Lineage', label=c('Others','B.1.351'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Botswana')

P_Botswana_lin


P_Zambia_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date2)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(subset(df_africa1, country=='Zambia'), !is.na(pangolin_africa)), mapping = aes(x = date2, fill=pangolin_africa=='B.1.351'))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 months")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('white','mediumorchid3'), name='Lineage', label=c('Others','B.1.351'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Zambia')

P_Zambia_lin



P_Zimbabwe_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date2)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(subset(df_africa1, country=='Zimbabwe'), !is.na(pangolin_africa)), mapping = aes(x = date2, fill=Zim_lin))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 months")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('darkorange2','mediumorchid3','cadetblue3'), name='Lineage', label=c('B.1.1.29','B.1.351','Others'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Zimbabwe')

P_Zimbabwe_lin


P_UG_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(df_africa1, country=='Uganda'), mapping = aes(x = date2, fill=UG_lin))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('tan4','peachpuff3'), name='Lineage', labels=c('A.23','A.23.1','Others'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Uganda')

P_UG_lin


P_Ghana_lin<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  #group_by(D614G_variant, Date_received_by_Sender) %>% #group
  #summarise(prop = sum(D614G_variant=="G")/n()) %>% #calculate proportion 
  ggplot(data=subset(df_africa1, country=='Ghana'), mapping = aes(x = date2, fill=Ghana_lin))+
  #ggplot(data=subset(data2, !is.na(strain)), mapping = aes(x = date,fill=District))+
  geom_bar(position='fill',width=10,color='black')+
  #geom_bar(width=5)+
  
  #geom_bar(width=5,color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month")+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=10, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  #scale_fill_manual(values=custom1,labels = c("A","B.1","B.1.1.1","B.1.1.119","B.1.1.56","B.1.1.54","B.1.351","B.1.416","C.1","H.1", "Others"), name='Lineage') +
  #scale_fill_brewer(palette="Spectral")+
  #scale_fill_manual(values = getPalette(24))+
  #scale_fill_manual(values = getPalette(24))+
  scale_fill_manual(values = c('grey','skyblue3'), name='Lineage', labels=c('B.1.1.7','B.1.525','Others'))+
  #scale_fill_manual(values=inauguration("inauguration_2021"))+
  #scale_fill_manual(values=c('white','gold2'),labels = c("Other lineages", "501Y.V2/B.1.351"), name='Lineage')+
  #scale_fill_manual(labels = c("B.1.1.54", "B.1.1.56","C.1","Others"))+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  #theme(axis.title.y = element_blank())+
  #theme(axis.title.x = element_blank())+
  #theme(legend.key.size = unit(0.2, "cm"))+
  xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')+
  ggtitle('Ghana')

P_Ghana_lin

plot_grid(P_SA_lin,P_Mozambique_lin,P_Zimbabwe_lin,P_Zambia_lin,P_Botswana_lin,P_Mayotte_lin,P_Nigeria_lin,P_Ghana_lin,P_UG_lin)

ggsave('Supp_lineage_progression_15April.pdf',width = 30, height = 20, units = "cm")



