library(ggplot2)
library("readxl")
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library(ggalluvial)
library("lubridate")

library(ggtree)
library(tidytree)
#library(ape)
library(treeio)


tree<-read.newick('data/timetree.nwk')
metadata_df <- read_excel("data/new_metadata.xlsx")


custom3<-c('antiquewhite2',"tan4",'peachpuff3','palegreen3','dodgerblue1','hotpink2','mediumorchid3','blue3','skyblue1','purple4',
           'mediumpurple','darkseagreen4','cadetblue3','thistle3','deeppink3',
           #custom2<-c("darkolivegreen",'darkseagreen3','darkseagreen2','antiquewhite2','peachpuff3','tan4',
           'red3','grey30','darkolivegreen','slategray1','plum2',
           'white','white')

p<-ggtree(tree, mrsd="2021-03-15",as.Date=TRUE, color='grey80',size=0.2) + theme_tree2()+
  scale_x_date(date_labels = "%b",date_breaks = "2 month")+
  expand_limits(y = 21000)+
  theme(axis.text.x = element_text(size=10,angle=90))

p

panelA <- p %<+% metadata_df + 
  geom_tippoint(aes(
    subset=(region=='Africa')),fill='white',size=2,align=F,stroke=0.2,color='grey60',shape=21)+
  geom_tippoint(aes(
    subset=(pango_lineage=='A.23' & region=='Africa')),fill='tan4',size=3, stroke=0.2,align=F, color='salmon4',shape=21)+
  geom_tippoint(aes(
    subset=(pango_lineage=='A.23.1' & region=='Africa')),fill='peachpuff3',size=3, stroke=0.2,align=F, color='peachpuff4',shape=21)+
  geom_tippoint(aes(
    subset=(pango_lineage=='B.1.1.7' & region=='Africa')),fill='grey',size=3,align=F, stroke=0.2,color='grey30',shape=21)+
  geom_tippoint(aes(
    subset=(pango_lineage=='B.1.351' & region=='Africa')),fill='mediumorchid3',size=3,align=F, stroke=0.2,color='purple4',shape=21)+
  geom_tippoint(aes(
    subset=(pango_lineage=='B.1.525' & region=='Africa')),fill='cadetblue3',size=3, stroke=0.2,align=F, color='cadetblue4',shape=21)+
  geom_tippoint(aes(
    subset=(pango_lineage=='C.1' & region=='Africa')),fill='skyblue4',size=3, stroke=0.2,align=F, color='darkslategrey',shape=21)+
  geom_tippoint(aes(
    subset=(pango_lineage=='C.1.1' & region=='Africa')),fill='skyblue1',size=3, stroke=0.2,align=F, color='skyblue4',shape=21)+
  scale_fill_manual(values=custom3, name='African genomes of 20 most\ncommon lineages')+
  theme(axis.text.x = element_text(vjust=0.5, hjust=1))

panelA

ggsave('Africa_tree_round3_test.pdf', width = 50, height = 100, units = "cm",limitsize = FALSE)



importexport<-read_excel('data/CountryAnnotatedImportExport.xlsx')

importexport$EventTime<-as.numeric(importexport$EventTime)


importexport$date <- as.Date(format(date_decimal(importexport$EventTime), "%Y-%m-%d"))

importexport$date

importexport_africa<-subset(importexport, link!='Other regions-Other regions')
importexport_africa<-subset(importexport_africa, link!='Africa-Other regions')

importexport_africa$days<-as.Date(cut(importexport_africa$date,
                        breaks = "day",
                        start.on.monday = FALSE))

importexport_africa$date<-as.Date(cut(importexport_africa$date,
                        breaks = "2 week",
                        start.on.monday = FALSE))

panelB<-importexport_africa %>%
  mutate(date=as.POSIXct(date)) %>%
  ggplot(data=importexport_africa, mapping = aes(x = date,fill=link))+
  geom_bar(width=10, color='black', position='fill')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16, angle=90))+
  theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=15))+
  scale_fill_manual(values=c('black','white'), labels=c('Other African countries','Rest of the world'))+
  scale_x_date(date_labels = "%b",date_breaks = "month", limits=as.Date(c('2020-01-01','2021-01-15')))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=14))+
  #theme(legend.position = c(0.2, 0.8))+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  theme(legend.position = "bottom")+
  labs(fill="Origin")+
  xlab('Date')+
  ylab('Proportion')+
  ggtitle("Sources of viral introductions into African countries")
panelB


data4<-read_excel('data/ImportExportRegionAnnotated.xlsx')

data4$EventTime<-as.numeric(data4$EventTime)


data4$date <- as.Date(format(date_decimal(data4$EventTime), "%Y-%m-%d"))

data4$date

data5<-data4
data5$days<-as.Date(cut(data5$date,
                        breaks = "day",
                        start.on.monday = FALSE))

data5$date<-as.Date(cut(data5$date,
                        breaks = "week",
                        start.on.monday = FALSE))

intro_df<-subset(data5, Destination=='Africa')

custom_pal2<-c('hotpink3','darkseagreen3','mediumpurple','lightskyblue1','royalblue3')


panelC<-intro_df %>%
  mutate(date=as.POSIXct(date)) %>%
  ggplot(data=subset(intro_df, !is.na(Origin)), mapping = aes(x = date,fill=Origin))+
  #geom_bar(position='fill',width=5)+
  #geom_bar(width=0.7,color='black',fill='darkorange')+
  geom_bar(width=5, color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16, angle=90))+
  theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=15))+
  scale_fill_manual(values=custom_pal2)+
  scale_x_date(date_labels = "%b",date_breaks = "month")+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=14))+
  #theme(legend.position = c(0.2, 0.8))+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  theme(legend.position = "bottom")+
  labs(fill="Origin")+
  xlab('Date')+
  ylab('Count')+
  ggtitle("International introductions into Africa")

panelC


panelD<-ggplot(subset(importexport,link!='Other regions-Other regions'),
               aes(axis1 = origin1, axis2 = destination1)) +
  geom_alluvium(width = 1/12, aes(fill=link),decreasing = TRUE, alpha=0.9) +
  geom_stratum(width = 1/12, fill='white', color='grey', decreasing = TRUE) +
  #geom_label_repel(size=3,stat = "stratum", min.y = 20,decreasing = TRUE,aes(label = after_stat(stratum)))+ ## Note: Can be uncommented to show country labels
  scale_fill_manual(values=c('deeppink3','dodgerblue3','grey90'), name='Connection', labels=c('Africa-Africa','Africa-World','World-Africa'))+
  scale_x_discrete(limits = c("Origin", "Destination"), expand = c(.05, .05)) +
  theme_minimal()+
  theme(legend.position='bottom')
panelD
# Note: This panel is labelled and further edited in Illustrator or Powerpoint





