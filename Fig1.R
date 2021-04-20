library(ggplot2)
library("readxl")
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)


library(sf)
library(raster)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
library(ggnewscale)

data2<-read_excel('data/new_metadata.xlsx')
data2<-subset(data2,region=='Africa')

df_africa<-data2


df_africa$date<-as.Date(df_africa$date)

df_africa$days<-as.Date(cut(df_africa$date,breaks = "day",start.on.monday = FALSE))
df_africa$date2<-as.Date(cut(df_africa$date,breaks = "2 weeks",start.on.monday = FALSE))
df_africa$date3<-as.Date(cut(df_africa$date,breaks = "1 month",start.on.monday = FALSE))

df_count <- df_africa %>% count(country)
names(df_count)[names(df_count) == "country"] <- "country"
names(df_count)[names(df_count) == "n"] <- "Count"
df_count

df_africa = df_africa %>% 
  left_join(df_count, by = c("country" = "country"))


africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(worldbank_df, by = "iso_a2") %>% 
  dplyr::select(name, subregion, gdpPercap, HDI, pop_growth) %>% 
  st_transform("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25")


df_count[df_count == "Republic of the Congo"] <- "Republic of Congo"
df_count[df_count == "Eswatini"] <- "Swaziland"
df_count[df_count == "Gambia"] <- "The Gambia"


africa = africa %>% 
  left_join(df_count, by = c("name" = "country"))



panelA<-ggplot(africa) +
  geom_sf(aes(geometry = geom, fill = Count))+
  #scale_x_continuous(trans = log2_trans())+
  scale_fill_distiller(palette = "BrBG", direction = 1,trans = "log",na.value = "white",breaks = c(3500,1000), labels = c(3500,1000)) + 
  # scale_fill_gradient(name = "count", trans = "log",
  #                     breaks = my_breaks, labels = my_breaks)
  geom_sf_label(aes(label = Count),alpha=0.5, size=3,label.padding = unit(0.1, "lines"))+
  #geom_sf_text(aes(label = name),alpha=0.5, size=4)+
  theme(legend.position = 'none')+
  theme_void()
panelA 
#### Important Note: Map is further edited in Illustrator to correct the border for Morocco as the shapefile is incorrect



epi_data<-read_excel('data/owid-covid-data_29March.xlsx')
epi_data<-subset(epi_data,continent=='Africa')
epi_data<-subset(epi_data, date==as.Date("2021/03/16"))

epi_data$location[epi_data$location == "Eswatini"] <- "Swaziland"
epi_data$location[epi_data$location == "Gambia"] <- "The Gambia"
epi_data$location[epi_data$location == "Congo"] <- "Republic of Congo"
epi_data$location[epi_data$location == "Democratic Republic of Congo"] <- "Democratic Republic of the Congo"
names(epi_data)[names(epi_data) == "location"] <- "country"

df_count = df_count %>% 
  left_join(epi_data, by = c("country" = "country"))

df_count=transform(df_count, total_cases = as.numeric(total_cases))
df_count=transform(df_count, Count = as.numeric(Count))

corr <- cor.test(x=df_count$total_cases, y=df_count$Count, method = 'pearson')
corr

df_count$country[df_count$country == "The Gambia"] <- "Gambia"


panelB<-ggplot(df_count, aes(x=total_cases/100000, y=Count)) + 
  theme_classic()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='deeppink3')+
  
  geom_point(color='#2980B9', aes(x=total_cases/100000, y=Count, size=Count), alpha=0.5) + 
  geom_text_repel(aes(label = ifelse(Count>400,as.character(country),'')),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50',
                  max.overlaps=15)+
  theme(axis.title.x = element_text(color="black", size=18, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=15))+
  theme(axis.title.y = element_text(color="black", size=18, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=15))+
  theme(legend.position=c(0.3,0.8), 
        #legend.justification='left',
        legend.direction='vertical')+
  scale_size(name   = "Number of Sequences",
             breaks = c(10,100,1000,4000))+
  annotate(geom="text", x=1000000/100000, y=1900, label="r = 0.91",
           color="hotpink3",size=5, fontface='bold')+
  ylab("Sequences")+
  xlab("Total Cases (x100,000)")


panelB


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

custom2<-c("tan4",'peachpuff3','antiquewhite2','darkseagreen2','darkseagreen4','darkolivegreen',
           #custom2<-c("darkolivegreen",'darkseagreen3','darkseagreen2','antiquewhite2','peachpuff3','tan4',
           'red3','deeppink3','hotpink2','plum2','grey',
           'thistle3','mediumpurple','mediumorchid3','purple4',
           'blue3','dodgerblue1','cadetblue3','skyblue1','slategray1', 'white')



custom1 <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
             "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352", "black","slategrey",
             "burlywood4","thistle4","lightsalmon2","coral2","lightpink3","indianred3","orange3","yellow3","yellowgreen",
             "khaki4","springgreen4","aquamarine3","lightblue2","lightblue4")


df_africa1 <- merge(df_africa,pangolin_count_top20,by="pango_lineage",all.x = TRUE)


panelC<-df_africa1 %>%
  mutate(date=as.POSIXct(date3)) %>%
  ggplot(data=subset(df_africa1, !is.na(strain)), mapping = aes(x = date2, fill=pangolin_africa))+
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
  scale_fill_manual(values = custom2, name='Lineage')+
   theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=14))+
  theme(legend.position = "top")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
   xlab('Date')+
  #ylab('Genome Count')
  ylab('Proportion of Genomes')
#ggtitle('Africa - Top 20 circulating lineages')

panelC


panelD<-ggplot(data=df_africa)  + theme_classic()+
  geom_segment(aes(x=min(df_africa1$date), y=reorder(country,Count), xend=max(df_africa1$date), yend=country, group=country), colour="grey80", size=5) +
  geom_point(aes(x=days, fill='Other Lineages',y=reorder(country,Count)),position = position_jitter(width=0.2, height=0.2), shape=21,stroke=0.05, col='grey70', size=4)+
  geom_point(data=subset(df_africa, pango_lineage=='B.1.351'),aes(fill='B.1.351',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=4, alpha=1)+
  geom_point(data=subset(df_africa, pango_lineage=='B.1.1.7'),aes(fill='B.1.1.7',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=4, alpha=1)+
  geom_point(data=subset(df_africa, pango_lineage=='B.1.525'),aes(fill='B.1.525',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=4, alpha=1)+
  geom_point(data=subset(df_africa, pango_lineage=='A.23.1'),aes(fill='A.23.1',x=days, y=reorder(country,Count)), position = position_jitter(width=0.2, height=0.2), stroke=0.2,shape=21, col='black', size=4, alpha=1)+
  ylab('')+ xlab('month')+ ggtitle('African countries with sequencing data') +
  scale_fill_manual(values=c('peachpuff3','grey','mediumorchid3','cadetblue3','white'), name='VOCs')+
  theme(legend.position="bottom") +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=11))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "1 month")+
  guides(fill = guide_legend(override.aes = list(size=5)))+
  xlab('Sampling Dates')

panelD



