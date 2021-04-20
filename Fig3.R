library("readxl")
library(lubridate)
library("scales")
library(ggplot2)


ggplotRegression <- function(fit){
  
  
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       #"Intercept =",signif(fit$coef[[1]],5 ),
                       "R = ",signif(summary(fit)$coef, 5),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

#Panel A - tempest plots
#Note: Same code for all 4 plots, just substitute the data file for the required cluster
#Further edits can be done in power point or illustrator

tempest_data<-read_excel('data/C.1_tempest.xlsx')
tempest_data$date2<-date_decimal(tempest_data$date)

tempest_data$date2<-as.Date(cut(tempest_data$date2,
                                breaks = "day",
                                start.on.monday = FALSE))
custom2<-c("tan4",'antiquewhite2','darkolivegreen',
           'red3','plum2',
           'purple3','dodgerblue1','slategray1', 'white')

custom3<-c("tan4",'antiquewhite2',
           'red3','plum2',
           'purple3','dodgerblue1','slategray1', 'white')

custom4<-c("tan4",'plum2',
           'dodgerblue1')

custom1<-c('plum2',
           'slategray1', 'white')


ggplotRegression(lm(distance ~ date2, data = tempest_data))


p_tempest<-ggplot(tempest_data, aes(date2,distance))+
  geom_point(color='black',size=2.5,aes(fill=location),shape=21, stroke=0.3)+
   geom_smooth(method=lm,se=T, color='black')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=10))+
  theme(axis.text.x =element_text(size=10))+
  theme(axis.title=element_text(size=10))+
  scale_x_date(date_labels = "%b\n%Y",breaks='3 month')+
  scale_fill_manual(values=custom4, name='Location')+
  ylab("Root-to-tip Distance")+
  xlab("Date")+
  scale_y_continuous(labels = scientific)+
  annotate(geom="text", x=as.Date("2020/07/01"),y=0.001,label="r = 0.82",
           color="black",size=4)+
  annotate(geom="text", x=as.Date("2020/07/01"),y=0.0009,label="r2 = 0.67",
           color="black",size=4)

p_tempest


#Panel B - Done with a JSON code and edited in Illustrator

