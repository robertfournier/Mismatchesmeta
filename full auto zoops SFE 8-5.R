library(zooper)
library(ggplot2)
library(dplyr)
library(lubridate)
library(sf)
library(deltafish)
library(tidyr)
library(padr)
library(MARSS)
library(imputeTS)
library(zoo)
library(forecast)
library(patchwork)
library(lsr)
library(cowplot)
library(tidyverse)
library(DescTools)

### Get the data
MyZoops <- Zoopsynther(Data_type = "Community", 
                       Sources = c("EMP"), 
                       Size_class = "Macro", 
                       Date_range = c("1972-01-01", "2020-01-01"))



MyZoops$Date<-as.Date(MyZoops$Date, format = "%Y-%m-%d")

MyZoops <- MyZoops[,-24]  
MyZoops<- mutate(MyZoops, Julian = yday(Date))
MyZoops<- mutate(MyZoops, Year = year(Date))
MyZoops<- mutate(MyZoops, Month = month(Date))

###Filter for completeness

Target <- MyZoops%>% 
  filter(Family =="Mysidae") 

Check <- Target %>% filter((CPUE>0)) %>% filter(Julian %in% (60:273))

Check1<-Check%>%              
  group_by(Station)%>%
  mutate(Yeartot=(n_distinct(Year)))%>%
  mutate(Monthtot=(n_distinct(Month)))

Checktot<-Check1%>%  
  group_by(Year, Station)%>%
  summarise(freq=n()
  )
Checktot1<-Checktot%>%filter(freq>4)

Checktot1<-Checktot1%>%  
  group_by(Station)%>% filter((n_distinct(Year)>10))
Checktot1

Stationlist<-unique(Checktot1[c("Station")])
Stationlist<-c(Stationlist)


##Impute the count data

Citsum <- Target %>% filter(Station %in% Stationlist$Station)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temperature, 3)))
Citsum <- Citsum %>% mutate(Citsum, rollsal = c(rep(NA, 3 - 1), rollmean(SalSurf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(CPUE)

Citsumall<-Citsummax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 


Citsumall<-Citsumall %>% distinct()

Citcounts<-select(Citsumall, c("Station", "Date_month", "CPUE"))

Citsumwide<-Citcounts %>% distinct() %>%
  pivot_wider(names_from = Station, values_from = CPUE, values_fill=NA)


Citimpdate1 <- Citsumwide %>% 
  arrange(Date_month) 

Citimpdate <- Citimpdate1[,-1]

Citimpdate <- as.ts(Citimpdate)


Time_series_inter <- matrix(NA,nrow(Citimpdate),(ncol(Citimpdate))) 
dim(Time_series_inter)
colnames(Time_series_inter) <- colnames(Citimpdate)
head(Time_series_inter)
Year <- lubridate::year(Citimpdate1$Date_month)
for (i in 1:ncol(Citimpdate)) {
  dat_imp <- data.frame(Year, y = Citimpdate[,i])
  y <- ts(dat_imp[2], start= c(1975,1), end = c(2020,12), frequency = 12) 
  
  Time_series_inter[,i] <- dat_imp[,c(-1)] 
  
  if(length(which(is.na(y))) > 0){
    
    fit <- auto.arima(y, seasonal = TRUE)  
    y_inter <- na_kalman(y,model=fit$model)
    id.na <- which(is.na(y))
    aa <- split(id.na, cumsum(c(1, diff(id.na) != 1)))
    last <- length(y)
    first <- 1
    is.first <- sapply(aa,function(x) length(which(x %in% first))) 
    is.last <- sapply(aa,function(x) length(which(x %in% last))) 
    
    if(sum(is.first) > 0 | sum(is.last) > 0){
      aa <- aa[-c(which(is.last == 1),which(is.first == 1))] 
      id.na <- unlist(aa)
    }
    Time_series_inter[id.na,i] <- y_inter[id.na]
  }
  print(i)
  
}

plot.ts(Time_series_inter)

Time_series_inter<-as.data.frame(Time_series_inter)

Time_series_inter <- as_tibble(Time_series_inter) %>% 
  mutate(Date_month = Citimpdate1$Date_month)


Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols=!"Date_month",names_to = "Station", values_to = "Count_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumall, by=c("Date_month", "Station"))

fit

forecast::checkresiduals(fit) 

##Impute Temp

Citsum <- Target %>% filter(Station %in% Stationlist$Station)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temperature, 3)))
Citsum <- Citsum %>% mutate(Citsum, rollsal = c(rep(NA, 3 - 1), rollmean(SalSurf, 3)))

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(rolltemp)

Citsumall<-Cittempmax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()


Cittemp<-select(Citsumall, c("Station", "Date_month", "rolltemp"))

#Cittemp<-select(Citsumall, c("Station", "Date_month", "Temperature"))

Citsumwide<-Cittemp %>% distinct() %>%
  pivot_wider(names_from = Station, values_from =rolltemp)

#Citsumwide<-Cittemp %>%
  #pivot_wider(names_from = Station, values_from =Temperature)


Citimpdate1 <- Citsumwide %>% 
  arrange(Date_month) 


Citimpdate <- Citimpdate1[,-1]

Citimpdate <- as.ts(Citimpdate)


Time_series_inter <- matrix(NA,nrow(Citimpdate),(ncol(Citimpdate))) 
dim(Time_series_inter)
colnames(Time_series_inter) <- colnames(Citimpdate)
head(Time_series_inter)

Year <- lubridate::year(Citimpdate1$Date_month)

for (i in 1:ncol(Citimpdate)) {
  dat_imp <- data.frame(Year, y = Citimpdate[,i])
  y <- ts(dat_imp[2], start= c(1975,1), end = c(2020,12), frequency = 12) 
  Time_series_inter[,i] <- dat_imp[,c(-1)]
  if(length(which(is.na(y))) > 0){
    
    fit <- auto.arima(y, seasonal = TRUE) 
    y_inter <- na_kalman(y,model=fit$model)
    id.na <- which(is.na(y))
    aa <- split(id.na, cumsum(c(1, diff(id.na) != 1))) 
    last <- length(y)
    first <- 1
    is.first <- sapply(aa,function(x) length(which(x %in% first))) 
    is.last <- sapply(aa,function(x) length(which(x %in% last))) 
    
    if(sum(is.first) > 0 | sum(is.last) > 0){
      aa <- aa[-c(which(is.last == 1),which(is.first == 1))] 
      id.na <- unlist(aa)
    }
    Time_series_inter[id.na,i] <- y_inter[id.na]
  }
  print(i)
  
}

plot.ts(Time_series_inter)

Time_series_inter<-as.data.frame(Time_series_inter)

Time_series_inter <- as_tibble(Time_series_inter) %>% 
  mutate(Date_month = Citimpdate1$Date_month)

Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols=!"Date_month",names_to = "Station", values_to = "Temp_interp")


Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))

##Salinity Impute

Citsum <- Target %>% filter(Station %in% Stationlist$Station)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temperature, 3)))
Citsum <- Citsum %>% mutate(Citsum, rollsal = c(rep(NA, 3 - 1), rollmean(SalSurf, 3)))

Citsalmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(rollsal)

Citsumall<-Citsalmax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()

#Cittemp<-select(Citsumall, c("Station", "Date_month", "Sal_surf"))

Cittemp<-select(Citsumall, c("Station", "Date_month", "rollsal"))


#Citsumwide<-Cittemp %>% distinct() %>%
#pivot_wider(names_from = Station, values_from =Sal_surf)

Citsumwide<-Cittemp %>% distinct() %>%
  pivot_wider(names_from = Station, values_from =rollsal)


Citimpdate1 <- Citsumwide %>% 
  arrange(Date_month) #%>%  


Citimpdate <- Citimpdate1[,-1]

Citimpdate <- as.ts(Citimpdate)


Time_series_inter <- matrix(NA,nrow(Citimpdate),(ncol(Citimpdate))) 
dim(Time_series_inter)
colnames(Time_series_inter) <- colnames(Citimpdate)
head(Time_series_inter)
Year <- lubridate::year(Citimpdate1$Date_month)
for (i in 1:ncol(Citimpdate)) {
  dat_imp <- data.frame(Year, y = Citimpdate[,i])
  y <- ts(dat_imp[2], start= c(1975,1), end = c(2020,12), frequency = 12)
  Time_series_inter[,i] <- dat_imp[,c(-1)] 
  if(length(which(is.na(y))) > 0){
    
    fit <- auto.arima(y, seasonal = TRUE) 
    y_inter <- na_kalman(y,model=fit$model)
    id.na <- which(is.na(y))
    aa <- split(id.na, cumsum(c(1, diff(id.na) != 1)))   
    last <- length(y)
    first <- 1
    is.first <- sapply(aa,function(x) length(which(x %in% first)))
    is.last <- sapply(aa,function(x) length(which(x %in% last))) 
    
    if(sum(is.first) > 0 | sum(is.last) > 0){
      aa <- aa[-c(which(is.last == 1),which(is.first == 1))]
      id.na <- unlist(aa)
    }
    Time_series_inter[id.na,i] <- y_inter[id.na]
  }
  print(i)
  
}

plot.ts(Time_series_inter)

Time_series_inter<-as.data.frame(Time_series_inter)

Time_series_inter <- as_tibble(Time_series_inter) %>% 
  mutate(Date_month = Citimpdate1$Date_month)

Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols=!"Date_month",names_to = "Station", values_to = "Sal_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))


###################################
Fullimpute<-Citsumfull


#####Correlation loop

Corrframe <- Fullimpute

Corrframe<-Corrframe[complete.cases(Corrframe$Count_interp),]

Corrframe$Year<-as.Date(Corrframe$Year, format="%Y")

Corrframe$Julian[is.na(Corrframe$Julian)] <- yday(Corrframe$Date_month[is.na(Corrframe$Julian)])

Corrframe$Date[is.na(Corrframe$Date)] <- Corrframe$Date_month[is.na(Corrframe$Date)]

Corrframe$Year[is.na(Corrframe$Year)] <- year(Corrframe$Date)[is.na(Corrframe$Date)]

Corrframe$Year<-year(Corrframe$Date)

Corrframe$Year<-as.Date(Corrframe$Year, format="%Y")

species.day = function(x, start.month = 12L, start.day=25L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, start.day)
  as.integer(difftime(x, start.date))
}

Corrframe<- mutate(Corrframe, Speciesday = species.day(Date))

Corrframe<- Corrframe %>% group_by(Station)%>%mutate(Mean.temp = mean(Temp_interp))

Corrframe<- Corrframe %>% group_by(Station)%>%mutate(Mean.sal = mean(Sal_interp))

Corrframemax <- Corrframe %>% group_by(Year, Station) %>% slice_max(Count_interp)

Corrframemax <-Corrframemax [complete.cases(Corrframemax$Year),]


Corrframemax$Year<-as.numeric(Corrframemax$Year)


st_name <- list() 
PT <- list()  
PC <- list()  
TC <- list()

Station<-unique(Corrframemax$Station)

for(i in unique(Station)) { 
  Newstation <- Corrframemax %>% filter(Station==i)
  st_name<- c(st_name, i)
  PT<- c(PT, print(cor(Newstation$Year, Newstation$Speciesday, method = "pearson")))
  PC<- c(PC, print(cor(Newstation$Temp_interp, Newstation$Speciesday, method = "pearson")))  
  TC<- c(TC, print(cor(Newstation$Year, Newstation$Temp_interp, method = "pearson")))
  results.test<-cbind(st_name,PT,PC,TC)
  
} 

results.test<-as.data.frame(results.test)


##Merge with remaining covariates

Corrframemax<-Corrframemax %>% group_by(Station) %>% mutate(Start = min(Date_month)) %>% mutate(End = max(Date_month))  %>% mutate(Length = (year(End))-(year(Start)))

Corrframelite<-Corrframemax %>% group_by(Station) %>% select(Station, Start, End, Length, Latitude, Longitude, Class, Order, Family, Mean.temp, Mean.sal)

Corrframelite<-Corrframelite %>% drop_na(Latitude)

Corrframelite<-Corrframelite[!duplicated(Corrframelite[,c("Station")]),]

Corrframelite<-as.data.frame(Corrframelite)

Corrcomplete <- merge(x=Corrframelite, y=results.test, by.x='Station', by.y='st_name', all.y = TRUE)

Corrcomplete<-Corrcomplete %>% gather("Test", "Value", c(PT, PC, TC))

Corrcomplete$Value<-as.numeric(Corrcomplete$Value)

Corrcomplete<- Corrcomplete %>% mutate(Z_cor = FisherZ(Value))





###
write.csv(Corrcomplete, "zoopstations.csv")




##Manual tests

Corrframe <- MYSIDAE ##%>% filter(Station=="NZ044")

Corrframe<-Corrframe[complete.cases(Corrframe$Count_interp),]

Corrframe$Year<-as.Date(Corrframe$Year, format="%Y")

Corrframe$Julian[is.na(Corrframe$Julian)] <- yday(Corrframe$Date_month[is.na(Corrframe$Julian)])

Corrframe$Date[is.na(Corrframe$Date)] <- Corrframe$Date_month[is.na(Corrframe$Date)]

Corrframe$Year[is.na(Corrframe$Year)] <- year(Corrframe$Date)[is.na(Corrframe$Date)]

Corrframe$Year<-year(Corrframe$Date)

Corrframe$Year<-as.Date(Corrframe$Year, format="%Y")

Corrframemax <- Corrframe %>% group_by(Year, Station) %>% slice_max(Count_interp)

Corrframemax <-Corrframemax [complete.cases(Corrframemax$Year),]


Corrframemax$Year<-as.numeric(Corrframemax$Year)

##Correlation test Temp

Corrframemax$Year<-as.numeric(Corrframemax$Year)

PT<-cor(Corrframemax$Year, Corrframemax$Julian, method = c("pearson"))
PC<-cor(Corrframemax$Temp_interp, Corrframemax$Julian, method = c("pearson"))
TC<-cor(Corrframemax$Year, Corrframemax$Temp_interp, method = c("pearson"))

PT
PC
TC

FisherZ(PT)
FisherZ(PC)
FisherZ(TC)

##Correlation test Sal

PT<-cor(Corrframemax$Year, Corrframemax$Julian, method = c("pearson"))
PC<-cor(Corrframemax$Sal_interp, Corrframemax$Julian, method = c("pearson"))
TC<-cor(Corrframemax$Year, Corrframemax$Sal_interp, method = c("pearson"))

PT
PC
TC

FisherZ(PT)
FisherZ(PC)
FisherZ(TC)

##for loop correlation


Station<-Corrframemax$Station


