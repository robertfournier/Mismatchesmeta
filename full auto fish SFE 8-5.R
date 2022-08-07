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

### Get the data

create_fish_db()
surv <- open_survey()
fish <- open_fish()

fish<-as.data.frame(fish)
surv<-as.data.frame(surv)
surv$Date<-as.Date(surv$Date, format = "%Y-%m-%d")

surv_Bay <- surv %>% 
  filter(Source == "Bay Study")

dfBay <- left_join(surv_Bay, fish) %>% 
  collect() 

dfBay <- mutate(dfBay, Julian = yday(Date))
dfBay<- mutate(dfBay, Year = format(dfBay[,"Date"], "%Y"))
dfBay<- mutate(dfBay, Month = format(dfBay[,"Date"], "%m-%Y"))


##Filter for completeness

Target <- dfBay %>% 
  filter(Taxa =="Clupea pallasii")

Check <- Target %>% filter((Count>0)) %>% filter(Julian %in% (60:273))

Check<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Check, sum)

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

Stationlist<-unique(Checktot1[c("Station")])
Stationlist


##Impute the count data

Cit <- Target %>% filter(Station %in% Stationlist$Station)


Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temp_surf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Citsumall<-Citsummax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()

Citcounts<-select(Citsumall, c("Station", "Date_month", "Count"))

Citsumwide<-Citcounts %>% distinct() %>%
  pivot_wider(names_from = Station, values_from = Count, values_fill=NA)


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
  y <- ts(dat_imp[2], start= c(1980,1), end = c(2020,12), frequency = 12) 
  
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

Cit <- Target %>% filter(Station %in% Stationlist$Station)

Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temp_surf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Temp_surf)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(rolltemp)

Citsumall<-Cittempmax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()


Cittemp<-select(Citsumall, c("Station", "Date_month", "rolltemp"))

Citsumwide<-Cittemp %>% distinct() %>%
  pivot_wider(names_from = Station, values_from =rolltemp)


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
  y <- ts(dat_imp[2], start= c(1980,1), end = c(2020,12), frequency = 12) 
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

Cit <- Target %>% filter(Station %in% Stationlist$Station)

Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum <- Citsum %>% mutate(Citsum, rollsal = c(rep(NA, 3 - 1), rollmean(Sal_surf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

#CitSalmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Sal_surf)

CitSalmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(rollsal)

Citsumall<-CitSalmax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

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
  y <- ts(dat_imp[2], start= c(1980,1), end = c(2020,12), frequency = 12)
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


###
Fullimpute <-Citsumfull

#####NEW

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


##merge with covariates

Corrframemax<-Corrframemax %>% group_by(Station) %>% mutate(Start = min(Date_month)) %>% mutate(End = max(Date_month))  %>% mutate(Length = (year(End))-(year(Start)))

Corrframelite<-Corrframemax %>% group_by(Station) %>% select(Station, Start, End, Length, Latitude, Longitude, Class, Order, Family, Mean.temp, Mean.sal)

Corrframelite<-Corrframelite %>% drop_na(Latitude)

Corrframelite<-Corrframelite[!duplicated(Corrframelite[,c("Station")]),]

Corrframelite<-as.data.frame(Corrframelite)

Corrcomplete <- merge(x=Corrframelite, y=results.test, by.x='Station', by.y='st_name', all.y = TRUE)

Corrcomplete<-Corrcomplete %>% gather("Test", "Value", c(PT, PC, TC))

Corrcomplete$Value<-as.numeric(Corrcomplete$Value)

Corrcomplete<- Corrcomplete %>% mutate(Z_cor = FisherZ(Value))
