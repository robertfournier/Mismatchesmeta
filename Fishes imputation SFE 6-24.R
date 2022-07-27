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


##Impute the count data

Cit <- dfBay %>% 
  filter(Taxa =="Spirinchus thaleichthys") %>% filter(Station %in% c("345", "736"))


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
  pivot_longer(cols="345":"736",names_to = "Station", values_to = "Count_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumall, by=c("Date_month", "Station"))

fit

forecast::checkresiduals(fit) 

##Impute Temp

Cit <- dfBay %>% 
  filter(Taxa =="Spirinchus thaleichthys") %>% filter(Station %in% c("345", "736"))

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
  pivot_longer(cols="345":"736",names_to = "Station", values_to = "Temp_interp")


Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))

##Salinity Impute

Cit <- dfBay %>% 
  filter(Taxa =="Spirinchus thaleichthys") %>% filter(Station %in% c("345", "736"))

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
  pivot_longer(cols="345":"736",names_to = "Station", values_to = "Sal_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))


###
SPECIESNAME<-Citsumfull

##Pull max with imputed data

Citsumfull <- Citsumfull %>% filter(Station==318)

Citsumfull$Year<-as.Date(Citsumfull$Year, format="%Y")

Citsumfull$Julian[is.na(Citsumfull$Julian)] <- yday(Citsumfull$Date_month[is.na(Citsumfull$Julian)])

Citsumfull$Date[is.na(Citsumfull$Date)] <- Citsumfull$Date_month[is.na(Citsumfull$Date)]

Citmax <- Citsumfull %>% group_by(Year) %>% slice_max(Count_interp)
