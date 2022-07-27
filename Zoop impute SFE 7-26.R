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



##Impute the count data

Citsum <- MyZoops%>% 
  filter(Family =="Mysidae") %>% filter(Station %in% c("NZ020",
                                                       "NZ022",
                                                       "NZ024",
                                                       "NZ028",
                                                       "NZ030",
                                                       "NZ032",
                                                       "NZ034",
                                                       "NZ036",
                                                       "NZ038",
                                                       "NZ040",
                                                       "NZ042",
                                                       "NZ044",
                                                       "NZ046",
                                                       "NZ048",
                                                       "NZ050",
                                                       "NZ052",
                                                       "NZ054",
                                                       "NZ056",
                                                       "NZ058",
                                                       "NZ060",
                                                       "NZ062",
                                                       "NZ064",
                                                       "NZ066",
                                                       "NZ068",
                                                       "NZ072",
                                                       "NZ074",
                                                       "NZ076",
                                                       "NZ078",
                                                       "NZ080",
                                                       "NZ082",
                                                       "NZ084",
                                                       "NZ086",
                                                       "NZ088",
                                                       "NZ090",
                                                       "NZ092",
                                                       "NZ098",
                                                       "NZ104",
                                                       "NZD11",
                                                       "NZD14",
                                                       "NZD15",
                                                       "NZD19",
                                                       "NZD28",
                                                       "NZC03",
                                                       "NZC09",
                                                       "NZM10",
                                                       "NZS42",
                                                       "NZ002",
                                                       "NZ004",
                                                       "NZD06",
                                                       "NZEZ2",
                                                       "NZEZ6",
                                                       "NZD16",
                                                       "NZD41",
                                                       "NZ41A"))


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
  pivot_longer(cols="NZ002":"NZS42",names_to = "Station", values_to = "Count_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumall, by=c("Date_month", "Station"))

fit

forecast::checkresiduals(fit) 

##Impute Temp

Citsum <- MyZoops%>% 
  filter(Family =="Mysidae") %>% filter(Station %in% c("NZ020",
                                                       "NZ022",
                                                       "NZ024",
                                                       "NZ028",
                                                       "NZ030",
                                                       "NZ032",
                                                       "NZ034",
                                                       "NZ036",
                                                       "NZ038",
                                                       "NZ040",
                                                       "NZ042",
                                                       "NZ044",
                                                       "NZ046",
                                                       "NZ048",
                                                       "NZ050",
                                                       "NZ052",
                                                       "NZ054",
                                                       "NZ056",
                                                       "NZ058",
                                                       "NZ060",
                                                       "NZ062",
                                                       "NZ064",
                                                       "NZ066",
                                                       "NZ068",
                                                       "NZ072",
                                                       "NZ074",
                                                       "NZ076",
                                                       "NZ078",
                                                       "NZ080",
                                                       "NZ082",
                                                       "NZ084",
                                                       "NZ086",
                                                       "NZ088",
                                                       "NZ090",
                                                       "NZ092",
                                                       "NZ098",
                                                       "NZ104",
                                                       "NZD11",
                                                       "NZD14",
                                                       "NZD15",
                                                       "NZD19",
                                                       "NZD28",
                                                       "NZC03",
                                                       "NZC09",
                                                       "NZM10",
                                                       "NZS42",
                                                       "NZ002",
                                                       "NZ004",
                                                       "NZD06",
                                                       "NZEZ2",
                                                       "NZEZ6",
                                                       "NZD16",
                                                       "NZD41",
                                                       "NZ41A"))



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
  pivot_longer(cols="NZ002":"NZS42",names_to = "Station", values_to = "Temp_interp")


Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))

##Salinity Impute

Citsum <- MyZoops%>% 
  filter(Family =="Mysidae") %>% filter(Station %in% c("NZ020",
                                                       "NZ022",
                                                       "NZ024",
                                                       "NZ028",
                                                       "NZ030",
                                                       "NZ032",
                                                       "NZ034",
                                                       "NZ036",
                                                       "NZ038",
                                                       "NZ040",
                                                       "NZ042",
                                                       "NZ044",
                                                       "NZ046",
                                                       "NZ048",
                                                       "NZ050",
                                                       "NZ052",
                                                       "NZ054",
                                                       "NZ056",
                                                       "NZ058",
                                                       "NZ060",
                                                       "NZ062",
                                                       "NZ064",
                                                       "NZ066",
                                                       "NZ068",
                                                       "NZ072",
                                                       "NZ074",
                                                       "NZ076",
                                                       "NZ078",
                                                       "NZ080",
                                                       "NZ082",
                                                       "NZ084",
                                                       "NZ086",
                                                       "NZ088",
                                                       "NZ090",
                                                       "NZ092",
                                                       "NZ098",
                                                       "NZ104",
                                                       "NZD11",
                                                       "NZD14",
                                                       "NZD15",
                                                       "NZD19",
                                                       "NZD28",
                                                       "NZC03",
                                                       "NZC09",
                                                       "NZM10",
                                                       "NZS42",
                                                       "NZ002",
                                                       "NZ004",
                                                       "NZD06",
                                                       "NZEZ2",
                                                       "NZEZ6",
                                                       "NZD16",
                                                       "NZD41",
                                                       "NZ41A"))


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
  pivot_longer(cols="NZ002":"NZS42",names_to = "Station", values_to = "Sal_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))


###################################
MYSIDAE<-Citsumfull

##Pull max with imputed data

Corrframe <- MYSIDAE %>% filter(Station=="NZ022")

Corrframe<-Corrframe[complete.cases(Corrframe$Count_interp),]

Corrframe$Year<-as.Date(Corrframe$Year, format="%Y")

Corrframe$Julian[is.na(Corrframe$Julian)] <- yday(Corrframe$Date_month[is.na(Corrframe$Julian)])

Corrframe$Date[is.na(Corrframe$Date)] <- Corrframe$Date_month[is.na(Corrframe$Date)]

Corrframe$Year[is.na(Corrframe$Year)] <- year(Corrframe$Date)[is.na(Corrframe$Date)]

Corrframe$Year<-year(Corrframe$Date)

Corrframe$Year<-as.Date(Corrframe$Year, format="%Y")

Corrframemax <- Corrframe %>% group_by(Year) %>% slice_max(Count_interp)

Corrframemax <-Corrframemax [complete.cases(Corrframemax$Year),]

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