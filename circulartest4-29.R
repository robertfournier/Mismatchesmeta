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
library(circular)
library(Directional)

### All stations get data

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

#dfFMWT <- left_join(surv_FMWT, fish) %>% 
  #collect() 

dfBay <- mutate(dfBay, Julian = yday(Date))
dfBay<- mutate(dfBay, Year = format(dfBay[,"Date"], "%Y"))
dfBay<- mutate(dfBay, Month = format(dfBay[,"Date"], "%m-%Y"))

##density plots


Cit1 <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))

Citsum1<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit1, sum)

Citsummax1 <- Citsum1 %>% group_by(Year, Month, Station) %>% slice_max(Count)


mean<-Citsummax1 %>% group_by(Taxa) %>%
  summarize(mean=mean(Julian))

p<-ggplot(Citsummax1, aes(Julian)) + geom_density() +geom_vline(data = mean, aes(xintercept = mean), size=0.5) +theme_classic()
p

mean<-Citsummax1 %>% group_by(Station) %>%
  summarize(mean=mean(Julian))

p<-ggplot(Citsummax1, aes(Julian)) + geom_density() +geom_vline(data = mean, aes(xintercept = mean), size=0.5) +facet_wrap(~Station)+theme_classic()
p


##Completeness check

Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211")) %>% filter(Julian %in% (60:273))

Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)


Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Temp_surf)

p<-ggplot(Citsum, aes(x = Month, y = Count, group = 1,)) +
  geom_smooth(method=lm)+ geom_point()+theme_classic() + facet_wrap(~Station)
p


##Full impute

Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))

Cit <- dfBay %>% 
  filter(Taxa =="Engraulis mordax") %>% filter(!(Station %in% c("324", "326", "427", "428", "429", "430","431","432","433","534")))

Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Temp_surf)

p<-ggplot(Citsum, aes(x = Month, y = Count, group = 1,)) +
 geom_smooth(method=lm)+ geom_point()+theme_classic() + facet_wrap(~Station)
p

Citsumall<-Citsummax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()


p<-ggplot(Citsummax1, aes(Julian)) + geom_density()
p

#all stations impute count

Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))

Cit <- dfBay %>% 
  filter(Taxa =="Clupea pallasii") %>% filter(Station %in% c("108", "109", "110", "211", "213", "214","215","216","317","318"))



Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temp_surf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Citsumall<-Citsummax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()

Citcounts<-select(Citsumall, c("Station", "Date_month", "Count"))

Citsumwide<-Citcounts %>% distinct() %>%
  pivot_wider(names_from = Station, values_from = Count, values_fill=NA)


Citimpdate1 <- Citsumwide %>% 
  arrange(Date_month) #%>%  # ensure it's in chronological order



Citimpdate <- Citimpdate1[,-1]

Citimpdate <- as.ts(Citimpdate)


Time_series_inter <- matrix(NA,nrow(Citimpdate),(ncol(Citimpdate))) #create empty matrix
dim(Time_series_inter)
colnames(Time_series_inter) <- colnames(Citimpdate)
head(Time_series_inter)
# Create a vector of years (for each month x year combo = 312)
Year <- lubridate::year(Citimpdate1$Date_month)
## start for-loop, make sure to refresh Time_series_inter above
for (i in 1:ncol(Citimpdate)) {
  #create time-series object
  dat_imp <- data.frame(Year, y = Citimpdate[,i])
  y <- ts(dat_imp[2], start= c(1980,1), end = c(2020,12), frequency = 12) # time series object
  
  #add raw values into empty matrix
  Time_series_inter[,i] <- dat_imp[,c(-1)] # removing the Year column c(-1)
  
  #interpolate missing values (if any)
  if(length(which(is.na(y))) > 0){
    
    #fit ARIMA and impute missing values
    fit <- auto.arima(y, seasonal = TRUE) # don't use lambda, already log10 transformed 
    y_inter <- na_kalman(y,model=fit$model)
    #identify missing values to impute (and replace in the matrix)
    id.na <- which(is.na(y))
    #remove missing values at the begining and end of time series
    aa <- split(id.na, cumsum(c(1, diff(id.na) != 1)))  #split sequences of missing values 
    last <- length(y)
    first <- 1
    is.first <- sapply(aa,function(x) length(which(x %in% first))) #identify series of NAs at the start of the time-series
    is.last <- sapply(aa,function(x) length(which(x %in% last))) #idem with last value
    
    if(sum(is.first) > 0 | sum(is.last) > 0){
      aa <- aa[-c(which(is.last == 1),which(is.first == 1))] #remove them from the list
      id.na <- unlist(aa)
    }
    #replace missing values by imputed values into matrix
    Time_series_inter[id.na,i] <- y_inter[id.na]
  }
  print(i)
  
}

plot.ts(Time_series_inter)

Time_series_inter<-as.data.frame(Time_series_inter)

Time_series_inter <- as_tibble(Time_series_inter) %>% 
  mutate(Date_month = Citimpdate1$Date_month)


Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols="211":"215",names_to = "Station", values_to = "Count_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumall, by=c("Date_month", "Station"))

fit

forecast::checkresiduals(fit) 

##all stations temp

Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))

Cit <- dfBay %>% 
  filter(Taxa =="Clupea pallasii") %>% filter(Station %in% c("108", "109", "110", "211", "213", "214","215","216","317","318"))

Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temp_surf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Temp_surf)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(rolltemp)

Citsumall<-Cittempmax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()



#Cittemp<-select(Citsumall, c("Station", "Date_month", "Temp_surf"))

Cittemp<-select(Citsumall, c("Station", "Date_month", "rolltemp"))

#Cittemp<-Cittemp%>%group_by(Date_month, Station)%>% distinct()

#Citsumwide<-Cittemp %>% distinct() %>%
  #pivot_wider(names_from = Station, values_from =Temp_surf)

Citsumwide<-Cittemp %>% distinct() %>%
  pivot_wider(names_from = Station, values_from =rolltemp)


Citimpdate1 <- Citsumwide %>% 
  arrange(Date_month) #%>%  # ensure it's in chronological order



Citimpdate <- Citimpdate1[,-1]

Citimpdate <- as.ts(Citimpdate)


Time_series_inter <- matrix(NA,nrow(Citimpdate),(ncol(Citimpdate))) #create empty matrix
dim(Time_series_inter)
colnames(Time_series_inter) <- colnames(Citimpdate)
head(Time_series_inter)
# Create a vector of years (for each month x year combo = 312)
Year <- lubridate::year(Citimpdate1$Date_month)
## start for-loop, make sure to refresh Time_series_inter above
for (i in 1:ncol(Citimpdate)) {
  #create time-series object
  dat_imp <- data.frame(Year, y = Citimpdate[,i])
  y <- ts(dat_imp[2], start= c(1980,1), end = c(2020,12), frequency = 12) # time series object
  
  #add raw values into empty matrix
  Time_series_inter[,i] <- dat_imp[,c(-1)] # removing the Year column c(-1)
  
  #interpolate missing values (if any)
  if(length(which(is.na(y))) > 0){
    
    #fit ARIMA and impute missing values
    fit <- auto.arima(y, seasonal = TRUE) # don't use lambda, already log10 transformed 
    y_inter <- na_kalman(y,model=fit$model)
    #identify missing values to impute (and replace in the matrix)
    id.na <- which(is.na(y))
    #remove missing values at the begining and end of time series
    aa <- split(id.na, cumsum(c(1, diff(id.na) != 1)))  #split sequences of missing values 
    last <- length(y)
    first <- 1
    is.first <- sapply(aa,function(x) length(which(x %in% first))) #identify series of NAs at the start of the time-series
    is.last <- sapply(aa,function(x) length(which(x %in% last))) #idem with last value
    
    if(sum(is.first) > 0 | sum(is.last) > 0){
      aa <- aa[-c(which(is.last == 1),which(is.first == 1))] #remove them from the list
      id.na <- unlist(aa)
    }
    #replace missing values by imputed values into matrix
    Time_series_inter[id.na,i] <- y_inter[id.na]
  }
  print(i)
  
}

plot.ts(Time_series_inter)

Time_series_inter<-as.data.frame(Time_series_inter)

Time_series_inter <- as_tibble(Time_series_inter) %>% 
  mutate(Date_month = Citimpdate1$Date_month)

Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols="211":"215",names_to = "Station", values_to = "Temp_interp")

#Time_series_inter<-Time_series_inter %>%
  #pivot_longer(cols="108":"318",names_to = "Station", values_to = "Temp_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))

##All stations sal


Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))

Cit <- dfBay %>% 
  filter(Taxa =="Clupea pallasii") %>% filter(Station %in% c("108", "109", "110", "211", "213", "214","215","216","317","318"))

Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum <- Citsum %>% mutate(Citsum, rollsal = c(rep(NA, 3 - 1), rollmean(Sal_surf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

#CitSalmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Sal_surf)

CitSalmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(rollsal)

Citsumall<-CitSalmax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()



#Cittemp<-select(Citsumall, c("Station", "Date_month", "Sal_surf"))

Cittemp<-select(Citsumall, c("Station", "Date_month", "rollsal"))

#Cittemp<-Cittemp%>%group_by(Date_month, Station)%>% distinct()

#Citsumwide<-Cittemp %>% distinct() %>%
  #pivot_wider(names_from = Station, values_from =Sal_surf)

Citsumwide<-Cittemp %>% distinct() %>%
  pivot_wider(names_from = Station, values_from =rollsal)


Citimpdate1 <- Citsumwide %>% 
  arrange(Date_month) #%>%  # ensure it's in chronological order



Citimpdate <- Citimpdate1[,-1]

Citimpdate <- as.ts(Citimpdate)


Time_series_inter <- matrix(NA,nrow(Citimpdate),(ncol(Citimpdate))) #create empty matrix
dim(Time_series_inter)
colnames(Time_series_inter) <- colnames(Citimpdate)
head(Time_series_inter)
# Create a vector of years (for each month x year combo = 312)
Year <- lubridate::year(Citimpdate1$Date_month)
## start for-loop, make sure to refresh Time_series_inter above
for (i in 1:ncol(Citimpdate)) {
  #create time-series object
  dat_imp <- data.frame(Year, y = Citimpdate[,i])
  y <- ts(dat_imp[2], start= c(1980,1), end = c(2020,12), frequency = 12) # time series object
  
  #add raw values into empty matrix
  Time_series_inter[,i] <- dat_imp[,c(-1)] # removing the Year column c(-1)
  
  #interpolate missing values (if any)
  if(length(which(is.na(y))) > 0){
    
    #fit ARIMA and impute missing values
    fit <- auto.arima(y, seasonal = TRUE) # don't use lambda, already log10 transformed 
    y_inter <- na_kalman(y,model=fit$model)
    #identify missing values to impute (and replace in the matrix)
    id.na <- which(is.na(y))
    #remove missing values at the begining and end of time series
    aa <- split(id.na, cumsum(c(1, diff(id.na) != 1)))  #split sequences of missing values 
    last <- length(y)
    first <- 1
    is.first <- sapply(aa,function(x) length(which(x %in% first))) #identify series of NAs at the start of the time-series
    is.last <- sapply(aa,function(x) length(which(x %in% last))) #idem with last value
    
    if(sum(is.first) > 0 | sum(is.last) > 0){
      aa <- aa[-c(which(is.last == 1),which(is.first == 1))] #remove them from the list
      id.na <- unlist(aa)
    }
    #replace missing values by imputed values into matrix
    Time_series_inter[id.na,i] <- y_inter[id.na]
  }
  print(i)
  
}

plot.ts(Time_series_inter)

Time_series_inter<-as.data.frame(Time_series_inter)

Time_series_inter <- as_tibble(Time_series_inter) %>% 
  mutate(Date_month = Citimpdate1$Date_month)

Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols="211":"215",names_to = "Station", values_to = "Sal_interp")

Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols="108":"318",names_to = "Station", values_to = "Sal_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))

##Pull max with imputed all stations

Citsumfull$Julian[is.na(Citsumfull$Julian)] <- yday(Citsumfull$Date_month[is.na(Citsumfull$Julian)])

Citsumfull$Date[is.na(Citsumfull$Date)] <- Citsumfull$Date_month[is.na(Citsumfull$Date)]

Citmax <- Citsumfull %>% group_by(Year, Taxa, Station) %>% slice_max(Count_interp)

#Citmax <- Citmax %>% filter(!(Count==0))

#Citmax <- Citmax %>% filter(Station==318)

Citmax$Year<-as.Date(Citmax$Year, format="%Y")

Citmax$Julian<-as.circular(Citmax$Julian, units="radians", zero=0)

Citmax$Julian <- circular(Citmax$Julian * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "clock", zero = 0)

Citmax$Julian<-as.numeric(Citmax$Julian)

f<-lm(Julian ~ Date, data = Citmax)
summary(f)

c<-lm.circular(x=Julian, y=Date, data=Citmax, order = 1, level = 0.05, type=('c-c'))
summary(c)

length(Citmax$Date)

lm.circular(y = Citmax$Julian, x = Citmax$Temp_interp)

c<-lm.circular(x=Julian, y=Date, data=Citmax, init = 0, type = "c-l", verbose = F)

cor.circular(y = Citmax$Julian, x = Citmax$Date)

x<-circlin.cor(Citmax$Julian, Citmax$Temp_interp, rads = FALSE)
summary(x)

cor(Citmax$Julian, Citmax$Temp_interp)

class(Citmax$Julian)

f<-lm(Julian ~ Sal_interp, data = Citmax)
summary(f)

f<-lm(Julian ~ Temp_interp, data = Citmax)
summary(f)

p1<-ggplot(Citmax, aes(x = Date, y = Julian, group = 1,)) +
  geom_smooth(method=lm)+ geom_point()+theme_classic() + facet_wrap(~Station)
p1

p2<-ggplot(Citmax, aes(x = Temp_interp, y = Julian, group = 1,)) +
  geom_smooth(method=lm)+ geom_point()+theme_classic() + facet_wrap(~Station)
p2

p3<-ggplot(Citmax, aes(x = Sal_interp, y = Julian, group = 1,)) +
  geom_smooth(method=lm)+ geom_point()+theme_classic() + facet_wrap(~Station)
p3


prow <- plot_grid(
  p1 + theme(text = element_text(size=20)),
  p2   + theme(text = element_text(size=20)),
  p3  + theme(text = element_text(size=20)),
  align = 'hv',
  nrow = 2,
  scale = 1
)
prow