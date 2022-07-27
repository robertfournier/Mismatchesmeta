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

#dfFMWT <- left_join(surv_FMWT, fish) %>% 
#collect() 

dfBay <- mutate(dfBay, Julian = yday(Date))
dfBay<- mutate(dfBay, Year = year(Date))
dfBay<- mutate(dfBay, Month = format(dfBay[,"Date"], "%m-%Y"))

##Zooplankton
MyZoops <- Zoopsynther(Data_type = "Community", 
                       Sources = c("EMP"), 
                       Size_class = "Micro", 
                       Date_range = c("1972-01-01", "2020-01-01"))

str(MyZoops)

MyZoops<-MyZoops%>%filter(!is.na(Family))


###Completeness check initial Fish

st<-unique(dfBay[c("Station")])
st

unique(dfBay[c("Taxa")])

Check <- dfBay %>% 
  filter(Taxa =="Engraulis mordax") %>% filter((Count>0))%>% filter(Julian %in% (60:273))

Check1<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Check, sum)

Checktot<-Check1%>%              
  group_by(Station)%>%
  mutate(Yeartot=(n_distinct(Year)))%>%
  mutate(Monthtot=(n_distinct(Month)))%>%
  mutate(Complete=Monthtot/(Yeartot*6))

Checktot1<-Checktot%>%filter((Yeartot >= 10))%>%filter((Complete >= 0.8))
Checktot1



##Zoop

st<-unique(MyZoops[c("Station")])
st
write.csv(st, "zoopstations.csv")

fa<-unique(MyZoops[c("Family")])
fa

Check <- MyZoops%>% 
  filter(Family =="Pseudodiaptomidae")

Check1 <- 
  Check %>%                   
  group_by(Station) %>%
  summarise(Year = n_distinct(Year))
Check1
write.csv(Check1, "zoopstations.csv")


##density plots to look at peaks


Cit1 <- dfBay %>% 
  filter(Taxa =="Cymatogaster aggregata") %>% filter(Station %in% c("102", "103", "106", "107", "108", "142", "211", "212", "214", "215", "317"))

Cit1 <- dfBay %>% 
  filter(Taxa =="Acanthogobius flavimanus") %>% filter(Station %in% c("430", "431", "447"))

Cit1 <- dfBay %>% 
  filter(Taxa =="Alosa sapidissima") %>% filter(Station %in% c("431", "534", "736", "750", "751", "752", "760", "761", "837", "853", "863", "864", "865"))

Cit1 <- dfBay %>% 
  filter(Taxa =="Atherinops affinis") %>% filter(Station %in% c("106"))
Cit1 <- dfBay %>% 
  filter(Taxa =="Atherinopsis californiensis") %>% filter(Station %in% c("103", "106", "108", "109", "110", "142", "211", "212", "214" ))

Cit1 <- dfBay %>% 
  filter(Taxa =="Spirinchus thaleichthys") %>% filter(Station %in% c("110", "211", "214", "215", "216", "320", "321", "322", "323" , "325" , "346" , "427" , "428" , "429" , "432"  ))

Cit1 <- dfBay %>% 
  filter(Taxa =="Morone saxatilis") %>% filter(Station %in% c("319", "427", "428", "429", "430", "431", "432", "447", "534" , "736" , "760" , "863" , "865" ))

Cit1 <- dfBay %>% 
  filter(Taxa =="Hyperprosopon argenteum") %>% filter(Station %in% c("212" ))

Cit1 <- dfBay %>% 
  filter(Taxa =="Genyonemus lineatus") %>% filter(Station %in% c("107", "110", "211", "214", "215"))

Cit1 <- dfBay %>% 
  filter(Taxa =="Clupea pallasii") %>% filter(Station %in% c("101", "102", "103", "105", "106", "107", "108", "109", "110" , "142" , "211" , "212" , "213" , "214" , "215" , "216" , "243", "244", "317" , "318" , "319" , "320" , "321" , "322", "323" , "325" , "345"))

Cit1 <- dfBay %>% 
  filter(Taxa =="Engraulis mordax") %>% filter(Station %in% c("101", "102", "103", "104", "105", "106", "107", "108", "109", "110" , "140" , "142" , "211" , "212" , "213" , "214" , "215" , "216" , "243", "244", "317" , "318" , "319" , "320" , "321" , "322", "323" , "325" , "345" , "346", "427", "428" , "432", "447"))



species.day = function(x, start.month = 12L, start.day=25L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, start.day)
  as.integer(difftime(x, start.date))
}

Cit1<- mutate(Cit1, Speciesday = species.day(Date))

Citsum1<-aggregate(Count ~ Taxa + Year + Speciesday + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit1, sum)

Citsummax1 <- Citsum1 %>% group_by(Year, Month, Station) %>% slice_max(Count)

mean<-Citsummax1 %>% group_by(Taxa) %>%
  summarize(mean=mean(Julian))
mean

p<-ggplot(Citsummax1, aes(Julian)) + geom_density() +geom_vline(data = mean, aes(xintercept = mean), size=0.5) +theme_classic()
p

mean<-Citsummax1 %>% group_by(Station) %>%
  summarize(mean=mean(Julian))
mean

p<-ggplot(Citsummax1, aes(Julian)) + geom_density() +geom_vline(data = mean, aes(xintercept = mean), size=0.5) +facet_wrap(~Station)+theme_classic()
p


#####NEW Completeness check by species

Cit <- dfBay %>% 
  filter(Taxa =="Hyperprosopon ellipticum") %>% filter((Count>0)) %>% filter(Julian %in% (60:273))


Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum<-Citsum%>%              
  group_by(Station)%>%
  mutate(Yeartot=(n_distinct(Year)))%>%
  mutate(Monthtot=(n_distinct(Month)))

Checktot<-Citsum%>%  
  group_by(Year, Station)%>%
  summarise(freq=n()
  )
Checktot1<-Checktot%>%filter(freq>4)

Checktot1<-Checktot1%>%  
  group_by(Station)%>% filter((n_distinct(Year)>10))
Checktot1

unique(Checktot1[c("Station")])

Checktot2<-Checktot1%>%group_by(Station)%>%summarise(Year = n_distinct(Year))

###


Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Temp_surf)

p<-ggplot(Citsum, aes(x = Month, y = Count, group = 1,)) +
  geom_smooth(method=lm)+ geom_point()+theme_classic() + facet_wrap(~Station)
p

##Zoops

MyZoops <- Zoopsynther(Data_type = "Community", 
                       Sources = c("EMP"), 
                       Size_class = "Micro", 
                       Date_range = c("1972-01-01", "2020-01-01"))

MyZoops<- mutate(MyZoops, Julian = yday(Date))

Check <- MyZoops%>% 
  filter(Family =="Centropagidae") %>% filter((CPUE>0)) %>% filter(Julian %in% (60:273))


Checktot<-Check%>%  
  group_by(Station)%>%
  summarise(freq=n()
  )
Checktot1<-mutate(Checktot, Years=(Check %>%                 
  group_by(Station) %>%
  summarise(Year = n_distinct(Year))))

Checktot1<-mutate(Checktot1, Expected=Years$Year * 7)
Checktot1<-mutate(Checktot1, Complete=freq/Expected)
Checktot1
write.csv(Checktot1, "zoopstations.csv")


###NEW ZOOPS FILTERING
MyZoops <- Zoopsynther(Data_type = "Community", 
                       Sources = c("EMP"), 
                       Size_class = "Micro", 
                       Date_range = c("1972-01-01", "2020-01-01"))

MyZoops<- mutate(MyZoops, Julian = yday(Date))
MyZoops<- mutate(MyZoops, Year = year(Date))
MyZoops<- mutate(MyZoops, Month = month(Date))
MyZoops<- mutate(MyZoops, Month = as.Date(MyZoops$Date, format ="%m-%Y"))



Check <- MyZoops%>% 
  filter(Family =="Centropagidae") %>% filter((CPUE>0)) %>% filter(Julian %in% (60:273))

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

Stations<-unique(Checktot1[c("Station")])
Stations

write.csv(Stations, "stationscheck.csv") 

Checktot2<-Checktot1%>%group_by(Station)%>%summarise(Year = n_distinct(Year))
Checktot2

##Impute the count data

Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))


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
  pivot_longer(cols="211":"215",names_to = "Station", values_to = "Count_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumall, by=c("Date_month", "Station"))

fit

forecast::checkresiduals(fit) 

##Impute Temp

Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))

Citsum<-aggregate(Count ~ Taxa + Year + Temp_surf + Station +  Julian + Latitude + Longitude +Sal_surf + Month + Date, Cit, sum)

Citsum <- Citsum %>% mutate(Citsum, rolltemp = c(rep(NA, 3 - 1), rollmean(Temp_surf, 3)))

Citsummax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Count)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(Temp_surf)

Cittempmax <- Citsum %>% group_by(Year, Month, Station) %>% slice_max(rolltemp)

Citsumall<-Cittempmax%>% group_by(Station)%>%thicken(interval = "month") %>% pad(by="Date_month") 

Citsumall<-Citsumall %>% distinct()


#Cittemp<-select(Citsumall, c("Station", "Date_month", "Temp_surf"))

Cittemp<-select(Citsumall, c("Station", "Date_month", "rolltemp"))

#Citsumwide<-Cittemp %>% distinct() %>%
#pivot_wider(names_from = Station, values_from =Temp_surf)

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
  pivot_longer(cols="211":"215",names_to = "Station", values_to = "Temp_interp")


Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))

##Salinity Impute


Cit <- dfBay %>% 
  filter(Taxa =="Citharichthys stigmaeus") %>% filter(Station %in% c("211", "214", "215"))

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
  pivot_longer(cols="211":"215",names_to = "Station", values_to = "Sal_interp")

Time_series_inter<-Time_series_inter %>%
  pivot_longer(cols="108":"318",names_to = "Station", values_to = "Sal_interp")

Time_series_inter$Date_month<-as.Date(Time_series_inter$Date_month, format = "%Y-%m-%d")

Citsumfull<- Time_series_inter %>% 
  left_join(Citsumfull, by=c("Date_month", "Station"))

##Pull max with imputed data

Citsumfull$Julian[is.na(Citsumfull$Julian)] <- yday(Citsumfull$Date_month[is.na(Citsumfull$Julian)])

Citsumfull$Date[is.na(Citsumfull$Date)] <- Citsumfull$Date_month[is.na(Citsumfull$Date)]

Citmax <- Citsumfull %>% group_by(Year, Taxa, Station) %>% slice_max(Count_interp)

#Citmax <- Citmax %>% filter(!(Count==0))

#Citmax <- Citmax %>% filter(Station==318)

Citmax$Year<-as.Date(Citmax$Year, format="%Y")

f<-lm(Julian ~ Date, data = Citmax)
summary(f)

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