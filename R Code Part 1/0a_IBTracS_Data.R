##########################################################################
##########################################################################
##########################################################################
# Title: IBTrACS Data Read and Processing Program 
# Author: Lindsey Dietz
# Updated by: Sakshi Arya
# Last Updated: 3/25/21
#
# Goal: Read the All Storms version of the IBTrACS
# database in netCDF format.  Extracts and modifies data
# to be in form necessary for modeling. 
#
# The data comes from the website:
# http://www.ncdc.noaa.gov/ibtracs/index.php?name=ibtracs-data
# The key for the netCDF file is located at:
# https://www.ncdc.noaa.gov/ibtracs/pdf/IBTrACS_v04_column_documentation.pdf
# Note that data file should be in working directory.
##########################################################################
##########################################################################
##########################################################################

#Clear out any R junk
rm(list=ls())

#Loads all necessary libraries
library(ncdf4)
library(fields)
library(chron)
library(plyr)
library(date)
library(geosphere)
##########################################################################
##########################################################################
##########################################################################
#Read in the IBTrACS netCDF

#Data from multiple sources-http://www.ncdc.noaa.gov/ibtracs/images/Fig_02_period_of_record.png
hurricanes= nc_open("Data/Original Data 21/IBTrACS.ALL.v04r00.nc")
print(hurricanes)

#Data only from WMO-sanctioned forecast agencies.
#Variable Definitions
#[1] "file name: Allstorms.ibtracs_all.v04r00.nc"
#[1] "file Allstorms.ibtracs_all.v04r00.nc has 12 dimensions:"
### Dimensions needed:
#[1] "storm   Size: 13476"
#[1] "time   Size: 180"
#[1] "charsn   Size: 13"
#[1] "date_time  Size: 360"
#[1] "char2  Size:2"
#[1] "char128 Size 128"
#[1] "niso Size 19"

#[1] "------------------------"
#[1] "file Allstorms.ibtracs_all.v04r00.nc has 150 variables (includes basin specific):"
#### Variables of interest:
##[1] "char name[char128,storm]  long_name: Name of system Missval:NA"
##[1] "short numobs[storm]  long_name: Number of observations per system units: 1_FillValue: -9999"
##[1] "short season[storm]  Longname:Year based on season Missval:NA"
##[1] "nature[char2,date_time,storm] long_name: Nature of the cyclone NR=Not_Reported DS=Disturbance TS=Tropical_System ET=Extratropical_System SS=Subtropical_System MX=MIXED_occurs_when_agencies_reported_inconsistent_types"
##[1] "short number[storm] long_name: Storm number (within season)"
##[1] "char basin[char2,date_time,storm]  long_name: Current basin: EP=East_Pacific NA=North_Atlantic NI=North_Indian SA=South_Atlantic SI=South_Indian SP=South_Pacific WP=Western_Pacific"
##[1] "short landfall[date_time,storm] long_name: Minimum distance to land between current location and next"
##[1] "char track_type[niso,storm]  long_name: Name of track type Missval:-127"
##[1] "char wmo_agency[niso,date_time,storm]  long_name: Official WMO agency Missval:NA"
##[1] "double time[date_time,storm] "long_name: time, units: days since 1858-11-17 00:00:00
##[1] "float lat[date_time,storm]  long_name: latitude, units: degrees_north"
##[1] "float lon[date_time,storm]  long_name: longitude, units: degrees_east"
##[1] "short wmo_wind[date_time,storm]  long_name:Maximum sustained wind speed from Official WMO agency units: kts Missval:-9999"
##[1] "short wmo_pres[date_time,storm] long_name: Minimum central preossure from  Official WMO agency units: mb #Missval:-9999"
##[1] "short dist2land[time,storm]  Longname:Distance to land at current location Missval:-9999"
##[1] "char iso_time[niso,date_time,storm] long_name: Time (ISO)"

dist2land= ncvar_get(nc= hurricanes,varid="dist2land") 
storm.number  =  ncvar_get(nc= hurricanes,varid="numobs") 
time.number= ncvar_get(nc= hurricanes,varid="time") 
name= ncvar_get(nc= hurricanes,varid="name")
nature_for_mapping= ncvar_get(nc= hurricanes,varid="nature")
lat_for_mapping = ncvar_get(nc= hurricanes,varid="lat")
lon_for_mapping = ncvar_get(nc= hurricanes,varid="lon")
source_wind= ncvar_get(nc= hurricanes,varid="wmo_wind")  
source_pres = ncvar_get(nc= hurricanes,varid="wmo_pres") 
basin<- ncvar_get(nc=hurricanes,varid="basin")
source_time<- ncvar_get(nc=hurricanes,varid="iso_time")

#According to the key, there are more possibilities
#However, the data only contains the following:
#basin[time,storm]
#0 = NA - North Atlantic
#1 = SA - South Atlantic
#2 = WP - West Pacific
#3 = EP - East Pacific
#4 = SP - South Pacific
#5 = NI - North Indian
#6 = SI - South Indian

##########################################################################
##########################################################################
##########################################################################
#Assumption 1: Basin of origin is defined as the first documented basin
# in the data set as the basin for the tropical cyclone
#Assumption 2: If multiple reports are found for minCP or maxWS, the 
# average of the reports are used

#Source time is stored as Modified Julian Day
#We subtract 36933.5 to put it into the usual time
# time_modifier<- 36933.5  Not needed in the new dataset as the date format is different

# empty space means data is missing

missing_val<- ""

storm.basin<-vector()
year<-vector()
data.storm<-list()

for(i in 1:length(storm.number)){
		if(sum(basin[,i]!= missing_val)==1){
			storm.basin[i]<-unique(basin[which(basin[,i]!= missing_val),][i])[1]
		} else{
			storm.basin[i]<-unique(basin[which(basin[,i]!= missing_val),][,i])[1]
		}
	
   
	year[i] <- na.omit(unique(as.integer(substring(source_time[,i],1,4))))[1]
}

## Some logitude values were > 180 deg, from the official google forum of IBTracs 
## these are values on the east side of the dateline. IBTrACS tries to keep the longitudes continuous for each individual storm, which results in the behavior you see. 
## The correct remedy is to subtract 360 from any values over 180. 
lon_over_180_cols <- which(apply(lon_for_mapping, 2, max, na.rm = TRUE) > 180)
for(i in lon_over_180_cols){
lon_for_mapping[which(lon_for_mapping[, i] > 180), i] <- lon_for_mapping[which(lon_for_mapping[, i] > 180), i] - 360
}

#Creating a data frame for each storm
for(storm in 1:length(storm.number)){
	
	data.storm[[storm]]<-data.frame(
	
		Name=rep(name[storm],sum(nature_for_mapping[, storm]!= missing_val)),
	
		Storm=rep(storm, sum(nature_for_mapping[, storm]!= missing_val)),
		
		Basin=rep(storm.basin[storm],sum(nature_for_mapping[, storm]!= missing_val)),
	
		Month= na.omit(as.integer(substring(source_time[,storm], 6,7))),
	
		Day=na.omit(as.integer(substring(source_time[,storm], 9,10))),
	
		Year=rep(year[storm],sum(nature_for_mapping[, storm]!= missing_val)),
	
		Classification=nature_for_mapping[nature_for_mapping[, storm]!= missing_val ,storm],
	
		Latitude=lat_for_mapping[nature_for_mapping[, storm]!= missing_val,storm],
		
		Longitude=lon_for_mapping[nature_for_mapping[, storm]!= missing_val,storm],
	
		# MaximumWS=colMeans(as.matrix(source_wind[, storm]),na.rm=T), 
		MaximumWS= source_wind[nature_for_mapping[, storm]!= missing_val,storm], 
		# MinimumCP =colMeans(as.matrix(source_pres[, storm]),na.rm=T),
		MinimumCP = source_pres[nature_for_mapping[,storm]!= missing_val, storm],
		
		
		Landfall= dist2land[nature_for_mapping[, storm]!= missing_val,storm]
		
	)
    
	if(storm==11881){
		data.storm[[storm]]$Basin<-0	
	} 
	obs.by.storm<-dim(data.storm[[storm]])[1] 
	dat<-data.storm[[storm]]
	dat$MaximumWS[is.na(dat$MaximumWS)] = mean(dat$MaximumWS, na.rm = TRUE)
	dat$MinimumCP[is.na(dat$MinimumCP)]= mean(dat$MinimumCP, na.rm = TRUE)
	Translation_Speed<-vector()
	Pressure_Time<-vector()
	if(dim(data.storm[[storm]])[1]>1){
	for(i in 1: obs.by.storm){
		if(i>1 & i<obs.by.storm){
		meters<-with(dat ,distHaversine(c(Longitude[i-1],Latitude[i-1]),c(Longitude[i+1],Latitude[i+1]) ))
		press_drop<-with(dat, MinimumCP[i+1]-MinimumCP[i-1])
		time_diff_sec<-12*60*60
		Translation_Speed[i]= meters/time_diff_sec
		Pressure_Time[i]= press_drop/12
		}
		if(i==1){
		meters<-with(dat ,distHaversine(c(Longitude[i],Latitude[i]),c(Longitude[i+1],Latitude[i+1]) ))
		press_drop<-with(dat, MinimumCP[i+1]-MinimumCP[i])
		time_diff_sec<-6*60*60
		Translation_Speed[i]= meters/time_diff_sec
		Pressure_Time[i]= press_drop/6
		}
		if( i== obs.by.storm){
		meters<-with(dat ,distHaversine(c(Longitude[i],Latitude[i]),c(Longitude[i-1],Latitude[i-1]) ))
		press_drop<-with(dat, MinimumCP[i]-MinimumCP[i-1])
		time_diff_sec<-6*60*60
		Translation_Speed[i]= meters/time_diff_sec
		Pressure_Time[i]= press_drop/6
		}
	
	}
		data.storm[[storm]]<-cbind(data.storm[[storm]],Translation_Speed, Pressure_Time)
	}else{
		data.storm[[storm]]<-cbind(data.storm[[storm]],Translation_Speed=NA, Pressure_Time=NA)
	}
}


#Combining all storms into one giant data frame
data.storms.all<-do.call("rbind", data.storm)

##########################################################################
##########################################################################
##########################################################################
#Separating the storms into their respective basins

North_Atlantic_Storms<- with(data.storms.all ,data.storms.all[which(Basin=="NA"),])
South_Atlantic_Storms<- with(data.storms.all ,data.storms.all[which(Basin=="SA"),])
West_Pacific_Storms<-with(data.storms.all , data.storms.all[which(Basin=="WP"),])
East_Pacific_Storms<- with(data.storms.all ,data.storms.all[which(Basin=="EP"),])
South_Pacific_Storms<-with(data.storms.all , data.storms.all[which(Basin=="SP"),])
North_Indian_Storms<- with(data.storms.all ,data.storms.all[which(Basin=="NI"),])
South_Indian_Storms<- with(data.storms.all ,data.storms.all[which(Basin=="SI"),])


## Rough
dim(North_Atlantic_Storms)
head(North_Atlantic_Storms)
hist(North_Atlantic_Storms$MaximumWS)

##########################################################################
##########################################################################
##########################################################################
#Keep only storms meeting Tropical Cyclone criteria from 1960-2014
qualifying_condition<-expression(Year>1959 & !is.na(MaximumWS) & MaximumWS >0 & !is.nan(MinimumCP) & MinimumCP >0)

North_Atlantic_Storms_1960<-with(North_Atlantic_Storms, North_Atlantic_Storms[which(eval(qualifying_condition)),])

South_Atlantic_Storms_1960<-with(South_Atlantic_Storms, South_Atlantic_Storms[which(eval(qualifying_condition)),])

West_Pacific_Storms_1960<-with(West_Pacific_Storms, West_Pacific_Storms[which(eval(qualifying_condition)),])

East_Pacific_Storms_1960<-with(East_Pacific_Storms, East_Pacific_Storms[which(eval(qualifying_condition)),])

South_Pacific_Storms_1960<-with(South_Pacific_Storms, South_Pacific_Storms[which(eval(qualifying_condition)),])

North_Indian_Storms_1960<-with(North_Indian_Storms, North_Indian_Storms[which(eval(qualifying_condition)),])

South_Indian_Storms_1960<-with(South_Indian_Storms, South_Indian_Storms[which(eval(qualifying_condition)),])

##########################################################################
##########################################################################
#Write out Basin Data Set files to working directory

write.table(North_Atlantic_Storms_1960,"North_Atlantic_Storms_1960_2019.txt")
write.table(South_Atlantic_Storms_1960,"South_Atlantic_Storms_1960_2019.txt")
write.table(West_Pacific_Storms_1960,"West_Pacific_Storms_1960_2019.txt")
write.table(East_Pacific_Storms_1960,"East_Pacific_Storms_1960_2019.txt")
write.table(South_Pacific_Storms_1960,"South_Pacific_Storms_1960_2019.txt")
write.table(North_Indian_Storms_1960,"North_Indian_Storms_1960_2019.txt")
write.table(South_Indian_Storms_1960,"South_Indian_Storms_1960_2019.txt")


dim(North_Atlantic_Storms_1960)
hist(North_Atlantic_Storms_1960$MaximumWS)
hist(North_Atlantic_Storms_1960$MinimumCP)
##########################################################################
##########################################################################
##########################################################################

##########################################################################
#Function to summarize data by tropical cyclone
#Assumption 1: The month associated with the TC is taken to be the first
# month listed in the data set for that storm
#Assumption 2: Latitude and Longitude are averaged to summarize

ddply_TC<-function(data){
	count_qualifying_condition<-expression(Year>1959)
	qualifying_condition<-expression(Year>1959 & !is.na(MaximumWS) & MaximumWS >0 & !is.na(MinimumCP) & MinimumCP >0)

	newdata<-with(data, data[which(eval(qualifying_condition)),])
	newdata$OnLand <- as.numeric(newdata$Landfall==0)
	newdata$Avg_Trans_Speed<-NA
	newdata$SeatoLand<-0
		for(Storm in unique(newdata$Storm)){
			temp.data<-newdata[which(newdata$Storm== Storm),]
			n<-dim(temp.data)[1]
			
			if(n>1){
				meters<-with(temp.data ,distHaversine(c(Longitude[1],Latitude[1]),c(Longitude[n],Latitude[n]) ))
				seconds<-6*(n-1)*60*60
				temp.data$Avg_Trans_Speed<-meters/seconds
			  for(i in 2:n){
				temp.data[i-1,]$SeatoLand<-temp.data[i,]$OnLand-
					temp.data[i-1,]$OnLand
			  }
			}
			newdata[which(newdata$Storm== Storm),]<-temp.data	
		}	
		
	newdata$Land_Speed <- 0
	newdata$Land_Pressure <-NA

	newdata[newdata$Landfall==0,]$Land_Speed<-
		newdata[newdata$Landfall==0 ,]$MaximumWS
	
	newdata[newdata$Landfall==0,]$Land_Pressure<-
		newdata[newdata$Landfall==0 ,]$MinimumCP	

	newdata[newdata$SeatoLand==1,]$Land_Speed<-
		newdata[ newdata$SeatoLand==1,]$MaximumWS
		
	newdata[newdata$SeatoLand==1,]$Land_Pressure<-
		newdata[ newdata$SeatoLand==1,]$MinimumCP	
		
		
	ddply(newdata,.(Storm,Year,Name), summarize, StartMonth= min(Month) , maxWS= max(MaximumWS,na.rm=T),minCP= min(MinimumCP,na.rm=T),avgLat= mean(Latitude,na.rm=T),avgLon= mean(Longitude,na.rm=T),madeLandfall=any(Landfall==0), Max_Landfall_Speed= max(Land_Speed,na.rm=T), Min_Landfall_Pressure= min(Land_Pressure,na.rm=T), Avg_Trans_Speed= mean(Avg_Trans_Speed))
}


ddply_TC2<-function(data){
	# count_qualifying_condition<-expression(Year>1959)
	 qualifying_condition<-expression(Year>1959) #& !is.na(MaximumWS) & MaximumWS >0 & !is.nan(MinimumCP) & MinimumCP >0)

	newdata<-with(data, data[which(eval(qualifying_condition)),])
	newdata$OnLand <- as.numeric(newdata$Landfall==0)
	newdata$Avg_Trans_Speed<-NA
	newdata$SeatoLand<-0
		for(Storm in unique(newdata$Storm)){
			temp.data<-newdata[which(newdata$Storm== Storm),]
			n<-dim(temp.data)[1]
			
			if(n>1){
				meters<-with(temp.data ,distHaversine(c(Longitude[1],Latitude[1]),c(Longitude[n],Latitude[n]) ))
				seconds<-6*(n-1)*60*60
				temp.data$Avg_Trans_Speed<-meters/seconds
			  for(i in 2:n){
				temp.data[i-1,]$SeatoLand<-temp.data[i,]$OnLand-
					temp.data[i-1,]$OnLand
			  }
			}
			newdata[which(newdata$Storm== Storm),]<-temp.data	
		}	
		
	newdata$Land_Speed <-0
	newdata$Land_Pressure <-NA

	newdata[newdata$Landfall==0,]$Land_Speed<-
		newdata[newdata$Landfall==0 ,]$MaximumWS
	
	newdata[newdata$Landfall==0,]$Land_Pressure<-
		newdata[newdata$Landfall==0 ,]$MinimumCP	

	newdata[newdata$SeatoLand==1,]$Land_Speed<-
		newdata[ newdata$SeatoLand==1,]$MaximumWS
		
	newdata[newdata$SeatoLand==1,]$Land_Pressure<-
		newdata[ newdata$SeatoLand==1,]$MinimumCP	
		
		
	ddply(newdata,.(Storm,Year,Name), summarize, StartMonth= min(Month) , maxWS= max(MaximumWS,na.rm=T),minCP= min(MinimumCP,na.rm=T),avgLat= mean(Latitude,na.rm=T),avgLon= mean(Longitude,na.rm=T),madeLandfall=any(Landfall==0))
}

##########################################################################
#Summarizing by tropical cyclone within Basin


summarized_North_Atlantic_Storms_1960<-ddply_TC(North_Atlantic_Storms)

sum_North_Atlantic_Storms_1960<-ddply_TC2(North_Atlantic_Storms)

summarized_South_Atlantic_Storms_1960 <-ddply_TC(South_Atlantic_Storms)

summarized_West_Pacific_Storms_1960 <-ddply_TC(West_Pacific_Storms)

summarized_East_Pacific_Storms_1960 <-ddply_TC(East_Pacific_Storms)

summarized_South_Pacific_Storms_1960 <-ddply_TC(South_Pacific_Storms)

summarized_North_Indian_Storms_1960 <-ddply_TC(North_Indian_Storms)

summarized_South_Indian_Storms_1960 <-ddply_TC(South_Indian_Storms)

##########################################################################
##########################################################################
#Write out summarized Basin Data Set files to working directory

write.table(summarized_North_Atlantic_Storms_1960,"summarized_North_Atlantic_Storms_1960_2019.txt")
write.table(sum_North_Atlantic_Storms_1960,"sum_North_Atlantic_Storms_1960_2019.txt")
write.table(summarized_South_Atlantic_Storms_1960,"summarized_South_Atlantic_Storms_1960_2019.txt")
write.table(summarized_West_Pacific_Storms_1960,"summarized_West_Pacific_Storms_1960_2019.txt")
write.table(summarized_East_Pacific_Storms_1960,"summarized_East_Pacific_Storms_1960_2019.txt")
write.table(summarized_South_Pacific_Storms_1960,"summarized_South_Pacific_Storms_1960_2019.txt")
write.table(summarized_North_Indian_Storms_1960,"summarized_North_Indian_Storms_1960_2019.txt")
write.table(summarized_South_Indian_Storms_1960,"summarized_South_Indian_Storms_1960_2019.txt")







