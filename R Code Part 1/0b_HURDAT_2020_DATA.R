#################################################################
#################################################################
#Title: HURDAT Processing for Bayesian Hierarchical Model
#Author: Lindsey Dietz (diet0146@umn.edu)
# Updated by: Sakshi Arya (aryax010@umn.edu)
#Objective: Process HURDAT to derive data set
#Last Updated: 3/25/21
#################################################################
#################################################################

#Loads all necessary libraries
library(stringr)

#################################################################
#Useful Functions
#################################################################

ddply_TC2<-function(data){
	qualifying_condition<-expression(substr(Name_date,1,4)>1959) #& !is.na(MaximumWS) & MaximumWS >0 & !is.nan(MinimumCP) & MinimumCP >0)

	newdata<-with(data, data[which(eval(qualifying_condition)),])		
		
	ddply(newdata,~ name, summarize, StartMonth= min(as.numeric(substr(Name_date,5,6))) , maxWS= max(Max_Wind_Speed,na.rm=T), minCP= min(Min_Central_Pressure,na.rm=T), avgLat= mean(as.numeric(substr(Latitude,1,4),na.rm=T)), avgLon= mean(as.numeric(substr(Longitude,1,4),na.rm=T)))
}

#############Creating Atlantic Tropical Cyclone Data#############
##https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2019-052520.txt

Atlantic_Hurricanes_1851_2019<-read.csv("hurdat2-1851-2019-052520.txt",header=F, na.strings=c(' -999',NA))

#Providing more useful names
names(Atlantic_Hurricanes_1851_2019)<-c("Name_date","Hour","Special","Classification","Latitude","Longitude","Max_Wind_Speed","Min_Central_Pressure","NE_34kt", "SE_34kt","SW_34kt", "NW_34kt", "NE_50kt", "SE_50kt", "SW_50kt", "NW_50kt", "NE_64kt", "SE_64kt", "SW_64kt", "NW_64kt")

#Dropping wind characteristics (not needed)
atlantic_hurricane_data<-Atlantic_Hurricanes_1851_2019

#Keeping only from 1960-2019 to eliminate data issues
atlantic_data_1960_2019<-atlantic_hurricane_data[substr(atlantic_hurricane_data$Name_date,1,3) %in% c('AL0','AL1','AL2','196','197','198','199','200','201') & !substr(atlantic_hurricane_data$Name_date,5,7) %in% c('185','186','187','188','189','190','191','192','193','194','195'),]

#Names of Atlantic storms
names_of_atlantic_storms<-paste(as.character(atlantic_data_1960_2019[substr(atlantic_data_1960_2019$Name_date,1,2)=='AL',]$Hour),substr(atlantic_data_1960_2019[substr(atlantic_data_1960_2019$Name_date,1,2)=='AL',]$Name_date,3,8),'')

#Atlantic storms are defined by the prefix 'AL'
rows_of_atlantic_names<-grep('AL',substr(atlantic_data_1960_2019$Name_date,1,2))

rows_of_atlantic_names[length(rows_of_atlantic_names)+1]<-dim(atlantic_data_1960_2019)[1]
atlantic_data_1960_2019$name<-NA

for(named_row in 1:(length(rows_of_atlantic_names)-1)){
	for(row in rows_of_atlantic_names[named_row]: rows_of_atlantic_names[named_row+1]){
		atlantic_data_1960_2019[row,]$name <-names_of_atlantic_storms[named_row]
	}
}		

atlantic_data_1960_2019_na_omitted<-atlantic_data_1960_2019

#Creating several new variables in the NA omitted data set
atlantic_data_1960_2019_na_omitted$year<-as.factor(substr(atlantic_data_1960_2019_na_omitted$Name_date,1,4))
atlantic_data_1960_2019_na_omitted$category<-NA
atlantic_data_1960_2019_na_omitted$V500max <-NA
atlantic_data_1960_2019_na_omitted$V500avg <-NA

atlantic_data_1960_2019_na_omitted <-atlantic_data_1960_2019_na_omitted[!is.na(atlantic_data_1960_2019_na_omitted $Max_Wind_Speed),]

#Giant if statement to assign categories based on Saffir-Simpson Scale
for(label in 1:dim(atlantic_data_1960_2019_na_omitted)[1]){	
	if(is.na(atlantic_data_1960_2019_na_omitted[label,]$Max_Wind_Speed)){
		atlantic_data_1960_2019_na_omitted[label,]$category<-NA
	}else if(atlantic_data_1960_2019_na_omitted[label,]$Max_Wind_Speed<39){atlantic_data_1960_2019_na_omitted[label,]$category<-"TD"} else 	
	  if(atlantic_data_1960_2019_na_omitted[label,]$Max_Wind_Speed<74){atlantic_data_1960_2019_na_omitted[label,]$category<-"TS"} else
	    if(atlantic_data_1960_2019_na_omitted[label,]$Max_Wind_Speed<96){atlantic_data_1960_2019_na_omitted[label,]$category<-'1'} else
		  if(atlantic_data_1960_2019_na_omitted[label,]$Max_Wind_Speed<111){atlantic_data_1960_2019_na_omitted[label,]$category<-'2'} else
			if(atlantic_data_1960_2019_na_omitted[label,]$Max_Wind_Speed<131){atlantic_data_1960_2019_na_omitted[label,]$category<-'3'} else
			  if(atlantic_data_1960_2019_na_omitted[label,]$Max_Wind_Speed<156){atlantic_data_1960_2019_na_omitted[label,]$category<-'4'}
				else {atlantic_data_1960_2019_na_omitted[label,]$category<-'5'}

atlantic_data_1960_2019_na_omitted[label,]$V500max<-	
	with(atlantic_data_1960_2019_na_omitted[label,] , 
	max(c(NE_34kt,SE_34kt,SW_34kt,NW_34kt),na.rm=T)/9-3)

atlantic_data_1960_2019_na_omitted[label,]$V500avg<-
	with(atlantic_data_1960_2019_na_omitted[label,] , 
	mean(c(NE_34kt,SE_34kt,SW_34kt,NW_34kt),na.rm=T)/9-3)
}
	
atlantic_data_1960_2019_na_omitted$lat<-
	with(atlantic_data_1960_2019_na_omitted ,as.numeric(substr(Latitude,1,5)))

atlantic_data_1960_2019_na_omitted$V500c<-with(atlantic_data_1960_2019_na_omitted ,
	Max_Wind_Speed*((66.785-0.0912* Max_Wind_Speed+1.0619*(lat-25))/500)^(0.1147+0.0055* 
	Max_Wind_Speed-0.001*(lat-25))
	)

atlantic_data_1960_2019_na_omitted$Smax<-atlantic_data_1960_2019_na_omitted$V500max/atlantic_data_1960_2019_na_omitted$V500c	
atlantic_data_1960_2019_na_omitted$Savg<-atlantic_data_1960_2019_na_omitted$V500avg/atlantic_data_1960_2019_na_omitted$V500c
	
#Adding basin name	
atlantic_data_1960_2019_na_omitted$basin<-"Atlantic"

write.table(atlantic_data_1960_2019_na_omitted,"atlantic_hurricane_data_1960_2019_new.txt")

sum_atl<-ddply_TC2(atlantic_data_1960_2019_na_omitted)

sum_atl$Year<-0
for(i in 1:dim(sum_atl)[1]){
  sum_atl[i,]$Year<-substr(strsplit(str_trim(sum_atl$name),' ')[[i]][2],3,6)
}

write.table(sum_atl,"summarized_atl_data_1960_2019.txt")