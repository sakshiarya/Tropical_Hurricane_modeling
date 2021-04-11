##########################################################################
##########################################################################
##########################################################################
# Title: Comparing WLS to Knaff-Zehr and Holland Models
# Author: Lindsey Dietz
# Updated by: Sakshi Arya
# Last Updated: 3/25/21
#
# Goal: Provides the code for all analyses presented within the paper
#
# Data has been created via 0_IBTracs_Data
#
#Lit. Source 1: Knaff, J. A. and Zehr, R. M. (2007). Reexamination of Tropical Cyclone 
# Wind– Pressure Relationships. Weather and Forecasting, 22(1):71–88. 
#Lit. Source 2: Holland, G. J. (2008). A Revised Hurricane Pressure–Wind Model. 
# Monthly Weather Review, 136(9):3432–3445 

##########################################################################
##########################################################################
##########################################################################

#Clear out any R junk
rm(list=ls())

#Loads all necessary libraries
library(ggplot2)
library(MASS)
library(reshape2)
library(multcomp)
library(car)
library(stringr)
library(lme4)

filep <-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure'

#################################################################
#################################################################
#####################Read in all data files######################
#################################################################
#################################################################

North_Atlantic_Storms_1960 <-read.table(file.path(filep, 
	'Data/Derived Data 21/North_Atlantic_Storms_1960_2019.txt'),header=T)

#Hurdat2 data
NA_Hurdat <-read.table('Data/Derived Data 21/atlantic_hurricane_data_1960_2019_new.txt')

#Removing data without minimum central pressure or v500
NA_HURDAT <-(NA_Hurdat[!is.na(NA_Hurdat$Min_Central_Pressure) & NA_Hurdat$V500max>0 & NA_Hurdat$Max_Wind_Speed>=64,])
NA_HURDAT$Name<-str_trim(NA_HURDAT$name)
NA_HURDAT$Year<-str_trim(NA_HURDAT$year)
NA_HURDAT$Latitude<-as.numeric(str_trim(NA_HURDAT$lat))
NA_HURDAT$MinimumCP <-str_trim(NA_HURDAT $Min_Central_Pressure)

#Data for V500 only exists from 2004 forward
NA_IBTrACs<-North_Atlantic_Storms_1960[North_Atlantic_Storms_1960$Year>=2004,]
NA_IBTrACs$Name<-str_trim(NA_IBTrACs$Name)
NA_IBTrACs$MinimumCP <-str_trim(round(NA_IBTrACs$MinimumCP,0))

atl_compare_data <-merge(NA_HURDAT, NA_IBTrACs,by=c('Name','Year','MinimumCP'))

#################################################################
#################################################################
##Creating all necessary variables following Sources 1 and 2  
#################################################################
#################################################################

#1.94384 is used to convert from m/s to knots
ms_to_knots<-1.94384 

#Source 1, Page 74
atl_compare_data$knaff_v_srm<-with(atl_compare_data,
	MaximumWS-1.5*(ms_to_knots*Translation_Speed)^0.63)

atl_compare_data$knaff_DP<-with(atl_compare_data,
	Min_Central_Pressure-1013)

#Source 1, Page 80, Equation 8
atl_compare_data$knaff_pred_vmax<-with(atl_compare_data,
	18.633 - 14.960*Savg - 0.755*lat - 0.518* knaff_DP + 
	9.738*sqrt(abs(knaff_DP))+
	1.5*(ms_to_knots*Translation_Speed)^0.63)

#Pulling out longitude
atl_compare_data$Longitude<-as.numeric(substr(atl_compare_data$Longitude.x,1,6))


#Source 2, 1st column of page 3435
atl_compare_data$holland_DP<-with(atl_compare_data,
	1015-Min_Central_Pressure)

#Source 2, Page 3435, above equation 9
atl_compare_data$p_rmw<-with(atl_compare_data, 
	Min_Central_Pressure +(holland_DP)/3.7)
	
#Source 2, Page 3435, equation 9
atl_compare_data$T_s<- with(atl_compare_data, 
	28-3*(abs(lat)-10)/20)
	
#Source 2, Page 3435, equation 9
atl_compare_data$q_m<-with(atl_compare_data, 
	(0.9*3.802/p_rmw)* exp(17.67*T_s/(243.5+T_s)))
	
#Source 2, Page 3435, equation 9
atl_compare_data$T_vs<-with(atl_compare_data, 
	(T_s+273.15)*(1+0.81*q_m))

#Gas constant is 287.058 J /(kg* K)= 287.058 newton*m/(kg* K)
#multiply p_rmw by 100 since it is in hPa and we need pascals 
#(equivalent to newtons/m^2)
#T_vs is in Kelvins (K)
#Surface Air density (rho) in kg/m^3
gas_constant<-287.058
atl_compare_data$holland_rho<- with(atl_compare_data, 
	100*p_rmw /(gas_constant * T_vs))

#Source 2, Page 3436, equation 11
atl_compare_data$holland_x<-with(atl_compare_data,
	0.6*(1-(holland_DP)/215))

#Source 2, Page 3436, equation 11
atl_compare_data$holland_bs<- with(atl_compare_data,
	-4.4e-5* holland_DP ^2 +
	0.01* holland_DP +
	0.03* Pressure_Time -
	0.014*abs(lat)+ 
	0.15* Translation_Speed ^ holland_x + 1)

atl_compare_data$holland_MWS_ms_predict<-with(atl_compare_data,
	(holland_bs/(holland_rho*exp(1))* holland_DP)^0.5)

#Source 2, Page 3433, equation 7
#v_m= ((b/(rho*e)*(press_diff in pascals))^0.5
#v_m= ((b/(rho*e)*100(press_diff in hPa))^0.5
#v_m= 10((b/(rho*e)*(press_diff in hPa))^0.5
atl_compare_data$holland_MWS_knots<-with(atl_compare_data, 
	ms_to_knots*holland_MWS_ms_predict*10)

#################################################################
#################################################################
#Fitting actual v. predicted models
#################################################################
#################################################################

#Our model
summary(lm3_ours<-lm(log(Max_Wind_Speed)~log(1013-Min_Central_Pressure)+I((log(1013-Min_Central_Pressure))^2)+ Latitude.y+ Translation_Speed, data= atl_compare_data, weights=log(1013-Min_Central_Pressure)))

atl_compare_data$pred_vmax_ours <-exp(predict(lm3_ours))

#Our mixed model
summary(lm3_ours_mm <-lmer(log(Max_Wind_Speed)~log(1013-Min_Central_Pressure)+I((log(1013-Min_Central_Pressure))^2)+ Latitude.y+ Translation_Speed+(1| Name), data= atl_compare_data,weights=log(1013-Min_Central_Pressure)))

atl_compare_data$pred_vmax_ours2 <-exp(predict(lm3_ours_mm))

#Seeing which model has the smallest squared error difference from 
#actual values of wind speed
with(atl_compare_data,sqrt(mean((Max_Wind_Speed-pred_vmax_ours)^2,na.rm=T)))
with(atl_compare_data,sqrt(mean((Max_Wind_Speed-pred_vmax_ours2)^2,na.rm=T)))
with(atl_compare_data,sqrt(mean((Max_Wind_Speed-knaff_pred_vmax)^2,na.rm=T)))
with(atl_compare_data,sqrt(mean((Max_Wind_Speed-holland_MWS_knots)^2,na.rm=T)))

#Looking at the estimated actual v. fitted lines for each model 
summary(pred.reg<-lm(Max_Wind_Speed~ pred_vmax_ours,data= atl_compare_data))
summary(pred.reg2<-lm(Max_Wind_Speed~ pred_vmax_ours2,data= atl_compare_data))
summary(mws.pred.knaff <-lm(Max_Wind_Speed ~knaff_pred_vmax ,data= atl_compare_data))
mws.coefs.knaff<-coef(mws.pred.knaff)
summary(mws.pred.holland <-lm(Max_Wind_Speed ~ holland_MWS_knots ,data= atl_compare_data))
mws.coefs.holland<-coef(mws.pred.holland)

#Testing the (0,1) hypothesis for the estimated actual v. fitted lines for each model 
linearHypothesis(pred.reg, c("(Intercept) = 0"))
linearHypothesis(pred.reg, c("pred_vmax_ours = 1"))
linearHypothesis(pred.reg, c("(Intercept) = 0","pred_vmax_ours = 1"))

linearHypothesis(pred.reg2, c("(Intercept) = 0"))
linearHypothesis(pred.reg2, c("pred_vmax_ours2 = 1"))
linearHypothesis(pred.reg2, c("(Intercept) = 0","pred_vmax_ours2 = 1"))

linearHypothesis(mws.pred.knaff, c("(Intercept) = 0"))
linearHypothesis(mws.pred.knaff, c("knaff_pred_vmax = 1"))
linearHypothesis(mws.pred.knaff, c("(Intercept) = 0","knaff_pred_vmax = 1"))

linearHypothesis(mws.pred.holland, c("(Intercept) = 0"))
linearHypothesis(mws.pred.holland, c("holland_MWS_knots = 1"))
linearHypothesis(mws.pred.holland, c("(Intercept) = 0","holland_MWS_knots = 1"))

#Plotting actual v. fitted regression results compared to perfect fit
plot1 <-ggplot(data= atl_compare_data, aes(knaff_pred_vmax,Max_Wind_Speed))+
	geom_point(size=0.75)+ 
	stat_function(fun=function(x){x},color='grey',lwd=1.5)+
	stat_function(fun=function(x){mws.coefs.knaff[1]+mws.coefs.knaff[2]*x}, color='blue',lwd=2)+
	xlim(50,150)+ylim(50,150)+
	xlab('Knaff/Zehr Predicted Max. Wind Speed')+
	ylab('Actual Maximum Wind Speed')+
	theme_bw()+
	theme(panel.grid.major = element_blank(), text = element_text(size=20), 
		panel.grid.minor = element_blank())
plot1

plot2 <-ggplot(data= atl_compare_data, aes(holland_MWS_knots,Max_Wind_Speed))+
	geom_point(size=0.75)+
	stat_function(fun=function(x){x},color='grey',lwd=1.5)+
	stat_function(fun=function(x){mws.coefs.holland[1]+ mws.coefs.holland[2]*x}, 
		color='darkgreen', lwd=2,lty=2)+
	xlim(50,150)+ylim(50,150)+
	ylab('Actual Maximum Wind Speed')+
	theme_bw()+xlab('Holland Predicted Max. Wind Speed')+
	theme(panel.grid.major = element_blank(), text = element_text(size=20), 
		panel.grid.minor = element_blank())
plot2

plot3<-ggplot(data= atl_compare_data, aes(pred_vmax_ours ,Max_Wind_Speed))+
	geom_point(size=0.75)+
	stat_function(fun=function(x){x},color='grey',lwd=1.5)+
	stat_function(fun=function(x){coef(pred.reg)[2]*x+ coef(pred.reg)[1]},color='red',lwd=2,lty=5)+
	xlim(50,150)+ylim(50,150)+
	xlab('WLS Predicted Max. Wind Speed')+
	ylab('Actual Maximum Wind Speed')+
	theme_bw()+
	theme(panel.grid.major = element_blank(), text = element_text(size=20), 
		panel.grid.minor = element_blank())
plot3

plot4<-ggplot(data= atl_compare_data, aes(pred_vmax_ours2 ,Max_Wind_Speed))+
	geom_point(size=0.75)+
	stat_function(fun=function(x){x},color='grey',lwd=1.5)+
	stat_function(fun=function(x){coef(pred.reg)[2]*x+ coef(pred.reg)[1]},
		color='orange',lwd=2,lty=6)+
	xlim(50,150)+ylim(50,150)+
	xlab('WLS MM Predicted Max. Wind Speed')+
	ylab('Actual Maximum Wind Speed')+
	theme_bw()+
	theme(panel.grid.major = element_blank(), text = element_text(size=20), 
		panel.grid.minor = element_blank())
plot4