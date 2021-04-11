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
# Note that data file should be in working directory.
##########################################################################
##########################################################################
##########################################################################

#Clear out any R junk
# rm(list = ls())

#Read in Dvorak/Koba Table info
dvorak_table<-read.csv('dvorak.csv')

#1 knot = 0.514444 m/s
KnotstoMetersPerSec <- 0.514444

#Create 10 minute conversions using C1=0.88, C2=0.93
dvorak_table$MWS_tenmin_C1_dv<-NA
dvorak_table$MWS_tenmin_C1_dv<-0.88*dvorak_table$MWS_onemin_dv

dvorak_table$MWS_tenmin_C2_dv<-NA
dvorak_table$MWS_tenmin_C2_dv<-0.93*dvorak_table$MWS_onemin_dv

#Create m/s wind speed
dvorak_table$MWS_ms<-NA
dvorak_table$MWS_ms<-KnotstoMetersPerSec *dvorak_table$MWS_onemin_dv

#NLS for discussion section of Cyclostrophic Methods
nls(MWS_ms ~C*(Pref-Atlantic_Pressure_dv)^n,
	data= dvorak_table,start=list(Pref=1013,C=4,n=0.5))

#West Pacific Pressure Graph
with(dvorak_table,plot(Wpacific_1984Pressure_dv, MWS_onemin_dv,type='b',pch=18,col='Purple',lty=4,
	ylim=c(0,200),xlim=c(850,1025),xlab='West Pacific Minimum Central Pressure (mb)',ylab='Maximum Wind Speed (knots)'))
abline(v= dvorak_table$Wpacific_1989Pressure_koba[seq(1,15,2)],col='lightgrey')
with(dvorak_table,lines(Wpacific_1984Pressure_dv, 
	MWS_tenmin_C2_dv,type='b',pch=17,col='darkgreen',lty=3))
with(dvorak_table,lines(Wpacific_1984Pressure_dv, 
	MWS_tenmin_C1_dv,type='b',pch=16,col='blue',lty=2))
with(dvorak_table,lines(Wpacific_1989Pressure_koba, 
	MWS_tenmin_koba,type='b',pch=16,col='red',lty=1))
text(dvorak_table$Wpacific_1989Pressure_koba[seq(1,15,2)],200,dvorak_table$CI[seq(1,15,2)],col='Black')
text(865,200,'Koba CI Value',col='Black')
legend(850,100,lty=c(4,3,2,1),legend=c("DVT","0.88 DVT","0.93 DVT","Koba"),col=c("Purple","darkgreen","blue","red"),cex=0.9)

#Atlantic Pressure Graph
with(dvorak_table,plot(Atlantic_Pressure_dv, MWS_onemin_dv,type='b',pch=18,col='Purple',lty=4,
	ylim=c(0,200),xlim=c(850,1025),xlab='North Atlantic Minimum Central Pressure (mb)',ylab='Maximum Wind Speed (knots)'))
abline(v= dvorak_table$Atlantic_Pressure_dv[seq(1,15,2)],col='lightgrey')
with(dvorak_table,lines(Atlantic_Pressure_dv, MWS_onemin_dv,type='b',pch=18,col='Purple',lty=4,
	ylim=c(0,200),xlim=c(850,1025),xlab='North Atlantic Minimum Central Pressure (mb)',ylab='Maximum Wind Speed (knots)'))
with(dvorak_table,lines(Atlantic_Pressure_dv, 
	MWS_tenmin_C2_dv,type='b',pch=17,col='darkgreen',lty=3))
with(dvorak_table,lines(Atlantic_Pressure_dv, 
	MWS_tenmin_C1_dv,type='b',pch=16,col='blue',lty=2))
text(dvorak_table$Atlantic_Pressure_dv[seq(1,15,2)],200,dvorak_table$CI[seq(1,15,2)],col='Black')
text(870,200,'DVT CI Value',col='Black')
legend(850,100,lty=c(4,3,2,1),legend=c("DVT","0.88 DVT","0.93 DVT"),col=c("Purple","darkgreen","blue"),cex=0.9)



