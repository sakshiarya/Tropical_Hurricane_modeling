##########################################################################
##########################################################################
##########################################################################
# Title: Linear and Weighted Linear Models
# Author: Lindsey Dietz
# Updated by: Sakshi Arya
# Last Updated: 3/25/21
#
# Goal: Provides the code for all analyses presented within the paper
#
# All data has been created via 0_IBTracs_Data
##########################################################################
##########################################################################
##########################################################################

#Clear out any R junk
rm(list=ls())

#Loads all necessary libraries
library(ggplot2)
library(MASS)
library(GGally)
library(reshape2)
library(multcomp)
library(DAAG)
library(car)
library(xtable)
library(mvnfast)
library(msm)
library(stringr)
library(gridExtra)
library(Matrix)
library(inline)
library(dplyr)
library(ggh4x)

filep <-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure'

#Loads all outside functions
source(file.path(filep,'Helper Functions/delta.method.wsp.r'))
source(file.path(filep,'Helper Functions/para.boot.lm.r'))
source(file.path(filep,'Helper Functions/boot.pred.r'))
source(file.path(filep,'Helper Functions/WLS_boot.r'))
source(file.path(filep,'Helper Functions/Basin_estimates.r'))

#Colors function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
#################################################################
#################################################################
#####################Read in all data files######################
#################################################################
#################################################################
NA_sum<-read.table(file.path(filep, 
	'Data/Derived Data 21/summarized_North_Atlantic_Storms_1960_2019.txt'),header=T)
SA_sum<-read.table(file.path(filep, 
	'Data/Derived Data 21/summarized_South_Atlantic_Storms_1960_2019.txt'),header=T)
WP_sum<-read.table(file.path(filep,
	'Data/Derived Data 21/summarized_West_Pacific_Storms_1960_2019.txt'),header=T)
EP_sum<-read.table(file.path(filep,
	'Data/Derived Data 21/summarized_East_Pacific_Storms_1960_2019.txt'),header=T)
SP_sum<-read.table(file.path(filep,
	'Data/Derived Data 21/summarized_South_Pacific_Storms_1960_2019.txt'),header=T)
NI_sum<-read.table(file.path(filep,
	'Data/Derived Data 21/summarized_North_Indian_Storms_1960_2019.txt'),header=T)
SI_sum<-read.table(file.path(filep,
	'Data/Derived Data 21/summarized_South_Indian_Storms_1960_2019.txt'),header=T)

SwI_sum <-SI_sum[SI_sum$avgLon<=100,]
SeI_sum<-SI_sum[SI_sum$avgLon>100,]



#Only keep data that reaches hurricane strength= 64 knots
all_usable_data<-data.frame(rbind(
	cbind(NA_sum[NA_sum$maxWS>=64,],Basin='NoAt'), 
	cbind(WP_sum[WP_sum$maxWS>=64,],Basin='WePa'),
	cbind(EP_sum[EP_sum$maxWS>=64,],Basin='EaPa'), 
	cbind(SP_sum[SP_sum$maxWS>=64,],Basin='SoPa'),
	cbind(NI_sum[NI_sum$maxWS>=64,],Basin='NoIn'), 
	cbind(SeI_sum[SeI_sum$maxWS>=64,],Basin='SoEIn'), 
	cbind(SwI_sum[SwI_sum$maxWS>=64,],Basin='SoWIn')))



#1 knot = 0.514444 m/s
KnotstoMetersPerSec <- 0.514444
modeling_data<-with(all_usable_data, 
	na.omit(data.frame(maxWS_knots= maxWS, minCP, Basin, avgLat, Avg_Trans_Speed, 
	maxWS_ms= maxWS*KnotstoMetersPerSec)))

# If already saved
# se_graph<-read.csv('AllSEGraphing.csv')

#################################################################
#################################################################
#Testing different Ambient Pressures in Cyclostrophic
#################################################################
#################################################################
ordered_basins= c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn')
modeling_data$Basin <- as.factor(modeling_data$Basin)
modeling_data$Basin <- factor(modeling_data$Basin,levels(modeling_data$Basin)[c(2,7,1,5,3,4,6)])

#1010 mb
summary(all_basins_1010lm<-lm(log(maxWS_knots)~log(1010-minCP)*Basin, data= modeling_data,weights=log(1010-minCP)))
summary(all_basins_1010lm_ms<-lm(log(maxWS_ms)~log(1010-minCP)*Basin, data= modeling_data,weights=log(1010-minCP)))
basin_est_1010<- Basin_estimates(ordered_basins= ordered_basins, 
	fit1= all_basins_1010lm,fit2= all_basins_1010lm_ms, Model_Type=1, exp.int=TRUE)
mb1010<- data.frame(Pref= 1010, Basin = ordered_basins, C_kt = round(exp(basin_est_1010[[1]][,c(5)]),3), SE_kt = round(basin_est_1010[[1]][,c(6)],3), C_ms= round(exp(basin_est_1010[[3]][,c(5)]),3), SE_ms= round(basin_est_1010[[3]][,c(6)],3), n = round(basin_est_1010[[2]][,c(5,6)],3))


#1013 mb
summary(all_basins_1013lm<-lm(log(maxWS_knots)~log(1013-minCP)*Basin, data= modeling_data,weights=log(1013-minCP)))
summary(all_basins_1013lm_ms<-lm(log(maxWS_ms)~log(1013-minCP)*Basin, data= modeling_data,weights=log(1013-minCP)))
basin_est_1013<- Basin_estimates(ordered_basins= ordered_basins, 
	fit1= all_basins_1013lm,fit2= all_basins_1013lm_ms, Model_Type=1, exp.int=TRUE)
mb1013<- data.frame(Pref= 1013, Basin = ordered_basins, C_kt = round(exp(basin_est_1013[[1]][,c(5)]),3), SE_kt = round(basin_est_1013[[1]][,c(6)],3), C_ms= round(exp(basin_est_1013[[3]][,c(5)]),3), SE_ms= round(basin_est_1013[[3]][,c(6)],3), n = round(basin_est_1013[[2]][,c(5,6)],3))

#1015 mb
summary(all_basins_1015lm<-lm(log(maxWS_knots)~log(1015-minCP)*Basin, data= modeling_data,weights=log(1015-minCP)))
summary(all_basins_1015lm_ms<-lm(log(maxWS_ms)~log(1015-minCP)*Basin, data= modeling_data,weights=log(1015-minCP)))
basin_est_1015<- Basin_estimates(ordered_basins= ordered_basins, 
	fit1= all_basins_1015lm,fit2= all_basins_1015lm_ms, Model_Type=1, exp.int=TRUE)
mb1015<- data.frame(Pref= 1015, Basin = ordered_basins, C_kt = round(exp(basin_est_1015[[1]][,c(5)]),3), SE_kt = round(basin_est_1015[[1]][,c(6)],3), C_ms= round(exp(basin_est_1015[[3]][,c(5)]),3), SE_ms= round(basin_est_1015[[3]][,c(6)],3), n = round(basin_est_1015[[2]][,c(5,6)],3))


#1017 mb
summary(all_basins_1017lm<-lm(log(maxWS_knots)~log(1017-minCP)*Basin, data= modeling_data,weights=log(1017-minCP)))
summary(all_basins_1017lm_ms<-lm(log(maxWS_ms)~log(1017-minCP)*Basin, data= modeling_data,weights=log(1017-minCP)))
basin_est_1017<- Basin_estimates(ordered_basins= ordered_basins, 
	fit1= all_basins_1017lm,fit2= all_basins_1017lm_ms, Model_Type=1, exp.int=TRUE)
mb1017<- data.frame(Pref= 1017, Basin = ordered_basins, C_kt = round(exp(basin_est_1017[[1]][,c(5)]),3), SE_kt = round(basin_est_1017[[1]][,c(6)],3), C_ms= round(exp(basin_est_1017[[3]][,c(5)]),3), SE_ms= round(basin_est_1017[[3]][,c(6)],3), n = round(basin_est_1017[[2]][,c(5,6)],3))

allPressTest<- data.frame(rbind(mb1010,mb1013,mb1015,mb1017))
print(xtable(x= cbind(allPressTest[,c(1:2)], 
		paste(format(allPressTest[,c(3)],nsmall=3),
		' (',format(allPressTest[,c(4)],nsmall=3),')',sep=''),
		paste(format(allPressTest[,c(5)],nsmall=3),
		' (',format(allPressTest[,c(6)],nsmall=3),')',sep=''),
		paste(format(allPressTest[,c(7)],nsmall=3),
		' (',format(allPressTest[,c(8)],nsmall=3),')',sep=''))),
		include.rownames =F)
		
#simultaneous tests of parameters
K<-rbind(
#NoAt B_0
 c(1,0,0,0,0,0,0,
   0,0,0,0,0,0,0),
#NoAt B_1
 c(0,1,0,0,0,0,0,
   0,0,0,0,0,0,0),
#WePa B_0
 c(1,0,1,0,0,0,0,
   0,0,0,0,0,0,0),
#WePa B_1
 c(0,1,0,0,0,0,0,
   0,1,0,0,0,0,0),
#EaPa B_0
 c(1,0,0,1,0,0,0,
   0,0,0,0,0,0,0),
#EaPa B_1
 c(0,1,0,0,0,0,0,
   0,0,1,0,0,0,0),
#NoIn B_0
 c(1,0,0,0,1,0,0,
   0,0,0,0,0,0,0),
#NoIn B_1
 c(0,1,0,0,0,0,0,
   0,0,0,1,0,0,0),
#SoPa B_0
 c(1,0,0,0,0,1,0,
   0,0,0,0,0,0,0),
#SoPa B_1
 c(0,1,0,0,0,0,0,
   0,0,0,0,1,0,0),
#SoEIn B_0
 c(1,0,0,0,0,0,1,
   0,0,0,0,0,0,0),
#SoEIn B_1
 c(0,1,0,0,0,0,0,
   0,0,0,0,0,1,0),
#SoWIn B_0
 c(1,0,0,0,0,0,0,
   1,0,0,0,0,0,0),
#SoWIn B_1
 c(0,1,0,0,0,0,0,
   0,0,0,0,0,0,1)
)

rownames(K)<-c('NoAt B_0','NoAt B_1' ,
'WePa B_0','WePa B_1',
'EaPa B_0','EaPa B_1',
'NoIn B_0','NoIn B_1',
'SoPa B_0','SoPa B_1',
'SoEIn B_0','SoEIn B_1', 
'SoWIn B_0','SoWIn B_1') 
testing_1010<- glht(all_basins_1010lm, linfct = K)
exp(confint(testing_1010)$confint[seq(1,13,2),])	
confint(testing_1010)$confint[seq(2,14,2),]	

testing_1010ms<- glht(all_basins_1010lm_ms, linfct = K)
exp(confint(testing_1010ms)$confint[seq(1,13,2),])	
confint(testing_1010ms)$confint[seq(2,14,2),]	

testing_1013<- glht(all_basins_1013lm, linfct = K)
exp(confint(testing_1013)$confint[seq(1,13,2),])	
confint(testing_1013)$confint[seq(2,14,2),]	

testing_1013ms<- glht(all_basins_1013lm_ms, linfct = K)
exp(confint(testing_1013ms)$confint[seq(1,13,2),])	
confint(testing_1013ms)$confint[seq(2,14,2),]	

testing_1015<- glht(all_basins_1015lm, linfct = K)
exp(confint(testing_1015)$confint[seq(1,13,2),])	
confint(testing_1015)$confint[seq(2,14,2),]	

testing_1015ms<- glht(all_basins_1015lm_ms, linfct = K)
exp(confint(testing_1015ms)$confint[seq(1,13,2),])	
confint(testing_1015ms)$confint[seq(2,14,2),]	

testing_1017<- glht(all_basins_1017lm, linfct = K)
exp(confint(testing_1017)$confint[seq(1,13,2),])	
confint(testing_1017)$confint[seq(2,14,2),]	

testing_1017ms<- glht(all_basins_1017lm_ms, linfct = K)
exp(confint(testing_1017ms)$confint[seq(1,13,2),])	
confint(testing_1017ms)$confint[seq(2,14,2),]	
		
#################################################################
#################################################################
########################Linear Models############################
#################################################################
#################################################################

#Fitting LM1
summary(lm1<-lm(log(maxWS_knots)~log(1013-minCP)*Basin, 
	data= modeling_data))
par(mfrow=c(2,2))
plot(lm1)

#Fitting LM2
summary(lm2<-lm(log(maxWS_knots)~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin, 
	data= modeling_data))
par(mfrow=c(2,2))
plot(lm2)

#Fitting LM3
summary(lm3<-lm(log(maxWS_knots)~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin+
	avgLat*Basin+ Avg_Trans_Speed*Basin, 
	data= modeling_data))
par(mfrow=c(2,2))
plot(lm3)

#################################################################
#################################################################
######################Weighted Linear Models#####################
#################################################################
#################################################################
modeling_data$Basin<-relevel(relevel(relevel(relevel(modeling_data$Basin,'NoIn'),'EaPa'),'WePa'),'NoAt')

modeling_data$Basinlong<-modeling_data$Basin 
levels(modeling_data$Basinlong)<-c("North Atlantic",  "West Pacific",  "East Pacific",
  "North Indian", "South Pacific",  "South East Indian" ,"South West Indian")
#################################################################
#Fitting WLM1
#################################################################

#Fitting the model
summary(wlm1<-lm(log(maxWS_knots)~log(1013-minCP)*Basin, 
	data= modeling_data, weights=log(1013-minCP)))

modeling_data$predwlsm1<-predict(wlm1)

#Residuals v. Fitted Plot WLM1
ggplot(data= modeling_data ,aes(predict(wlm1),resid(wlm1,type='pearson'),color= Basinlong))+
	geom_point(size=.9)+
	geom_abline(intercept=0,slope=0,lty=2,color='grey')+
	facet_wrap(~ Basinlong,nrow=2)+
	xlab('Fitted Values')+
	ylab("Pearson Residuals")+ 
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")+ 
	scale_x_continuous(breaks=c(4,4.5,5))


#Q-Q Plot WLM1
ggplot(data= modeling_data ,aes(sample=scale(resid(wlm1,type='pearson')),color= Basinlong))+ 
	geom_abline(intercept=quantile(scale(resid(wlm1,type='pearson')),0.25)-
		qnorm(0.25)*diff(quantile(scale(resid(wlm1,type='pearson')),
			c(0.25,0.75)))/diff(qnorm(c(0.25, 0.75))), 
			slope=diff(quantile(scale(resid(wlm1,type='pearson')),
			c(0.25,0.75)))/diff(qnorm(c(0.25, 0.75))), lty=2,color='black')+
	geom_point(stat='qq',size=.9)+
	facet_wrap(~ Basinlong,nrow=2)+
	scale_x_continuous(limits=c(-4.1,4.1))+
	scale_y_continuous(limits=c(-4.1,4.1))+
	xlab('Standardized Pearson Residual Quantiles')+
	ylab('Theoretical Quantiles')+ 
	theme_bw()+ coord_fixed()+ coord_flip()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")


#Log Scale Data and Predictions from WLS1
ggplot(modeling_data,aes(log(1013-minCP),log(maxWS_knots),color= Basinlong))+ 
	geom_point(size=.9)+ 
	geom_line(aes(log(1013-minCP),predwlsm1),color='black',alpha=.5)+ 
	facet_wrap(~ Basinlong,nrow=2)+  
	xlab('Minimum Central Pressure (mb)')+ 
	ylab('Maximum Wind Speed (knots)')+ 
	theme_bw()+ 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")

#Original Scale Data and Predictions from WLS1
ggplot(modeling_data,aes(minCP,maxWS_knots,color= Basinlong))+
	geom_point(size=.9)+
	geom_line(aes(minCP,exp(predwlsm1)),color='black',alpha=.5)+
	facet_wrap(~ Basinlong,nrow=2)+
	xlab('Minimum Central Pressure (mb)')+ 
	ylab('Maximum Wind Speed (knots)')+ 
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")
#################################################################
#Fitting WLM2
#################################################################

#Fitting the model
summary(wlm2<-lm(log(maxWS_knots)~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin, 
	data= modeling_data, weights=log(1013-minCP)))

modeling_data$predwlsm2<-predict(wlm2)

#Residuals v. Fitted Plot WLM2
ggplot(data= modeling_data ,aes(predict(wlm2),resid(wlm2,type='pearson'),color= Basinlong))+
	geom_point(size=.9)+
	geom_abline(intercept=0,slope=0,lty=2,color='grey')+
	facet_wrap(~ Basinlong,nrow=2)+
	xlab('Fitted Values')+
	ylab("Pearson Residuals")+ 
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")+ 
	scale_x_continuous(breaks=c(4,4.5,5))


#Q-Q Plot WLM2	
ggplot(data= modeling_data ,aes(sample=scale(resid(wlm2,type='pearson')),color= Basinlong))+ 
	geom_abline(intercept=quantile(scale(resid(wlm2,type='pearson')),0.25)-
		qnorm(0.25)*diff(quantile(scale(resid(wlm2,type='pearson')),
			c(0.25,0.75)))/diff(qnorm(c(0.25, 0.75))), 
			slope=diff(quantile(scale(resid(wlm2,type='pearson')),
			c(0.25,0.75)))/diff(qnorm(c(0.25, 0.75))), lty=2,color='black')+
	geom_point(stat='qq',size=.9)+
	facet_wrap(~ Basinlong,nrow=2)+
	scale_x_continuous(limits=c(-4.1,4.1))+
	scale_y_continuous(limits=c(-4.1,4.1))+
	xlab('Standardized Pearson Residual Quantiles')+
	ylab('Theoretical Quantiles')+ 
	theme_bw()+ coord_fixed()+ coord_flip()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")

#Log Scale Data and Predictions from WLS2
ggplot(modeling_data,aes(log(1013-minCP),log(maxWS_knots),color= Basin))+
	geom_point(size=.9)+
	geom_line(aes(log(1013-minCP),(predwlsm2)),color='black',alpha=.5)+ 
	facet_wrap(~Basin,nrow=2)+
	xlab('Minimum Central Pressure (mb)')+
	ylab('Maximum Wind Speed (knots)')+ 
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")

#Original Scale Data and Predictions from WLS2
ggplot(modeling_data,aes(minCP,maxWS_knots,color= Basin))+
	geom_point(size=.9)+
	geom_line(aes(minCP,exp(predwlsm2)),color='black',alpha=.5)+ 
	facet_wrap(~Basin,nrow=2)+xlab('Minimum Central Pressure (mb)')+
	ylab('Maximum Wind Speed (knots)')+ 
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")
#################################################################
#Fitting WLM3
#################################################################

#Fitting the model
summary(wlm3<-lm(log(maxWS_knots)~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin+
	avgLat*Basin+ Avg_Trans_Speed*Basin, 
	data= modeling_data, weights=log(1013-minCP)))

modeling_data$predwlsm3<-predict(wlm3)

#Residuals v. Fitted Plot WLM3
ggplot(data= modeling_data ,aes(predict(wlm3),resid(wlm3,type='pearson'),color= Basinlong))+
	geom_point(size=.9)+
	geom_abline(intercept=0,slope=0,lty=2,color='grey')+
	facet_wrap(~ Basinlong,nrow=2)+
	xlab('Fitted Values')+
	ylab("Pearson Residuals")+ 
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")+ 
	scale_x_continuous(breaks=c(4,4.5,5))

#Q-Q Plot WLM3
ggplot(data= modeling_data ,aes(sample=scale(resid(wlm3,type='pearson')),color= Basinlong))+ 
	geom_abline(intercept=quantile(scale(resid(wlm3,type='pearson')),0.25)-
		qnorm(0.25)*diff(quantile(scale(resid(wlm3,type='pearson')),
			c(0.25,0.75)))/diff(qnorm(c(0.25, 0.75))), 
			slope=diff(quantile(scale(resid(wlm3,type='pearson')),
			c(0.25,0.75)))/diff(qnorm(c(0.25, 0.75))), lty=2,color='black')+
	geom_point(stat='qq',size=.9)+
	facet_wrap(~ Basinlong,nrow=2)+
	scale_x_continuous(limits=c(-4.1,4.1))+
	scale_y_continuous(limits=c(-4.1,4.1))+
	xlab('Standardized Pearson Residual Quantiles')+
	ylab('Theoretical Quantiles')+ 
	theme_bw()+ coord_fixed()+ coord_flip()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")

#Log Scale Actual vs. Predictions from WLS3
ggplot(modeling_data,aes(predwlsm3,log(maxWS_knots),color= Basinlong))+
	geom_point(size=.9)+
	facet_wrap(~ Basinlong,nrow=2)+geom_abline(intercept = 0, slope = 1,color='grey')+
	xlab('Predicted log(Maximum Wind Speed)')+
	ylab('Actual log(Maximum Wind Speed)')+
	theme_bw()+
	theme(legend.position="none")+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")

#Original Scale Actual vs. Predictions from WLS3
options(repr.plot.width = 4, repr.plot.height = 2)
ggplot(modeling_data,aes(exp(predwlsm3),maxWS_knots,color= Basinlong))+
	geom_point(size=.9)+
	facet_wrap(~ Basinlong,nrow=2)+
	geom_abline(intercept = 0, slope = 1,color='grey')+
	xlab('Predicted Maximum Wind Speed (knots)')+
	ylab('Actual Maximum Wind Speed (knots)')+
	theme_bw()+
	theme(legend.position="none")+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.background = element_blank(), 
		axis.line = element_line(colour = "black"),legend.position="none")

#################################################################
#WLS1/WLS2/WLS3 Testing
#################################################################

#LRT for WLS1 v. WLS2
anova(wlm1,wlm2)
#LRT for WLS2 v. WLS3
anova(wlm2,wlm3)
#LRT for WLS1 v. WLS3
anova(wlm1,wlm3)

#AIC
AIC(wlm1)
AIC(wlm2)
AIC(wlm3)

#BIC
BIC(wlm1)
BIC(wlm2)
BIC(wlm3)

#Shapiro-Wilk Test
shapiro.test((resid(wlm1)))
shapiro.test((resid(wlm2)))
shapiro.test((resid(wlm3)))

#Asymptotic Standard Errors for WLSs
delta.method.wsp(model.object=wlm1,model.type='lm',
	model.num=1,num.par=c(13,20,34)[1])
delta.method.wsp(model.object=wlm2,model.type='lm',
	model.num=2,num.par=c(13,20,34)[2])
delta.method.wsp(model.object=wlm3,model.type='lm',
	model.num=3,num.par=c(13,20,34)[3])

#################################################################
#################################################################
#Bootstrapping SE calculations and plotting
#################################################################
#################################################################
file.boot<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/Data/Derived Data 21/Model Output Data'

boots<- 1000
model.type1=c('1','2','3')[1]
model.type2=c('1','2','3')[2]
model.type3=c('1','2','3')[3]
m.percent1 = c(0.6666, 0.75, 1.00)[1]
m.percent2 = c(0.6666, 0.75, 1.00)[2]
m.percent3 = c(0.6666, 0.75, 1.00)[3]
#################################################################
#Bootstrapping WLS 1
#1000 bootstraps then outputs to 3 csv files
#################################################################

#m/n = 0.6666, MooN parametric subsampling

if(file.exists(file.path(file.boot,'wlm1.66.boot.all.csv'))){
wlm1.66.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'wlm1.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type1, ordered_basins= ordered_basins)
} else {
	set.seed(1313123)
	pbwlm1<-para.boot.lm(boot.num= boots,model.object=wlm1,
		model.type= model.type1 ,m.percent= m.percent1)
	write.csv(rbind(pbwlm1$boot.est.basins,pbwlm1$boot.se.basins),
		file.path(file.boot,'wlm1.66.boot.basin.csv'))
	write.csv(rbind(pbwlm1$boot.est.orig,pbwlm1$boot.se.orig),
		file.path(file.boot,'wlm1.66.boot.orig.csv'))
	temp1<-data.frame(pbwlm1$all.coeff.est)
	names(temp1)<-names(coef(wlm1))
	write.csv(temp1,file.path(file.boot,'wlm1.66.boot.all.csv'))
}

#m/n = 0.6666, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm1.66.boot.all.csv'))){
repwlm1.66.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repwlm1.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type1, ordered_basins= ordered_basins)
} else {
	set.seed(1313123)
	pbwlm1<-para.boot.lm(boot.num=boots,model.object=wlm1,
		model.type= model.type1,m.percent= m.percent1,rep=TRUE)
	write.csv(rbind(pbwlm1$boot.est.basins,pbwlm1$boot.se.basins),
		file.path(file.boot,'repwlm1.66.boot.basin.csv'))
	write.csv(rbind(pbwlm1$boot.est.orig,pbwlm1$boot.se.orig),
		file.path(file.boot,'reprewlm1.66.boot.orig.csv'))
	temp1<-data.frame(pbwlm1$all.coeff.est)
	names(temp1)<-names(coef(wlm1))
	write.csv(temp1,	file.path(file.boot,'repwlm1.66.boot.all.csv'))
}

#m/n = 0.75, MooN parametric subsampling

if(file.exists(file.path(file.boot,'wlm1.75.boot.all.csv'))){
wlm1.75.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'wlm1.75.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent= m.percent2, model.type= model.type1, ordered_basins= ordered_basins)
} else {
  set.seed(1391039)
  pbwlm1<-para.boot.lm(boot.num=boots,model.object=wlm1,
  	model.type= model.type1,m.percent= m.percent2)
  write.csv(rbind(pbwlm1$boot.est.basins,pbwlm1$boot.se.basins),
  	file.path(file.boot,'wlm1.75.boot.basin.csv'))
  write.csv(rbind(pbwlm1$boot.est.orig,pbwlm1$boot.se.orig),
  	file.path(file.boot,'wlm1.75.boot.orig.csv'))
  temp1<-data.frame(pbwlm1$all.coeff.est)
  names(temp1)<-names(coef(wlm1))
  write.csv(temp1, file.path(file.boot,'wlm1.75.boot.all.csv'))
}


#m/n = 0.75, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm1.75.boot.all.csv'))){
repwlm1.75.boot<- WLS_Boot(data= modeling_data, 
  	data.file=read.csv(file.path(file.boot,'repwlm1.75.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
  	m.percent= m.percent2,model.type= model.type1, ordered_basins= ordered_basins)
} else {
  set.seed(1391039)
  pbwlm1<-para.boot.lm(boot.num= boots,model.object=wlm1,
  	model.type= model.type1 ,m.percent= m.percent2,rep=TRUE)
  write.csv(rbind(pbwlm1$boot.est.basins,pbwlm1$boot.se.basins),
  	file.path(file.boot,'repwlm1.75.boot.basin.csv'))
  write.csv(rbind(pbwlm1$boot.est.orig,pbwlm1$boot.se.orig),
  	file.path(file.boot,'repwlm1.75.boot.orig.csv'))
  temp1<-data.frame(pbwlm1$all.coeff.est)
  names(temp1)<-names(coef(wlm1))
  write.csv(temp1, file.path(file.boot,'repwlm1.75.boot.all.csv'))
}

#m/n = 1.00, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm1.boot.all.csv'))){
repwlm1.boot<-  WLS_Boot(data= modeling_data, 
  	data.file=read.csv(file.path(file.boot,'repwlm1.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
  	m.percent = m.percent3, model.type= model.type1, ordered_basins= ordered_basins)
} else {
  set.seed(131313)
  pbwlm1.wrep<-para.boot.lm(boot.num=boots,model.object=wlm1,
  	model.type= model.type1, m.percent = m.percent3,rep=TRUE)
  write.csv(rbind(pbwlm1.wrep$boot.est.basins,pbwlm1.wrep$boot.se.basins),
  	file.path(file.boot,'repwlm1.boot.basin.csv'))
  write.csv(rbind(pbwlm1.wrep$boot.est.orig, pbwlm1.wrep$boot.se.orig),
  	file.path(file.boot,'repwlm1.boot.orig.csv'))
  temp1.rep <-data.frame(pbwlm1.wrep$all.coeff.est)
  names(temp1.rep)<-names(coef(wlm1))
  write.csv(temp1.rep,file.path(file.boot,'repwlm1.boot.all.csv'))
}

#################################################################
#Bootstrapping WLS 2
#1000 bootstraps then outputs to 3 csv files
#################################################################

#m/n = 0.6666, MooN parametric subsampling
if(file.exists(file.path(file.boot,'wlm2.66.boot.all.csv'))){
wlm2.66.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'wlm2.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type2, ordered_basins= ordered_basins)
} else {
	set.seed(131314323)
	pbwlm2<-para.boot.lm(boot.num= boots,model.object=wlm2,
		model.type= model.type2 ,m.percent= m.percent1)
	write.csv(rbind(pbwlm2$boot.est.basins,pbwlm2$boot.se.basins),
		file.path(file.boot,'wlm2.66.boot.basin.csv'))
	write.csv(rbind(pbwlm2$boot.est.orig,pbwlm2$boot.se.orig),
		file.path(file.boot,'wlm2.66.boot.orig.csv'))
	temp1<-data.frame(pbwlm2$all.coeff.est)
	names(temp1)<-names(coef(wlm2))
	write.csv(temp1,file.path(file.boot,'wlm2.66.boot.all.csv'))
}

#m/n = 0.6666, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm2.66.boot.all.csv'))){
repwlm2.66.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repwlm2.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type2, ordered_basins= ordered_basins)
} else {
	set.seed(131314323)
	pbwlm2rep<-para.boot.lm(boot.num=boots,model.object=wlm2,
		model.type= model.type2,m.percent= m.percent1,rep=TRUE)
	write.csv(rbind(pbwlm2rep$boot.est.basins, pbwlm2rep$boot.se.basins),
		file.path(file.boot,'repwlm2.66.boot.basin.csv'))
	write.csv(rbind(pbwlm2 $boot.est.orig, pbwlm2 $boot.se.orig),
		file.path(file.boot,'repwlm2.66.boot.orig.csv'))
	temp2rep <-data.frame(pbwlm2rep$all.coeff.est)
	names(temp2rep)<-names(coef(wlm2))
	write.csv(temp2rep,file.path(file.boot,'repwlm2.66.boot.all.csv'))
}


#m/n = 0.75, MooN parametric subsampling
if(file.exists(file.path(file.boot,'wlm2.75.boot.all.csv'))){
wlm2.75.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'wlm2.75.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent2, model.type= model.type2, ordered_basins= ordered_basins)
}else {
	set.seed(1391354039)
	pbwlm2 <-para.boot.lm(boot.num=boots ,model.object= wlm2,
		model.type= model.type2,m.percent= m.percent2)
	write.csv(rbind(pbwlm2$boot.est.basins, pbwlm2$boot.se.basins),
		file.path(file.boot,'wlm2.75.boot.basin.csv'))
	write.csv(rbind(pbwlm2$boot.est.orig, pbwlm2$boot.se.orig),
	    file.path(file.boot,'wlm2.75.boot.orig.csv'))
	temp2 <-data.frame(pbwlm2$all.coeff.est)
	names(temp2)<-names(coef(wlm2))
	write.csv(temp2,file.path(file.boot,'wlm2.75.boot.all.csv'))
}

#m/n = 0.75, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm2.75.boot.all.csv'))){
repwlm2.75.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repwlm2.75.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent2, model.type= model.type2, ordered_basins= ordered_basins)
}else {
	set.seed(1391354039)
	pbwlm2rep<-para.boot.lm(boot.num=boots,model.object=wlm2,
		model.type= model.type2,m.percent= m.percent2,rep=TRUE)
	write.csv(rbind(pbwlm2rep$boot.est.basins, pbwlm2rep$boot.se.basins),
		file.path(file.boot,'repwlm2.75.boot.basin.csv'))
	write.csv(rbind(pbwlm2rep$boot.est.orig, pbwlm2rep $boot.se.orig),
		file.path(file.boot,'repwlm2.75.boot.orig.csv'))
	temp2rep <-data.frame(pbwlm2rep$all.coeff.est)
	names(temp2rep)<-names(coef(wlm2))
	write.csv(temp2rep,file.path(file.boot,'repwlm2.75.boot.all.csv'))
}

#m/n = 1.00, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm2.boot.all.csv'))){
repwlm2.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repwlm2.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent3, model.type= model.type2, ordered_basins= ordered_basins)
}else {
	set.seed(3452345)
	pbwlm2rep<-para.boot.lm(boot.num=boots,model.object=wlm2,
		model.type= model.type2,m.percent = m.percent3,rep=TRUE)
	write.csv(rbind(pbwlm2rep$boot.est.basins, pbwlm2rep$boot.se.basins),
		file.path(file.boot,'repwlm2.boot.basin.csv'))
	write.csv(rbind(pbwlm2rep$boot.est.orig, pbwlm2 $boot.se.orig),
		file.path(file.boot,'repwlm2.boot.orig.csv'))
	temp2rep <-data.frame(pbwlm2rep$all.coeff.est)
	names(temp2rep)<-names(coef(wlm2))
	write.csv(temp2rep,file.path(file.boot,'repwlm2.boot.all.csv'))
}

#################################################################
#Bootstrapping WLS 3
#1000 bootstraps then outputs to 3 csv files
#################################################################

#m/n = 0.6666, MooN parametric subsampling
if(file.exists(file.path(file.boot,'wlm3.66.boot.all.csv'))){
wlm3.66.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'wlm3.66.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent1, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(4323)
	pbwlm3<-para.boot.lm(boot.num=boots,model.object=wlm3,
		model.type= model.type3,m.percent= m.percent1)
	write.csv(rbind(pbwlm3$boot.est.basins, pbwlm3$boot.se.basins),
		file.path(file.boot,'wlm3.66.boot.basin.csv'))
	write.csv(rbind(pbwlm3$boot.est.orig, pbwlm3$boot.se.orig),
		file.path(file.boot,'wlm3.66.boot.orig.csv'))
	temp3<-data.frame(pbwlm3$all.coeff.est)
	names(temp3)<-names(coef(wlm3))
	write.csv(temp3,file.path(file.boot,'wlm3.66.boot.all.csv'))
}

#m/n = 0.6666, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm3.66.boot.all.csv'))){
repwlm3.66.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repwlm3.66.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent1, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(4323)
	pbwlm3rep<-para.boot.lm(boot.num=boots,model.object=wlm3,
		model.type= model.type3,m.percent= m.percent1,rep=TRUE)
	write.csv(rbind(pbwlm3rep$boot.est.basins, pbwlm3rep$boot.se.basins),
		file.path(file.boot,'repwlm3.66.boot.basin.csv'))
	write.csv(rbind(pbwlm3rep$boot.est.orig, pbwlm3rep $boot.se.orig),
		file.path(file.boot,'repwlm3.66.boot.orig.csv'))
	temp3rep <-data.frame(pbwlm3rep$all.coeff.est)
	names(temp3rep)<-names(coef(wlm3))
	write.csv(temp3rep,file.path(file.boot,'repwlm3.66.boot.all.csv'))
}

#m/n = 0.75, MooN parametric subsampling
if(file.exists(file.path(file.boot,'wlm3.75.boot.all.csv'))){
wlm3.75.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'wlm3.75.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent2, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(54039)
	pbwlm3 <-para.boot.lm(boot.num=boots, model.object= wlm3,
		model.type= model.type3,m.percent= m.percent2)
	write.csv(rbind(pbwlm3$boot.est.basins, pbwlm3$boot.se.basins),
		file.path(file.boot,'wlm3.75.boot.basin.csv'))
	write.csv(rbind(pbwlm3$boot.est.orig, pbwlm3$boot.se.orig),
		file.path(file.boot,'wlm3.75.boot.orig.csv'))
	temp3 <-data.frame(pbwlm3$all.coeff.est)
	names(temp3)<-names(coef(wlm3))
	write.csv(temp3,file.path(file.boot,'wlm3.75.boot.all.csv'))
}

#m/n = 0.75, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm3.75.boot.all.csv'))){
repwlm3.75.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repwlm3.75.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent2, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(54039)
	pbwlm3rep<-para.boot.lm(boot.num=boots,model.object=wlm3,
		model.type= model.type3,m.percent= m.percent2,rep=TRUE)
	write.csv(rbind(pbwlm3rep$boot.est.basins, pbwlm3rep$boot.se.basins),
		file.path(file.boot,'repwlm3.75.boot.basin.csv'))
	write.csv(rbind(pbwlm3rep$boot.est.orig, pbwlm3rep $boot.se.orig),
		file.path(file.boot,'repwlm3.75.boot.orig.csv'))
	temp3rep <-data.frame(pbwlm3rep$all.coeff.est)
	names(temp3rep)<-names(coef(wlm3))
	write.csv(temp3rep,file.path(file.boot,'repwlm3.75.boot.all.csv'))
}

#m/n = 1.00, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repwlm3.boot.all.csv'))){
repwlm3.boot<- WLS_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repwlm3.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent3, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(324)
	pbwlm3rep<-para.boot.lm(boot.num=boots,model.object=wlm3,
		model.type= model.type3, m.percent = m.percent3,rep=TRUE)
	write.csv(rbind(pbwlm3rep$boot.est.basins, pbwlm3rep$boot.se.basins),
		file.path(file.boot,'repwlm3.boot.basin.csv'))
	write.csv(rbind(pbwlm3rep$boot.est.orig, pbwlm3rep $boot.se.orig),
		file.path(file.boot,'repwlm3.boot.orig.csv'))
	temp3rep <-data.frame(pbwlm3rep$all.coeff.est)
	names(temp3rep)<-names(coef(wlm3))
	write.csv(temp3rep,file.path(file.boot,'repwlm3.boot.all.csv'))
}

SE_boots<-rbind(
	data.frame(SE= wlm1.66.boot$Boot.SE, Model_Type=1, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='No'),
	 
	data.frame(SE=repwlm1.66.boot$Boot.SE, Model_Type=1, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	 
	data.frame(SE=wlm1.75.boot$Boot.SE, Model_Type=1, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='No'),
	 
	data.frame(SE=repwlm1.75.boot$Boot.SE, Model_Type=1, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	 
	data.frame(SE=repwlm1.boot$Boot.SE, Model_Type=1, Method='Par_1.00', 
	 Par=c(rep('beta0',7),rep('beta1',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	
	
	data.frame(SE= wlm2.66.boot$Boot.SE, Model_Type=2, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='No'),
	 
	data.frame(SE=repwlm2.66.boot$Boot.SE, Model_Type=2, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	 
	data.frame(SE=wlm2.75.boot$Boot.SE, Model_Type=2, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='No'),
	 
	data.frame(SE=repwlm2.75.boot$Boot.SE, Model_Type=2, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	 
	data.frame(SE=repwlm2.boot$Boot.SE, Model_Type=2, Method='Par_1.00', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	
	
	data.frame(SE= wlm3.66.boot$Boot.SE, Model_Type=3, Method='Par_0.66',
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='No'),
	 
	data.frame(SE=repwlm3.66.boot$Boot.SE, Model_Type=3, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	 
	data.frame(SE=wlm3.75.boot$Boot.SE, Model_Type=3, Method='Par_0.75', 
	  Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	  rep('beta3',7),rep('beta4',7)),
	  Basin=ordered_basins, Model='WLS', Replacement='No'),
	  
	data.frame(SE=repwlm3.75.boot$Boot.SE, Model_Type=3, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes'),
	 
	data.frame(SE=repwlm3.boot$Boot.SE, Model_Type=3, Method='Par_1.00', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7)),
	 Basin=ordered_basins, Model='WLS', Replacement='Yes')
)
#################################################################
#################################################################
#Calculating the Basin estimates
#################################################################
#################################################################
ordered_basins<-c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn')

#################################################################
#WLS Model 1 Basin Coefficients
#################################################################

wlm1.est<-Basin_estimates(ordered_basins= ordered_basins, 
	fit1= wlm1,fit2=NULL, Model_Type=1, exp.int=FALSE)

Intercepts_M1<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=1,Par='beta0',Estimate=wlm1.est$Intercepts_kt$Estimate,
				SE=wlm1.est$Intercepts_kt$SE, Method='delta')

Slopes_M1<- data.frame(Model='WLS',Basin=ordered_basins, 
			Model_Type=1,Par='beta1',Estimate=wlm1.est$Slopes_kt$Estimate,
				SE=wlm1.est$Slopes_kt$SE, Method='delta')

mod1<-rbind(Intercepts_M1, Slopes_M1)

#################################################################
#WLS Model 2 Basin Coefficients
#################################################################

wlm2.est<-Basin_estimates(ordered_basins= ordered_basins, 
	fit1= wlm2,fit2=NULL, Model_Type=2, exp.int=FALSE)

Intercepts_M2<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=2,Par='beta0',Estimate=wlm2.est$Intercepts_kt$Estimate,
				SE=wlm2.est$Intercepts_kt$SE, Method='delta')

slopes_M2<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=2,Par='beta1',Estimate=wlm2.est$Slopes_kt$Estimate,
				SE=wlm2.est$Slopes_kt$SE, Method='delta')

slopesquared_M2<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=2,Par='beta2',Estimate=wlm2.est$Slopes_sq_kt$Estimate,
				SE=wlm2.est$Slopes_sq_kt$SE, Method='delta')

mod2<-rbind(Intercepts_M2, slopes_M2, slopesquared_M2)

#################################################################
#WLS Model 3 Basin Coefficients
#################################################################

wlm3.est<-Basin_estimates(ordered_basins= ordered_basins, 
	fit1= wlm3,fit2=NULL, Model_Type=3, exp.int=FALSE)

Intercepts_M3<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=3,Par='beta0',Estimate=wlm3.est$Intercepts_kt$Estimate,
				SE=wlm3.est$Intercepts_kt$SE, Method='delta')

slopes_M3<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=3,Par='beta1',Estimate=wlm3.est$Slopes_kt$Estimate,
				SE=wlm3.est$Slopes_kt $SE, Method='delta')

slopesquared_M3<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=3,Par='beta2',Estimate=wlm3.est$Slopes_sq_kt$Estimate,
				SE=wlm3.est$Slopes_sq_kt $SE, Method='delta')


lat_M3<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=3,Par='beta3',Estimate=wlm3.est$lat_kt$Estimate,
				SE=wlm3.est$lat_kt $SE, Method='delta')

ts_M3<- data.frame(Model='WLS',Basin=ordered_basins, 
				Model_Type=3,Par='beta4',Estimate=wlm3.est$ts_kt$Estimate,
				SE=wlm3.est$ts_kt $SE, Method='delta')

mod3<-rbind(Intercepts_M3, slopes_M3, slopesquared_M3, lat_M3, ts_M3)


#Table of WLS Estimates
wlmmods<-rbind(mod1,mod2,mod3)

beta0 <- wlmmods[wlmmods$Par=='beta0',][,c(2:3,5)]
names(beta0)[3]<-'beta0'
beta1 <- wlmmods[wlmmods$Par=='beta1',][,c(2:3,5)]
names(beta1)[3]<-'beta1'
beta2 <- wlmmods[wlmmods$Par=='beta2',][,c(2:3,5)]
names(beta2)[3]<-'beta2'
beta3 <- wlmmods[wlmmods$Par=='beta3',][,c(2:3,5)]
names(beta3)[3]<-'beta3'
beta4 <- wlmmods[wlmmods$Par=='beta4',][,c(2:3,5)]
names(beta4)[3]<-'beta4'

wlstable<-merge(beta0,
  merge(beta1,
	merge(beta2,
	  merge(beta3, beta4,by=c('Basin','Model_Type'),all.x=T),
	by=c('Basin','Model_Type'),all.x=T),
  by=c('Basin','Model_Type'),all.x=T),
by=c('Basin','Model_Type'),all.x=T)

wlstable <- wlstable[order(wlstable$Model_Type),]
wlstable <- rbind(
wlstable[wlstable$Basin=='NoAt' & wlstable$Model_Type==1,],
wlstable[wlstable$Basin=='EaPa' & wlstable$Model_Type==1,],
wlstable[wlstable$Basin=='WePa' & wlstable$Model_Type==1,],
wlstable[wlstable$Basin=='NoIn' & wlstable$Model_Type==1,],
wlstable[wlstable$Basin=='SoPa' & wlstable$Model_Type==1,],
wlstable[wlstable$Basin=='SoEIn' & wlstable$Model_Type==1,],
wlstable[wlstable$Basin=='SoWIn' & wlstable$Model_Type==1,],
wlstable[wlstable$Basin=='NoAt' & wlstable$Model_Type==2,],
wlstable[wlstable$Basin=='EaPa' & wlstable$Model_Type==2,],
wlstable[wlstable$Basin=='WePa' & wlstable$Model_Type==2,],
wlstable[wlstable$Basin=='NoIn' & wlstable$Model_Type==2,],
wlstable[wlstable$Basin=='SoPa' & wlstable$Model_Type==2,],
wlstable[wlstable$Basin=='SoEIn' & wlstable$Model_Type==2,],
wlstable[wlstable$Basin=='SoWIn' & wlstable$Model_Type==2,],
wlstable[wlstable$Basin=='NoAt' & wlstable$Model_Type==3,],
wlstable[wlstable$Basin=='EaPa' & wlstable$Model_Type==3,],
wlstable[wlstable$Basin=='WePa' & wlstable$Model_Type==3,],
wlstable[wlstable$Basin=='NoIn' & wlstable$Model_Type==3,],
wlstable[wlstable$Basin=='SoPa' & wlstable$Model_Type==3,],
wlstable[wlstable$Basin=='SoEIn' & wlstable$Model_Type==3,],
wlstable[wlstable$Basin=='SoWIn' & wlstable$Model_Type==3,]
)
print(xtable(x=cbind(wlstable[,c(1:2)],
	format(round(wlstable[,c(3:7)],3),nsmall=3))), include.rownames=F)
	
#################################################################
#################################################################
#Plotting Model Coefficients with Standard Errors
#################################################################
#################################################################


SE_boot<-merge(SE_boots, wlmmods[,-c(6,7)],by=c('Basin','Model_Type','Par','Model'),all=T)
SE_WLS<-rbind(data.frame(wlmmods, Replacement="NotAPP"), SE_boot)
write.csv(SE_WLS,file.path(file.boot,'AllSEGraphingWLS.csv'), row.names = FALSE)

se_graph <-SE_WLS
se_graph$lo95<-se_graph$Estimate-2* se_graph$SE
se_graph$hi95<-se_graph$Estimate+2* se_graph$SE

se_graph$Estimate <-as.numeric(se_graph$Estimate)
se_graph$Replacement <- as.factor(se_graph$Replacement)
se_graph$Basin <- as.factor(se_graph$Basin)
se_graph$lo95<-se_graph$Estimate-2* se_graph$SE
se_graph$hi95<-se_graph$Estimate+2* se_graph$SE

se_graph<-rbind(se_graph[which(se_graph$Method!='Par_1.00'),],
se_graph[(se_graph$Method=='Par_1.00' & se_graph$Replacement!='No'),])
se_graph$Replacement<-relevel(se_graph$Replacement,'Yes')
se_graph$Method<-factor(se_graph $Method)

se_graph <-se_graph[order(se_graph$Method),]
se_graph$Method<-relevel(relevel(relevel(relevel(se_graph$Method,'delta'),'Par_1.00'),'Par_0.75'),'Par_0.66')
se_graph$Basin<-relevel(relevel(relevel(relevel(relevel(relevel(relevel(se_graph$Basin,'NoAt'),'EaPa'),'WePa'),'NoIn'),'SoPa'),'SoEIn'),'SoWIn')

#################################################################
#################################################################
#####WLS Model 1 Coefficients +2SE Plotting
#################################################################
#################################################################

#Beta0
pdf("beta0mod1wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==1 & se_graph$Model=='WLS' & se_graph$Par =='beta0',],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method, 
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	ylim(c(0,4))+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()

#Beta1
pdf("beta1mod1wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==1 &se_graph$Model=='WLS' & se_graph$Par =='beta1',],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14), axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()
#################################################################
#################################################################
#####WLS Model 2 Coefficients+2SE Plotting
#################################################################
#################################################################

#Beta0
pdf("beta0mod2wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==2 &se_graph$Model=='WLS' & se_graph$Par =='beta0' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()

#Beta1
pdf("beta1mod2wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==2 &se_graph$Model=='WLS' & se_graph$Par =='beta1' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()

#Beta2
pdf("beta2mod2wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==2 &se_graph$Model=='WLS' & se_graph$Par =='beta2' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()
#################################################################
#################################################################
#####WLS Model 3 Coefficients +2SE Plotting
#################################################################
#################################################################

#Beta0
pdf("beta0mod3wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='WLS' & se_graph$Par =='beta0' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()
	  
#Beta1
pdf("beta1mod3wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='WLS' & se_graph$Par =='beta1' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()

#Beta2
pdf("beta2mod3wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='WLS' & se_graph$Par =='beta2' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()

#Beta3
pdf("beta3mod3wls_upd.pdf")
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='WLS' & se_graph$Par =='beta3' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()

#Beta4
pdf("beta4mod3wls_upd.pdf")
options(repr.plot.width = 14, repr.plot.height = 8)
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='WLS' & se_graph$Par =='beta4' ,],
  aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,
  group= Method:Replacement,lty=Replacement))+ 
	geom_linerange(size=1,position= position_dodge(.5))+
	xlab('Basin')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	coord_flip()+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = 'none',axis.text.y = element_text(size=14),
	  axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
	  axis.title.y = element_text(size=16))
dev.off()
#theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))

#Running the MSPE/MAPE for bootstrap predictions

#WLM1
wlm1.66<- boot.pred(file=file.path(file.boot,'wlm1.66.boot.all.csv'),modeldat=modeling_data, model.type='1')
wlm1.75<- boot.pred(file=file.path(file.boot,'wlm1.75.boot.all.csv'),modeldat=modeling_data, model.type='1')
repwlm1.66<- boot.pred(file=file.path(file.boot,'repwlm1.66.boot.all.csv'),modeldat=modeling_data, model.type='1')
repwlm1.75<- boot.pred(file=file.path(file.boot,'repwlm1.75.boot.all.csv'),modeldat=modeling_data, model.type='1')
repwlm1<- boot.pred(file=file.path(file.boot,'repwlm1.boot.all.csv'),modeldat=modeling_data, model.type='1')

#WLM2
wlm2.66 <- boot.pred(file=file.path(file.boot,'wlm2.66.boot.all.csv'),modeldat=modeling_data, model.type='2')
wlm2.75 <- boot.pred(file=file.path(file.boot,'wlm2.75.boot.all.csv'),modeldat=modeling_data, model.type='2')
repwlm2.66 <- boot.pred(file=file.path(file.boot,'repwlm2.66.boot.all.csv'),modeldat=modeling_data, model.type='2')
repwlm2.75 <- boot.pred(file=file.path(file.boot,'repwlm2.75.boot.all.csv'),modeldat=modeling_data, model.type='2')
repwlm2 <- boot.pred(file=file.path(file.boot,'repwlm2.boot.all.csv'),modeldat=modeling_data, model.type='2')

#WLM3
wlm3.66 <- boot.pred(file=file.path(file.boot,'wlm3.66.boot.all.csv'),modeldat=modeling_data, model.type='3')
wlm3.75 <- boot.pred(file=file.path(file.boot,'wlm3.75.boot.all.csv'),modeldat=modeling_data, model.type='3')
repwlm3.66 <- boot.pred(file=file.path(file.boot,'repwlm3.66.boot.all.csv'),modeldat=modeling_data, model.type='3')
repwlm3.75 <- boot.pred(file=file.path(file.boot,'repwlm3.75.boot.all.csv'),modeldat=modeling_data, model.type='3')
repwlm3 <- boot.pred(file=file.path(file.boot,'repwlm3.boot.all.csv'),modeldat=modeling_data, model.type='3')

xtable(format(round(rbind(
  c(wlm1.66$Avg_MSPE_ws, repwlm1.66$Avg_MSPE_ws, wlm2.66$Avg_MSPE_ws, repwlm2.66$Avg_MSPE_ws, wlm3.66$Avg_MSPE_ws, repwlm3.66$Avg_MSPE_ws),
  c(wlm1.75$Avg_MSPE_ws, repwlm1.75$Avg_MSPE_ws, wlm2.75$Avg_MSPE_ws, repwlm2.75$Avg_MSPE_ws, wlm3.75$Avg_MSPE_ws, repwlm3.75$Avg_MSPE_ws),
  c(NA, repwlm1$Avg_MSPE_ws, NA, repwlm2$Avg_MSPE_ws, NA, repwlm3$Avg_MSPE_ws)
),5),4))


xtable(format(round(rbind(
  c(wlm1.66$Avg_MAPE_ws, repwlm1.66$Avg_MAPE_ws, wlm2.66$Avg_MAPE_ws, repwlm2.66$Avg_MAPE_ws, wlm3.66$Avg_MAPE_ws, repwlm3.66$Avg_MAPE_ws),
  c(wlm1.75$Avg_MAPE_ws, repwlm1.75$Avg_MAPE_ws, wlm2.75$Avg_MAPE_ws, repwlm2.75$Avg_MAPE_ws, wlm3.75$Avg_MAPE_ws, repwlm3.75$Avg_MAPE_ws),
  c(NA, repwlm1$Avg_MAPE_ws, NA, repwlm2$Avg_MAPE_ws, NA, repwlm3$Avg_MAPE_ws)
),5),4))






