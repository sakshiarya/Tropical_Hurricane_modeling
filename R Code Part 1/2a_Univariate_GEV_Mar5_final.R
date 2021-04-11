##########################################################################
##########################################################################
##########################################################################
# Title: Univariate GEV Modeling  
# Author: Lindsey Dietz
# Last Updated: 3/12/21 (By Sakshi Arya)
#
# Goal: Fits GEV models for log(maxWS) as a function of log(1013.25-minCP) 
#
# Note data files should be in working directory.
##########################################################################
##########################################################################
##########################################################################

rm(list=ls())

library(ggplot2)
library(ismev)
library(MASS)
library(GGally)
library(extRemes)	
library(evd)
library(reshape2)
library(multcomp)
library(car)
library(msm)
library(xtable)
library(parallel)
library(doParallel)
filep<- '/home/aryax010/private/LindseyThesis/Chapter1'

#Loads all outside functions
source(file.path(filep,'Helper Functions/boot.pred.R'))
source(file.path(filep,'Helper Functions/Basin_estimates_GEV.R'))
source(file.path(filep,'Helper Functions/GEV_Boot.R'))
source(file.path(filep,'Helper Functions/para.boot.gev.R'))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

##########################################################################
##########################################################################
##########################################################################
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

##########################################################################
##########################################################################
##########################################################################

all_usable_data<-data.frame(rbind(cbind(NA_sum[NA_sum$maxWS>=64,],Basin='NoAt'), 
	cbind(WP_sum[WP_sum$maxWS>=64,],Basin='WePa'),
	cbind(EP_sum[EP_sum$maxWS>=64,],Basin='EaPa'), 
	cbind(SP_sum[SP_sum$maxWS>=64,],Basin='SoPa'),
	cbind(NI_sum[NI_sum$maxWS>=64,],Basin='NoIn'), 
	cbind(SeI_sum[SeI_sum$maxWS>=64,],Basin='SoEIn'), 
	cbind(SwI_sum[SwI_sum$maxWS>=64,],Basin='SoWIn')))

#1 knot = 0.514444 m/s
KnotstoMetersPerSec <- 0.514444
modeling_data<-with(all_usable_data, 
	na.omit(data.frame(maxWS_knots= maxWS, minCP,Basin, avgLat, Avg_Trans_Speed,
	maxWS_ms= maxWS* KnotstoMetersPerSec)))

modeling_data$cat<-1
modeling_data[modeling_data$maxWS_knots>=83,]$cat<-2
modeling_data[modeling_data$maxWS_knots>=96,]$cat<-3
modeling_data[modeling_data$maxWS_knots>=113,]$cat<-4
modeling_data[modeling_data$maxWS_knots>=137,]$cat<-5

modeling_data$min1013cp<-1013.25-modeling_data$minCP
modeling_data$logmin1013cp<-log(1013.25-modeling_data$minCP)
modeling_data$scalelogmin1013cp<-scale(log(1013.25-modeling_data$minCP))

ordered_basins=c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn')

##########################################################################
##########################################################################
#######################Log Data Maximum Wind Speed########################
##########################################################################
##########################################################################
#Scaling is recommended in extRemes docs but it doesn't seem necessary
gum_mod0 <-fevd(log(modeling_data$maxWS_knots),data= modeling_data,
	location.fun=~1,
	scale.fun=~1,
	shape.fun=~1,type='Gumbel')
	
gev_mod0 <-fevd(log(modeling_data$maxWS_knots),data= modeling_data,
	location.fun=~1,
	scale.fun=~1,
	shape.fun=~1,type='GEV')
				
gev_mod1 <-fevd(log(modeling_data$maxWS_knots),data= modeling_data,
	location.fun=~log(1013.25-minCP)*Basin,
	scale.fun=~1,
	shape.fun=~1,type='GEV', period.basis='Tropical Cyclone')
		
gev_mod2 <-fevd(log(modeling_data$maxWS_knots),data= modeling_data,
	location.fun=~log(1013.25-minCP)*Basin+I((log(1013.25-minCP))^2)* Basin,
	scale.fun=~1,
	shape.fun=~1,type='GEV', period.basis='Tropical Cyclone')	

# Added this line in the code for grlevd function: z[which(z<0)] <- abs(z[which(z<0)] )
#trace("grlevd",edit = TRUE)	
#untrace("grlevd")	
gev_mod3 <-fevd(log(modeling_data$maxWS_knots),data= modeling_data,
	location.fun=~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin+
		avgLat*Basin+ Avg_Trans_Speed*Basin,
	scale.fun=~1,
	shape.fun=~1,type='GEV',method = "MLE", period.basis='Tropical Cyclone')

gev_mod4 <-fevd(log(modeling_data$maxWS_knots),data= modeling_data,
                location.fun=~log(1013-minCP)*Basin+
                  avgLat*Basin+ Avg_Trans_Speed*Basin,
                scale.fun=~1,
                shape.fun=~1,type='GEV',method = "MLE", period.basis='Tropical Cyclone')


#Stationary GEV against a Gumbel distribution using a LRT 
lr.test(gum_mod0 , gev_mod0)
zerov0<- lr.test(gum_mod0 , gev_mod0)$p.value
#GEV Model 1 v. Model 2
lr.test(gev_mod1 , gev_mod2)
Onev2<- lr.test(gev_mod1 , gev_mod2)$p.value
#GEV Model 2 v. Model 3
lr.test(gev_mod2 , gev_mod3)
Twov3<- lr.test(gev_mod2 , gev_mod3)$p.value
#GEV Model 1 v. Model 3
lr.test(gev_mod1 , gev_mod3)
Onev3<- lr.test(gev_mod1 , gev_mod3)$p.value
#GEV Model 3 bs Model 4
lr.test(gev_mod2,gev_mod4)

gev_mod1_AIC<- summary(gev_mod1)$AIC; gev_mod1_BIC<-summary(gev_mod1)$BIC
gev_mod2_AIC<- summary(gev_mod2)$AIC; gev_mod2_BIC<-summary(gev_mod2)$BIC
gev_mod3_AIC<- summary(gev_mod3)$AIC; gev_mod3_BIC<-summary(gev_mod3)$BIC


xtable(rbind(c(1, format(round(gev_mod1$results$par['scale'],4),4), format(round(gev_mod1$results$par['shape'],4),4), 
	round(gev_mod1_AIC,0), round(gev_mod1_BIC,0)),
c(2, format(round(gev_mod2$results$par['scale'],4),4), format(round(gev_mod2$results$par['shape'],4),4), 
	round(gev_mod2_AIC,0), round(gev_mod2_BIC,0)),
c(3, format(round(gev_mod3$results$par['scale'],4),4), format(round(gev_mod3$results$par['shape'],4),4), round(gev_mod3_AIC,0), round(gev_mod3_BIC,0))), include.rownames=F)



#################################################################
#GEV Model 1 Basin Coefficients
#################################################################
gevest1<-Basin_estimates_GEV(ordered_basins,fit1=gev_mod1,data=modeling_data,Model_Type=1)

Intercepts_M1<- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=1,Par='beta0',Estimate= gevest1$Intercepts_kt$Estimate,
				SE= gevest1$Intercepts_kt$SE, Method='delta')

Slopes_M1 <- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=1,Par='beta1',Estimate= gevest1$Slopes_kt$Estimate,
				SE= gevest1$Slopes_kt$SE, Method='delta')

scale_M1 <- data.frame(Model='GEV',Basin=NA, 
				Model_Type=1,Par='scale',Estimate= gev_mod1$results$par['scale'],
				SE= summary(gev_mod1)$se.theta['scale'], Method='delta')

shape_M1 <- data.frame(Model='GEV',Basin=NA, 
				Model_Type=1,Par='shape',Estimate= gev_mod1$results$par['shape'],
				SE= summary(gev_mod1)$se.theta['shape'], Method='delta')

mod1<-rbind(Intercepts_M1, Slopes_M1, scale_M1, shape_M1)

#################################################################
#GEV Model 2 Basin Coefficients
#################################################################
gevest2<-Basin_estimates_GEV(ordered_basins,fit1=gev_mod2,data=modeling_data,Model_Type=2)

Intercepts_M2<- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=2,Par='beta0',Estimate= gevest2$Intercepts_kt$Estimate,
				SE= rep(gevest2$Intercepts_kt$SE,length(gevest2$Intercepts_kt$Estimate)), Method='delta')

Slopes_M2 <- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=2,Par='beta1',Estimate= gevest2$Slopes_kt$Estimate,
				SE= gevest2$Slopes_kt$SE, Method='delta')

Slopes_sq_M2 <- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=2,Par='beta2',Estimate= gevest2$Slopes_sq_kt$Estimate,
				SE= gevest2$Slopes_sq_kt$SE, Method='delta')

scale_M2 <- data.frame(Model='GEV',Basin=NA, 
				Model_Type=2,Par='scale',Estimate= gev_mod2$results$par['scale'],
				SE= NA, Method='delta')

shape_M2 <- data.frame(Model='GEV',Basin=NA, 
				Model_Type=2,Par='shape',Estimate= gev_mod2$results$par['shape'],
				SE= NA, Method='delta')

mod2<-rbind(Intercepts_M2, Slopes_M2, Slopes_sq_M2, scale_M2, shape_M2)

#################################################################
#GEV Model 3 Basin Coefficients
#################################################################
gevest3<-Basin_estimates_GEV(ordered_basins,fit1=gev_mod3,data=modeling_data,Model_Type=3)

Intercepts_M3<- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=3,Par='beta0',Estimate= gevest3$Intercepts_kt$Estimate,
				SE= gevest3$Intercepts_kt$SE, Method='delta')

Slopes_M3 <- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=3,Par='beta1',Estimate= gevest3$Slopes_kt$Estimate,
				SE= gevest3$Slopes_kt$SE, Method='delta')

Slopes_sq_M3 <- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=3,Par='beta2',Estimate= gevest3$Slopes_sq_kt$Estimate,
				SE= gevest3$Slopes_sq_kt$SE, Method='delta')

lat_M3 <- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=3,Par='beta3',Estimate= gevest3$lat_kt$Estimate,
				SE= gevest3$lat_kt$SE, Method='delta')

ts_M3 <- data.frame(Model='GEV',Basin=ordered_basins, 
				Model_Type=3,Par='beta4',Estimate= gevest3$ts_kt$Estimate,
				SE= gevest3$ts_kt$SE, Method='delta')

scale_M3 <- data.frame(Model='GEV',Basin=NA, 
				Model_Type=3,Par='scale',Estimate= gev_mod3$results$par['scale'],
				SE= NA, Method='delta')

shape_M3 <- data.frame(Model='GEV',Basin=NA, 
				Model_Type=3,Par='shape',Estimate= gev_mod3$results$par['shape'],
				SE=NA, Method='delta')

mod3<-rbind(Intercepts_M3, Slopes_M3, Slopes_sq_M3, lat_M3, ts_M3, scale_M3, shape_M3)

#Table of GEV Estimates
gevmods<-rbind(mod1,mod2,mod3)

beta0 <- gevmods[gevmods$Par=='beta0',][,c(2:3,5)]
names(beta0)[3]<-'beta0'
beta1 <- gevmods[gevmods$Par=='beta1',][,c(2:3,5)]
names(beta1)[3]<-'beta1'
beta2 <- gevmods[gevmods$Par=='beta2',][,c(2:3,5)]
names(beta2)[3]<-'beta2'
beta3 <- gevmods[gevmods$Par=='beta3',][,c(2:3,5)]
names(beta3)[3]<-'beta3'
beta4 <- gevmods[gevmods$Par=='beta4',][,c(2:3,5)]
names(beta4)[3]<-'beta4'

gevtable<-merge(beta0,
                merge(beta1,
                      merge(beta2,
                            merge(beta3, beta4,by=c('Basin','Model_Type'),all.x=T),
                            by=c('Basin','Model_Type'),all.x=T),
                      by=c('Basin','Model_Type'),all.x=T),
                by=c('Basin','Model_Type'),all.x=T)

gevtable <- gevtable[order(gevtable$Model_Type),]
gevtable <- rbind(
gevtable[gevtable$Basin=='NoAt' & gevtable$Model_Type==1,],
gevtable[gevtable$Basin=='EaPa' & gevtable$Model_Type==1,],
gevtable[gevtable$Basin=='WePa' & gevtable$Model_Type==1,],
gevtable[gevtable$Basin=='NoIn' & gevtable$Model_Type==1,],
gevtable[gevtable$Basin=='SoPa' & gevtable$Model_Type==1,],
gevtable[gevtable$Basin=='SoEIn' & gevtable$Model_Type==1,],
gevtable[gevtable$Basin=='SoWIn' & gevtable$Model_Type==1,],
gevtable[gevtable$Basin=='NoAt' & gevtable$Model_Type==2,],
gevtable[gevtable$Basin=='EaPa' & gevtable$Model_Type==2,],
gevtable[gevtable$Basin=='WePa' & gevtable$Model_Type==2,],
gevtable[gevtable$Basin=='NoIn' & gevtable$Model_Type==2,],
gevtable[gevtable$Basin=='SoPa' & gevtable$Model_Type==2,],
gevtable[gevtable$Basin=='SoEIn' & gevtable$Model_Type==2,],
gevtable[gevtable$Basin=='SoWIn' & gevtable$Model_Type==2,],
gevtable[gevtable$Basin=='NoAt' & gevtable$Model_Type==3,],
gevtable[gevtable$Basin=='EaPa' & gevtable$Model_Type==3,],
gevtable[gevtable$Basin=='WePa' & gevtable$Model_Type==3,],
gevtable[gevtable$Basin=='NoIn' & gevtable$Model_Type==3,],
gevtable[gevtable$Basin=='SoPa' & gevtable$Model_Type==3,],
gevtable[gevtable$Basin=='SoEIn' & gevtable$Model_Type==3,],
gevtable[gevtable$Basin=='SoWIn' & gevtable$Model_Type==3,]
)
print(xtable(x=cbind(gevtable[,c(1:2)],
	format(round(gevtable[,c(3:7)],3),nsmall=3))), include.rownames=F)


#εt = [1/ξ(t)] log {1 + ξ(t) [Xt − μ(t)] / σ(t)}
par(mar=c(3,3,0,0)+0.1, mgp =c(2,1,0))#sets margins of plotting area
pdf("GEVdens1_new.pdf")
plot(gev_mod1, type = "density",main='',xlab='Log(Max. Wind Speed)', cex.axis=1.25, cex.lab=1.25)
dev.off()
pdf("GEVdens2_new.pdf")
plot(gev_mod2, type = "density",main='',xlab='Log(Max. Wind Speed)', cex.axis=1.25, cex.lab=1.25)
dev.off()
pdf("GEVdens3_new.pdf")
plot(gev_mod3, type = "density",main='',xlab='Log(Max. Wind Speed)', cex.axis=1.25, cex.lab=1.25)	
dev.off()
#################################################################
#Return Level Plotting Model 1
#################################################################
return.levels.gev1<-data.frame(cbind(data.frame(matrix(
	exp(return.level(gev_mod1,c(2,5,20,100))),ncol=4),
	Basin=modeling_data$Basin),Storm =rownames(modeling_data)))
names(return.levels.gev1)[1:4]<-c('2 TC','5 TC','20 TC','100 TC')
plotting_returns_gev1 <-melt(id.var=c('Basin','Storm'), return.levels.gev1)
names(plotting_returns_gev1)[3:4]<-c('Return_Level','Maximum_Wind_Speed')
levels(plotting_returns_gev1$Basin)<-c("North Atlantic", "West Pacific", "East Pacific",
	"South Pacific", "North Indian", "South East Indian","South West Indian")
#plotting_returns_gev1$Basin<-relevel(relevel(relevel(relevel(plotting_returns_gev1$Basin,
#	'North Indian'),"East Pacific"),"West Pacific"),"North Atlantic")
pdf("Returns1_new.pdf")
ggplot(plotting_returns_gev1,aes(Return_Level, Maximum_Wind_Speed,color= Return_Level,group=Return_Level))+
	geom_boxplot(width=0.5,position=position_dodge(.6))+
	facet_wrap(~Basin,nrow=2)+xlab('Return Period')+ylab('Predicted Maximum Wind Speed (knots)')+
	theme_bw()+ 
	geom_hline(yintercept=c(65,84,96,114,135),color='grey',alpha=.5,lty=3, show.legend =TRUE)+ 
	annotate('text',4.5,65,label="1",size=3)+ 
	annotate('text',4.5,84,label="2",size=3)+ 
	annotate('text',4.5,96,label="3",size=3)+ 
	annotate('text',4.5,114,label="4",size=3)+
	annotate('text',4.5,135,label="5",size=3)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		legend.position="none")
dev.off()

#################################################################
#Return Level Plotting Model 2
#################################################################
return.levels.gev2<-data.frame(cbind(data.frame(matrix(
	exp(return.level(gev_mod2,c(2,5,20,100))),ncol=4),
	Basin=modeling_data$Basin),Storm =rownames(modeling_data)))
names(return.levels.gev2)[1:4]<-c('2 TC','5 TC','20 TC','100 TC')
plotting_returns_gev2 <-melt(id.var=c('Basin','Storm'), return.levels.gev2)
names(plotting_returns_gev2)[3:4]<-c('Return_Level','Maximum_Wind_Speed')
levels(plotting_returns_gev2$Basin)<-c("North Atlantic", "West Pacific", "East Pacific",
	 "South Pacific", "North Indian", "South East Indian","South West Indian")
#plotting_returns_gev2$Basin<-relevel(relevel(relevel(relevel(plotting_returns_gev1$Basin,
#	'North Indian'),"East Pacific"),"West Pacific"),"North Atlantic")

pdf("Returns2_new.pdf")
ggplot(plotting_returns_gev2,aes(Return_Level, Maximum_Wind_Speed,color= Return_Level,group=Return_Level))+
	geom_boxplot(width=0.5,position=position_dodge(.6))+
	facet_wrap(~Basin,nrow=2)+xlab('Return Period')+
	ylab('Predicted Maximum Wind Speed (knots)')+
	theme_bw()+ geom_hline(yintercept=c(65,84,96,114,135),color='grey',
		alpha=.5,lty=3, show.legend =TRUE)+ 
	annotate('text',4.5,65,label="1",size=3)+ 
	annotate('text',4.5,84,label="2",size=3)+ 
	annotate('text',4.5,96,label="3",size=3)+ 
	annotate('text',4.5,114,label="4",size=3)+
	annotate('text',4.5,135,label="5",size=3)+
	theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
		legend.position="none")
dev.off()


#################################################################
#Return Level Plotting Model 3
#################################################################
return.levels.gev3<-data.frame(cbind(data.frame(matrix(
	exp(return.level(gev_mod3,c(2,5,20,100))),ncol=4),
	Basin=modeling_data$Basin),Storm =rownames(modeling_data)))
names(return.levels.gev3)[1:4]<-c('2 TC','5 TC','20 TC','100 TC')
plotting_returns_gev3<-melt(id.var=c('Basin','Storm'), return.levels.gev3)
names(plotting_returns_gev3)[3:4]<-c('Return_Level','Maximum_Wind_Speed')
levels(plotting_returns_gev3$Basin)<-c("North Atlantic", "West Pacific", "East Pacific", 
	"South Pacific", "North Indian", "South East Indian","South West Indian")
#plotting_returns_gev3$Basin<-relevel(relevel(relevel(relevel(plotting_returns_gev1$Basin,
#	'North Indian'),"East Pacific"),"West Pacific"),"North Atlantic")
pdf("Returns3_new.pdf")
ggplot(plotting_returns_gev3,aes(Return_Level, Maximum_Wind_Speed,color= Return_Level,group=Return_Level))+
	geom_boxplot(width=0.5,position=position_dodge(.6))+
	facet_wrap(~Basin,nrow=2)+xlab('Return Period')+
	ylab('Predicted Maximum Wind Speed (knots)')+
	theme_bw()+ 
	geom_hline(yintercept=c(65,84,96,114,135),color='grey',alpha=.5,lty=3, show.legend =TRUE)+ 
	annotate('text',4.5,65,label="1",size=3)+ 
	annotate('text',4.5,84,label="2",size=3)+ 
	annotate('text',4.5,96,label="3",size=3)+ 
	annotate('text',4.5,114,label="4",size=3)+
	annotate('text',4.5,135,label="5",size=3)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		legend.position="none")
dev.off()
#################################################################
#################################################################
#Bootstrapping SE calculations and plotting
#################################################################
#################################################################
file.boot <- '/home/aryax010/private/LindseyThesis/Chapter1/Data/Derived Data 21/Model Output Data'

boots<-1000
model.type1=c('1','2','3')[1]
model.type2=c('1','2','3')[2]
model.type3=c('1','2','3')[3]
m.percent1 = c(0.6666, 0.75, 1.00)[1]
m.percent2 = c(0.6666, 0.75, 1.00)[2]
m.percent3 = c(0.6666, 0.75, 1.00)[3]
#################################################################
#Bootstrapping GEV 1
#1000 bootstraps then outputs to 3 csv files
#################################################################

#m/n = 0.6666, MooN parametric subsampling

if(file.exists(file.path(file.boot,'gev1.66.boot.all.csv'))){
gev1.66.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'gev1.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type1, ordered_basins= ordered_basins)
} else {
	set.seed(10093)
	pbgev1<-para.boot.gev.Mar11(boot.num= boots,model.object= gev_mod1,
		model.type= model.type1 ,m.percent= m.percent1, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev1$boot.est.basins,pbgev1$boot.se.basins),
		file.path(file.boot,'gev1.66.boot.basin.csv'))
	temp1<-pbgev1$all.coeff.est
	write.csv(temp1,file.path(file.boot,'gev1.66.boot.all.csv'))
}

#m/n = 0.6666, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev1.66.boot.all.csv'))){
repgev1.66.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repgev1.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type1, ordered_basins= ordered_basins)
} else {
 	set.seed(10093)
	pbgev1<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod1,
		model.type= model.type1,m.percent= m.percent1,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev1$boot.est.basins,pbgev1$boot.se.basins),
		file.path(file.boot,'repgev1.66.boot.basin.csv'))
	temp1<-pbgev1$all.coeff.est
	write.csv(temp1,file.path(file.boot,'repgev1.66.boot.all.csv'))
}


#m/n = 0.75, MooN parametric subsampling

if(file.exists(file.path(file.boot,'gev1.75.boot.all.csv'))){
gev1.75.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'gev1.75.boot.all.csv'),
	header=T, check.names = FALSE)[,-1], m.percent= m.percent2, 
	model.type= model.type1, ordered_basins= ordered_basins)
} else {
  set.seed(1023423121)
  pbgev1<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod1,
  	model.type= model.type1,m.percent= m.percent2, 
		data= modeling_data, ordered_basins= ordered_basins)
  write.csv(rbind(pbgev1$boot.est.basins,pbgev1$boot.se.basins),
  	file.path(file.boot,'gev1.75.boot.basin.csv'))
  temp1<-pbgev1$all.coeff.est
  write.csv(temp1, file.path(file.boot,'gev1.75.boot.all.csv'))
}


#m/n = 0.75, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev1.75.boot.all.csv'))){
repgev1.75.boot<- GEV_Boot(data= modeling_data, 
  	data.file=read.csv(file.path(file.boot,'repgev1.75.boot.all.csv'),
  	header=T, check.names = FALSE)[,-1], m.percent= m.percent2,
  	model.type= model.type1, ordered_basins= ordered_basins)
} else {
  set.seed(1023423121)
  pbgev1<-para.boot.gev.Mar11(boot.num= boots,model.object= gev_mod1,
  	model.type= model.type1 ,m.percent= m.percent2,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
  write.csv(rbind(pbgev1$boot.est.basins,pbgev1$boot.se.basins),
  	file.path(file.boot,'repgev1.75.boot.basin.csv'))
  temp1<-pbgev1$all.coeff.est
  write.csv(temp1, file.path(file.boot,'repgev1.75.boot.all.csv'))
}

#m/n = 1.00, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev1.boot.all.csv'))){
repgev1.boot<-  GEV_Boot(data= modeling_data, 
  	data.file=read.csv(file.path(file.boot,'repgev1.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
  	m.percent = m.percent3, model.type= model.type1, ordered_basins= ordered_basins)
} else {
  set.seed(131313)
  pbgev1.wrep<-para.boot.gev.Mar11(boot.num=boots, model.object= gev_mod1,
  	model.type= model.type1,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
  write.csv(rbind(pbgev1.wrep$boot.est.basins,pbgev1.wrep$boot.se.basins),file.path(file.boot,'repgev1.boot.basin.csv'))
  temp1.rep <-pbgev1.wrep$all.coeff.est
  write.csv(temp1.rep,file.path(file.boot,'repgev1.boot.all.csv'))
}
#################################################################
#Bootstrapping GEV 2
#1000 bootstraps then outputs to 3 csv files
#################################################################

#m/n = 0.6666, MooN parametric subsampling
if(file.exists(file.path(file.boot,'gev2.66.boot.all.csv'))){
gev2.66.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'gev2.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type2, ordered_basins= ordered_basins)
} else {
	set.seed(75699)
	pbgev2<-para.boot.gev.Mar11(boot.num= boots,model.object= gev_mod2,
		model.type= model.type2 ,m.percent= m.percent1, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev2$boot.est.basins,pbgev2$boot.se.basins),
		file.path(file.boot,'gev2.66.boot.basin.csv'))
	temp1<-pbgev2$all.coeff.est
	write.csv(temp1,file.path(file.boot,'gev2.66.boot.all.csv'))
}

#m/n = 0.6666, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev2.66.boot.all.csv'))){
repgev2.66.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repgev2.66.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent1, model.type= model.type2, ordered_basins= ordered_basins)
} else {
	set.seed(75699)
	pbgev2rep<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod2,
		model.type= model.type2,m.percent= m.percent1,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev2rep$boot.est.basins, pbgev2rep$boot.se.basins),
		file.path(file.boot,'repgev2.66.boot.basin.csv'))
	temp2rep <-pbgev2rep$all.coeff.est
	write.csv(temp2rep,file.path(file.boot,'repgev2.66.boot.all.csv'))
}

#m/n = 0.75, MooN parametric subsampling
if(file.exists(file.path(file.boot,'gev2.75.boot.all.csv'))){
gev2.75.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'gev2.75.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent2, model.type= model.type2, ordered_basins= ordered_basins)
}else {
	set.seed(3547898)
	pbgev2 <-para.boot.gev.Mar11(boot.num=boots ,model.object= gev_mod2,
		model.type= model.type2,m.percent= m.percent2, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev2$boot.est.basins, pbgev2$boot.se.basins),
		file.path(file.boot,'gev2.75.boot.basin.csv'))
	temp2 <-pbgev2$all.coeff.est
	write.csv(temp2,file.path(file.boot,'gev2.75.boot.all.csv'))
}

#m/n = 0.75, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev2.75.boot.all.csv'))){
repgev2.75.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repgev2.75.boot.all.csv'),header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent2, model.type= model.type2, ordered_basins= ordered_basins)
}else {
	set.seed(3547898)
	pbgev2rep<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod2,
		model.type= model.type2,m.percent= m.percent2,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev2rep$boot.est.basins, pbgev2rep$boot.se.basins),
		file.path(file.boot,'repgev2.75.boot.basin.csv'))
	temp2rep <-pbgev2rep$all.coeff.est
	write.csv(temp2rep,file.path(file.boot,'repgev2.75.boot.all.csv'))
}

#m/n = 1.00, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev2.boot.all.csv'))){
repgev2.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repgev2.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent3, model.type= model.type2, ordered_basins= ordered_basins)
}else {
	set.seed(42342)
	pbgev2rep<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod2,
		model.type= model.type2,m.percent = m.percent3,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev2rep$boot.est.basins, pbgev2rep$boot.se.basins),
		file.path(file.boot,'repgev2.boot.basin.csv'))
	temp2rep <-pbgev2rep$all.coeff.est
	write.csv(temp2rep,file.path(file.boot,'repgev2.boot.all.csv'))
}

#################################################################
#Bootstrapping GEV 3
#1000 bootstraps then outputs to 3 csv files
#################################################################

#m/n = 0.6666, MooN parametric subsampling
if(file.exists(file.path(file.boot,'gev3.66.boot.all.csv'))){
gev3.66.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'gev3.66.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent1, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(4323)
	pbgev3<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod3,
		model.type= model.type3,m.percent= m.percent1, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev3$boot.est.basins, pbgev3$boot.se.basins),
		file.path(file.boot,'gev3.66.boot.basin.csv'))
	temp3<-pbgev3$all.coeff.est
	write.csv(temp3,file.path(file.boot,'gev3.66.boot.all.csv'))
}

#m/n = 0.6666, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev3.66.boot.all.csv'))){
repgev3.66.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repgev3.66.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent1, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(4323)
	pbgev3rep<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod3,
		model.type= model.type3,m.percent= m.percent1,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev3rep$boot.est.basins, pbgev3rep$boot.se.basins),
		file.path(file.boot,'repgev3.66.boot.basin.csv'))
	temp3rep <-cbind(pbgev3rep$all.coeff.est,pbgev3rep$convergence)
	write.csv(temp3rep,file.path(file.boot,'repgev3.66.boot.all.csv'))
}

#m/n = 0.75, MooN parametric subsampling
if(file.exists(file.path(file.boot,'gev3.75.boot.all.csv'))){
gev3.75.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'gev3.75.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent2, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(54039)
	pbgev3 <-para.boot.gev.Mar11(boot.num=boots, model.object= gev_mod3,
		model.type= model.type3,m.percent= m.percent2, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev3$boot.est.basins, pbgev3$boot.se.basins),
		file.path(file.boot,'gev3.75.boot.basin.csv'))
	temp3 <-cbind(pbgev3$all.coeff.est,pbgev3$convergence)
	write.csv(temp3,file.path(file.boot,'gev3.75.boot.all.csv'))
}

#m/n = 0.75, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev3.75.boot.all.csv'))){
repgev3.75.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repgev3.75.boot.all.csv'), header=T, check.names = FALSE)[,-1], m.percent = m.percent2, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(54039)
	pbgev3rep<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod3,
		model.type= model.type3,m.percent= m.percent2,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev3rep$boot.est.basins, pbgev3rep$boot.se.basins),
		file.path(file.boot,'repgev3.75.boot.basin.csv'))
	temp3rep <-cbind(pbgev3rep$all.coeff.est,pbgev3rep$convergence)
	write.csv(temp3rep,file.path(file.boot,'repgev3.75.boot.all.csv'))
}

#m/n = 1.00, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repgev3.boot.all.csv'))){
repgev3.boot<- GEV_Boot(data= modeling_data, 
	data.file=read.csv(file.path(file.boot,'repgev3.boot.all.csv'), header=T, check.names = FALSE)[,-1], 
	m.percent = m.percent3, model.type= model.type3, ordered_basins= ordered_basins)
}else {
	set.seed(324)
	pbgev3rep<-para.boot.gev.Mar11(boot.num=boots,model.object= gev_mod3,
		model.type= model.type3, m.percent = m.percent3,rep=TRUE, 
		data= modeling_data, ordered_basins= ordered_basins)
	write.csv(rbind(pbgev3rep$boot.est.basins, pbgev3rep$boot.se.basins),
		file.path(file.boot,'repgev3.boot.basin.csv'))
	temp3rep <-cbind(pbgev3rep$all.coeff.est, pbgev3rep$convergence)
	write.csv(temp3rep,file.path(file.boot,'repgev3.boot.all.csv'))
}

SE_boots<-rbind(
	data.frame(SE= gev1.66.boot$Boot.SE, Model_Type=1, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,2),NA,NA), Model='GEV', Replacement='No'),
	 
	data.frame(SE=repgev1.66.boot$Boot.SE, Model_Type=1, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,2),NA,NA), Model='GEV', Replacement='Yes'),
	 
	data.frame(SE=gev1.75.boot$Boot.SE, Model_Type=1, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,2),NA,NA), Model='GEV', Replacement='No'),
	 
	data.frame(SE=repgev1.75.boot$Boot.SE, Model_Type=1, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,2),NA,NA), Model='GEV', Replacement='Yes'),
	 
	data.frame(SE=repgev1.boot$Boot.SE, Model_Type=1, Method='Par_1.00', 
	 Par=c(rep('beta0',7),rep('beta1',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,2),NA,NA), Model='GEV', Replacement='Yes'),
	
	
	data.frame(SE= gev2.66.boot$Boot.SE, Model_Type=2, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,3),NA,NA), Model='GEV', Replacement='No'),
	 
	data.frame(SE=repgev2.66.boot$Boot.SE, Model_Type=2, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,3),NA,NA), Model='GEV', Replacement='Yes'),
	 
	data.frame(SE=gev2.75.boot$Boot.SE, Model_Type=2, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,3),NA,NA), Model='GEV', Replacement='No'),
	 
	data.frame(SE=repgev2.75.boot$Boot.SE, Model_Type=2, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,3),NA,NA), Model='GEV', Replacement='Yes'),
	 
	data.frame(SE=repgev2.boot$Boot.SE, Model_Type=2, Method='Par_1.00', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,3),NA,NA), Model='GEV', Replacement='Yes'),
	
	
	data.frame(SE= gev3.66.boot$Boot.SE, Model_Type=3, Method='Par_0.66',
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,5),NA,NA), Model='GEV', Replacement='No'),
	 
	data.frame(SE=repgev3.66.boot$Boot.SE, Model_Type=3, Method='Par_0.66', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,5),NA,NA), Model='GEV', Replacement='Yes'),
	 
	data.frame(SE=gev3.75.boot$Boot.SE, Model_Type=3, Method='Par_0.75', 
	  Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	  rep('beta3',7),rep('beta4',7),'scale','shape'),
	  Basin=c(rep(ordered_basins,5),NA,NA), Model='GEV', Replacement='No'),
	  
	data.frame(SE=repgev3.75.boot$Boot.SE, Model_Type=3, Method='Par_0.75', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,5),NA,NA), Model='GEV', Replacement='Yes'),
	 
	data.frame(SE=repgev3.boot$Boot.SE, Model_Type=3, Method='Par_1.00', 
	 Par=c(rep('beta0',7),rep('beta1',7),rep('beta2',7),
	 rep('beta3',7),rep('beta4',7),'scale','shape'),
	 Basin=c(rep(ordered_basins,5),NA,NA), Model='GEV', Replacement='Yes')
)


#################################################################
#################################################################
#Plotting Model Coefficients with Standard Errors
#################################################################
#################################################################


SE_boot<-merge(SE_boots, gevmods[,-c(6,7)],by=c('Basin','Model_Type','Par','Model'),all=T)
SE_GEV<-rbind(data.frame(gevmods, Replacement="NotAPP"), SE_boot)
write.csv(SE_GEV,file.path(file.boot,'AllSEGraphingGEV.csv'), row.names = FALSE)

se_graph <-SE_GEV
se_graph$lo95<-se_graph$Estimate-2*se_graph$SE
se_graph$hi95<-se_graph$Estimate+2*se_graph$SE

se_graph$Estimate <-as.numeric(se_graph$Estimate)
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
#####GEV Model 1 Coefficients +2SE Plotting
#################################################################
#################################################################

#Beta0
ggplot(se_graph[se_graph$Model_Type==1 & se_graph$Model=='GEV' & se_graph$Par =='beta0',],
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

#Beta1
ggplot(se_graph[se_graph$Model_Type==1 &se_graph$Model=='GEV' & se_graph$Par =='beta1',],
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

#Scale/Shape
ggplot(se_graph[se_graph$Model_Type==1 &se_graph$Model=='GEV' & 
		se_graph$Par %in%c('scale','shape'),],
		aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,
			group= Method:Replacement,lty= Replacement))+ 
	geom_linerange(lwd=1,position= position_dodge(.5))+
	xlab('Parameter')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	coord_flip()+
	scale_x_discrete(label=c(expression(paste(sigma)),expression(paste(xi))))+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),
	axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),
	axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))
#################################################################
#################################################################
#####GEV Model 2 Coefficients+2SE Plotting
#################################################################
#################################################################

#Beta0
ggplot(se_graph[se_graph$Model_Type==2 &se_graph$Model=='GEV' & se_graph$Par =='beta0' ,],
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


#Beta1
ggplot(se_graph[se_graph$Model_Type==2 &se_graph$Model=='GEV' & se_graph$Par =='beta1' ,],
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
	  
#Beta2
ggplot(se_graph[se_graph$Model_Type==2 &se_graph$Model=='GEV' & se_graph$Par =='beta2' ,],
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

#Scale/Shape
ggplot(se_graph[se_graph$Model_Type==2 &se_graph$Model=='GEV' & 
		se_graph$Par %in%c('scale','shape'),],
		aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,
			group= Method:Replacement,lty= Replacement))+ 
	geom_linerange(lwd=1,position= position_dodge(.5))+
	xlab('Parameter')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	coord_flip()+
	scale_x_discrete(label=c(expression(paste(sigma)),expression(paste(xi))))+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),
	axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),
	axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))
#################################################################
#################################################################
#####GEV Model 3 Coefficients +2SE Plotting
#################################################################
#################################################################

#Beta0
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='GEV' & se_graph$Par =='beta0' ,],
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
	  

#Beta1
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='GEV' & se_graph$Par =='beta1' ,],
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


#Beta2
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='GEV' & se_graph$Par =='beta2' ,],
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

#Beta3
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='GEV' & se_graph$Par =='beta3' ,],
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

#Beta4
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='GEV' & se_graph$Par =='beta4' ,],
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

#Scale/Shape
ggplot(se_graph[se_graph$Model_Type==3 &se_graph$Model=='GEV' & 
		se_graph$Par %in%c('scale','shape'),],
		aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,
			group= Method:Replacement,lty= Replacement))+ 
	geom_linerange(lwd=1,position= position_dodge(.5))+
	xlab('Parameter')+ylab('Estimate +/- 2 SE')+
	theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	coord_flip()+
	scale_x_discrete(label=c(expression(paste(sigma)),expression(paste(xi))))+
	geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),
	axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),
	axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))



## MSPE-MAPE 

## Model 1	
# bootstrap 0.66 
Mod1_66_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod1, model.type = '1', iterEach = 200, data = modeling_data, m.percent = 0.66, rep = FALSE)

# subsampling 0.66 
Mod1_66_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod1, model.type = '1', iterEach = 200, data = modeling_data, m.percent = 0.66, rep = TRUE)

# bootstrap 0.75 
Mod1_75_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod1, model.type = '1', iterEach = 200, data = modeling_data, m.percent = 0.75, rep = FALSE)

# subsampling 0.75
Mod1_75_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod1, model.type = '1', iterEach = 200, data = modeling_data, m.percent = 0.75, rep = TRUE)

# subsampling 1
Mod1_1_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod1, model.type = '1', iterEach = 200, data = modeling_data, m.percent = 1, rep = TRUE)


## Model 2
# bootstrap 0.66 
Mod2_66_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod2, model.type = '2', iterEach = 200, data = modeling_data, m.percent = 0.66, rep = FALSE)


# subsampling 0.66 
Mod2_66_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod2, model.type = '2', iterEach = 200, data = modeling_data, m.percent = 0.66, rep = TRUE)


# bootstrap 0.75 
Mod2_75_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod2, model.type = '2', iterEach = 200, data = modeling_data, m.percent = 0.75, rep = FALSE)


# subsampling 0.75
Mod2_75_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod2, model.type = '2', iterEach = 200, data = modeling_data, m.percent = 0.75, rep = TRUE)

# subsampling 1
Mod2_1_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod2, model.type = '2', iterEach = 200, data = modeling_data, m.percent = 1, rep = TRUE)


## Model 3
# bootstrap 0.66 
Mod3_66_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod3, model.type = '3', iterEach = 200, data = modeling_data, m.percent = 0.66, rep = FALSE)

# subsampling 0.66 
Mod3_66_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod3, model.type = '3', iterEach = 200, data = modeling_data, m.percent = 0.66, rep = TRUE)

# bootstrap 0.75 
Mod3_75_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod3, model.type = '3', iterEach = 200, data = modeling_data, m.percent = 0.75, rep = FALSE)

# subsampling 0.75
Mod3_75_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod3, model.type = '3', iterEach = 200, data = modeling_data, m.percent = 0.75, rep = TRUE)

# subsampling 1
Mod3_1_rep_gev_boot <- para.gev.pred(boot.num = 1000, model.object = gev_mod3, model.type = '3', iterEach = 200, data = modeling_data, m.percent = 1, rep = TRUE)



MSPE_Mod123 <- rbind(cbind(Mod1_66_gev_boot$Avg_MSPE, Mod1_66_rep_gev_boot$Avg_MSPE, Mod2_66_gev_boot$Avg_MSPE, Mod2_66_rep_gev_boot$Avg_MSPE, Mod3_66_gev_boot$Avg_MSPE, Mod3_66_rep_gev_boot$Avg_MSPE),
                     cbind(Mod1_75_gev_boot$Avg_MSPE, Mod1_75_rep_gev_boot$Avg_MSPE, Mod2_75_gev_boot$Avg_MSPE, Mod2_75_rep_gev_boot$Avg_MSPE, Mod3_75_gev_boot$Avg_MSPE, Mod3_75_rep_gev_boot$Avg_MSPE),
                     cbind(NA, Mod1_1_rep_gev_boot$Avg_MSPE, NA, Mod2_1_rep_gev_boot$Avg_MSPE, NA, Mod3_1_rep_gev_boot$Avg_MSPE))
MAPE_Mod123 <- rbind(cbind(Mod1_66_gev_boot$Avg_MAPE, Mod1_66_rep_gev_boot$Avg_MAPE, Mod2_66_gev_boot$Avg_MAPE, Mod2_66_rep_gev_boot$Avg_MAPE, Mod3_66_gev_boot$Avg_MAPE, Mod3_66_rep_gev_boot$Avg_MAPE),
                     cbind(Mod1_75_gev_boot$Avg_MAPE, Mod1_75_rep_gev_boot$Avg_MAPE, Mod2_75_gev_boot$Avg_MAPE, Mod2_75_rep_gev_boot$Avg_MAPE, Mod3_75_gev_boot$Avg_MAPE, Mod3_75_rep_gev_boot$Avg_MAPE),
                     cbind(NA, Mod1_1_rep_gev_boot$Avg_MAPE, NA, Mod2_1_rep_gev_boot$Avg_MAPE, NA, Mod3_1_rep_gev_boot$Avg_MAPE))


xtable(rbind(MSPE_Mod123, MAPE_Mod123), digits = 5)

