##########################################################################
##########################################################################
##########################################################################
# Title: Bivariate GEV Modeling  
# Author: Lindsey Dietz
# Edited by: Sakshi Arya
# Created by Lindsey: 3/1/16
# Last updated by Sakshi: 3/4/2021
#
# Goal: Fits bivariate GEV models for log(maxWS), log(1013-minCP) 
#
# Note data files should be in working directory.
##########################################################################
##########################################################################
##########################################################################

rm(list=ls())

library(xtable)
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
library(doParallel)
library(parallel)

filep<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure'

#Loads all outside functions
source(file.path(filep,'Helper Functions/BVEVD_Boot.R'))
source(file.path(filep,'Helper Functions/para.boot.bvevd.2021.R'))
#################################################################
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

#################################################################
#################################################################
#################################################################

all_usable_data<-data.frame(rbind(cbind(NA_sum[NA_sum$maxWS>=64,],Basin='NoAt'), 
	cbind(WP_sum[WP_sum$maxWS>=64,],Basin='WePa'),
	cbind(EP_sum[EP_sum$maxWS>=64,],Basin='EaPa'), 
	cbind(SP_sum[SP_sum$maxWS>=64,],Basin='SoPa'),
	cbind(NI_sum[NI_sum$maxWS>=64,],Basin='NoIn'), 
	cbind(SeI_sum[SeI_sum$maxWS>=64,],Basin='SoEIn'), 
	cbind(SwI_sum[SwI_sum$maxWS>=64,],Basin='SoWIn')))

#################################################################
#################################################################
#1 knot = 0.514444 m/s
KnotstoMetersPerSec <- 0.514444
modeling_data<-with(all_usable_data, 
	na.omit(data.frame(maxWS_knots= maxWS, minCP,Basin, avgLat, Avg_Trans_Speed,
	maxWS_ms= maxWS* KnotstoMetersPerSec)))


#modeling_data$Basin <- as.factor(modeling_data$Basin)
#levels(modeling_data$Basin)
#ordered_basins=c('EaPa','NoAt','NoIn','SoEIn','SoPa','SoWIn','WePa')

data.for.bivariate<-with(modeling_data, cbind(logmaxWS=log(maxWS_knots), logminCP=log(1013-minCP), Basin))
data.for.loc<-with(modeling_data, data.frame(Int=1,avgLat= avgLat, Avg_Trans_Speed = Avg_Trans_Speed))

nloc<-data.frame(with(modeling_data ,model.matrix(log(maxWS_knots) ~ avgLat*Basin+ Avg_Trans_Speed* Basin)))
data.for.loc1<-with(modeling_data, data.frame(avgLat))
data.for.loc2<-with(modeling_data, data.frame(Avg_Trans_Speed))
nloc.basin<-data.frame(with(modeling_data ,model.matrix(~Basin)))

data.for.bivariate <- as.data.frame(data.for.bivariate)
data.for.bivariate[,1] <- as.numeric(data.for.bivariate[,1])
data.for.bivariate[,2] <- as.numeric(data.for.bivariate[,2])

#Excluded via LRT
m1<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "log")
#Excluded via Qualitative Comparison
m2<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')], 
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "alog")
#Final Model
m3<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "bilog")
anova(m2,m1)
anova(m3,m1)

#Excluded via LRT
m4<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "neglog", std.err = FALSE)
#Excluded via AIC
m5<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "aneglog", std.err = FALSE)
#Excluded via AIC
m6<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "negbilog")
anova(m5,m4)

#Excluded via AIC
m7<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "hr")
#Excluded via AIC
m8<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "ct")
#Excluded via Nonconvergence
m9<-fbvevd(data.for.bivariate[,c('logmaxWS','logminCP')],
	nsloc1= nloc.basin[,-1], nsloc2= nloc.basin[,-1], model = "amix")

AIC(m1);AIC(m2);AIC(m3);
AIC(m4);AIC(m5);AIC(m6);
AIC(m7);AIC(m8);AIC(m9)

#Plot of dependence function for alog
plot(m2, which = 4, nplty = 1, lty = 3,lwd=2,col='blue',main='Dependence: alog')
legend('bottomleft', legend=c('Nonparametric','alog'),
	col=c('blue','black'),lwd=c(2,1),lty=c(1,3))

#Plot of dependence function for bilog
plot(m3, which = 4, nplty = 1, lty = 4,lwd=2,col='blue',main='Dependence: bilog')
legend('bottomleft',legend =c('Nonparametric','bilog'),
	col=c('blue','black'),lwd=c(2,1),lty=c(1,4))

final.model<-m3
bivariate.estimates<-final.model$estimate

loc1.EP<-bivariate.estimates['loc1']
loc2.EP<-bivariate.estimates['loc2']
loc1.NA<-loc1.EP+ bivariate.estimates['loc1BasinNoAt']
loc2.NA<-loc2.EP+ bivariate.estimates['loc2BasinNoAt']
loc1.NI<-loc1.EP+ bivariate.estimates['loc1BasinNoIn']
loc2.NI<-loc2.EP+ bivariate.estimates['loc2BasinNoIn']
loc1.SEI<-loc1.EP+ bivariate.estimates['loc1BasinSoEIn']
loc2.SEI<-loc2.EP+ bivariate.estimates['loc2BasinSoEIn']
loc1.SP<-loc1.EP+ bivariate.estimates['loc1BasinSoPa']
loc2.SP<-loc2.EP+ bivariate.estimates['loc2BasinSoPa']
loc1.SWI<-loc1.EP+ bivariate.estimates['loc1BasinSoWIn']
loc2.SWI<-loc2.EP+ bivariate.estimates['loc2BasinSoWIn']
loc1.WP<-loc1.EP+ bivariate.estimates['loc1BasinWePa']
loc2.WP<- loc2.EP+ bivariate.estimates['loc2BasinWePa']

scale1<-bivariate.estimates['scale1']
scale2<-bivariate.estimates['scale2']

shape1<-bivariate.estimates['shape1']
shape2<-bivariate.estimates['shape2']

alpha<-bivariate.estimates['alpha']
beta<-bivariate.estimates['beta']

cov.mat<-vcov(final.model)

Ints_M1<- c(loc1.EP,
				loc1.NA,
				loc1.NI,
				loc1.SEI,
				loc1.SP,
				loc1.SWI,
				loc1.WP)

Ints_SE_M1<-c(
	deltamethod(~x1, bivariate.estimates, cov.mat, ses=TRUE), 
	deltamethod(~x1+x2, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x1+x3, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x1+x4, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x1+x5, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x1+x6, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x1+x7, bivariate.estimates, cov.mat, ses=TRUE))

Ints_M2<- c(loc2.EP,
				loc2.NA,
				loc2.NI,
				loc2.SEI,
				loc2.SP,
				loc2.SWI,
				loc2.WP)

Ints_SE_M2<-c(
	deltamethod(~x10, bivariate.estimates, cov.mat, ses=TRUE), 
	deltamethod(~x10+x11, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x10+x12, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x10+x13, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x10+x14, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x10+x15, bivariate.estimates, cov.mat, ses=TRUE),
	deltamethod(~x10+x16, bivariate.estimates, cov.mat, ses=TRUE))

B<-c('EP','NA','NI','SEI','SP','SWI','WP')

xtable(data.frame(B, paste(format(round(c(Ints_M1),3),3),
			' (',format(round(unlist(Ints_SE_M1),3)),')',sep=''), 
			paste(format(round(c(Ints_M2),3),3),
			' (',format(round(unlist(Ints_SE_M2),3)),')',sep='')))

ests<-c(loc1.EP, loc1.NA, loc1.NI, loc1.SEI, loc1.SP, loc1.SWI, loc1.WP, scale1, shape1, alpha,loc2.EP, loc2.NA, loc2.NI, loc2.SEI, loc2.SP, loc2.SWI, loc2.WP, scale2, shape2, beta )


set.seed(611984)
sim.EP=rbvevd(n=1000,alpha= alpha,beta= beta,model='bilog',mar1=c(loc1.EP, scale1, shape1),mar2=c(loc2.EP, scale2, shape2))
sim.NA=rbvevd(n=1000,alpha= alpha,beta= beta,model='bilog',mar1=c(loc1.NA, scale1, shape1),mar2=c(loc2.NA, scale2, shape2))
sim.NI=rbvevd(n=1000,alpha= alpha,beta= beta,model='bilog',mar1=c(loc1.NI, scale1, shape1),mar2=c(loc2.NI, scale2, shape2))
sim.SEI=rbvevd(n=1000,alpha= alpha,beta= beta,model='bilog',mar1=c(loc1.SEI, scale1, shape1),mar2=c(loc2.SEI, scale2, shape2))
sim.SP=rbvevd(n=1000,alpha= alpha,beta= beta,model='bilog',mar1=c(loc1.SP, scale1, shape1),mar2=c(loc2.SP, scale2, shape2))
sim.SWI=rbvevd(n=1000,alpha= alpha,beta= beta,model='bilog',mar1=c(loc1.SWI, scale1, shape1),mar2=c(loc2.SWI, scale2, shape2))
sim.WP=rbvevd(n=1000,alpha= alpha,beta= beta,model='bilog',mar1=c(loc1.WP, scale1, shape1),mar2=c(loc2.WP, scale2, shape2))

par(mar = c(4,4,1,1))
plot(sim.EP,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(sim.EP[,1], sim.EP[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
plot(data.for.bivariate,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(data.for.bivariate[data.for.bivariate[,'Basin']==1,][,1], data.for.bivariate[data.for.bivariate[,'Basin']==1,][,2], n=50),add=T,labcex = 1, col='blue',font=2)

plot(sim.NA,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(sim.NA[,1], sim.NA[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
plot(data.for.bivariate,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(data.for.bivariate[data.for.bivariate[,'Basin']==2,][,1], data.for.bivariate[data.for.bivariate[,'Basin']==2,][,2], n=50),add=T,labcex = 1, col='blue',font=2)

plot(sim.NI,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(sim.NI[,1], sim.NI[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
plot(data.for.bivariate,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(data.for.bivariate[data.for.bivariate[,'Basin']==3,][,1], data.for.bivariate[data.for.bivariate[,'Basin']==3,][,2], n=50),add=T,labcex = 1, col='blue',font=2)

plot(sim.SEI,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(sim.SEI[,1], sim.SEI[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
plot(data.for.bivariate,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(data.for.bivariate[data.for.bivariate[,'Basin']==4,][,1], data.for.bivariate[data.for.bivariate[,'Basin']==4,][,2], n=50),add=T,labcex = 1, col='blue',font=2)

plot(sim.SP,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(sim.SP[,1], sim.SP[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
plot(data.for.bivariate,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(data.for.bivariate[data.for.bivariate[,'Basin']==5,][,1], data.for.bivariate[data.for.bivariate[,'Basin']==5,][,2], n=50),add=T,labcex = 1, col='blue',font=2)

plot(sim.SWI,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(sim.SWI[,1], sim.SWI[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
plot(data.for.bivariate,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(data.for.bivariate[data.for.bivariate[,'Basin']==6,][,1], data.for.bivariate[data.for.bivariate[,'Basin']==6,][,2], n=50),add=T,labcex = 1, col='blue',font=2)

plot(sim.WP,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(sim.WP[,1], sim.WP[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
plot(data.for.bivariate,pch=15,cex=.25,col='grey',xlab='log(1013-MinCP)',ylab='log(MaxWS)',xlim=c(4,5.4),ylim=c(2.5,5.2))
contour(kde2d(data.for.bivariate[data.for.bivariate[,'Basin']==7,][,1], data.for.bivariate[data.for.bivariate[,'Basin']==7,][,2], n=50),add=T,labcex = 1, col='blue',font=2)



delta.ses<-c(Ints_SE_M1, final.model$std.err[c('scale1', 'shape1', 'alpha')],
	Ints_SE_M2,final.model$std.err[c('scale2', 'shape2', 'beta')])

bvevdmods<-data.frame(Estimate=ests,SE=delta.ses,Model_Type=1, Method='delta', 
	 Par=c(rep('beta01',7),'scale1','shape1','alpha',
	 	rep('beta02',7),'scale2','shape2','beta'),
	 Basin=c(ordered_basins, rep(NA,3), ordered_basins,rep(NA,3)), 
	 Model='BVEVD', Replacement='NotApp')
#################################################################
#################################################################
#Bootstrapping SE calculations and plotting
#################################################################
#################################################################
file.boot<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/Data/Derived Data/Model Output Data 2021 3b'

boots<- 1000
model.type1=c('1','2','3')[1]
m.percent1 = c(0.6666, 0.75, 1.00)[1]
m.percent2 = c(0.6666, 0.75, 1.00)[2]
m.percent3 = c(0.6666, 0.75, 1.00)[3]
#################################################################
#Bootstrapping
#1000 bootstraps then outputs to 3 csv files
#################################################################
ordered_basins= c('EaPa','NoAt','NoIn','SoEIn','SoPa','SoWIn','WePa')

#m/n = 0.6666, MooN parametric subsampling
system.time({
if(file.exists(file.path(file.boot,'bvevd.66.boot.all.csv'))){
dat.bv.66<-	read.csv(file.path(file.boot,'bvevd.66.boot.all.csv'), header=F, check.names = FALSE)[,-1]
colnames(dat.bv.66)<-names(bivariate.estimates)
bvevd.66.boot<- BVEVD_Boot(data= modeling_data, 
	data.file= dat.bv.66, m.percent = m.percent1, ordered_basins= ordered_basins)
} else {
	set.seed(1313123)
	pb66bvevd <-para.boot.bvevd.2021(boot.num= boots,model.object=final.model, 
	data= data.for.bivariate, cov= nloc.basin,m.percent= m.percent1)
	write.csv(rbind(pb66bvevd$boot.est.basins, pb66bvevd$boot.se.basins),
		file.path(file.boot,'bvevd.66.boot.basin.csv'))
	write.csv(rbind(pb66bvevd$boot.est.orig, pb66bvevd$boot.se.orig),
		file.path(file.boot,'bvevd.66.boot.orig.csv'))
	temp1<-data.frame(pb66bvevd$all.coeff.est)
	names(temp1)<-names(coef(pb66bvevd))
	write.csv(temp1,file.path(file.boot,'bvevd.66.boot.all.csv'))
}})

#m/n = 0.6666, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repbvevd.66.boot.all.csv'))){
dat.repbv.66<-read.csv(file.path(file.boot,'repbvevd.66.boot.all.csv'), header=F, check.names = FALSE)[,-1]
colnames(dat.repbv.66)<-names(bivariate.estimates)
repbvevd.66.boot<- BVEVD_Boot(data= modeling_data, 
	data.file= dat.repbv.66, m.percent = m.percent1, ordered_basins= ordered_basins)
} else {
	set.seed(1313123)
	rep66bvevd <-para.boot.bvevd.2021(boot.num= boots,model.object=final.model, 
	data= data.for.bivariate, cov= nloc.basin,m.percent= m.percent1, rep=TRUE)
	write.csv(rbind(rep66bvevd$boot.est.basins, rep66bvevd$boot.se.basins),
		file.path(file.boot,'repbvevd.66.boot.basin.csv'))
	write.csv(rbind(rep66bvevd$boot.est.orig, rep66bvevd$boot.se.orig),
		file.path(file.boot,'reprebvevd.66.boot.orig.csv'))
	temp1<-data.frame(rep66bvevd$all.coeff.est)
	names(temp1)<-names(coef(rep66bvevd))
	write.csv(temp1, file.path(file.boot,'repbvevd.66.boot.all.csv'))
}

#m/n = 0.75, MooN parametric subsampling

if(file.exists(file.path(file.boot,'bvevd.75.boot.all.csv'))){
dat.bv.75<-read.csv(file.path(file.boot,'bvevd.75.boot.all.csv'), header=F, check.names = FALSE)[,-1]
colnames(dat.bv.75)<-names(bivariate.estimates)
bvevd.75.boot<- BVEVD_Boot(data= modeling_data, 
	data.file= dat.bv.75, m.percent= m.percent2, ordered_basins= ordered_basins)
} else {
  set.seed(139140)
  pb75bvevd<-para.boot.bvevd(boot.num= boots,model.object=final.model, 
	data= data.for.bivariate, cov= nloc.basin,m.percent= m.percent2)
  write.csv(rbind(pb75bvevd$boot.est.basins, pb75bvevd$boot.se.basins),
  	file.path(file.boot,'bvevd.75.boot.basin.csv'))
  write.csv(rbind(pb75bvevd$boot.est.orig, pb75bvevd$boot.se.orig),
  	file.path(file.boot,'bvevd.75.boot.orig.csv'))
  temp1<-data.frame(pb75bvevd$all.coeff.est)
  names(temp1)<-names(coef(pb75bvevd))
  write.csv(temp1, file.path(file.boot,'bvevd.75.boot.all.csv'))
}


#m/n = 0.75, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repbvevd.75.boot.all.csv'))){
dat.repbv.75<-read.csv(file.path(file.boot,'repbvevd.75.boot.all.csv'), header=F, check.names = FALSE)[,-1]
colnames(dat.repbv.75)<-names(bivariate.estimates)
repbvevd.75.boot<- BVEVD_Boot(data= modeling_data, 
  	data.file= dat.repbv.75, m.percent= m.percent2, ordered_basins= ordered_basins)
} else {
  set.seed(139140)
  reppb75bvevd<-para.boot.bvevd(boot.num= boots,model.object=final.model, 
	data= data.for.bivariate, cov= nloc.basin,m.percent= m.percent2,rep=TRUE)
  write.csv(rbind(reppb75bvevd$boot.est.basins, reppb75bvevd$boot.se.basins),
  	file.path(file.boot,'repbvevd.75.boot.basin.csv'))
  write.csv(rbind(reppb75bvevd$boot.est.orig, reppb75bvevd$boot.se.orig),
  	file.path(file.boot,'repbvevd.75.boot.orig.csv'))
  temp1<-data.frame(reppb75bvevd$all.coeff.est)
  names(temp1)<-names(coef(reppb75bvevd))
  write.csv(temp1, file.path(file.boot,'repbvevd.75.boot.all.csv'))
}

#m/n = 1.00, MooN parametric bootstrapping
if(file.exists(file.path(file.boot,'repbvevd.boot.all.csv'))){
dat.repbv<-read.csv(file.path(file.boot,'repbvevd.boot.all.csv'), header=F, check.names = FALSE)[,-1]
colnames(dat.repbv)<-names(bivariate.estimates)
repbvevd.boot<-  BVEVD_Boot(data= modeling_data, 
  	data.file= dat.repbv, m.percent = m.percent3, ordered_basins= ordered_basins)
} else {
  set.seed(131313)
  pbbvevd.wrep<-para.boot.bvevd(boot.num= boots,model.object=final.model, 
	data= data.for.bivariate, cov= nloc.basin,m.percent= m.percent3,rep=TRUE)
  write.csv(rbind(pbbvevd.wrep$boot.est.basins,pbbvevd.wrep$boot.se.basins),file.path(file.boot,'repbvevd.boot.basin.csv'))
  write.csv(rbind(pbbvevd.wrep$boot.est.orig, pbbvevd.wrep$boot.se.orig),file.path(file.boot,'repbvevd.boot.orig.csv'))
  temp1.rep <-data.frame(pbbvevd.wrep$all.coeff.est)
  names(temp1.rep)<-names(coef(pbbvevd.wrep))
  write.csv(temp1.rep,file.path(file.boot,'repbvevd.boot.all.csv'))
}

SE_boots<-rbind(
	data.frame(SE= bvevd.66.boot$Boot.SE, Model_Type=1, Method='Par_0.66', 
	 Par=c(rep('beta01',7),'scale1','shape1','alpha',
	 	rep('beta02',7),'scale2','shape2','beta'),
	 Basin=c(ordered_basins, rep(NA,3), ordered_basins,rep(NA,3)), 
	 Model='BVEVD', Replacement='No'),
	 
	data.frame(SE=repbvevd.66.boot$Boot.SE, Model_Type=1, Method='Par_0.66', 
	 Par=c(rep('beta01',7),'scale1','shape1','alpha',
	 	rep('beta02',7),'scale2','shape2','beta'),	 
	 Basin=c(ordered_basins, rep(NA,3), ordered_basins,rep(NA,3)), 
	 Model='BVEVD', Replacement='Yes'),
	 	 
	data.frame(SE=bvevd.75.boot$Boot.SE, Model_Type=1, Method='Par_0.75', 
	 Par=c(rep('beta01',7),'scale1','shape1','alpha',
	 	rep('beta02',7),'scale2','shape2','beta'),
	 Basin=c(ordered_basins, rep(NA,3), ordered_basins,rep(NA,3)), 
	 Model='BVEVD', Replacement='No'),
	 	 
	data.frame(SE=repbvevd.75.boot$Boot.SE, Model_Type=1, Method='Par_0.75', 
	 Par=c(rep('beta01',7),'scale1','shape1','alpha',
	 	rep('beta02',7),'scale2','shape2','beta'),
	 Basin=c(ordered_basins, rep(NA,3), ordered_basins,rep(NA,3)), 
	 Model='BVEVD', Replacement='Yes'),
	 
	data.frame(SE=repbvevd.boot$Boot.SE, Model_Type=1, Method='Par_1.00', 
	 Par=c(rep('beta01',7),'scale1','shape1','alpha',
	 	rep('beta02',7),'scale2','shape2','beta'),
	 Basin=c(ordered_basins, rep(NA,3), ordered_basins,rep(NA,3)), 
	 Model='BVEVD', Replacement='Yes')
)
#################################################################
#################################################################
#Calculating the Basin estimates
#################################################################
#################################################################
ordered_basins<-c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn')

#################################################################
#################################################################
#Plotting Model Coefficients with Standard Errors
#################################################################
#################################################################


SE_boot<-merge(SE_boots, bvevdmods[,c('Estimate','Basin','Model_Type','Par','Model')],
	by=c('Basin','Model_Type','Par','Model'),all=T)
SE_BVEVD<-rbind(bvevdmods, SE_boot)
write.csv(SE_BVEVD,file.path(file.boot,'AllSEGraphingBVEVD.csv'), row.names = FALSE)

se_graph <-SE_BVEVD
se_graph$lo95<-se_graph$Estimate-2* se_graph$SE
se_graph$hi95<-se_graph$Estimate+2* se_graph$SE

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
#####BVEVD Model 1 Coefficients +2SE Plotting
#################################################################
#################################################################

#Beta01
ggplot(se_graph[se_graph$Model_Type==1 &se_graph$Model=='BVEVD' & se_graph$Par =='beta01',],aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty=Replacement))+ geom_linerange(size=1,position= position_dodge(.5))+xlab('Basin')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+#+geom_hline(yintercept =0,lty=5, alpha=0.5)+ 
theme(legend.position = 'none',axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))

#Beta02
ggplot(se_graph[se_graph$Model_Type==1 &se_graph$Model=='BVEVD' & se_graph$Par =='beta02',],aes(x=Basin, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty=Replacement))+ geom_linerange(size=1,position= position_dodge(.5))+xlab('Basin')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+#+geom_hline(yintercept =0,lty=5, alpha=0.5)+ 
theme(legend.position = 'none',axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))

#Margin1: scale, shape, alpha 
ggplot(se_graph[se_graph$Model_Type==1 &se_graph$Model=='BVEVD' & se_graph$Par %in%c('scale1','shape1','alpha'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
scale_x_discrete(label=c(expression(paste(alpha)),expression(paste(sigma[1])),expression(paste(xi[1]))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))

#Margin2: scale, shape, alpha 
ggplot(se_graph[se_graph$Model_Type==1 &se_graph$Model=='BVEVD' & se_graph$Par %in%c('scale2','shape2','beta'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
scale_x_discrete(label=c(expression(paste(beta)),expression(paste(sigma[2])),expression(paste(xi[2]))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))