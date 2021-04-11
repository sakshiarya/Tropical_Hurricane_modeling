##########################################################################
##########################################################################
##########################################################################
# Title: Bivariate GEV Modeling  
# Author: Lindsey Dietz
# Updated/Edited by: Sakshi Arya
# Last Updated: 3/8/2021
#
# Goal: Fits bivariate GEV models for damage and (log(maxWS), log(1013-minCP)) 
#
# Note data files should be in working directory.
##########################################################################
##########################################################################
##########################################################################

rm(list=ls())

library(evd)
library(MASS)
library(xtable)
library(car)
#library(KWPQuant)
#install.packages("KWPQuant")

filep<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/'
source(file.path(filep,'Helper Functions/para.boot.bvevd.damage.R'))
source(file.path(filep,'Helper Functions/dam.boot.R'))
storms<-read.csv('Data/Derived Data 21/Categorized_Storm_1960_2019.csv')

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

names(storms)
storm_damage<-storms[storms$minCP< Inf & storms$CURRENT.DAMAGE.2021>0 & !is.na(storms$CURRENT.DAMAGE.2021),]

tri_var<-storm_damage[,c('minCP','maxWS','CURRENT.DAMAGE.2021')]
tri_var$logdamage<-log(tri_var$CURRENT.DAMAGE.2021)
tri_var$logmaxWS<-log(tri_var$maxWS)
tri_var$logminCP<-log(1013-tri_var$minCP)
tri_var$avgLat <- storm_damage$avgLat


dam.ws<-fbvevd(tri_var[,c(5,4)], model = "bilog")
dam.ws.log<-fbvevd(tri_var[,c(5,4)], model = "log")
anova(dam.ws,dam.ws.log) # No need for bilog

dam.cp<-fbvevd(tri_var[,c(6,4)], model = "bilog")
dam.cp.b<-fbvevd(tri_var[,c(6,4)], model = "log")
anova(dam.cp, dam.cp.b)

ws.cp<-fbvevd(tri_var[,c(5,6)], model = "bilog")
ws.cp.b <-fbvevd(tri_var[,c(5,6)], model = "log")
anova(ws.cp, ws.cp.b)

damagesmodel1<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6]), model = "bilog",std.err=F)
damagesmodel1.log<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6]), model = "log",std.err=F)
anova(damagesmodel1, damagesmodel1.log)

anova(damagesmodel1, dam.ws)
anova(damagesmodel1.log, dam.ws.log)

damagesmodel2<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6],(tri_var[,6])^2), model = "bilog",std.err=F)
damagesmodel2.log<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6],(tri_var[,6])^2), model = "log",std.err=F)
anova(damagesmodel2, damagesmodel2.log)
anova(damagesmodel2, dam.ws)
anova(damagesmodel2, damagesmodel1)

damagesmodel3<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), model = "bilog",std.err=F)
damagesmodel3.log<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), model = "log",std.err=F)
anova(damagesmodel3, damagesmodel3.log)
anova(damagesmodel3, dam.ws)

damagesmodel4<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6], tri_var[,7]), model = "bilog",std.err=F)
damagesmodel4.log<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6], tri_var[,7]), model = "log",std.err=F)
anova(damagesmodel4, damagesmodel4.log)
anova(damagesmodel4, dam.ws)

damagesmodel5<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6]), nsloc2= data.frame(tri_var[,6]), model = "bilog",std.err=F)
damagesmodel5.log<-fbvevd(tri_var[,c(5,4)],nsloc1= data.frame(tri_var[,6]),nsloc2= data.frame(tri_var[,6]), model = "log",std.err=F)
anova(damagesmodel5, damagesmodel5.log)
anova(damagesmodel5, dam.ws)



AIC(dam.ws);AIC(damagesmodel1); AIC(damagesmodel2); AIC(damagesmodel3); AIC(damagesmodel4)
AIC(dam.ws,k=log(dim(tri_var)[1])); AIC(damagesmodel1,k=log(dim(tri_var)[1])); AIC(damagesmodel2,k=log(dim(tri_var)[1])); AIC(damagesmodel3,k=log(dim(tri_var)[1])); AIC(damagesmodel4,k=log(dim(tri_var)[1]))

#damagesmodel3 <- fbvevd(tri_var[,c(6,4)], nsloc1 = data.frame(tri_var[,5],(tri_var[,5])^2), model = "bilog",std.err=F)
###########################################
###########################################
storm_damage.2008<-storms[storms$minCP< Inf & storms$CURRENT.DAMAGE.2021>0 & !is.na(storms$CURRENT.DAMAGE.2021) & storms$Year<2009,]
test.set<-storms[storms$minCP< Inf & storms$CURRENT.DAMAGE.2021>0 & !is.na(storms$CURRENT.DAMAGE.2021) & storms$Year>=2009,]
tri_var_2008<-storm_damage.2008[,c('minCP','maxWS','CURRENT.DAMAGE.2021')]
tri_var_2008 $logdamage<-log(tri_var_2008 $CURRENT.DAMAGE.2021)
tri_var_2008 $logmaxWS<-log(tri_var_2008 $maxWS)
tri_var_2008 $logminCP<-log(1013-tri_var_2008 $minCP)
tri_var_2008$avgLat <- storm_damage.2008$avgLat

damagesmodel1.2008<-fbvevd(tri_var_2008[,c(5,4)],nsloc1= data.frame(tri_var_2008[,6]), model = "bilog",std.err=F)

damagesmodel2.2008<-fbvevd(tri_var_2008[,c(5,4)],nsloc1= data.frame(tri_var_2008[,6],(tri_var_2008[,6])^2), model = "bilog",std.err=F)

damagesmodel3.2008<-fbvevd(tri_var_2008[,c(5,4)],nsloc1= data.frame(tri_var_2008[,6],(tri_var_2008[,6])^2, tri_var_2008[,7]), model = "bilog",std.err=F)

damagesmodel4.2008<-fbvevd(tri_var_2008[,c(5,4)],nsloc1= data.frame(tri_var_2008[,6], tri_var_2008[,7]), model = "bilog",std.err=F)

mahal.pval<-function(model.type=1,dat.cp=test.set[1,] ,num.p.val= 500, num.sim.bvevd=999){
  pvalue <-vector()
  for(i in 1:num.p.val){
    
    if(model.type==1){
      sim.dam.ws<-rbvevd(num.sim.bvevd,
                         alpha=damagesmodel1.2008$estimate['alpha'], 
                         beta= damagesmodel1.2008$estimate['beta'],model='bilog',
                         mar1=c(damagesmodel1.2008$estimate[1]+ 
                                  damagesmodel1.2008$estimate[2]*log(1013-dat.cp$minCP), 
                                damagesmodel1.2008$estimate[3], damagesmodel1.2008$estimate[4]),
                         mar2= damagesmodel1.2008$estimate[5:7])
    } else if (model.type==2){
      sim.dam.ws<-rbvevd(num.sim.bvevd,
                         alpha=damagesmodel2.2008$estimate['alpha'], 
                         beta= damagesmodel2.2008$estimate['beta'],model='bilog',
                         mar1=c(damagesmodel2.2008$estimate[1]+ 
                                  damagesmodel2.2008$estimate[2]*log(1013-dat.cp$minCP)+ 
                                  damagesmodel2.2008$estimate[3]*(log(1013-dat.cp$minCP))^2, 
                                damagesmodel2.2008$estimate[4], damagesmodel2.2008$estimate[5]),
                         mar2= damagesmodel2.2008$estimate[6:8])
    }else if (model.type==3){
      sim.dam.ws<-rbvevd(num.sim.bvevd,
                         alpha=damagesmodel3.2008$estimate['alpha'], 
                         beta= damagesmodel3.2008$estimate['beta'],model='bilog',
                         mar1=c(damagesmodel3.2008$estimate[1]+ 
                                  damagesmodel3.2008$estimate[2]*log(1013-dat.cp$minCP)+ 
                                  damagesmodel3.2008$estimate[3]*(log(1013-dat.cp$minCP))^2 +
                                  damagesmodel3.2008$estimate[4]*dat.cp$avgLat, 
                                damagesmodel3.2008$estimate[5], damagesmodel3.2008$estimate[6]),
                         mar2= damagesmodel3.2008$estimate[7:9])
    }else if (model.type==4){
      sim.dam.ws<-rbvevd(num.sim.bvevd,
                         alpha=damagesmodel4.2008$estimate['alpha'], 
                         beta= damagesmodel4.2008$estimate['beta'],model='bilog',
                         mar1=c(damagesmodel4.2008$estimate[1]+ 
                                  damagesmodel4.2008$estimate[2]*log(1013-dat.cp$minCP)+ 
                                  damagesmodel4.2008$estimate[3]*dat.cp$avgLat, 
                                damagesmodel4.2008$estimate[4], damagesmodel4.2008$estimate[5]),
                         mar2= damagesmodel4.2008$estimate[6:8])
    }
    cov.dam.ws<-var(sim.dam.ws)
    mean.dam.ws<-colMeans(sim.dam.ws)	
    median.dam.ws<-apply(sim.dam.ws,2,median)	
    #		plot(sim.dam.ws,pch=16,cex=0.4,xlab='Simulated log(MaxWS)',ylab='Simulated log(Damage)',col='grey')
    plot(sim.dam.ws,pch=16,cex=0.9,xlab='Simulated log(MaxWS)',ylab='Simulated log(Damage)',col='grey',cex.axis=1.8, cex.lab= 1.6)
    contour(kde2d(x=sim.dam.ws[,1], y=sim.dam.ws[,2], n=50),add=T,labcex = 1, col='darkgreen',font=2)
    points(mean.dam.ws[1], mean.dam.ws[2],col='blue',pch=16,cex=1.9)
    points(log(dat.cp$maxWS),log(dat.cp$CURRENT.DAMAGE.2021),col='red',pch=15,cex=1.9)
    
    location.in.mah<-order(c(mahalanobis(x=c(log(dat.cp$maxWS), 
                                             log(dat.cp$CURRENT.DAMAGE.2021)),center= mean.dam.ws, 
                                         cov= cov.dam.ws), 
                             mahalanobis(x= sim.dam.ws,center= mean.dam.ws,cov= cov.dam.ws)))[1]
    pvalue[i]<-2*min(location.in.mah/(num.sim.bvevd+1), (num.sim.bvevd+1-location.in.mah)/(num.sim.bvevd+1))
  }
  return(pvalue)
}

#par(mar=c(2,3,0,0)+0.1, mgp =c(2,1,0))#sets margins of plotting area
par(mar=c(5,6,1,3)+0.1)#sets margins of plotting area

set.seed(1093728)
summary(alex1<-mahal.pval(model.type=1, dat.cp=test.set[1,] ,num.p.val= 200, num.sim.bvevd=999))
summary(alex2<- mahal.pval(model.type=2, dat.cp=test.set[1,] ,num.p.val= 200, num.sim.bvevd=999))
summary(alex3<- mahal.pval(model.type=3, dat.cp=test.set[1,] ,num.p.val= 200, num.sim.bvevd=999))
summary(alex4<- mahal.pval(model.type=4, dat.cp=test.set[1,] ,num.p.val= 200, num.sim.bvevd=999))



set.seed(347285)
summary(earl1<-mahal.pval(model.type=1, dat.cp=test.set[2,] ,num.p.val= 200, num.sim.bvevd=999))
summary(earl2<-mahal.pval(model.type=2, dat.cp=test.set[2,] ,num.p.val= 200, num.sim.bvevd=999))
summary(earl3<-mahal.pval(model.type=3, dat.cp=test.set[2,] ,num.p.val= 200, num.sim.bvevd=999))
summary(earl4<-mahal.pval(model.type=4, dat.cp=test.set[2,] ,num.p.val= 200, num.sim.bvevd=999))


set.seed(342075987)
summary(hermine1<-mahal.pval(model.type=1, dat.cp=test.set[3,] ,num.p.val= 200, num.sim.bvevd=999))
summary(hermine2<- mahal.pval(model.type=2, dat.cp=test.set[3,] ,num.p.val= 200, num.sim.bvevd=999))
summary(hermine3<-mahal.pval(model.type=3, dat.cp=test.set[3,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=0.9)
summary(hermine4<- mahal.pval(model.type=4, dat.cp=test.set[3,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=0.9)


set.seed(43687)
summary(irene1<-mahal.pval(model.type=1, dat.cp=test.set[4,] ,num.p.val= 1, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F, cex = 2, pt.cex=2.5)
summary(irene2<- mahal.pval(model.type=2, dat.cp=test.set[4,] ,num.p.val= 200, num.sim.bvevd=999))
summary(irene3<- mahal.pval(model.type=3, dat.cp=test.set[4,] ,num.p.val= 200, num.sim.bvevd=999))
summary(irene4<- mahal.pval(model.type=4, dat.cp=test.set[4,] ,num.p.val= 200, num.sim.bvevd=999))


set.seed(687)
summary(isaac1<-mahal.pval(model.type=1, dat.cp=test.set[5,] ,num.p.val= 1, num.sim.bvevd=999))
#legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=0.9)
summary(isaac2<- mahal.pval(model.type=2, dat.cp=test.set[5,] ,num.p.val= 1, num.sim.bvevd=999))
#legend("bottomleft", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=0.9)
summary(isaac3<- mahal.pval(model.type=3, dat.cp=test.set[5,] ,num.p.val= 200, num.sim.bvevd=999))
#legend("bottomleft", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=0.9)
summary(isaac4<- mahal.pval(model.type=4, dat.cp=test.set[5,] ,num.p.val= 200, num.sim.bvevd=999))
#legend("bottomleft", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=0.9)


set.seed(67)
summary(sandy1<- mahal.pval(model.type=1, dat.cp=test.set[6,] ,num.p.val= 200, num.sim.bvevd=999))
summary(sandy2<- mahal.pval(model.type=2, dat.cp=test.set[6,] ,num.p.val= 200, num.sim.bvevd=999))


set.seed(67)
summary(sandy3<- mahal.pval(model.type=3, dat.cp=test.set[6,] ,num.p.val= 200, num.sim.bvevd=999))
#legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1.5)
summary(sandy4<- mahal.pval(model.type=4, dat.cp=test.set[6,] ,num.p.val= 200, num.sim.bvevd=999))
#legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1.5)


set.seed(687)
summary(arthur1<-mahal.pval(model.type=1, dat.cp=test.set[7,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)
summary(arthur2<- mahal.pval(model.type=2, dat.cp=test.set[7,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)



set.seed(687)
summary(Hermine21<-mahal.pval(model.type=1, dat.cp=test.set[8,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)
summary(Hermine22<- mahal.pval(model.type=2, dat.cp=test.set[8,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)



set.seed(687)
summary(Matthew1<-mahal.pval(model.type=1, dat.cp=test.set[9,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)
summary(Matthew2<- mahal.pval(model.type=2, dat.cp=test.set[9,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)


set.seed(687)
summary(Harvey1<-mahal.pval(model.type=1, dat.cp=test.set[10,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)
summary(Harvey2<- mahal.pval(model.type=2, dat.cp=test.set[10,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=0.8)


set.seed(687)
summary(Irma1<-mahal.pval(model.type=1, dat.cp=test.set[11,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)
summary(Irma2<- mahal.pval(model.type=2, dat.cp=test.set[11,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)


set.seed(687)
summary(Nate1<-mahal.pval(model.type=1, dat.cp=test.set[12,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)
summary(Nate2<- mahal.pval(model.type=2, dat.cp=test.set[12,] ,num.p.val= 200, num.sim.bvevd=999))
legend("bottomright", inset=.01,  c("Simulated Mean","Actual"),pch=c(16,15),col=c('blue','red'), horiz=F,cex=1)



xtable(
  matrix(c(quantile(alex1,0.025),quantile(alex1,0.975),mean(alex1<=0.05)*100,
           quantile(alex2,0.025),quantile(alex2,0.975),mean(alex2<=0.05)*100,
           quantile(earl1,0.025),quantile(earl1,0.975),mean(earl1<=0.05)*100,
           quantile(earl2,0.025),quantile(earl2,0.975),mean(earl2<=0.05)*100,
           quantile(hermine1,0.025),quantile(hermine1,0.975),mean(hermine1 <=0.05)*100,
           quantile(hermine2,0.025),quantile(hermine2,0.975),mean(hermine2 <=0.05)*100,
           quantile(irene1,0.025),quantile(irene1,0.975),mean(irene1 <=0.05)*100,
           quantile(irene2,0.025),quantile(irene2,0.975),mean(irene2 <=0.05)*100,
           quantile(isaac1,0.025),quantile(isaac1,0.975),mean(isaac1 <=0.05)*100,
           quantile(isaac2,0.025),quantile(isaac2,0.975),mean(isaac2 <=0.05)*100,
           quantile(sandy1,0.025),quantile(sandy1,0.975),mean(sandy1 <=0.05)*100,
           quantile(sandy2,0.025),quantile(sandy2,0.975),mean(sandy2 <=0.05)*100,
           quantile(arthur1,0.025),quantile(arthur1,0.975),mean(arthur1 <=0.05)*100,
           quantile(arthur2,0.025),quantile(arthur2,0.975),mean(arthur2 <=0.05)*100,
           quantile(Hermine21,0.025),quantile(Hermine21,0.975),mean(Hermine21 <=0.05)*100,
           quantile(Hermine22,0.025),quantile(Hermine22,0.975),mean(Hermine22 <=0.05)*100,
           quantile(Matthew1,0.025),quantile(Matthew1,0.975),mean(Matthew1 <=0.05)*100,
           quantile(Matthew2,0.025),quantile(Matthew2,0.975),mean(Matthew2 <=0.05)*100,
           quantile(Harvey1,0.025),quantile(Harvey1,0.975),mean(Harvey1 <=0.05)*100,
           quantile(Harvey2,0.025),quantile(Harvey2,0.975),mean(Harvey2 <=0.05)*100,
           quantile(Irma1,0.025),quantile(Irma1,0.975),mean(Irma1 <=0.05)*100,
           quantile(Irma2,0.025),quantile(Irma2,0.975),mean(Irma2 <=0.05)*100,
           quantile(Nate1,0.025),quantile(Nate1,0.975),mean(Nate1 <=0.05)*100,
           quantile(Nate2,0.025),quantile(Nate2,0.975),mean(Nate2 <=0.05)*100),ncol=6,byrow=T),digits=3)

###########################################
###########################################
dam.ws$estimate
damagesmodel1$estimate
damagesmodel2$estimate
damagesmodel3$estimate
damagesmodel4$estimate


results <- matrix(NA, ncol = length(names(damagesmodel2$estimate)), nrow = 3)
colnames(results) <- names(damagesmodel2$estimate)
results[1, names(dam.ws$estimate)] <- dam.ws$estimate
results[2, names(damagesmodel1$estimate)] <- damagesmodel1$estimate
results[3, names(damagesmodel2$estimate)] <- damagesmodel2$estimate
xtable(results, digits = 3)

anova(damagesmodel2, damagesmodel1)
anova(damagesmodel2, dam.ws)

par(mfrow=c(2,1))

mu1<-damagesmodel1$estimate['loc1']+damagesmodel1$estimate['loc1tri_var...6.']* tri_var$logminCP

scale1<-damagesmodel1$estimate['scale1']
shape1<-damagesmodel1$estimate['shape1']
trans1<-(1+(tri_var$logmaxWS- mu1)*(shape1/scale1))^(-1/shape1)
trans2<-(tri_var$logmaxWS- mu1)*(shape1/scale1)
trans3<-(tri_var$logmaxWS- mu1+ damagesmodel1$estimate['loc1'])

damagesmodeltest<-fbvevd(cbind(trans3,tri_var$logdamage), model = "bilog",std.err=F)
plot(density(rbvevd(n=10000,model='bilog',alpha= damagesmodeltest$estimate['alpha'], beta= damagesmodeltest$estimate['beta'],mar1=c(damagesmodeltest$estimate['loc1'], damagesmodeltest$estimate['scale1'], damagesmodeltest$estimate['shape1']),mar2=c(damagesmodeltest$estimate['loc2'], damagesmodeltest$estimate['scale2'], damagesmodeltest$estimate['shape2']))[,1]),main='',ylim=c(0,3))
lines(density(trans3),lty=2)


par(mar=c(3,3,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
pdf("margws1_new.pdf")
plot(damagesmodel1,mar=1,which=3,ylim=c(0,3),xlab= expression(paste(t[1],'(MaxWS)',sep='')),main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()

pdf("margdam1_new.pdf")
plot(damagesmodel1,mar=2,which=3,ylim=c(0,0.15),xlab='log(Damage)',main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()

AIC(damagesmodel1)
AIC(damagesmodel1,k=log(dim(tri_var)[1]))

par(mfrow=c(2,1))
pdf("margws0_new.pdf")
plot(dam.ws,mar=1,which=3,ylim=c(0,1.5),xlab='log(MaxWS)',main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()

pdf("margdam0_new.pdf")
plot(dam.ws,mar=2,which=3,ylim=c(0,0.15),xlab='log(Damage)',main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()
AIC(dam.ws)
AIC(dam.ws,k=log(dim(tri_var)[1]))


par(mfrow=c(2,1))
pdf("margws2_new.pdf")
plot(damagesmodel2,mar=1,which=3,ylim=c(0,4),xlab= expression(paste(t[2],'(MaxWS)',sep='')),main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()

pdf("margdam2_new.pdf")
plot(damagesmodel2,mar=2,which=3,ylim=c(0,0.15),xlab='log(Damage)',main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()
AIC(damagesmodel2)
AIC(damagesmodel2,k=log(dim(tri_var)[1]))

par(mfrow=c(2,1))
pdf("margws3_new.pdf")
plot(damagesmodel3,mar=1,which=3,ylim=c(0,4),xlab= expression(paste(t[3],'(MaxWS)',sep='')),main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()

pdf("margdam3_new.pdf")
plot(damagesmodel3,mar=2,which=3,ylim=c(0,0.15),xlab='log(Damage)',main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()


par(mfrow=c(2,1))
pdf("margws4_new.pdf")
plot(damagesmodel4,mar=1,which=3,ylim=c(0,4),xlab= expression(paste(t[4],'(MaxWS)',sep='')),main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()


pdf("margdam4_new.pdf")
plot(damagesmodel4,mar=2,which=3,ylim=c(0,0.15),xlab='log(Damage)',main='')
legend('topright',legend =c('Model','Empirical'),lty=c(1,2))
dev.off()


## Bootstrapping and subsamling SEs
file.boot <- "/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/Data/Derived Data 21/Model Output Data"


set.seed(987987)
a1rep<-para.boot.bvevd.damage(boot.num=1000,model.type=0,model.object= dam.ws,data= tri_var[,c(5,4)], cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.6666,rep=TRUE)
write.csv(a1rep$all.boots,file.path(file.boot,'damage.bvevd.66.rep.boots.csv'))
a1rep<-read.csv(file.path(file.boot,'damage.bvevd.66.rep.boots.csv'),header=TRUE)[-1]
a1rep_SE <- dam.boot(file = a1rep,boot.num = 1000, model.type = c(0,1,2)[1], data= tri_var[,c(5,4)], m.percent=0.66)

set.seed(987987)
a1<-para.boot.bvevd.damage(boot.num=1000,model.type=0,model.object= dam.ws,data= tri_var[,c(5,4)], cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.6666,rep=FALSE)
write.csv(a1$all.boots,'damage.bvevd.66.boots.csv')
a1<-read.csv(file.path(file.boot,'damage.bvevd.66.boots.csv'),header=TRUE)[-1]
a1_SE <- dam.boot(file = a1,boot.num = 1000, model.type = c(0,1,2)[1], data= tri_var[,c(5,4)], m.percent=0.6666)

set.seed(1432412)
a2rep<-para.boot.bvevd.damage(boot.num=1000,model.type=0,model.object= dam.ws,data= tri_var[,c(5,4)], cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.75,rep=TRUE)
write.csv(a2rep$all.boots,'damage.bvevd.75.rep.boots.csv')
a2rep<-read.csv(file.path(file.boot,'damage.bvevd.75.rep.boots.csv'),header=TRUE)[-1]
a2rep_SE <- dam.boot(file = a2rep,boot.num = 1000, model.type = c(0,1,2)[1], data= tri_var[,c(5,4)], m.percent=0.75)


set.seed(1432412)
a2<-para.boot.bvevd.damage(boot.num=1000,model.type=0,model.object= dam.ws,data= tri_var[,c(5,4)], cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.75,rep=FALSE)
write.csv(a2$all.boots,'damage.bvevd.75.boots.csv')
a2<-read.csv(file.path(file.boot,'damage.bvevd.75.boots.csv'),header=TRUE)[-1]
a2_SE <- dam.boot(file = a2,boot.num = 1000, model.type = c(0,1,2)[1], data= tri_var[,c(5,4)], m.percent=0.75)

set.seed(102123)
a3rep<-para.boot.bvevd.damage(boot.num=1000,model.type=0,model.object= dam.ws,data= tri_var[,c(5,4)], cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=1,rep=TRUE)
write.csv(a3rep$all.boots,'damage.bvevd.1.rep.boots.csv')
a3rep<-read.csv(file.path(file.boot,'damage.bvevd.1.rep.boots.csv'),header=TRUE)[-1]
a3rep_SE <- dam.boot(file = a3rep,boot.num = 1000, model.type = c(0,1,2)[1], data= tri_var[,c(5,4)], m.percent=1)



set.seed(23476)
b0rep<-para.boot.bvevd.damage(boot.num=1000,model.type=1,model.object= damagesmodel1,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.66,rep=TRUE)
write.csv(b0rep$all.boots,'mod1.damage.bvevd.66.rep.boots.csv')
b0rep<-read.csv(file.path(file.boot,'mod1.damage.bvevd.66.rep.boots.csv'),header=TRUE)[-1]
b0rep_SE <- dam.boot(file = b0rep,boot.num = 1000, model.type = c(0,1,2)[2], data= tri_var[,c(5,4)], m.percent=0.66)
b0rep_SE$boot.se.basins

set.seed(234827439)
b0<-para.boot.bvevd.damage(boot.num=1000,model.type=1,model.object= damagesmodel1,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.66,rep=FALSE)
write.csv(b0$all.boots,'mod1.damage.bvevd.66.boots.csv')
b0<-read.csv(file.path(file.boot,'mod1.damage.bvevd.66.boots.csv'),header=TRUE)[-1]
b0_SE <- dam.boot(file = b0,boot.num = 1000, model.type = c(0,1,2)[2], data= tri_var[,c(5,4)], m.percent=0.66)
b0_SE$boot.se.basins


set.seed(432987)
 b1rep<-para.boot.bvevd.damage(boot.num=1000,model.type=1,model.object= damagesmodel1,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.75,rep=TRUE)
 write.csv(b1rep$all.boots,'mod1.damage.bvevd.75.rep.boots.csv')
 b1rep<-read.csv(file.path(file.boot,'mod1.damage.bvevd.75.rep.boots.csv'),header=TRUE)[-1]
 b1rep_SE <- dam.boot(file = b1rep,boot.num = 1000, model.type = c(0,1,2)[2], data= tri_var[,c(5,4)], m.percent=0.75)
 b1rep_SE$boot.se.basins


 set.seed(243908)
 b1<-para.boot.bvevd.damage(boot.num=1000,model.type=1,model.object= damagesmodel1,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.75,rep= FALSE)
 write.csv(b1$all.boots,'mod1.damage.bvevd.75.boots.csv')
 b1<-read.csv(file.path(file.boot,'mod1.damage.bvevd.75.boots.csv'),header=TRUE)[-1]
 b1_SE <- dam.boot(file = b1,boot.num = 1000, model.type = c(0,1,2)[2], data= tri_var[,c(5,4)], m.percent=0.75)
 b1_SE$boot.se.basins
 

 set.seed(453981)
 b2rep<-para.boot.bvevd.damage(boot.num=1000,model.type=1,model.object= damagesmodel1,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=1,rep=TRUE)
 write.csv(b2rep$all.boots,'mod1.damage.bvevd.1.rep.boots.csv')
 b2rep<-read.csv(file.path(file.boot,'mod1.damage.bvevd.1.rep.boots.csv'),header=TRUE)[-1]
 b2rep_SE <- dam.boot(file = b2rep,boot.num = 1000, model.type = c(0,1,2)[2], data= tri_var[,c(5,4)], m.percent=1)
 b2rep_SE$boot.se.basins


 set.seed(2343476)
 c0rep<-para.boot.bvevd.damage(boot.num=1000,model.type=2,model.object= damagesmodel2,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.66,rep=TRUE)
 write.csv(c0rep$all.boots,'mod2.damage.bvevd.66.rep.boots.csv')
 c0rep<-read.csv(file.path(file.boot,'mod2.damage.bvevd.66.rep.boots.csv'),header=TRUE)[-1]
 c0rep_SE <- dam.boot(file = c0rep,boot.num = 1000, model.type = c(0,1,2)[3], data= tri_var[,c(5,4)], m.percent=0.66)
 c0rep_SE$boot.se.basins


 set.seed(2342339)
 c0<-para.boot.bvevd.damage(boot.num=1000,model.type=2,model.object= damagesmodel2,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.66,rep=FALSE)
 write.csv(c0$all.boots,'mod2.damage.bvevd.66.boots.csv')
 c0<-read.csv(file.path(file.boot,'mod2.damage.bvevd.66.boots.csv'),header=TRUE)[-1]
 c0_SE <- dam.boot(file = c0,boot.num = 1000, model.type = c(0,1,2)[3], data= tri_var[,c(5,4)], m.percent=0.66)
 c0_SE$boot.se.basins

 set.seed(243234)
 c1rep<-para.boot.bvevd.damage(boot.num=1000,model.type=2,model.object= damagesmodel2,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.75,rep=TRUE)
 write.csv(c1rep$all.boots,'mod2.damage.bvevd.75.rep.boots.csv')
 c1rep<-read.csv(file.path(file.boot,'mod2.damage.bvevd.75.rep.boots.csv'),header=TRUE)[-1]
 c1rep_SE <- dam.boot(file = c1rep,boot.num = 1000, model.type = c(0,1,2)[3], data= tri_var[,c(5,4)], m.percent=0.75)
 c1rep_SE$boot.se.basins

 set.seed(435665)
 c1<-para.boot.bvevd.damage(boot.num=1000,model.type=2,model.object= damagesmodel2,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=0.75,rep= FALSE)
 write.csv(c1$all.boots,'mod2.damage.bvevd.75.boots.csv')
 c1<-read.csv(file.path(file.boot,'mod2.damage.bvevd.75.boots.csv'),header=TRUE)[-1]
 c1_SE <- dam.boot(file = c1,boot.num = 1000, model.type = c(0,1,2)[3], data= tri_var[,c(5,4)], m.percent=0.75)
 c1_SE$boot.se.basins

 set.seed(4564879)
 c2rep<-para.boot.bvevd.damage(boot.num=1000,model.type=2,model.object= damagesmodel2,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2), m.percent=1,rep=TRUE)
 write.csv(c2rep$all.boots,'mod2.damage.bvevd.1.rep.boots.csv')
 c2rep<-read.csv(file.path(file.boot,'mod2.damage.bvevd.1.rep.boots.csv'),header=TRUE)[-1]
 c2rep_SE <- dam.boot(file = c2rep,boot.num = 1000, model.type = c(0,1,2)[3], data= tri_var[,c(5,4)], m.percent=1)
 c2rep_SE$boot.se.basins

 
## Model 3
 set.seed(2342339)
 d0<-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=3,model.object= damagesmodel3,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.66,rep=FALSE)
 write.csv(d0$all.boots,'mod3.damage.bvevd.66.boots.csv')
 d0$boot.se.basins
 
 set.seed(2342339)
 d0_rep<-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=3,model.object= damagesmodel3,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.66,rep=TRUE)
 write.csv(d0_rep$all.boots,'mod3.damage.bvevd.66.rep.boots.csv')
 d0_rep$boot.se.basins
 
 set.seed(2342339)
 d1<-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=3,model.object= damagesmodel3,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.75,rep=FALSE)
 write.csv(d1$all.boots,'mod3.damage.bvevd.75.boots.csv')
 d1$boot.se.basins
 
 set.seed(2342339)
 d1_rep <-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=3,model.object= damagesmodel3,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.75,rep=TRUE)
 write.csv(d1_rep$all.boots,'mod3.damage.bvevd.75.rep.boots.csv')
 d1_rep$boot.se.basins
 
 set.seed(2342339)
 d2_rep <-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=3,model.object= damagesmodel3,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=1,rep=TRUE)
 write.csv(d2_rep$all.boots,'mod3.damage.bvevd.1.rep.boots.csv')
 d2_rep$boot.se.basins
 #d0<-read.csv('mod3.damage.bvevd.66.boots.csv',header=TRUE)[-1]
 #d0_SE <- dam.boot(file = d0,boot.num = 1000, model.type = c(0,1,2)[3], data= tri_var[,c(5,4)], m.percent=0.66)
 #c0_SE$boot.se.basins
# Model 4
 ## Model 3
 set.seed(2342339)
 e0<-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=4,model.object= damagesmodel4,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.66,rep=FALSE)
 write.csv(e0$all.boots,'mod3.damage.bvevd.66.boots.csv')
 e0$boot.se.basins
 
 
 set.seed(2342339)
 e0_rep<-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=4,model.object= damagesmodel4,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.66,rep=TRUE)
 write.csv(e0_rep$all.boots,'mod3.damage.bvevd.66.rep.boots.csv')
 e0_rep$boot.se.basins
 
 set.seed(2342339)
 e1<-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=4,model.object= damagesmodel4,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.75,rep=FALSE)
 write.csv(e1$all.boots,'mod3.damage.bvevd.75.boots.csv')
 e1$boot.se.basins
 
 set.seed(2342339)
 e1_rep <-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=4,model.object= damagesmodel4,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=0.75,rep=TRUE)
 write.csv(e1_rep$all.boots,'mod3.damage.bvevd.75.rep.boots.csv')
 e1_rep$boot.se.basins
 
 set.seed(2342339)
 e2_rep <-para.boot.bvevd.damage.mod34(boot.num=1000,model.type=4,model.object= damagesmodel4,data= data.frame(tri_var[,c(5,4)]), cov = data.frame(tri_var[,6],(tri_var[,6])^2, tri_var[,7]), m.percent=1,rep=TRUE)
 write.csv(e2_rep$all.boots,'mod3.damage.bvevd.1.rep.boots.csv')
 e2_rep$boot.se.basins

damagemods<-data.frame(Par=c('beta01','beta11','beta21','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta'),
                    rbind(cbind(Model_Type=0, Estimate =c(dam.ws$estimate[1],NA,NA,NA,dam.ws$estimate[c(4,2,5,3,6,7,8)])), 
                            cbind(Model_Type=1,Estimate =c(damagesmodel1$estimate[1:2],NA,NA, damagesmodel1$estimate[c(5,3,6,4,7,8,9)])), 
                            cbind(Model_Type=2,Estimate =c(damagesmodel2$estimate[c(1:3)],NA,damagesmodel2$estimate[c(6,4,7,5,8,9,10)] )),
                            cbind(Model_Type=3,Estimate =c(damagesmodel3$estimate[c(1:4)],damagesmodel3$estimate[c(7,5,8,6,9,10,11)])),
                            cbind(Model_Type=4,Estimate =c(damagesmodel4$estimate[c(1:2)],NA,damagesmodel4$estimate[c(3,6,4,7,5,8,9,10)]))
                          ))
damagemods

SE_boots<-rbind(
  data.frame(SE= c(a1rep_SE$boot.se.basins[1], NA, NA,a1rep_SE$boot.se.basins[2:8]), Model_Type=0, 
             Method='Par_0.66',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= c(a1_SE$boot.se.basins[1], NA, NA,a1_SE$boot.se.basins[2:8]), Model_Type=0, 
             Method='Par_0.66',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c(a2rep_SE$boot.se.basins[1], NA, NA,a2rep_SE$boot.se.basins[2:8]), Model_Type=0, 
             Method='Par_0.75',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c(a2_SE$boot.se.basins[1], NA, NA,a2_SE$boot.se.basins[2:8]), Model_Type=0, 
             Method='Par_0.75',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c(a3rep_SE$boot.se.basins[1], NA, NA,a3rep_SE$boot.se.basins[2:8]), Model_Type=0, 
             Method='Par_1',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  
  data.frame(SE= c(b0rep_SE$boot.se.basins[1:2], NA,b0rep_SE$boot.se.basins[3:9]), Model_Type=1, 
             Method='Par_0.66',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c(b0_SE$boot.se.basins[1:2], NA,b0_SE$boot.se.basins[3:9]), Model_Type=1, 
             Method='Par_0.66',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c(b1rep_SE$boot.se.basins[1:2], NA,b1rep_SE$boot.se.basins[3:9]), Model_Type=1, 
             Method='Par_0.75',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c(b1_SE$boot.se.basins[1:2], NA,b1_SE$boot.se.basins[3:9]), Model_Type=1, 
             Method='Par_0.75',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c(b2rep_SE$boot.se.basins[1:2], NA,b2rep_SE$boot.se.basins[3:9]), Model_Type=1, 
             Method='Par_1',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  
  data.frame(SE= c0rep_SE$boot.se.basins, Model_Type=2, 
             Method='Par_0.66',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c0_SE$boot.se.basins, Model_Type=2, 
             Method='Par_0.66',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c1rep_SE$boot.se.basins, Model_Type=2, 
             Method='Par_0.75',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c1_SE$boot.se.basins, Model_Type=2, 
             Method='Par_0.75',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= c2rep_SE$boot.se.basins, Model_Type=2, 
             Method='Par_1',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= d0$boot.se.basins, Model_Type=3, 
             Method='Par_0.66',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= d0_rep$boot.se.basins, Model_Type=3, 
             Method='Par_0.66',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= d1$boot.se.basins, Model_Type=3, 
             Method='Par_0.75',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta21','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= d1_rep$boot.se.basins, Model_Type=3, 
             Method='Par_0.75',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= d2_rep$boot.se.basins, Model_Type=3, 
             Method='Par_1',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta21','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  
  data.frame(SE= e0$boot.se.basins, Model_Type=4, 
             Method='Par_0.66',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= e0_rep$boot.se.basins, Model_Type=4, 
             Method='Par_0.66',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= e1$boot.se.basins, Model_Type=4, 
             Method='Par_0.75',
             Replacement = 'No', 
             Par=c('beta01','beta11','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= e1_rep$boot.se.basins, Model_Type=4, 
             Method='Par_0.75',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta')),
  data.frame(SE= e2_rep$boot.se.basins, Model_Type=4, 
             Method='Par_1',
             Replacement = 'Yes', 
             Par=c('beta01','beta11','beta31','beta02','sigma1','sigma2','xi1','xi2','alpha','beta'))
)

#se_graph<-read.csv('damage.bvevd.se.csv')
se_graph <-merge(damagemods, SE_boots, by=c('Model_Type', 'Par'),all=T)
se_graph$lo95<-se_graph$Estimate-2* se_graph$SE
se_graph$hi95<-se_graph$Estimate+2* se_graph$SE

se_graph$Estimate <-as.numeric(se_graph$Estimate)
se_graph$lo95<-se_graph$Estimate-2* se_graph$SE
se_graph$hi95<-se_graph$Estimate+2* se_graph$SE

se_graph$Replacement <- as.factor(se_graph$Replacement)
se_graph$Replacement<-relevel(se_graph$Replacement,'Yes')
se_graph$Method<-factor(se_graph $Method)

se_graph <-se_graph[order(se_graph$Method),]


#####Model 1 est+2SE Plotting#####
se_graph$Method<-relevel(relevel(relevel(se_graph$Method,'Par_1'),'Par_0.75'),'Par_0.66')

se_graph$Par <- as.factor(se_graph$Par)
se_graph$Par<-relevel(relevel(relevel(relevel(relevel(relevel(se_graph$Par,'beta01'),'beta11'),'beta21'),'sigma1'),'xi1'),'alpha')


ggplot(se_graph[se_graph$Model_Type==0 & se_graph$Par %in%c('xi1','sigma1','beta01','alpha'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(alpha)),expression(paste(xi[1])),expression(paste(sigma[1])),expression(paste(gamma['01']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))



ggplot(se_graph[se_graph$Model_Type==1 & se_graph$Par %in%c('xi1','sigma1','beta01','beta11','alpha'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(alpha)),expression(paste(xi[1])),expression(paste(sigma[1])),expression(paste(gamma['11'])),expression(paste(gamma['01']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


ggplot(se_graph[se_graph$Model_Type==2 & se_graph$Par %in%c('xi1','sigma1','beta01','beta11','beta21','alpha'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(alpha)),expression(paste(xi[1])),expression(paste(sigma[1])),expression(paste(gamma['21'])),expression(paste(gamma['11'])),expression(paste(gamma['01']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


ggplot(se_graph[se_graph$Model_Type==3 & se_graph$Par %in%c('xi1','sigma1','beta01','beta11','beta21','beta31','alpha'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(alpha)),expression(paste(xi[1])),expression(paste(sigma[1])),expression(paste(gamma['21'])),expression(paste(gamma['11'])),expression(paste(gamma['01'])),expression(paste(gamma['31']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


ggplot(se_graph[se_graph$Model_Type==4 & se_graph$Par %in%c('xi1','sigma1','beta01','beta11','beta31','alpha'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(alpha)),expression(paste(xi[1])),expression(paste(sigma[1])),expression(paste(gamma['11'])),expression(paste(gamma['01'])),expression(paste(gamma['31']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


se_graph$Par<-relevel(relevel(relevel(relevel(se_graph$Par,'beta02'),'sigma2'),'xi2'),'beta')


ggplot(se_graph[se_graph$Model_Type==0 & se_graph$Par %in%c('xi2','sigma2','beta02','beta'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(beta)),expression(paste(xi[2])),expression(paste(sigma[2])),expression(paste(gamma['02']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


ggplot(se_graph[se_graph$Model_Type==1 & se_graph$Par %in%c('xi2','sigma2','beta02','beta'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(beta)),expression(paste(xi[2])),expression(paste(sigma[2])),expression(paste(gamma['02']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


ggplot(se_graph[se_graph$Model_Type==2 & se_graph$Par %in%c('xi2','sigma2','beta02','beta'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(beta)),expression(paste(xi[2])),expression(paste(sigma[2])),expression(paste(gamma['02']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


ggplot(se_graph[se_graph$Model_Type==3 & se_graph$Par %in%c('xi2','sigma2','beta02','beta'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(beta)),expression(paste(xi[2])),expression(paste(sigma[2])),expression(paste(gamma['02']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


ggplot(se_graph[se_graph$Model_Type==4 & se_graph$Par %in%c('xi2','sigma2','beta02','beta'),],aes(x=Par, ymin= lo95,ymax= hi95,pch= Method,color= Method,group= Method:Replacement,lty= Replacement))+ geom_linerange(lwd=1,position= position_dodge(.5))+xlab('Parameter')+ylab('Estimate +/- 2 SE')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +coord_flip()+
  scale_x_discrete(label=c(expression(paste(beta)),expression(paste(xi[2])),expression(paste(sigma[2])),expression(paste(gamma['02']))))+geom_hline(yintercept =0,lty=5, alpha=0.5) +theme(legend.position = 'none', legend.background = element_rect(colour = 'black', fill = "white"),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))


