#################################################################
#################################################################
#Parametric Bootstrap wLM Function

para.boot.lm<-function(boot.num=100,model.object=wlm3,model.type='3',data= modeling_data, m.percent=1,rep= FALSE){
	
	coeffs<-matrix(nrow= boot.num,ncol=length(coef(model.object))) 
	Intercepts<-matrix(nrow= boot.num,ncol=7)
	slopes<-matrix(nrow= boot.num,ncol=7)
	slopesquared <-matrix(nrow= boot.num,ncol=7)
	lat <-matrix(nrow= boot.num,ncol=7)
	ts <-matrix(nrow= boot.num,ncol=7)
	sigma<-vector(length= boot.num)
	
	data$mean.par<-predict(model.object)
	cov.par<-Diagonal(length(model.object$model$'(weights)'),1/sqrt(model.object$model$'(weights)')* summary(model.object)$sigma)
	data$covs<-diag(cov.par)
	
	data$boot.y<-NA
	
	NA_set<-data[which(data$Basin=='NoAt'),]; NA_len<-dim(NA_set)[1]
	WP_set<-data[which(data$Basin=='WePa'),]; WP_len<-dim(WP_set)[1]
	EP_set<-data[which(data$Basin=='EaPa'),]; EP_len<-dim(EP_set)[1]
	SP_set<-data[which(data$Basin=='SoPa'),];	 SP_len<-dim(SP_set)[1]
	NI_set<-data[which(data$Basin=='NoIn'),];	 NI_len<-dim(NI_set)[1]
	SEI_set<-data[which(data$Basin=='SoEIn'),]; SEI_len<-dim(SEI_set)[1]
	SWI_set<-data[which(data$Basin=='SoWIn'),]; SWI_len<-dim(SWI_set)[1]
	
	sample.lengths<- floor(m.percent*
		c(NA_len, WP_len, EP_len, SP_len, NI_len, SEI_len, SWI_len))
	
	sample.size.m<-sum(sample.lengths)

	
	for(i in 1:boot.num){
		
		NA_samp<-sample(NA_len,sample.lengths[1],replace= rep)
		WP_samp<-sample(WP_len,sample.lengths[2],replace= rep)
		EP_samp<-sample(EP_len,sample.lengths[3],replace= rep)
		SP_samp<-sample(SP_len,sample.lengths[4],replace= rep)
		NI_samp<-sample(NI_len,sample.lengths[5],replace= rep)
		SEI_samp<-sample(SEI_len,sample.lengths[6],replace= rep)
		SWI_samp<-sample(SWI_len,sample.lengths[7],replace= rep)
		
		condition<-c(NA_samp, NA_len+WP_samp, WP_len+ NA_len +EP_samp, 
		EP_len +WP_len+ NA_len +SP_samp, SP_len +EP_len +WP_len+ NA_len+ NI_samp, 
		NI_len + SP_len +EP_len +WP_len+ NA_len+ SEI_samp,
		SEI_len + NI_len + SP_len +EP_len +WP_len+ NA_len+ SWI_samp)
		
		data[eval(condition),]$boot.y<-c(rmvn(1, data[eval(condition),]$mean.par, 
				Diagonal(length(data[eval(condition),]$covs),data[eval(condition),]$covs)))
		
		data.boot<-na.omit(data)

		if(model.type=='1'){
			temp<-lm(boot.y ~log(1013-minCP)*Basin,weights=log(1013-minCP), data.boot)
			coeffs[i,]<-coef(temp)
		
		Intercepts[i,]<- c(coeffs[i,1],
			coeffs[i,1]+coeffs[i,1+2],
			coeffs[i,1]+coeffs[i,1+3],
			coeffs[i,1]+coeffs[i,1+4],
			coeffs[i,1]+coeffs[i,1+5],
			coeffs[i,1]+coeffs[i,1+6],
			coeffs[i,1]+coeffs[i,1+7])
			
		slopes[i,]<- c(coeffs[i,2],
			coeffs[i,2]+coeffs[i,7+2],
			coeffs[i,2]+coeffs[i,7+3],
			coeffs[i,2]+coeffs[i,7+4],
			coeffs[i,2]+coeffs[i,7+5],
			coeffs[i,2]+coeffs[i,7+6],
			coeffs[i,2]+coeffs[i,7+7])
					}
		
		if(model.type=='2'){
			temp<-lm(boot.y ~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin,weights=log(1013-minCP), data.boot)
			coeffs[i,]<-coef(temp)
		
		Intercepts[i,]<- c(coeffs[i,1],
			coeffs[i,1]+coeffs[i,1+2],
			coeffs[i,1]+coeffs[i,1+3],
			coeffs[i,1]+coeffs[i,1+4],
			coeffs[i,1]+coeffs[i,1+5],
			coeffs[i,1]+coeffs[i,1+6],
			coeffs[i,1]+coeffs[i,1+7])
			
		slopes[i,]<- c(coeffs[i,2],
			coeffs[i,2]+coeffs[i,1+7+2],
			coeffs[i,2]+coeffs[i,1+7+3],
			coeffs[i,2]+coeffs[i,1+7+4],
			coeffs[i,2]+coeffs[i,1+7+5],
			coeffs[i,2]+coeffs[i,1+7+6],
			coeffs[i,2]+coeffs[i,1+7+7])
				
		slopesquared[i,]<- c(coeffs[i,9],
			coeffs[i,9]+ coeffs[i,9+5+2],
			coeffs[i,9]+ coeffs[i,9+5+3],
			coeffs[i,9]+ coeffs[i,9+5+4],
			coeffs[i,9]+ coeffs[i,9+5+5],
			coeffs[i,9]+ coeffs[i,9+5+6],
			coeffs[i,9]+ coeffs[i,9+5+7])
		}
			
		if(model.type=='3'){
			temp<-lm(boot.y ~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin+avgLat*Basin+ Avg_Trans_Speed*Basin,weights=log(1013-minCP), data.boot)
			coeffs[i,]<-coef(temp)
		
		Intercepts[i,]<- c(coeffs[i,1],
			coeffs[i,1]+coeffs[i,1+2],
			coeffs[i,1]+coeffs[i,1+3],
			coeffs[i,1]+coeffs[i,1+4],
			coeffs[i,1]+coeffs[i,1+5],
			coeffs[i,1]+coeffs[i,1+6],
			coeffs[i,1]+coeffs[i,1+7])
			
		slopes[i,]<- c(coeffs[i,2],
			coeffs[i,2]+coeffs[i,1+9+2],
			coeffs[i,2]+coeffs[i,1+9+3],
			coeffs[i,2]+coeffs[i,1+9+4],
			coeffs[i,2]+coeffs[i,1+9+5],
			coeffs[i,2]+coeffs[i,1+9+6],
			coeffs[i,2]+coeffs[i,1+9+7])
			
		slopesquared[i,]<- c(coeffs[i,9],
			coeffs[i,9]+ coeffs[i,9+7+2],
			coeffs[i,9]+ coeffs[i,9+7+3],
			coeffs[i,9]+ coeffs[i,9+7+4],
			coeffs[i,9]+ coeffs[i,9+7+5],
			coeffs[i,9]+ coeffs[i,9+7+6],
			coeffs[i,9]+ coeffs[i,9+7+7])
			
		lat[i,]<- c(coeffs[i,10],
			coeffs[i,10]+coeffs[i,10+12+2],
			coeffs[i,10]+coeffs[i,10+12+3],
			coeffs[i,10]+coeffs[i,10+12+4],
			coeffs[i,10]+coeffs[i,10+12+5],
			coeffs[i,10]+coeffs[i,10+12+6],
			coeffs[i,10]+coeffs[i,10+12+7])
			
			
		ts[i,]<- c(coeffs[i,11],
			coeffs[i,11]+coeffs[i,11+17+2],
			coeffs[i,11]+coeffs[i,11+17+3],
			coeffs[i,11]+coeffs[i,11+17+4],
			coeffs[i,11]+coeffs[i,11+17+5],
			coeffs[i,11]+coeffs[i,11+17+6],
			coeffs[i,11]+coeffs[i,11+17+7])
		}
		
		sigma[i]<-summary(temp)$sigma
	
	}
	
	sample.size.n<-dim(data)[1]
		
	boot.est.sigma<-mean(sigma)
	boot.est.sigma.se<-sd(sigma)

	boot.est.basin.Int<-colMeans(Intercepts)
	boot.se.basin.Int<-apply(Intercepts, 2, sd)
	names(boot.est.basin.Int)<-c('NoAt Int.','WePa Int.' ,'EaPa Int.' ,'SoPa Int.' ,'NoIn Int.' ,'SoEIn Int.','SoWIn Int.')
	names(boot.se.basin.Int)<-c('NoAt Int.','WePa Int.' ,'EaPa Int.' ,'SoPa Int.' ,'NoIn Int.' ,'SoEIn Int.','SoWIn Int.')

	boot.est.basin.slope<-colMeans(slopes)
	boot.se.basin.slope<-apply(slopes, 2, sd)
	names(boot.est.basin.slope)<-c('NoAt log(1013 - minCP)','WePa log(1013 - minCP)' ,'EaPa log(1013 - minCP)' ,'SoPa log(1013 - minCP)' ,'NoIn log(1013 - minCP)' ,'SoEIn log(1013 - minCP)','SoWIn log(1013 - minCP)')
	names(boot.se.basin.slope)<-c('NoAt log(1013 - minCP)','WePa log(1013 - minCP)' ,'EaPa log(1013 - minCP)' ,'SoPa log(1013 - minCP)' ,'NoIn log(1013 - minCP)' ,'SoEIn log(1013 - minCP)','SoWIn log(1013 - minCP)')
	
	if(model.type%in%c('2','3')){
		boot.est.basin.slopesquared<-colMeans(slopesquared)
		boot.se.basin.slopesquared<-apply(slopesquared, 2, sd)
		names(boot.est.basin.slopesquared)<-c('NoAt (log(1013 - minCP))^2','WePa log(1013 - minCP))^2' ,'EaPa (log(1013 - minCP))^2' ,'SoPa (log(1013 - minCP))^2' ,'NoIn (log(1013 - minCP))^2' ,'SoEIn (log(1013 - minCP))^2','SoWIn (log(1013 - minCP))^2')
		names(boot.se.basin.slopesquared)<-c('NoAt (log(1013 - minCP))^2','WePa log(1013 - minCP))^2' ,'EaPa (log(1013 - minCP))^2' ,'SoPa (log(1013 - minCP))^2' ,'NoIn (log(1013 - minCP))^2' ,'SoEIn (log(1013 - minCP))^2','SoWIn (log(1013 - minCP))^2')

		if(model.type=='3'){
			boot.est.basin.avgLat<-colMeans(lat)
			boot.se.basin.avgLat<-apply(lat, 2, sd)
			names(boot.est.basin.avgLat)<-c('NoAt avgLat','WePa avgLat' ,'EaPa avgLat' ,'SoPa avgLat' ,'NoIn avgLat' ,'SoEIn avgLat','SoWIn avgLat')
			names(boot.se.basin.avgLat)<-c('NoAt avgLat','WePa avgLat' ,'EaPa avgLat' ,'SoPa avgLat' ,'NoIn avgLat' ,'SoEIn avgLat','SoWIn avgLat')
			
			boot.est.basin.avgTS<-colMeans(ts)
			boot.se.basin.avgTS <-apply(ts, 2, sd)
			names(boot.est.basin.avgTS)<-c('NoAt avgTS','WePa avgTS' ,'EaPa avgTS' ,'SoPa avgTS' ,'NoIn avgTS' ,'SoEIn avgTS','SoWIn avgTS')
			names(boot.se.basin.avgTS)<-c('NoAt avgTS','WePa avgTS' ,'EaPa avgTS' ,'SoPa avgTS' ,'NoIn avgTS' ,'SoEIn avgTS','SoWIn avgTS')

			}
			
		}	


	if(model.type=='1'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, sigma=boot.est.sigma)
		boot.se.basins<-c(boot.se.basin.Int, boot.se.basin.slope, sigma=boot.est.sigma.se)*sqrt(sample.size.n/sample.size.m)
	}
	if(model.type=='2'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, boot.est.basin.slopesquared,sigma=boot.est.sigma)
		boot.se.basins<-c(boot.se.basin.Int, boot.se.basin.slope, boot.se.basin.slopesquared, sigma=boot.est.sigma.se)*sqrt(sample.size.n/sample.size.m)
	}
	if(model.type=='3'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, boot.est.basin.slopesquared, boot.est.basin.avgLat, boot.est.basin.avgTS, sigma=boot.est.sigma)
		boot.se.basins<-c(boot.se.basin.Int, boot.se.basin.slope, boot.se.basin.slopesquared, boot.se.basin.avgLat, boot.se.basin.avgTS, sigma=boot.est.sigma.se)*sqrt(sample.size.n/sample.size.m)
	}
	
	boot.est.orig<-c(colMeans(coeffs), boot.est.sigma)
	names(boot.est.orig)<-c(names(coef(temp)),'sigma')

	boot.se.orig<-c(apply(coeffs, 2, sd), boot.est.sigma.se)*sqrt(sample.size.n/sample.size.m)
	names(boot.se.orig)<-c(names(coef(temp)),'sigma')
	
	return(list(
		boot.est.basins = boot.est.basins, boot.se.basins = boot.se.basins,
		boot.est.orig = boot.est.orig , boot.se.orig = boot.se.orig ,all.coeff.est= coeffs))
}

