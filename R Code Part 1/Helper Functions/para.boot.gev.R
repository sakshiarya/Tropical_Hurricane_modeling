#################################################################
#################################################################
#Parametric Bootstrap Function

para.boot.gev.Mar11<-function(boot.num=100,model.object=gev3,model.type='3',data= modeling_data, m.percent=1, ordered_basins=c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn'),rep=FALSE){
	
	coeffs<-matrix(nrow= boot.num,ncol=length(model.object$results$par)) 
	Intercepts<-matrix(nrow= boot.num,ncol=7)
	slopes<-matrix(nrow= boot.num,ncol=7)
	slopesquared <-matrix(nrow= boot.num,ncol=7)
	lat <-matrix(nrow= boot.num,ncol=7)
	ts <-matrix(nrow= boot.num,ncol=7)
	scale.vec<-vector(length= boot.num)
	shape.vec<-vector(length= boot.num)
	conv.vec <- vector(length = boot.num)
	
	total.parameters<- c(16,23,37)[as.numeric(model.type)]
	
	colnames_fit1<- c(colnames(model.matrix(model.object$par.models$location, 
		data)),'shape','scale')
	coef1<-model.object$results$par
	names(coef1)<-colnames_fit1
	coef.cov<-summary(model.object)$cov.theta
	colnames(coef.cov)<-colnames_fit1
	rownames(coef.cov)<-colnames_fit1

	WePa_coefs <- coef1[grepl(ordered_basins[2], colnames_fit1)]
	WePa_coefs_loc <- grep(ordered_basins[2], colnames_fit1)
	
	EaPa_coefs <- coef1[grepl(ordered_basins[3], colnames_fit1)]
	EaPa_coefs_loc <- grep(ordered_basins[3], colnames_fit1)
	
	SoPa_coefs <- coef1[grepl(ordered_basins[4], colnames_fit1)]
	SoPa_coefs_loc <- grep(ordered_basins[4], colnames_fit1)
	
	NoIn_coefs <- coef1[grepl(ordered_basins[5], colnames_fit1)]
	NoIn_coefs_loc <- grep(ordered_basins[5], colnames_fit1)
	
	SoEIn_coefs <- coef1[grepl(ordered_basins[6], colnames_fit1)]
	SoEIn_coefs_loc <- grep(ordered_basins[6], colnames_fit1)
	
	SoWIn_coefs <- coef1[grepl(ordered_basins[7], colnames_fit1)]
	SoWIn_coefs_loc <- grep(ordered_basins[7], colnames_fit1)
	
	
	if(model.type=='1'){			
		loc.par<-model.object$results$par[1:14]
		covariates<-model.matrix(~log(1013-minCP)*Basin,modeling_data)	
	}else if(model.type=='2'){
		loc.par<-model.object$results$par[1:21]
		covariates<-model.matrix(~log(1013-minCP)*Basin+
			I((log(1013-minCP))^2)*Basin,modeling_data)	
	}else {
		loc.par<-model.object$results$par[1:35]
		covariates<-model.matrix(~log(1013-minCP)*Basin+
			I((log(1013-minCP))^2)*Basin+ 
			avgLat*Basin+ Avg_Trans_Speed*Basin,modeling_data)	
		}
		 
	data$location<-covariates%*%loc.par
	scale.par<-model.object$results$par['scale']
	shape.par<-model.object$results$par['shape']
	
	NA_set<-data[which(data$Basin=='NoAt'),'location']; NA_len<-length(NA_set)
	WP_set<-data[which(data$Basin=='WePa'),'location']; WP_len<-length(WP_set)
	EP_set<-data[which(data$Basin=='EaPa'),'location']; EP_len<-length(EP_set)
	SP_set<-data[which(data$Basin=='SoPa'),'location'];	 SP_len<-length(SP_set)
	NI_set<-data[which(data$Basin=='NoIn'),'location'];	 NI_len<-length(NI_set)
	SEI_set<-data[which(data$Basin=='SoEIn'),'location']; SEI_len<-length(SEI_set)
	SWI_set<-data[which(data$Basin=='SoWIn'),'location']; SWI_len<-length(SWI_set)
	
	sample.lengths<- floor(m.percent*
		c(NA_len, WP_len, EP_len, SP_len, NI_len, SEI_len, SWI_len))
	
	sample.size.m<-sum(sample.lengths)
	
	data$boot.y<-NA
	
	for(i in 1:boot.num){
		NA_samp<-sample(NA_len,sample.lengths[1],replace=rep)
		WP_samp<-sample(WP_len,sample.lengths[2],replace=rep)
		EP_samp<-sample(EP_len,sample.lengths[3],replace=rep)
		SP_samp<-sample(SP_len,sample.lengths[4],replace=rep)
		NI_samp<-sample(NI_len,sample.lengths[5],replace=rep)
		SEI_samp<-sample(SEI_len,sample.lengths[6],replace=rep)
		SWI_samp<-sample(SWI_len,sample.lengths[7],replace=rep)
		
			boot.locs<-c(NA_set[NA_samp],
			WP_set[WP_samp],
			EP_set[EP_samp],
			SP_set[SP_samp],
			NI_set[NI_samp],
			SEI_set[SEI_samp],
			SWI_set[SWI_samp])
		
		data[c(NA_samp, NA_len+WP_samp, WP_len+ NA_len +EP_samp, 
		EP_len +WP_len+ NA_len +SP_samp, SP_len +EP_len +WP_len+ NA_len+ NI_samp, 
		NI_len + SP_len +EP_len +WP_len+ NA_len+ SEI_samp,
		SEI_len + NI_len + SP_len +EP_len +WP_len+ NA_len+ SWI_samp),]$boot.y<-revd(sample.size.m,loc= boot.locs,scale= scale.par, shape= shape.par,type='GEV')
		
		data.boot<-na.omit(data)

		if(model.type=='1'){
		
		coeffs<-matrix(nrow= boot.num, ncol=length(colnames_fit1))
		colnames(coeffs)<-colnames_fit1
			
		temp <-with(data.boot,fevd(boot.y, data.boot,
			location.fun=~log(1013-minCP)*Basin,
			scale.fun=~1,
			shape.fun=~1,type='GEV',method = "MLE", period.basis='Tropical Cyclone'))
		
		coeffs[i,]<-temp$results$par
				
			WePa_coefs <- coeffs[i,grepl(ordered_basins[2], colnames_fit1)]
			WePa_coefs_loc <- grep(ordered_basins[2], colnames_fit1)
	
			EaPa_coefs <- coeffs[i,grepl(ordered_basins[3], colnames_fit1)]
			EaPa_coefs_loc <- grep(ordered_basins[3], colnames_fit1)
	
			SoPa_coefs <- coeffs[i,grepl(ordered_basins[4], colnames_fit1)]
			SoPa_coefs_loc <- grep(ordered_basins[4], colnames_fit1)
	
			NoIn_coefs <- coeffs[i,grepl(ordered_basins[5], colnames_fit1)]
			NoIn_coefs_loc <- grep(ordered_basins[5], colnames_fit1)
	
			SoEIn_coefs <- coeffs[i,grepl(ordered_basins[6], colnames_fit1)]
			SoEIn_coefs_loc <- grep(ordered_basins[6], colnames_fit1)
	
			SoWIn_coefs <- coeffs[i,grepl(ordered_basins[7], colnames_fit1)]
			SoWIn_coefs_loc <- grep(ordered_basins[7], colnames_fit1)
		
		Intercepts[i,]<- c(coeffs[i,1],
			 coeffs[i,1]+ WePa_coefs[1],
			 coeffs[i,1]+ EaPa_coefs[1],
			 coeffs[i,1]+ SoPa_coefs[1],
			 coeffs[i,1]+ NoIn_coefs[1],
			 coeffs[i,1]+ SoEIn_coefs[1],
			 coeffs[i,1]+ SoWIn_coefs[1])
			 
		#	 c(coeffs[i,1],
		#	coeffs[i,1]+coeffs[i,1+2],
		#	coeffs[i,1]+coeffs[i,1+3],
		#	coeffs[i,1]+coeffs[i,1+4],
		#	coeffs[i,1]+coeffs[i,1+5],
		#	coeffs[i,1]+coeffs[i,1+6],
		#	coeffs[i,1]+coeffs[i,1+7])
			
		slopes[i,]<- c(coeffs[i,2],
		 	 coeffs[i,2]+ WePa_coefs[2],
	 		 coeffs[i,2]+ EaPa_coefs[2],
			 coeffs[i,2]+ SoPa_coefs[2],
			 coeffs[i,2]+ NoIn_coefs[2],
			 coeffs[i,2]+ SoEIn_coefs[2],
			 coeffs[i,2]+ SoWIn_coefs[2])
		#c(coeffs[i,2],
		#	coeffs[i,2]+coeffs[i,7+2],
		#	coeffs[i,2]+coeffs[i,7+3],
		#	coeffs[i,2]+coeffs[i,7+4],
		#	coeffs[i,2]+coeffs[i,7+5],
		#	coeffs[i,2]+coeffs[i,7+6],
		#	coeffs[i,2]+coeffs[i,7+7])
			
		scale.vec[i]<-temp$results$par['scale']	
		shape.vec[i]<-temp$results$par['shape']

		}
		
		if(model.type=='2'){
			
			coeffs<-matrix(nrow= boot.num, ncol=length(colnames_fit1))
			colnames(coeffs)<-colnames_fit1

			temp <-with(data.boot,fevd(boot.y, data.boot,
			location.fun=~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin,
			scale.fun=~1,
			shape.fun=~1,type='GEV',method = "MLE", period.basis='Tropical Cyclone'))
			
			coeffs[i,]<-temp$results$par
			WePa_coefs <- coeffs[i,grepl(ordered_basins[2], colnames_fit1)]
			WePa_coefs_loc <- grep(ordered_basins[2], colnames_fit1)
	
			EaPa_coefs <- coeffs[i,grepl(ordered_basins[3], colnames_fit1)]
			EaPa_coefs_loc <- grep(ordered_basins[3], colnames_fit1)
	
			SoPa_coefs <- coeffs[i,grepl(ordered_basins[4], colnames_fit1)]
			SoPa_coefs_loc <- grep(ordered_basins[4], colnames_fit1)
	
			NoIn_coefs <- coeffs[i,grepl(ordered_basins[5], colnames_fit1)]
			NoIn_coefs_loc <- grep(ordered_basins[5], colnames_fit1)
	
			SoEIn_coefs <- coeffs[i,grepl(ordered_basins[6], colnames_fit1)]
			SoEIn_coefs_loc <- grep(ordered_basins[6], colnames_fit1)
	
			SoWIn_coefs <- coeffs[i,grepl(ordered_basins[7], colnames_fit1)]
			SoWIn_coefs_loc <- grep(ordered_basins[7], colnames_fit1)
		
		Intercepts[i,]<- c(coeffs[i,1],
			 coeffs[i,1]+ WePa_coefs[1],
			 coeffs[i,1]+ EaPa_coefs[1],
			 coeffs[i,1]+ SoPa_coefs[1],
			 coeffs[i,1]+ NoIn_coefs[1],
			 coeffs[i,1]+ SoEIn_coefs[1],
			 coeffs[i,1]+ SoWIn_coefs[1])

		# Intercepts[i,]<- c(coeffs[i,1],
			# coeffs[i,1]+coeffs[i,1+2],
			# coeffs[i,1]+coeffs[i,1+3],
			# coeffs[i,1]+coeffs[i,1+4],
			# coeffs[i,1]+coeffs[i,1+5],
			# coeffs[i,1]+coeffs[i,1+6],
			# coeffs[i,1]+coeffs[i,1+7])
			
		slopes[i,]<- c(coeffs[i,2],
		 	 coeffs[i,2]+ WePa_coefs[2],
	 		 coeffs[i,2]+ EaPa_coefs[2],
			 coeffs[i,2]+ SoPa_coefs[2],
			 coeffs[i,2]+ NoIn_coefs[2],
			 coeffs[i,2]+ SoEIn_coefs[2],
			 coeffs[i,2]+ SoWIn_coefs[2])

		# slopes[i,]<- c(coeffs[i,2],
			# coeffs[i,2]+coeffs[i,1+7+2],
			# coeffs[i,2]+coeffs[i,1+7+3],
			# coeffs[i,2]+coeffs[i,1+7+4],
			# coeffs[i,2]+coeffs[i,1+7+5],
			# coeffs[i,2]+coeffs[i,1+7+6],
			# coeffs[i,2]+coeffs[i,1+7+7])
				
		slopesquared[i,]<- c(coeffs[i, 9],
		 	 coeffs[i, 9]+ WePa_coefs[3],
	 		 coeffs[i, 9]+ EaPa_coefs[3],
			 coeffs[i, 9]+ SoPa_coefs[3],
			 coeffs[i, 9]+ NoIn_coefs[3],
			 coeffs[i, 9]+ SoEIn_coefs[3],
			 coeffs[i, 9]+ SoWIn_coefs[3])
		
		# c(coeffs[i,9],
			# coeffs[i,9]+ coeffs[i,9+5+2],
			# coeffs[i,9]+ coeffs[i,9+5+3],
			# coeffs[i,9]+ coeffs[i,9+5+4],
			# coeffs[i,9]+ coeffs[i,9+5+5],
			# coeffs[i,9]+ coeffs[i,9+5+6],
			# coeffs[i,9]+ coeffs[i,9+5+7])
			
		scale.vec[i]<-temp$results$par['scale']	
		shape.vec[i]<-temp$results$par['shape']

		}
			
		if(model.type=='3'){

		coeffs<-matrix(nrow= boot.num, ncol=length(colnames_fit1))
		colnames(coeffs)<-colnames_fit1

		temp <-with(data.boot,fevd(boot.y, data.boot,
			location.fun=~log(1013-minCP)*Basin+I((log(1013-minCP))^2)*Basin+ 
				avgLat*Basin+ Avg_Trans_Speed*Basin,
			scale.fun=~1,
			shape.fun=~1,type='GEV',method = "MLE", period.basis='Tropical Cyclone'))

			
			coeffs[i,]<-temp$results$par
		WePa_coefs <- coeffs[i,grepl(ordered_basins[2], colnames_fit1)]
			WePa_coefs_loc <- grep(ordered_basins[2], colnames_fit1)
	
			EaPa_coefs <- coeffs[i,grepl(ordered_basins[3], colnames_fit1)]
			EaPa_coefs_loc <- grep(ordered_basins[3], colnames_fit1)
	
			SoPa_coefs <- coeffs[i,grepl(ordered_basins[4], colnames_fit1)]
			SoPa_coefs_loc <- grep(ordered_basins[4], colnames_fit1)
	
			NoIn_coefs <- coeffs[i,grepl(ordered_basins[5], colnames_fit1)]
			NoIn_coefs_loc <- grep(ordered_basins[5], colnames_fit1)
	
			SoEIn_coefs <- coeffs[i,grepl(ordered_basins[6], colnames_fit1)]
			SoEIn_coefs_loc <- grep(ordered_basins[6], colnames_fit1)
	
			SoWIn_coefs <- coeffs[i,grepl(ordered_basins[7], colnames_fit1)]
			SoWIn_coefs_loc <- grep(ordered_basins[7], colnames_fit1)
		
		Intercepts[i,]<- c(coeffs[i,1],
			 coeffs[i,1]+ WePa_coefs[1],
			 coeffs[i,1]+ EaPa_coefs[1],
			 coeffs[i,1]+ SoPa_coefs[1],
			 coeffs[i,1]+ NoIn_coefs[1],
			 coeffs[i,1]+ SoEIn_coefs[1],
			 coeffs[i,1]+ SoWIn_coefs[1])

		# Intercepts[i,]<- c(coeffs[i,1],
			# coeffs[i,1]+coeffs[i,1+2],
			# coeffs[i,1]+coeffs[i,1+3],
			# coeffs[i,1]+coeffs[i,1+4],
			# coeffs[i,1]+coeffs[i,1+5],
			# coeffs[i,1]+coeffs[i,1+6],
			# coeffs[i,1]+coeffs[i,1+7])
			
		slopes[i,]<- c(coeffs[i,2],
		 	 coeffs[i,2]+ WePa_coefs[2],
	 		 coeffs[i,2]+ EaPa_coefs[2],
			 coeffs[i,2]+ SoPa_coefs[2],
			 coeffs[i,2]+ NoIn_coefs[2],
			 coeffs[i,2]+ SoEIn_coefs[2],
			 coeffs[i,2]+ SoWIn_coefs[2])
			 
#  	   slopes[i,]<- c(coeffs[i,2],
			# coeffs[i,2]+coeffs[i,1+9+2],
			# coeffs[i,2]+coeffs[i,1+9+3],
			# coeffs[i,2]+coeffs[i,1+9+4],
			# coeffs[i,2]+coeffs[i,1+9+5],
			# coeffs[i,2]+coeffs[i,1+9+6],
			# coeffs[i,2]+coeffs[i,1+9+7])
			
		slopesquared[i,]<- c(coeffs[i, 9],
		 	 coeffs[i, 9]+ WePa_coefs[3],
	 		 coeffs[i, 9]+ EaPa_coefs[3],
			 coeffs[i, 9]+ SoPa_coefs[3],
			 coeffs[i, 9]+ NoIn_coefs[3],
			 coeffs[i, 9]+ SoEIn_coefs[3],
			 coeffs[i, 9]+ SoWIn_coefs[3])

		# slopesquared[i,]<- c(coeffs[i,9],
			# coeffs[i,9]+ coeffs[i,9+7+2],
			# coeffs[i,9]+ coeffs[i,9+7+3],
			# coeffs[i,9]+ coeffs[i,9+7+4],
			# coeffs[i,9]+ coeffs[i,9+7+5],
			# coeffs[i,9]+ coeffs[i,9+7+6],
			# coeffs[i,9]+ coeffs[i,9+7+7])
			
		lat[i,]<- c(coeffs[i, 10],
		 	 coeffs[i, 10]+ WePa_coefs[4],
	 		 coeffs[i, 10]+ EaPa_coefs[4],
			 coeffs[i, 10]+ SoPa_coefs[4],
			 coeffs[i, 10]+ NoIn_coefs[4],
			 coeffs[i, 10]+ SoEIn_coefs[4],
			 coeffs[i, 10]+ SoWIn_coefs[4])
			 
		# c(coeffs[i,10],
			# coeffs[i,10]+coeffs[i,10+12+2],
			# coeffs[i,10]+coeffs[i,10+12+3],
			# coeffs[i,10]+coeffs[i,10+12+4],
			# coeffs[i,10]+coeffs[i,10+12+5],
			# coeffs[i,10]+coeffs[i,10+12+6],
			# coeffs[i,10]+coeffs[i,10+12+7])
			
			
		ts[i,]<- c(coeffs[i, 11],
		 	 coeffs[i, 11]+ WePa_coefs[5],
	 		 coeffs[i, 11]+ EaPa_coefs[5],
			 coeffs[i, 11]+ SoPa_coefs[5],
			 coeffs[i, 11]+ NoIn_coefs[5],
			 coeffs[i, 11]+ SoEIn_coefs[5],
			 coeffs[i, 11]+ SoWIn_coefs[5])
			 
			# c(coeffs[i,11],
			# coeffs[i,11]+coeffs[i,11+17+2],
			# coeffs[i,11]+coeffs[i,11+17+3],
			# coeffs[i,11]+coeffs[i,11+17+4],
			# coeffs[i,11]+coeffs[i,11+17+5],
			# coeffs[i,11]+coeffs[i,11+17+6],
			# coeffs[i,11]+coeffs[i,11+17+7])
			
		scale.vec[i]<-temp$results$par['scale']	
		shape.vec[i]<-temp$results$par['shape']
		}
	conv.vec[i] <- temp$results$convergence	
	}
	
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

	sample.size.n<-dim(data)[1]
	
	boot.est.scale<-mean(scale.vec)
	boot.se.scale<-sd(scale.vec)* sqrt(sample.size.n/sample.size.m)

	boot.est.shape<-mean(shape.vec)
	boot.se.shape<-sd(shape.vec)* sqrt(sample.size.n/sample.size.m)

	if(model.type=='1'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, scale=boot.est.scale, shape=boot.est.shape)
		boot.se.basins<-c(boot.se.basin.Int, boot.se.basin.slope, scale=boot.se.scale, shape=boot.se.shape)* sqrt(sample.size.n/sample.size.m)
		all.est<-data.frame(Intercepts, slopes, scale.vec, shape.vec)	
		names(all.est)<-c(names(boot.est.basin.Int), names(boot.est.basin.slope),
			'scale','shape')

	}
	if(model.type=='2'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, boot.est.basin.slopesquared, scale=boot.est.scale, shape=boot.est.shape)
		boot.se.basins<-c(boot.se.basin.Int, boot.se.basin.slope, boot.se.basin.slopesquared, scale=boot.se.scale, shape=boot.se.shape)* sqrt(sample.size.n/sample.size.m)
		
		all.est<-data.frame(Intercepts, slopes, slopesquared, scale.vec, shape.vec)	
		names(all.est)<-c(names(boot.est.basin.Int), names(boot.est.basin.slope),
		names(boot.est.basin.slopesquared),'scale','shape')

	}
	if(model.type=='3'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, boot.est.basin.slopesquared, boot.est.basin.avgLat, boot.est.basin.avgTS, scale=boot.est.scale, shape=boot.est.shape)
		boot.se.basins<-c(boot.se.basin.Int, boot.se.basin.slope, boot.se.basin.slopesquared, boot.se.basin.avgLat, boot.se.basin.avgTS, scale=boot.se.scale, shape=boot.se.shape)* sqrt(sample.size.n/sample.size.m)
		
	all.est<-data.frame(Intercepts, slopes, slopesquared, lat, ts, scale.vec, shape.vec)	
	names(all.est)<-c(names(boot.est.basin.Int), names(boot.est.basin.slope),
		names(boot.est.basin.slopesquared), names(boot.est.basin.avgLat),
		names(boot.est.basin.avgTS),'scale','shape')
	}
	
	#boot.est.orig<-colMeans(coeffs)
	#names(boot.est.orig)<-c(colnames(covariates),'scale','shape')
	#boot.se.orig<-apply(coeffs, 2, sd)* sqrt(sample.size.n/sample.size.m)
	#names(boot.se.orig)<-c(colnames(covariates),'scale','shape')
	
	#colnames(coeffs)<-c(colnames(covariates),'scale','shape')
	
	return(list(
		boot.est.basins = boot.est.basins, boot.se.basins = boot.se.basins,
		all.coeff.est= all.est, convergence = conv.vec))
}

