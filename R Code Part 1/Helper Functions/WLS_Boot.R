WLS_Boot<-function(data= modeling_data, data.file=read.csv('repwlm3.75.boot.all.csv',header=T, check.names = FALSE)[,-1], m.percent =0.75,model.type='3', ordered_basins= ordered_basins){

	boot.num =dim(data.file)[1]
	
	Intercepts<-matrix(nrow= boot.num,ncol=7)
	slopes<-matrix(nrow= boot.num,ncol=7)
	slopesquared <-matrix(nrow= boot.num,ncol=7)
	lat <-matrix(nrow= boot.num,ncol=7)
	ts <-matrix(nrow= boot.num,ncol=7)
	
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

	colnames_boot<-names(data.file)
	
	coeffs <- data.file
	
	WePa_coefs <- coeffs[,grepl(ordered_basins[2], colnames_boot)]
	EaPa_coefs <- coeffs[,grepl(ordered_basins[3], colnames_boot)]
	SoPa_coefs <- coeffs[,grepl(ordered_basins[4], colnames_boot)]
	NoIn_coefs <- coeffs[,grepl(ordered_basins[5], colnames_boot)]
	SoEIn_coefs <- coeffs[,grepl(ordered_basins[6], colnames_boot)]
	SoWIn_coefs <- coeffs[,grepl(ordered_basins[7], colnames_boot)]

	
		if(model.type=='1'){
		
		Intercepts<- cbind(coeffs[,1],
			coeffs[,1]+ WePa_coefs[,1],
			coeffs[,1]+ EaPa_coefs[,1],
			coeffs[,1]+ SoPa_coefs[,1],
			coeffs[,1]+ NoIn_coefs[,1],
			coeffs[,1]+ SoEIn_coefs[,1],
			coeffs[,1]+ SoWIn_coefs[,1])
			
		slopes<- cbind(coeffs[,2],
			coeffs[,2]+ WePa_coefs[,2],
			coeffs[,2]+ EaPa_coefs[,2],
			coeffs[,2]+ SoPa_coefs[,2],
			coeffs[,2]+ NoIn_coefs[,2],
			coeffs[,2]+ SoEIn_coefs[,2],
			coeffs[,2]+ SoWIn_coefs[,2])
			
		}
		
		if(model.type=='2'){
					
		Intercepts<- cbind(coeffs[,1],
			coeffs[,1]+ WePa_coefs[,1],
			coeffs[,1]+ EaPa_coefs[,1],
			coeffs[,1]+ SoPa_coefs[,1],
			coeffs[,1]+ NoIn_coefs[,1],
			coeffs[,1]+ SoEIn_coefs[,1],
			coeffs[,1]+ SoWIn_coefs[,1])
			
		slopes<- cbind(coeffs[,2],
			coeffs[,2]+ WePa_coefs[,2],
			coeffs[,2]+ EaPa_coefs[,2],
			coeffs[,2]+ SoPa_coefs[,2],
			coeffs[,2]+ NoIn_coefs[,2],
			coeffs[,2]+ SoEIn_coefs[,2],
			coeffs[,2]+ SoWIn_coefs[,2])
				
		slopesquared<- cbind(coeffs[,9],
			coeffs[,9]+ WePa_coefs[,3],
			coeffs[,9]+ EaPa_coefs[,3],
			coeffs[,9]+ SoPa_coefs[,3],
			coeffs[,9]+ NoIn_coefs[,3],
			coeffs[,9]+ SoEIn_coefs[,3],
			coeffs[,9]+ SoWIn_coefs[,3])
			
		}
			
		if(model.type=='3'){
					
		Intercepts<- cbind(coeffs[,1],
			coeffs[,1]+ WePa_coefs[,1],
			coeffs[,1]+ EaPa_coefs[,1],
			coeffs[,1]+ SoPa_coefs[,1],
			coeffs[,1]+ NoIn_coefs[,1],
			coeffs[,1]+ SoEIn_coefs[,1],
			coeffs[,1]+ SoWIn_coefs[,1])
			
		slopes<- cbind(coeffs[,2],
			coeffs[,2]+ WePa_coefs[,2],
			coeffs[,2]+ EaPa_coefs[,2],
			coeffs[,2]+ SoPa_coefs[,2],
			coeffs[,2]+ NoIn_coefs[,2],
			coeffs[,2]+ SoEIn_coefs[,2],
			coeffs[,2]+ SoWIn_coefs[,2])
				
		slopesquared<- cbind(coeffs[,9],
			coeffs[,9]+ WePa_coefs[,3],
			coeffs[,9]+ EaPa_coefs[,3],
			coeffs[,9]+ SoPa_coefs[,3],
			coeffs[,9]+ NoIn_coefs[,3],
			coeffs[,9]+ SoEIn_coefs[,3],
			coeffs[,9]+ SoWIn_coefs[,3])
						
		lat<-  cbind(coeffs[,10],
			coeffs[,10]+ WePa_coefs[,4],
			coeffs[,10]+ EaPa_coefs[,4],
			coeffs[,10]+ SoPa_coefs[,4],
			coeffs[,10]+ NoIn_coefs[,4],
			coeffs[,10]+ SoEIn_coefs[,4],
			coeffs[,10]+ SoWIn_coefs[,4])
			
			
		ts<- cbind(coeffs[,11],
			coeffs[,11]+ WePa_coefs[,5],
			coeffs[,11]+ EaPa_coefs[,5],
			coeffs[,11]+ SoPa_coefs[,5],
			coeffs[,11]+ NoIn_coefs[,5],
			coeffs[,11]+ SoEIn_coefs[,5],
			coeffs[,11]+ SoWIn_coefs[,5])
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
		
	if(model.type=='1'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope)
		boot.se.basins<-sqrt(sample.size.n/sample.size.m)* c(boot.se.basin.Int, boot.se.basin.slope)
	}
	if(model.type=='2'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, boot.est.basin.slopesquared)
		boot.se.basins<-sqrt(sample.size.n/sample.size.m)*c(boot.se.basin.Int, boot.se.basin.slope, boot.se.basin.slopesquared)
	}
	if(model.type=='3'){
		boot.est.basins<-c(boot.est.basin.Int, boot.est.basin.slope, boot.est.basin.slopesquared, boot.est.basin.avgLat, boot.est.basin.avgTS)
		boot.se.basins<-sqrt(sample.size.n/sample.size.m)*c(boot.se.basin.Int, boot.se.basin.slope, boot.se.basin.slopesquared, boot.se.basin.avgLat, boot.se.basin.avgTS)
	}
	
	return(list(Boot.est=boot.est.basins, Boot.SE=boot.se.basins))
}