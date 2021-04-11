BVEVD_Boot<-function(data= modeling_data, 
	data.file=read.csv('bvevd.66.boot.all.csv',header=T, check.names = FALSE)[,-1], 
	m.percent =0.66, ordered_basins= c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn')){

	boot.num =dim(data.file)[1]
	
	Intercepts_M1<-matrix(nrow= boot.num,ncol=7)
	Intercepts_M2<-matrix(nrow= boot.num,ncol=7)
	scale1.vec<-vector(length= boot.num)
	scale2.vec<-vector(length= boot.num)
	shape1.vec<-vector(length= boot.num)
	shape2.vec<-vector(length= boot.num)
	alpha1.vec<-vector(length= boot.num)
	beta2.vec<-vector(length= boot.num)
	
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
	
	noats<-!(grepl(ordered_basins[2], colnames_boot)+
	  grepl(ordered_basins[3], colnames_boot)+
	  grepl(ordered_basins[4], colnames_boot)+
	  grepl(ordered_basins[5], colnames_boot)+
	  grepl(ordered_basins[6], colnames_boot)+
	  grepl(ordered_basins[7], colnames_boot)+
	  grepl('scale', colnames_boot)+
	  grepl('shape', colnames_boot)+
	  grepl('alpha', colnames_boot)+
	  grepl('beta', colnames_boot))
	
	NoAt_coefs <- coeffs[, noats]
	WePa_coefs <- coeffs[,grepl(ordered_basins[2], colnames_boot)]
	EaPa_coefs <- coeffs[,grepl(ordered_basins[3], colnames_boot)]
	SoPa_coefs <- coeffs[,grepl(ordered_basins[4], colnames_boot)]
	NoIn_coefs <- coeffs[,grepl(ordered_basins[5], colnames_boot)]
	SoEIn_coefs <- coeffs[,grepl(ordered_basins[6], colnames_boot)]
	SoWIn_coefs <- coeffs[,grepl(ordered_basins[7], colnames_boot)]

	Intercepts_M1 <- cbind(NoAt_coefs[,1],
			NoAt_coefs[,1]+WePa_coefs[,1],
			NoAt_coefs[,1]+EaPa_coefs[,1],
			NoAt_coefs[,1]+SoPa_coefs[,1],
			NoAt_coefs[,1]+NoIn_coefs[,1],
			NoAt_coefs[,1]+SoEIn_coefs[,1],
			NoAt_coefs[,1]+SoWIn_coefs[,1])
			
	Intercepts_M2 <- cbind(NoAt_coefs[,2],
			NoAt_coefs[,2] + WePa_coefs[,2],
			NoAt_coefs[,2] + EaPa_coefs[,2],
			NoAt_coefs[,2] + SoPa_coefs[,2],
			NoAt_coefs[,2] + NoIn_coefs[,2],
			NoAt_coefs[,2] + SoEIn_coefs[,2],
			NoAt_coefs[,2] + SoWIn_coefs[,2])
					
	scale1.vec<- coeffs[,'scale1']
	shape1.vec<- coeffs[,'shape1']
	alpha1.vec<- coeffs[,'alpha']
	
	scale2.vec<- coeffs[,'scale2']
	shape2.vec<- coeffs[,'shape2']
	beta2.vec<- coeffs[,'beta']


	boot.est.basin.IntM1<-colMeans(Intercepts_M1)
	boot.se.basin.IntM1<-apply(Intercepts_M1, 2, sd)
	names(boot.est.basin.IntM1)<-c('NoAt Int. M1','WePa Int. M1' ,
	'EaPa Int. M1' ,'SoPa Int. M1' ,'NoIn Int. M1' ,'SoEIn Int. M1','SoWIn Int. M1')
	names(boot.se.basin.IntM1)<-c('NoAt Int. M1','WePa Int. M1' ,
	'EaPa Int. M1' ,'SoPa Int. M1' ,'NoIn Int. M1' ,'SoEIn Int. M1','SoWIn Int. M1')

	boot.est.basin.IntM2<-colMeans(Intercepts_M2)
	boot.se.basin.IntM2<-apply(Intercepts_M2, 2, sd)
	names(boot.est.basin.IntM2)<-c('NoAt Int. M2','WePa Int. M2' ,
	'EaPa Int. M2' ,'SoPa Int. M2' ,'NoIn Int. M2' ,'SoEIn Int. M2','SoWIn Int. M2')
	names(boot.se.basin.IntM2)<-c('NoAt Int. M2','WePa Int. M2' ,
	'EaPa Int. M2' ,'SoPa Int. M2' ,'NoIn Int. M2' ,'SoEIn Int. M2','SoWIn Int. M2')

	
	boot.est.scale1<-mean(scale1.vec)
	boot.se.scale1<-sd(scale1.vec)
	
	boot.est.shape1<-mean(shape1.vec)
	boot.se.shape1<-sd(shape1.vec)
		
	boot.est.alpha1<-mean(alpha1.vec)
	boot.se.alpha1<-sd(alpha1.vec)

	boot.est.scale2<-mean(scale2.vec)
	boot.se.scale2<-sd(scale2.vec)
	
	boot.est.shape2<-mean(shape2.vec)
	boot.se.shape2<-sd(shape2.vec)
	
	boot.est.beta2<-mean(beta2.vec)
	boot.se.beta2<-sd(beta2.vec)
	
		
	sample.size.n<-dim(data)[1]
		
	boot.est.basins<-c(boot.est.basin.IntM1, scale1= boot.est.scale1,
			shape1= boot.est.shape1, alpha1= boot.est.alpha1, 
		boot.est.basin.IntM2, scale2= boot.est.scale2, 
			shape2= boot.est.shape2, beta2= boot.est.beta2)
		
	boot.se.basins<- sqrt(sample.size.n/sample.size.m)* 
		c(boot.se.basin.IntM1, scale1= boot.se.scale1, 
			shape1= boot.se.shape1, alpha1= boot.se.alpha1, 
		boot.se.basin.IntM2, scale2= boot.se.scale2, 
			shape2= boot.se.shape2, beta2= boot.se.beta2)
	
	return(list(Boot.est=boot.est.basins, Boot.SE=boot.se.basins))
}