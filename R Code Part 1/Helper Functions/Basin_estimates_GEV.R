
Basin_estimates_GEV <-function(ordered_basins=c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn'), 
	fit1= gev_mod3, data= modeling_data, Model_Type=3){

	Intercepts <-  NULL 
	Slopes <-  NULL  
	Slopes_sq<-  NULL 
	lat <-  NULL 	
	ts <-  NULL 
	shape <- NULL
	scale <- NULL
	
	total.parameters<- c(16,23,37)[Model_Type]
	
	colnames_fit1<- c(colnames(model.matrix(fit1$par.models$location, data)), 'shape','scale')
	coef1<-fit1$results$par
	names(coef1)<-colnames_fit1
	coef.cov<-summary(fit1)$ cov.theta
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
	
	Intercept_SE_M1<-c(
	deltamethod(as.formula(paste('~x1')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x1+x', WePa_coefs_loc[1],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x1+x', EaPa_coefs_loc[1],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x1+x', SoPa_coefs_loc[1],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x1+x', NoIn_coefs_loc[1],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x1+x', SoEIn_coefs_loc[1],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x1+x', SoWIn_coefs_loc[1],sep='')), coef1, coef.cov, ses=TRUE))
		
	Intercepts<- data.frame(Model='GEV',Basin=ordered_basins, 
			Model_Type= Model_Type, Par='beta0',
			Estimate=c(coef1[1],
			 coef1[1]+ WePa_coefs[1],
			 coef1[1]+ EaPa_coefs[1],
			 coef1[1]+ SoPa_coefs[1],
			 coef1[1]+ NoIn_coefs[1],
			 coef1[1]+ SoEIn_coefs[1],
			 coef1[1]+ SoWIn_coefs[1]),
			SE= Intercept_SE_M1)
			 	
	Slope_SE_M1<-c(
	deltamethod(as.formula(paste('~x2')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x2+x', WePa_coefs_loc[2],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x2+x', EaPa_coefs_loc[2],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x2+x', SoPa_coefs_loc[2],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x2+x', NoIn_coefs_loc[2],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x2+x', SoEIn_coefs_loc[2],sep='')), coef1, coef.cov, ses=TRUE), 
	deltamethod(as.formula(paste('~x2+x', SoWIn_coefs_loc[2],sep='')), coef1, coef.cov, ses=TRUE))
		 	
	Slopes<- data.frame(Model='GEV',Basin=ordered_basins, 
			Model_Type= Model_Type,Par='beta1',
			Estimate=c(coef1[2],
		 	 coef1[2]+ WePa_coefs[2],
	 		 coef1[2]+ EaPa_coefs[2],
			 coef1[2]+ SoPa_coefs[2],
			 coef1[2]+ NoIn_coefs[2],
			 coef1[2]+ SoEIn_coefs[2],
			 coef1[2]+ SoWIn_coefs[2]),
			 SE= Slope_SE_M1)
	
	if(Model_Type > 1){

		Slopes_sq_SE_M1<-c(
		deltamethod(as.formula(paste('~x9')), coef1, coef.cov, ses=TRUE), 
		deltamethod(as.formula(paste('~x9+x', WePa_coefs_loc[3],sep='')), 
			coef1, coef.cov, ses=TRUE),
		deltamethod(as.formula(paste('~x9+x', EaPa_coefs_loc[3],sep='')), 
			coef1, coef.cov, ses=TRUE), 
		deltamethod(as.formula(paste('~x9+x', SoPa_coefs_loc[3],sep='')), 
			coef1, coef.cov, ses=TRUE), 
		deltamethod(as.formula(paste('~x9+x', NoIn_coefs_loc[3],sep='')), 
			coef1, coef.cov, ses=TRUE), 
		deltamethod(as.formula(paste('~x9+x', SoEIn_coefs_loc[3],sep='')), 
			coef1, coef.cov, ses=TRUE), 
		deltamethod(as.formula(paste('~x9+x', SoWIn_coefs_loc[3],sep='')), 
			coef1, coef.cov, ses=TRUE))

		Slopes_sq<- data.frame(Model='GEV',Basin=ordered_basins, 
			Model_Type= Model_Type,Par='beta2',
			Estimate=c(coef1[9],
		 	 coef1[9]+ WePa_coefs[3],
	 		 coef1[9]+ EaPa_coefs[3],
			 coef1[9]+ SoPa_coefs[3],
			 coef1[9]+ NoIn_coefs[3],
			 coef1[9]+ SoEIn_coefs[3],
			 coef1[9]+ SoWIn_coefs[3]),
			 SE= Slopes_sq_SE_M1)			 
			 
			 if(Model_Type > 2){
			 				
				lat_SE_M1<-c(
				deltamethod(as.formula(paste('~x10')), coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~x10+x', WePa_coefs_loc[4],sep='')), 
					coef1, coef.cov, ses=TRUE),
				deltamethod(as.formula(paste('~ x10 +x', EaPa_coefs_loc[4],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x10 +x', SoPa_coefs_loc[4],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x10 +x', NoIn_coefs_loc[4],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x10 +x', SoEIn_coefs_loc[4],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x10 +x', SoWIn_coefs_loc[4],sep='')), 
					coef1, coef.cov, ses=TRUE))

				lat<- data.frame(Model='GEV',Basin=ordered_basins, 
					Model_Type= Model_Type,Par='beta3',
					Estimate=c(coef1[10],
		 			coef1[10]+ WePa_coefs[4],
	 		 		coef1[10]+ EaPa_coefs[4],
			 		coef1[10]+ SoPa_coefs[4],
			 		coef1[10]+ NoIn_coefs[4],
			 		coef1[10]+ SoEIn_coefs[4],
			 		coef1[10]+ SoWIn_coefs[4]),
			 		SE= lat_SE_M1)
			 		
			 	ts_SE_M1<-c(
				deltamethod(as.formula(paste('~x11')), coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x11 +x', WePa_coefs_loc[5],sep='')), 
					coef1, coef.cov, ses=TRUE),
				deltamethod(as.formula(paste('~ x11 +x', EaPa_coefs_loc[5],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x11 +x', SoPa_coefs_loc[5],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x11 +x', NoIn_coefs_loc[5],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x11 +x', SoEIn_coefs_loc[5],sep='')), 
					coef1, coef.cov, ses=TRUE), 
				deltamethod(as.formula(paste('~ x11 +x', SoWIn_coefs_loc[5],sep='')), 
					coef1, coef.cov, ses=TRUE))

			 	ts<- data.frame(Model='GEV',Basin=ordered_basins, 
					Model_Type= Model_Type,Par='beta4',
					Estimate=c(coef1[11],
		 			coef1[11]+ WePa_coefs[5],
	 		 		coef1[11]+ EaPa_coefs[5],
			 		coef1[11]+ SoPa_coefs[5],
			 		coef1[11]+ NoIn_coefs[5],
			 		coef1[11]+ SoEIn_coefs[5],
			 		coef1[11]+ SoWIn_coefs[5]),
			 		SE= ts_SE_M1)
			}
	}		 
			 
	 	return(list(Intercepts_kt = Intercepts, Slopes_kt = Slopes ,
	 	Slopes_sq_kt= Slopes_sq, lat_kt= lat, ts_kt= ts))		 
}