
Basin_estimates<-function(ordered_basins= c('NoAt','WePa','EaPa','SoPa','NoIn','SoEIn','SoWIn'), fit1= all_basins_1010lm,fit2=NULL, Model_Type=1, exp.int=FALSE){

	Intercepts <-  NULL 
	Slopes <-  NULL  
	Slopes_sq<-  NULL 
	lat <-  NULL 	
	ts <-  NULL 
	
	total.parameters<- c(13,20,34)[Model_Type]
	
	colnames_fit1<-colnames(model.matrix(fit1))
	coef1<-coef(fit1)
	
	WePa_coefs <- coef1[grepl(ordered_basins[2], colnames_fit1)]
	EaPa_coefs <- coef1[grepl(ordered_basins[3], colnames_fit1)]
	SoPa_coefs <- coef1[grepl(ordered_basins[4], colnames_fit1)]
	NoIn_coefs <- coef1[grepl(ordered_basins[5], colnames_fit1)]
	SoEIn_coefs <- coef1[grepl(ordered_basins[6], colnames_fit1)]
	SoWIn_coefs <- coef1[grepl(ordered_basins[7], colnames_fit1)]
	
	if(exp.int==FALSE){
	Intercept_exp_SE_M1<-c(
	deltaMethod(fit1, "b0", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
	deltaMethod(fit1, "b0+b2", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b0+b3", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b0+b5", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b0+b4", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b0+b6", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b0+b7", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
	}else{
		Intercept_exp_SE_M1<-c(
	deltaMethod(fit1, "exp(b0)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
	deltaMethod(fit1, "exp(b0+b2)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "exp(b0+b3)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "exp(b0+b5)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "exp(b0+b4)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "exp(b0+b6)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "exp(b0+b7)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
	}
	
	Intercepts<- data.frame(Model='WLS',Basin=ordered_basins, 
			Model_Type= Model_Type, Par='beta0',
			Estimate=c(coef1[1],
			 coef1[1]+ WePa_coefs[1],
			 coef1[1]+ EaPa_coefs[1],
			 coef1[1]+ SoPa_coefs[1],
			 coef1[1]+ NoIn_coefs[1],
			 coef1[1]+ SoEIn_coefs[1],
			 coef1[1]+ SoWIn_coefs[1]),
			SE= Intercept_exp_SE_M1)
			 	
	Slope_SE_M1<-c(
	deltaMethod(fit1, "b1", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
	deltaMethod(fit1, "b1+b8", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b1+b9", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b1+b11", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b1+b10", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b1+b12", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
	deltaMethod(fit1, "b1+b13", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
		 	
	Slopes<- data.frame(Model='WLS',Basin=ordered_basins, 
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
		
		Slope_SE_M1<-c(
		deltaMethod(fit1, "b1", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
		deltaMethod(fit1, "b1+b9", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b1+b10", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b1+b12", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b1+b11", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b1+b13", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b1+b14", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])

		Slopes$SE <- Slope_SE_M1
		
		Slopes_sq_SE_M1<-c(
		deltaMethod(fit1, "b8", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
		deltaMethod(fit1, "b8+b15", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b8+b16", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b8+b18", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b8+b17", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b8+b19", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit1, "b8+b20", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
		
		Slopes_sq<- data.frame(Model='WLS',Basin=ordered_basins, 
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
			 			
				Slope_SE_M1<-c(
				deltaMethod(fit1, "b1", 
					parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
				deltaMethod(fit1, "b1+b11", 
					parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
				deltaMethod(fit1, "b1+b12", 
					parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
				deltaMethod(fit1, "b1+b14", 
					parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
				deltaMethod(fit1, "b1+b13", 
					parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
				deltaMethod(fit1, "b1+b15", 
					parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
				deltaMethod(fit1, "b1+b16", 
					parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])

				Slopes$SE <- Slope_SE_M1
				
				Slopes_sq_SE_M1<-c(
					deltaMethod(fit1, "b8", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
					deltaMethod(fit1, "b8+b17", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b8+b18", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b8+b20", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b8+b19", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b8+b21", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b8+b22", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
		
				Slopes_sq$SE <- Slopes_sq_SE_M1
				
				lat_SE_M1<-c(
					deltaMethod(fit1, "b9", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
					deltaMethod(fit1, "b9+b23", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b9+b24", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b9+b26", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b9+b25", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b9+b27", 
						parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b9+b28", 
						parameterNames = paste("b", 0: total.parameters, sep=""))[[2]])

				lat<- data.frame(Model='WLS',Basin=ordered_basins, 
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
					deltaMethod(fit1, "b10", 
						parameterNames= paste("b", 0: total.parameters, sep=""))[[2]], 
					deltaMethod(fit1, "b10+b29", 
						parameterNames= paste("b", 0: total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b10+b30", 
						parameterNames= paste("b", 0: total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b10+b32", 
						parameterNames= paste("b", 0: total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b10+b31", 
						parameterNames= paste("b", 0: total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b10+b33", 
						parameterNames= paste("b", 0: total.parameters, sep=""))[[2]],
					deltaMethod(fit1, "b10+b34", 
						parameterNames= paste("b", 0: total.parameters, sep=""))[[2]])

			 	ts<- data.frame(Model='WLS',Basin=ordered_basins, 
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
			 
	if(!is.null(fit2)){
		colnames_fit2 <-colnames(model.matrix(fit2))
		coef2<-coef(fit2)
		
		WePa_coefs2 <- coef2[grepl(ordered_basins[2], colnames_fit2)]
		EaPa_coefs2 <- coef2[grepl(ordered_basins[3], colnames_fit2)]
		SoPa_coefs2 <- coef2[grepl(ordered_basins[4], colnames_fit2)]
		NoIn_coefs2 <- coef2[grepl(ordered_basins[5], colnames_fit2)]
		SoEIn_coefs2 <- coef2[grepl(ordered_basins[6], colnames_fit2)]
		SoWIn_coefs2 <- coef2[grepl(ordered_basins[7], colnames_fit2)]
		
		if(exp.int==FALSE){
			Intercept_exp_SE_M1_ms <-c(
			deltaMethod(fit1, "b0", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
			deltaMethod(fit1, "b0+b2", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit1, "b0+b3", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit1, "b0+b5", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit1, "b0+b4", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit1, "b0+b6", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit1, "b0+b7", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
		}else{
		Intercept_exp_SE_M1_ms<-c(
		deltaMethod(fit2,"exp(b0)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
		deltaMethod(fit2,"exp(b0+b2)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2,"exp(b0+b3)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2,"exp(b0+b5)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2,"exp(b0+b4)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2,"exp(b0+b6)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2,"exp(b0+b7)", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
		}
		
		Intercepts_ms<- data.frame(Model='WLS',Basin=ordered_basins, 
			Model_Type= Model_Type, Par='beta0',
			Estimate=c(coef2[1],
			 coef2[1]+ WePa_coefs2[1],
			 coef2[1]+ EaPa_coefs2[1],
			 coef2[1]+ SoPa_coefs2[1],
			 coef2[1]+ NoIn_coefs2[1],
			 coef2[1]+ SoEIn_coefs2[1],
			 coef2[1]+ SoWIn_coefs2[1]),
			 SE= Intercept_exp_SE_M1_ms)
			 			 
		Slope_SE_ms<-c(
		deltaMethod(fit2, "b1", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
		deltaMethod(fit2, "b1+b8", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2, "b1+b9", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2, "b1+b11", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2, "b1+b10", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2, "b1+b12", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
		deltaMethod(fit2, "b1+b13", parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])
		 	
		Slopes_ms<- data.frame(Model='WLS',Basin=ordered_basins, 
			Model_Type= Model_Type,Par='beta1',
			Estimate=c(coef2[2],
		 	 coef2[2]+ WePa_coefs2[2],
	 		 coef2[2]+ EaPa_coefs2[2],
			 coef2[2]+ SoPa_coefs2[2],
			 coef2[2]+ NoIn_coefs2[2],
			 coef2[2]+ SoEIn_coefs2[2],
			 coef2[2]+ SoWIn_coefs2[2]),
			 SE= Slope_SE_ms)
			 
		if(Model_Type > 1){
			
			Slopes_sq_SE_M1_ms<-c(
			deltaMethod(fit2, "b8", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]], 
			deltaMethod(fit2, "b8+b15", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit2, "b8+b16", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit2, "b8+b18", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit2, "b8+b17", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit2, "b8+b19", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]],
			deltaMethod(fit2, "b8+b20", 
			parameterNames= paste("b", 0:total.parameters, sep=""))[[2]])

			
			Slopes_sq_ms<- data.frame(Model='WLS',Basin=ordered_basins, 
			Model_Type= Model_Type,Par='beta2',
			Estimate=c(coef2[9],
		 	 coef2[9]+ WePa_coefs2[3],
	 		 coef2[9]+ EaPa_coefs2[3],
			 coef2[9]+ SoPa_coefs2[3],
			 coef2[9]+ NoIn_coefs2[3],
			 coef2[9]+ SoEIn_coefs2[3],
			 coef2[9]+ SoWIn_coefs2[3]),
			 SE= Slopes_sq_SE_M1_ms)
			 
			
			 if(Model_Type > 2){
				lat_ms<- data.frame(Model='WLS',Basin=ordered_basins, 
					Model_Type= Model_Type,Par='beta3',
					Estimate=c(coef2[10],
		 			coef2[10]+ WePa_coefs2[4],
	 		 		coef2[10]+ EaPa_coefs2[4],
			 		coef2[10]+ SoPa_coefs2[4],
			 		coef2[10]+ NoIn_coefs2[4],
			 		coef2[10]+ SoEIn_coefs2[4],
			 		coef2[10]+ SoWIn_coefs2[4]))
			 		
			 	ts_ms<- data.frame(Model='WLS',Basin=ordered_basins, 
					Model_Type= Model_Type,Par='beta4',
					Estimate=c(coef2[11],
		 			coef2[11]+ WePa_coefs2[5],
	 		 		coef2[11]+ EaPa_coefs2[5],
			 		coef2[11]+ SoPa_coefs2[5],
			 		coef2[11]+ NoIn_coefs2[5],
			 		coef2[11]+ SoEIn_coefs2[5],
			 		coef2[11]+ SoWIn_coefs2[5]))
			 	
			 	return(list(Intercepts_kt = Intercepts, 
			 	Slopes_kt = Slopes ,
			 	Slopes_sq_kt= Slopes_sq, lat_kt= lat, ts_kt= ts,
				Intercepts_ms = Intercepts_ms, Slopes_ms = Slopes_ms,
				Slopes_sq_ms= Slopes_sq_ms, lat_ms= lat_ms, ts_ms= ts_ms))		 
			}
			return(list(Intercepts_kt = Intercepts, 
				Slopes_kt = Slopes ,
			 	Slopes_sq_kt= Slopes_sq,
				Intercepts_ms = Intercepts_ms, Slopes_ms = Slopes_ms,
				Slopes_sq_ms= Slopes_sq_ms))		 
		}		
		return(list(Intercepts_kt = Intercepts, 
				Slopes_kt = Slopes ,
				Intercepts_ms = Intercepts_ms, Slopes_ms = Slopes_ms))		  
			 
	}else{
	 	return(list(Intercepts_kt = Intercepts, 
	 	Slopes_kt = Slopes ,
	 	Slopes_sq_kt= Slopes_sq, lat_kt= lat, ts_kt= ts))		 

	}
}