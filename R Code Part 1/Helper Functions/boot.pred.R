boot.pred<-function(file='wlm3.boot.all.csv',modeldat=modeling_data, model.type='3'){

	data<-read.csv(file, check.names = FALSE)[,-1]
	MSPE_ws<-vector()
	MAPE_ws<-vector()
	MSPE_cp<-vector()
	MAPE_cp<-vector()

	for(i in 1:dim(data)[1]){

		if(model.type=='1'){			
		  covariates<-model.matrix(~log(1013-minCP)*Basin, modeldat)	
   		  loc.par<-data[i,1:14]
		} else if(model.type=='2'){
		  covariates<-model.matrix(~log(1013-minCP)*Basin+
			I((log(1013-minCP))^2)*Basin, modeldat)	
			loc.par<-data[i,1:21]
		} else if(model.type=='BVEVD'){
		  covariates<-model.matrix(~Basin, modeldat)	
			loc.par1<-data[i,c(1:7)]
			loc.par2<-data[i,c(10:16)]
		} else {
		  covariates<-model.matrix(~log(1013-minCP)*Basin+
			I((log(1013-minCP))^2)*Basin+ 
			avgLat*Basin+ Avg_Trans_Speed*Basin, modeldat)	
			loc.par<-data[i,1:35]
		}
		
		if(model.type!='BVEVD'){	 	
			modeldat$loc.pred<- covariates %*% t(as.matrix(loc.par))
			MSPE_ws[i]<-mean(with(modeldat, (log(maxWS_knots)-loc.pred)^2))
			MAPE_ws[i]<-mean(with(modeldat, abs(log(maxWS_knots)-loc.pred)))
			
		} else{
			modeldat$loc.pred1<- covariates %*% t(as.matrix(loc.par1))
			modeldat$loc.pred2<- covariates %*% t(as.matrix(loc.par2))
			
			MSPE_ws[i]<-mean(with(modeldat, (log(maxWS_knots)-loc.pred1)^2))
			MAPE_ws[i]<-mean(with(modeldat, abs(log(maxWS_knots)-loc.pred1)))

			MSPE_cp[i]<-mean(with(modeldat, (log(1013-minCP)-loc.pred2)^2))
			MAPE_cp[i]<-mean(with(modeldat, abs(log(1013-minCP)-loc.pred2)))
			
		}	
	}
		
		
		if(model.type!='BVEVD'){	 	
		return(list(Avg_MSPE_ws= mean(MSPE_ws), SD_n_MSPE_ws =sd(MSPE_ws)/length(MSPE_ws),
				Avg_MAPE_ws = mean(MAPE_ws),SD_n_MAPE_ws =sd(MAPE_ws)/length(MAPE_ws) 
				#MSPE_ws = MSPE_ws, MAPE_ws = MAPE_ws
				))
		} else{
		return(list(Avg_MSPE_ws=mean(MSPE_ws), SD_n_MSPE_ws =sd(MSPE_ws)/length(MSPE_ws),
				Avg_MAPE_ws =mean(MAPE_ws),SD_n_MAPE_ws =sd(MAPE_ws)/length(MAPE_ws), 
				#MSPE_ws = MSPE_ws, MAPE_ws = MAPE_ws, 
				Avg_MSPE_cp=mean(MSPE_cp), SD_n_MSPE_cp =sd(MSPE_cp)/length(MSPE_cp),
				Avg_MAPE_cp =mean(MAPE_cp),SD_n_MAPE_cp =sd(MAPE_cp)/length(MAPE_cp)
				#MSPE_cp = MSPE_cp, MAPE_cp = MAPE_cp
				))
		}
}		