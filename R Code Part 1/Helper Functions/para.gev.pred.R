#################################################################
#################################################################
#Parametric Bootstrap Function for prediction (MSPE and MAPE) for GEV univariate models

para.gev.pred<-function(boot.num=100,model.object=gev3,model.type='3',iterEach = 100,data= modeling_data, m.percent=1,rep=FALSE){
	MSPE <- vector(length = boot.num)
	MAPE <- vector(length = boot.num)
	
	
	if(model.type=='1'){			
		loc.par<-model.object$results$par[1:14]
		covariates<-model.matrix(~log(1013-minCP)*Basin,data)	
	}else if(model.type=='2'){
		loc.par<-model.object$results$par[1:21]
		covariates<-model.matrix(~log(1013-minCP)*Basin+
			I((log(1013-minCP))^2)*Basin,data)	
	}else {
		loc.par<-model.object$results$par[1:35]
		covariates<-model.matrix(~log(1013-minCP)*Basin+
			I((log(1013-minCP))^2)*Basin+ 
			avgLat*Basin+ Avg_Trans_Speed*Basin,data)	
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
	boot.y_iter <- matrix(nrow = iterEach, ncol = sample.size.m)
	data$boot.y<-NA
	
	for(i in 1:boot.num){
		print(i)
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
		
		for(j in 1:iterEach){
        boot.y_iter[j,] <- revd(sample.size.m,loc= boot.locs,scale= scale.par, shape= shape.par,type='GEV')
        boot.y_avg <- colMeans(boot.y_iter)
		data[c(NA_samp, NA_len+WP_samp, WP_len+ NA_len +EP_samp, 
		EP_len +WP_len+ NA_len +SP_samp, SP_len +EP_len +WP_len+ NA_len+ NI_samp, 
		NI_len + SP_len +EP_len +WP_len+ NA_len+ SEI_samp,
		SEI_len + NI_len + SP_len +EP_len +WP_len+ NA_len+ SWI_samp),]$boot.y<- boot.y_avg
		}


		data.boot<-na.omit(data)
    	
    	MSPE[i] <- mean(with(data.boot, (log(data.boot$maxWS_knots)-data.boot$boot.y)^2))
    	MAPE[i] <- mean(with(data.boot, abs(log(data.boot$maxWS_knots)-data.boot$boot.y)))

		}
	return(list(Avg_MSPE=mean(MSPE), SD_MSPE =sd(MSPE)/length(MSPE),
				Avg_MAPE=mean(MAPE),SD_MAPE =sd(MAPE)/length(MAPE)))
}


# para.gev.pred(boot.num = 100, model.object = gev_mod3, model.type = '3', iterEach = 100, data = modeling_data, m.percent = 0.66, rep = TRUE)

