#################################################################
#################################################################
#Delta Method Function

delta.method.wsp<-function(model.object=lm3,model.type='lm',model.num=3,num.par=c(13,20,34)[3]){
	B<-c('NA','WP','EP','SP','NI','SEI','SWI')

  if(model.type=='lm'){
  	coefs<-coef(model.object)
  	
	Intercepts<- c(coefs[1],
			coefs[1]+coefs[1+2],
			coefs[1]+coefs[1+3],
			coefs[1]+coefs[1+4],
			coefs[1]+coefs[1+5],
			coefs[1]+coefs[1+6],
			coefs[1]+coefs[1+7])

	Intercept_SE<-c(
	deltaMethod(model.object, "b0", parameterNames= paste("b", 0:num.par, sep=""))[2], 
	deltaMethod(model.object, "b0+b2", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b0+b3", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b0+b4", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b0+b5", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b0+b6", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b0+b7", parameterNames= paste("b", 0: num.par, sep=""))[2])
	
	 if(model.num==1){
	  slopes<- c(coefs[2],
			coefs[2]+coefs[7+2],
			coefs[2]+coefs[7+3],
			coefs[2]+coefs[7+4],
			coefs[2]+coefs[7+5],
			coefs[2]+coefs[7+6],
			coefs[2]+coefs[7+7])

	slopes_SE<-c(
	deltaMethod(model.object, "b1", parameterNames= paste("b", 0: num.par, sep=""))[2], 
	deltaMethod(model.object, "b1+b8", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b9", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b10", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b11", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b12", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b13", parameterNames= paste("b", 0: num.par, sep=""))[2])

  	  	table.print<-data.frame('Basin'=B,
  	"$\\hat{\\beta_0}$"=paste(format(round(c(Intercepts),3),3), 
  	' (',format(round(unlist(Intercept_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_1}$"=paste(format(round(c(slopes),3),3), 
  	' (', format(round(unlist(slopes_SE),3)),')',sep=''), 
	check.names = FALSE)
     
     print(xtable(table.print),sanitize.text.function = function(x){x})

  	}
	
	 if(model.num==2){
	  slopes<- c(coefs[2],
			coefs[2]+coefs[1+7+2],
			coefs[2]+coefs[1+7+3],
			coefs[2]+coefs[1+7+4],
			coefs[2]+coefs[1+7+5],
			coefs[2]+coefs[1+7+6],
			coefs[2]+coefs[1+7+7])

	slopes_SE<-c(
	deltaMethod(model.object, "b1", parameterNames= paste("b", 0: num.par, sep=""))[2], 
	deltaMethod(model.object, "b1+b9", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b10", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b11", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b12", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b13", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b14", parameterNames= paste("b", 0: num.par, sep=""))[2])

	slopesquared<- c(coefs[9],
			coefs[9]+coefs[9+5+2],
			coefs[9]+coefs[9+5+3],
			coefs[9]+coefs[9+5+4],
			coefs[9]+coefs[9+5+5],
			coefs[9]+coefs[9+5+6],
			coefs[9]+coefs[9+5+7])

	slopesquared_SE<-c(
	deltaMethod(model.object, "b8", parameterNames= paste("b", 0: num.par, sep=""))[2], 
	deltaMethod(model.object, "b8+b15", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b16", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b17", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b18", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b19", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b20", parameterNames= paste("b", 0: num.par, sep=""))[2])

  	  	table.print<-data.frame('Basin'=B,
  	"$\\hat{\\beta_0}$"=paste(format(round(c(Intercepts),3),3), 
  	' (',format(round(unlist(Intercept_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_1}$"=paste(format(round(c(slopes),3),3), 
  	' (', format(round(unlist(slopes_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_2}$"=paste(format(round(c(slopesquared),3),3), 
  	' (',format(round(unlist(slopesquared_SE),3)),')',sep=''), 
	check.names = FALSE)
     
     print(xtable(table.print),sanitize.text.function = function(x){x})

  	}
  	
  	if(model.num==3){
	  slopes<- c(coefs[2],
			coefs[2]+coefs[1+9+2],
			coefs[2]+coefs[1+9+3],
			coefs[2]+coefs[1+9+4],
			coefs[2]+coefs[1+9+5],
			coefs[2]+coefs[1+9+6],
			coefs[2]+coefs[1+9+7])

	slopes_SE<-c(
	deltaMethod(model.object, "b1", parameterNames= paste("b", 0: num.par, sep=""))[2], 
	deltaMethod(model.object, "b1+b11", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b12", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b13", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b14", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b15", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b1+b16", parameterNames= paste("b", 0: num.par, sep=""))[2])

	slopesquared<- c(coefs[9],
			coefs[9]+coefs[9+7+2],
			coefs[9]+coefs[9+7+3],
			coefs[9]+coefs[9+7+4],
			coefs[9]+coefs[9+7+5],
			coefs[9]+coefs[9+7+6],
			coefs[9]+coefs[9+7+7])

	slopesquared_SE<-c(
	deltaMethod(model.object, "b8", parameterNames= paste("b", 0: num.par, sep=""))[2], 
	deltaMethod(model.object, "b8+b17", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b18", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b19", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b20", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b21", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b8+b22", parameterNames= paste("b", 0: num.par, sep=""))[2])

	lat<- c(coefs[10],
			coefs[10]+coefs[10+12+2],
			coefs[10]+coefs[10+12+3],
			coefs[10]+coefs[10+12+4],
			coefs[10]+coefs[10+12+5],
			coefs[10]+coefs[10+12+6],
			coefs[10]+coefs[10+12+7])

	lat_SE<-c(
	deltaMethod(model.object, "b9", parameterNames= paste("b", 0: num.par, sep=""))[2], 
	deltaMethod(model.object, "b9+b23", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b9+b24", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b9+b25", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b9+b26", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b9+b27", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b9+b28", parameterNames= paste("b", 0: num.par, sep=""))[2])

	ts<- c(coefs[11],
			coefs[11]+coefs[11+17+2],
			coefs[11]+coefs[11+17+3],
			coefs[11]+coefs[11+17+4],
			coefs[11]+coefs[11+17+5],
			coefs[11]+coefs[11+17+6],
			coefs[11]+coefs[11+17+7])

	ts_SE<-c(
	deltaMethod(model.object, "b10", parameterNames= paste("b", 0: num.par, sep=""))[2], 
	deltaMethod(model.object, "b10+b29", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b10+b30", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b10+b31", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b10+b32", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b10+b33", parameterNames= paste("b", 0: num.par, sep=""))[2],
	deltaMethod(model.object, "b10+b34", parameterNames= paste("b", 0: num.par, sep=""))[2])
  	
  	
  	table.print<-data.frame('Basin'=B,
  	"$\\hat{\\beta_0}$"=paste(format(round(c(Intercepts),3),3), 
  	' (',format(round(unlist(Intercept_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_1}$"=paste(format(round(c(slopes),3),3), 
  	' (', format(round(unlist(slopes_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_2}$"=paste(format(round(c(slopesquared),3),3), 
  	' (',format(round(unlist(slopesquared_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_3}$"=paste(format(round(c(lat),3),3),  
  	' (',format(round(unlist(lat_SE),3)),')',sep=''),
  	"$\\hat{\\beta_4}$"=paste(format(round(c(ts),3),3), 
  	' (',format(round(unlist(ts_SE),3)),')',sep=''),
     check.names = FALSE)
     
     print(xtable(table.print),sanitize.text.function = function(x){x})
  	
  	}
  }
  
  if(model.type=='ev'){
  	coefs<-model.object$results$par
  	covs<-summary(model.object)$cov.theta

  Intercepts<- c(coefs[1],
			coefs[1]+coefs[1+2],
			coefs[1]+coefs[1+3],
			coefs[1]+coefs[1+4],
			coefs[1]+coefs[1+5],
			coefs[1]+coefs[1+6],
			coefs[1]+coefs[1+7])

	Intercept_SE<-c(
	deltamethod(~x1, coefs, covs, ses=TRUE), 
	deltamethod(~x1+x3, coefs, covs, ses=TRUE),
	deltamethod(~x1+x4, coefs, covs, ses=TRUE),
	deltamethod(~x1+x5, coefs, covs, ses=TRUE),
	deltamethod(~x1+x6, coefs, covs, ses=TRUE),
	deltamethod(~x1+x7, coefs, covs, ses=TRUE),
	deltamethod(~x1+x8, coefs, covs, ses=TRUE))
	
	 if(model.num==1){
	  slopes<- c(coefs[2],
			coefs[2]+coefs[7+2],
			coefs[2]+coefs[7+3],
			coefs[2]+coefs[7+4],
			coefs[2]+coefs[7+5],
			coefs[2]+coefs[7+6],
			coefs[2]+coefs[7+7])

	slopes_SE<-c(
	deltamethod(~x2, coefs, covs, ses=TRUE), 
	deltamethod(~x2+x9, coefs, covs, ses=TRUE),
	deltamethod(~x2+x10, coefs, covs, ses=TRUE),
	deltamethod(~x2+x11, coefs, covs, ses=TRUE),
	deltamethod(~x2+x12, coefs, covs, ses=TRUE),
	deltamethod(~x2+x13, coefs, covs, ses=TRUE),
	deltamethod(~x2+x14, coefs, covs, ses=TRUE))

  	  	table.print<-data.frame('Basin'=B,
  	"$\\hat{\\beta_0}$"=paste(format(round(c(Intercepts),3),3), 
  	' (',format(round(unlist(Intercept_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_1}$"=paste(format(round(c(slopes),3),3), 
  	' (', format(round(unlist(slopes_SE),3)),')',sep=''), 
	check.names = FALSE)
     
     print(xtable(table.print),sanitize.text.function = function(x){x})
  			return(pars=c(Intercept_SE,slopes_SE))

  	}
	
	 if(model.num==2){
	  slopes<- c(coefs[2],
			coefs[2]+coefs[1+7+2],
			coefs[2]+coefs[1+7+3],
			coefs[2]+coefs[1+7+4],
			coefs[2]+coefs[1+7+5],
			coefs[2]+coefs[1+7+6],
			coefs[2]+coefs[1+7+7])

	slopes_SE<-c(
	deltamethod(~x2, coefs, covs, ses=TRUE), 
	deltamethod(~x2+x10, coefs, covs, ses=TRUE),
	deltamethod(~x2+x11, coefs, covs, ses=TRUE),
	deltamethod(~x2+x12, coefs, covs, ses=TRUE),
	deltamethod(~x2+x13, coefs, covs, ses=TRUE),
	deltamethod(~x2+x14, coefs, covs, ses=TRUE),
	deltamethod(~x2+x15, coefs, covs, ses=TRUE))

	slopesquared<- c(coefs[9],
			coefs[9]+coefs[9+5+2],
			coefs[9]+coefs[9+5+3],
			coefs[9]+coefs[9+5+4],
			coefs[9]+coefs[9+5+5],
			coefs[9]+coefs[9+5+6],
			coefs[9]+coefs[9+5+7])

	slopesquared_SE<-c(
	deltamethod(~x9, coefs, covs, ses=TRUE), 
	deltamethod(~x9+x16, coefs, covs, ses=TRUE),
	deltamethod(~x9+x17, coefs, covs, ses=TRUE),
	deltamethod(~x9+x18, coefs, covs, ses=TRUE),
	deltamethod(~x9+x19, coefs, covs, ses=TRUE),
	deltamethod(~x9+x20, coefs, covs, ses=TRUE),
	deltamethod(~x9+x21, coefs, covs, ses=TRUE))

  	  	table.print<-data.frame('Basin'=B,
  	"$\\hat{\\beta_0}$"=paste(format(round(c(Intercepts),3),3), 
  	' (',format(round(unlist(Intercept_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_1}$"=paste(format(round(c(slopes),3),3), 
  	' (', format(round(unlist(slopes_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_2}$"=paste(format(round(c(slopesquared),3),3), 
  	' (',format(round(unlist(slopesquared_SE),3)),')',sep=''), 
	check.names = FALSE)
     
     print(xtable(table.print),sanitize.text.function = function(x){x})
		return(pars=c(Intercept_SE,slopes_SE,slopesquared_SE))
  	}
  	
  	if(model.num==3){
	  slopes<- c(coefs[2],
			coefs[2]+coefs[1+9+2],
			coefs[2]+coefs[1+9+3],
			coefs[2]+coefs[1+9+4],
			coefs[2]+coefs[1+9+5],
			coefs[2]+coefs[1+9+6],
			coefs[2]+coefs[1+9+7])

	slopes_SE<-c(
	deltamethod(~x2, coefs, covs, ses=TRUE), 
	deltamethod(~x2+x12, coefs, covs, ses=TRUE),
	deltamethod(~x2+x13, coefs, covs, ses=TRUE),
	deltamethod(~x2+x14, coefs, covs, ses=TRUE),
	deltamethod(~x2+x15, coefs, covs, ses=TRUE),
	deltamethod(~x2+x16, coefs, covs, ses=TRUE),
	deltamethod(~x2+x17, coefs, covs, ses=TRUE))

	slopesquared<- c(coefs[9],
			coefs[9]+coefs[9+7+2],
			coefs[9]+coefs[9+7+3],
			coefs[9]+coefs[9+7+4],
			coefs[9]+coefs[9+7+5],
			coefs[9]+coefs[9+7+6],
			coefs[9]+coefs[9+7+7])

	slopesquared_SE<-c(
	deltamethod(~x9, coefs, covs, ses=TRUE), 
	deltamethod(~x9+x18, coefs, covs, ses=TRUE),
	deltamethod(~x9+x19, coefs, covs, ses=TRUE),
	deltamethod(~x9+x20, coefs, covs, ses=TRUE),
	deltamethod(~x9+x21, coefs, covs, ses=TRUE),
	deltamethod(~x9+x22, coefs, covs, ses=TRUE),
	deltamethod(~x9+x23, coefs, covs, ses=TRUE))

	lat<- c(coefs[10],
			coefs[10]+coefs[10+12+2],
			coefs[10]+coefs[10+12+3],
			coefs[10]+coefs[10+12+4],
			coefs[10]+coefs[10+12+5],
			coefs[10]+coefs[10+12+6],
			coefs[10]+coefs[10+12+7])

	lat_SE<-c(
	deltamethod(~x10, coefs, covs, ses=TRUE), 
	deltamethod(~x10+x24, coefs, covs, ses=TRUE),
	deltamethod(~x10+x25, coefs, covs, ses=TRUE),
	deltamethod(~x10+x26, coefs, covs, ses=TRUE),
	deltamethod(~x10+x27, coefs, covs, ses=TRUE),
	deltamethod(~x10+x28, coefs, covs, ses=TRUE),
	deltamethod(~x10+x29, coefs, covs, ses=TRUE))

	ts<- c(coefs[11],
			coefs[11]+coefs[11+17+2],
			coefs[11]+coefs[11+17+3],
			coefs[11]+coefs[11+17+4],
			coefs[11]+coefs[11+17+5],
			coefs[11]+coefs[11+17+6],
			coefs[11]+coefs[11+17+7])

	ts_SE<-c(
	deltamethod(~x11, coefs, covs, ses=TRUE), 
	deltamethod(~x11+x30, coefs, covs, ses=TRUE),
	deltamethod(~x11+x31, coefs, covs, ses=TRUE),
	deltamethod(~x11+x32, coefs, covs, ses=TRUE),
	deltamethod(~x11+x33, coefs, covs, ses=TRUE),
	deltamethod(~x11+x34, coefs, covs, ses=TRUE),
	deltamethod(~x11+x35, coefs, covs, ses=TRUE))
  	
  	
  	table.print<-data.frame('Basin'=B,
  	"$\\hat{\\beta_0}$"=paste(format(round(c(Intercepts),3),3), 
  	' (',format(round(unlist(Intercept_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_1}$"=paste(format(round(c(slopes),3),3), 
  	' (', format(round(unlist(slopes_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_2}$"=paste(format(round(c(slopesquared),3),3), 
  	' (',format(round(unlist(slopesquared_SE),3)),')',sep=''), 
  	"$\\hat{\\beta_3}$"=paste(format(round(c(lat),3),3),  
  	' (',format(round(unlist(lat_SE),3)),')',sep=''),
  	"$\\hat{\\beta_4}$"=paste(format(round(c(ts),3),3), 
  	' (',format(round(unlist(ts_SE),3)),')',sep=''),
     check.names = FALSE)
     
     print(xtable(table.print),sanitize.text.function = function(x){x})
  			return(pars=c(Intercept_SE,slopes_SE,slopesquared_SE, lat_SE, ts_SE))

  	}
  }
  
}