#################################################################
#################################################################
#Parametric Bootstrap Function BVEVD Damages

para.boot.bvevd.damage<-function(boot.num=10,model.type=c(0,1,2)[1],model.object= dam.ws,data= tri_var[,c(5,4)], cov = nloc.basin, m.percent=1,rep=FALSE){
	
	sample.size.n<-dim(data)[1]
	data2<-matrix(NA,ncol=4,nrow= sample.size.n)
	data2[,1:2]<-as.matrix(data)

	Int_m1<-matrix(nrow= boot.num,ncol=(model.type+1))
	Int_m2<-matrix(nrow= boot.num,ncol=1)
	scale<-matrix(nrow= boot.num,ncol=2)
	shape<-matrix(nrow= boot.num,ncol=2)
	alpha<-vector(length=boot.num)
	beta<-vector(length=boot.num)
			
	loc_1<-model.object$estimate['loc1']
	if(model.type==1){
		loc_11<-model.object$estimate[2]
	}
	if(model.type==2){
		loc_11<-model.object$estimate[2]
		loc_12<-model.object$estimate[3]
	}
		
	scale_1<-model.object$estimate['scale1']
	shape_1<-model.object$estimate['shape1']
	
	loc_2<-model.object$estimate['loc2']

	scale_2<-model.object$estimate['scale2']
	shape_2<-model.object$estimate['shape2']

	alpha_1<-model.object$estimate['alpha']
	beta_2<-model.object$estimate['beta']
			
	sample.lengths<- floor(m.percent* sample.size.n)
	
	sample.size.m<-sum(sample.lengths)
	
	for(i in 1:boot.num){
    print(i)
		NA_samp<-sample(sample.size.n,sample.lengths,replace= rep)
		
		sampled<-c(NA_samp)
	
		if(model.type==1){
			loc1<-loc_1+loc_11*cov[sampled,1]
		}
		if(model.type==2){
			loc1<-loc_1+loc_11*cov[sampled,1]+loc_12*cov[sampled,2]
		}
	
		for(j in 1:length(NA_samp)){
			data2[NA_samp[j],3:4]<-
			rbvevd(1,
			alpha=alpha_1,beta= beta_2, 
			model='bilog',mar1=c(loc1[j],scale_1,shape_1), 
			mar2=c(loc_2,scale_2,shape_2))
		}
					
		if(model.type==0){			
			temp <-fbvevd(data2[sampled,3:4],model = "bilog",std.err=F)
		} else if (model.type==1){	
			temp <-fbvevd(data2[sampled,3:4],nsloc1=data.frame(cov[sampled,1]),model = "bilog",std.err=F)
		}  else{
			temp <-fbvevd(data2[sampled,3:4],nsloc1=data.frame(cov[sampled,1:2]),model = "bilog",std.err=F)
		}	
			
		Int_m1[i,1]<- c(temp$estimate['loc1'])
		if(model.type==1){
		Int_m1[i,2]<-temp $estimate[2]
			}
		if(model.type==2){
			Int_m1[i,2]<-temp$estimate[2]
			Int_m1[i,3]<-temp$estimate[3]
		}
	
		Int_m2[i,]<- c(temp$estimate['loc2'])
			
		scale[i,]<- c(temp$estimate['scale1'],temp$estimate['scale2'])
		shape[i,]<- c(temp$estimate['shape1'],temp$estimate['shape2'])
			
		alpha[i]<-temp$estimate['alpha']		
		beta[i]<-temp$estimate['beta']	

		}
			
	boot.est.basin.Int1<-colMeans(Int_m1)
	boot.se.basin.Int1<-apply(Int_m1, 2, sd)
	
	boot.est.basin.Int2<-colMeans(Int_m2)
	boot.se.basin.Int2<-apply(Int_m2, 2, sd)

	boot.est.scale<-colMeans(scale)
	names(boot.est.scale)<-c('scale1','scale2')
	boot.se.scale<-apply(scale, 2, sd)
	names(boot.se.scale)<-c('scale1','scale2')

	boot.est.shape<-colMeans(shape)		
	names(boot.est.shape)<-c('shape1','shape2')
	boot.se.shape<-apply(shape, 2, sd)
	names(boot.se.shape)<-c('shape1','shape2')

	boot.est.alpha<-mean(alpha)
	boot.se.alpha<-sd(alpha)

	boot.est.beta<-mean(beta)
	boot.se.beta<-sd(beta)

	boot.est.basins<-c(loc1=boot.est.basin.Int1,
		loc2= boot.est.basin.Int2,
		boot.est.scale,
		boot.est.shape,
		alpha=boot.est.alpha,
		beta=boot.est.beta)
		
	se.factor= sqrt(sample.size.n/sample.size.m)
		
	boot.se.basins<-c(loc1=boot.se.basin.Int1,
		loc2=boot.se.basin.Int2,
		boot.se.scale,
		boot.se.shape,
		alpha=boot.se.alpha,
		beta=boot.se.beta)* se.factor
		
	all.boots<-data.frame(Int_m1, Int_m2, scale, shape, alpha, beta)	
	if(model.type==0){
		names(all.boots)<-c('loc1','loc2','scale1','scale2','shape1','shape2','alpha','beta')
	}else if(model.type==1){
		names(all.boots)<-c('loc1','loc11','loc2','scale1','scale2','shape1','shape2','alpha','beta')
	} else{
		names(all.boots)<-c('loc1','loc11','loc12','loc2','scale1','scale2','shape1','shape2','alpha','beta')}


	return(list(
		boot.est.basins = boot.est.basins, boot.se.basins = boot.se.basins, se.factor= se.factor, all.boots= all.boots))
}

