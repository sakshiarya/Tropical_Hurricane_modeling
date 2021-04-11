#################################################################
#################################################################
##  Individual iteration function for Bootstrap
bvevd_iteri <- function(boot.i, boot.num=100,model.object=bvevd3,data= data.for.bivariate, cov = nloc.basin, m.percent=1,rep=FALSE){
  sample.size.n<-dim(data)[1]
  
  data2<- data.frame(cbind(data, rep(NA, sample.size.n), 
                           rep(NA, sample.size.n)), ncol = 5)
  # /data2 <- as.matrix(data2)	
  coeffs<-rep(NA, length(model.object$estimate)) 
  Int_m1<- rep(NA,ncol=7)
  Int_m2<-rep(NA,ncol=7)
  scale<-rep(NA,ncol=2)
  shape<-rep(NA,ncol=2)
  alpha<-NA
  beta<-NA
  
  loc_1<-c(model.object$estimate[1], 
           sum(model.object$estimate[1:2]),
           sum(model.object$estimate[c(1,3)]),
           sum(model.object$estimate[c(1,4)]), 
           sum(model.object$estimate[c(1,5)]),
           sum(model.object$estimate[c(1,6)]), 
           sum(model.object$estimate[c(1,7)]))
  
  scale_1<-model.object$estimate['scale1']
  shape_1<-model.object$estimate['shape1']
  
  loc_2<-c(model.object$estimate[10], 
           sum(model.object$estimate[10:11]),
           sum(model.object$estimate[c(10,12)]),
           sum(model.object$estimate[c(10,13)]), 
           sum(model.object$estimate[c(10,14)]),
           sum(model.object$estimate[c(10,15)]), 
           sum(model.object$estimate[c(10,16)]))
  
  scale_2<-model.object$estimate['scale2']
  shape_2<-model.object$estimate['shape2']
  alpha_1<-model.object$estimate['alpha']
  beta_2<-model.object$estimate['beta']
  
  
  EP_set<-data[,3]=="EaPa"; EP_len<-sum(EP_set)
  NA_set<-data[,3]=="NoAt"; NA_len<-sum(NA_set)
  NI_set<-data[,3]=="NoIn"; NI_len<-sum(NI_set)
  SEI_set<-data[,3]=="SoEIn"; SEI_len<-sum(SEI_set)
  SP_set<-data[,3]=="SoPa"; SP_len<-sum(SP_set)
  SWI_set<-data[,3]=="SoWIn"; SWI_len<-sum(SWI_set)
  WP_set<-data[,3]=="WePa"; WP_len<-sum(WP_set)
  
  sample.lengths<- floor(m.percent*
                           c(EP_len, NA_len, NI_len, SEI_len, SP_len, SWI_len, WP_len))
  
  sample.size.m<-sum(sample.lengths)

  
EP_samp<-sample(EP_len,sample.lengths[1],replace= rep)
NA_samp<-sample(NA_len,sample.lengths[2],replace= rep)
NI_samp<-sample(NI_len,sample.lengths[3],replace= rep)
SEI_samp<-sample(SEI_len,sample.lengths[4],replace= rep)
SP_samp<-sample(SP_len,sample.lengths[5],replace= rep)
SWI_samp<-sample(SWI_len,sample.lengths[6],replace= rep)
WP_samp<-sample(WP_len,sample.lengths[7],replace= rep)

sampled<-c(EP_samp,
           EP_len+NA_samp,
           EP_len+NA_len+NI_samp,
           EP_len+NA_len+NI_len+SEI_samp,
           EP_len+NA_len+NI_len+SEI_len+SP_samp,
           EP_len+NA_len+NI_len+SEI_len+SP_len+SWI_samp,
           EP_len+NA_len+NI_len+SEI_len+SP_len+SWI_len+ WP_samp)

data2[EP_samp,4:5]<-evd::rbvevd(sample.lengths[1],
                                alpha=alpha_1, beta= beta_2, 
                                model='bilog', mar1=c(loc_1[1],scale_1,shape_1), 
                                mar2=c(loc_2[1],scale_2,shape_2))

data2[EP_len+NA_samp,4:5]<-rbvevd(sample.lengths[2], 
                                  alpha=alpha_1, beta= beta_2, 
                                  model='bilog', mar1=c(loc_1[2],scale_1,shape_1), 
                                  mar2=c(loc_2[2],scale_2,shape_2))

data2[EP_len+NA_len+NI_samp,4:5]<-rbvevd(sample.lengths[3],
                                         alpha=alpha_1, beta= beta_2, 
                                         model='bilog', mar1=c(loc_1[3],scale_1,shape_1), 
                                         mar2=c(loc_2[3],scale_2,shape_2))

data2[EP_len+NA_len+NI_len+SEI_samp,4:5]<-rbvevd(sample.lengths[4],
                                                 alpha=alpha_1, beta= beta_2, 
                                                 model='bilog', mar1=c(loc_1[4],scale_1,shape_1), 
                                                 mar2=c(loc_2[4],scale_2,shape_2))

data2[EP_len+NA_len+NI_len+SEI_len+SP_samp,4:5]<-
  rbvevd(sample.lengths[5],
         alpha=alpha_1, beta= beta_2,
         model='bilog', mar1=c(loc_1[5],scale_1,shape_1), 
         mar2=c(loc_2[5],scale_2,shape_2))

data2[EP_len+NA_len+NI_len+SEI_len+SP_len+SWI_samp,4:5]<-
  rbvevd(sample.lengths[6],
         alpha=alpha_1, beta= beta_2, 
         model='bilog', mar1=c(loc_1[6],scale_1,shape_1), 
         mar2=c(loc_2[6],scale_2,shape_2))

data2[EP_len+NA_len+NI_len+SEI_len+SP_len+SWI_len+WP_samp,4:5]<-
  rbvevd(sample.lengths[7],
         alpha=alpha_1, beta= beta_2,
         model='bilog', mar1=c(loc_1[7],scale_1,shape_1), 
         mar2=c(loc_2[7],scale_2,shape_2))

temp <-fbvevd(data2[sampled,4:5],
              nsloc1= cov[c(sampled),-1],
              nsloc2= cov[c(sampled),-1], 
              model = "bilog",std.err=F)

coeffs <-temp$estimate

Int_m1<- c(coeffs[1],
           coeffs[1]+coeffs[1+1],
           coeffs[1]+coeffs[1+2],
           coeffs[1]+coeffs[1+3],
           coeffs[1]+coeffs[1+4],
           coeffs[1]+coeffs[1+5],
           coeffs[1]+coeffs[1+6])

Int_m2<- c(coeffs[10],
           coeffs[10]+coeffs[10+1],
           coeffs[10]+coeffs[10+2],
           coeffs[10]+coeffs[10+3],
           coeffs[10]+coeffs[10+4],
           coeffs[10]+coeffs[10+5],
           coeffs[10]+coeffs[10+6])

scale<- c(temp$estimate['scale1'],temp$estimate['scale2'])
shape<- c(temp$estimate['shape1'],temp$estimate['shape2'])

alpha<-temp$estimate['alpha']		
beta<-temp$estimate['beta']	
return(c(Int_m1, Int_m2, scale, shape, alpha, beta, coeffs))
}

#Parametric Bootstrap Function

para.boot.bvevd.2021<-function(boot.num=100,model.object=bvevd3,data= data.for.bivariate, cov = nloc.basin, m.percent=1,rep=FALSE){
  sample.size.n<-dim(data)[1]
  EP_set<-data[,3]=="EaPa"; EP_len<-sum(EP_set)
  NA_set<-data[,3]=="NoAt"; NA_len<-sum(NA_set)
  NI_set<-data[,3]=="NoIn"; NI_len<-sum(NI_set)
  SEI_set<-data[,3]=="SoEIn"; SEI_len<-sum(SEI_set)
  SP_set<-data[,3]=="SoPa"; SP_len<-sum(SP_set)
  SWI_set<-data[,3]=="SoWIn"; SWI_len<-sum(SWI_set)
  WP_set<-data[,3]=="WePa"; WP_len<-sum(WP_set)
  
  sample.lengths<- floor(m.percent*
                           c(EP_len, NA_len, NI_len, SEI_len, SP_len, SWI_len, WP_len))
  
  sample.size.m<-sum(sample.lengths)
  results <- list()
  results <- mclapply(1:boot.num, bvevd_iteri,model.object=model.object,data= data, cov = cov, m.percent=m.percent,rep=rep, mc.cores = 10)
  results <- matrix(unlist(results), ncol = 2*length(model.object$estimate), byrow = TRUE)
  Int_m1 <- results[,1:7]
  Int_m2 <- results[,8:14]
  scale <- results[,15:16]
  shape <- results[,17:18]
  alpha <- results[,19]
  beta <- results[,20]
  coeffs <- results[,21:dim(results)[2]]
  boot.est.basin.Int1<-colMeans(Int_m1)
  boot.se.basin.Int1<-apply(Int_m1, 2, sd)
  names(boot.est.basin.Int1)<-c('EaPa Int.','NoAt Int.' ,'NoIn Int.' ,'SoEIn Int.' ,'SoPa Int.' ,'SoWIn Int.','WePa Int.')
  names(boot.se.basin.Int1)<-c('EaPa Int.','NoAt Int.' ,'NoIn Int.' ,'SoEIn Int.' ,'SoPa Int.' ,'SoWIn Int.','WePa Int.')
  
  boot.est.basin.Int2<-colMeans(Int_m2)
  boot.se.basin.Int2<-apply(Int_m2, 2, sd)
  names(boot.est.basin.Int2)<-c('EaPa Int.','NoAt Int.' ,'NoIn Int.' ,'SoEIn Int.' ,'SoPa Int.' ,'SoWIn Int.','WePa Int.')
  names(boot.se.basin.Int2)<-c('EaPa Int.','NoAt Int.' ,'NoIn Int.' ,'SoEIn Int.' ,'SoPa Int.' ,'SoWIn Int.','WePa Int.')
  
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
  
  boot.est.basins<-c(boot.est.basin.Int1,
                     boot.est.basin.Int2,
                     boot.est.scale,
                     boot.est.shape,
                     alpha=boot.est.alpha,
                     beta=boot.est.beta)
  
  boot.se.basins<-c(boot.se.basin.Int1,
                    boot.se.basin.Int2,
                    boot.se.scale,
                    boot.se.shape,
                    alpha= boot.se.alpha,
                    beta= boot.se.beta)*sqrt(sample.size.n/sample.size.m)
  
  boot.est.orig<-colMeans(coeffs)
  names(boot.est.orig)<-names(model.object$estimate)
  boot.se.orig<-apply(coeffs, 2, sd)*sqrt(sample.size.n/sample.size.m)
  names(boot.se.orig)<-names(model.object$estimate)		
  return(list(
    boot.est.basins = boot.est.basins, boot.se.basins = boot.se.basins,
    boot.est.orig = boot.est.orig , boot.se.orig = boot.se.orig ,all.coeff.est= coeffs))
}

