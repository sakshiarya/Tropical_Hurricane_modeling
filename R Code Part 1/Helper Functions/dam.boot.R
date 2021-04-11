dam.boot <- function(file = 'a1', boot.num = 1000, model.type = c(0,1,2)[1], data= tri_var[,c(5,4)], m.percent=0.66){
  sample.size.n<-dim(data)[1]
  sample.lengths<- floor(m.percent* sample.size.n)
  sample.size.m<-sum(sample.lengths)
  Int_m1<-matrix(nrow= boot.num,ncol=(model.type+1))
  Int_m2<-matrix(nrow= boot.num,ncol=1)
  scale<-matrix(nrow= boot.num,ncol=2)
  shape<-matrix(nrow= boot.num,ncol=2)
  alpha<-rep(NA,boot.num)
  beta<- rep(NA,boot.num)
if(model.type == 0) {
  Int_m1 <- file$loc1
  Int_m2 <- file$loc2
  scale <- file[,c('scale1','scale2')]
  shape <- file[,c('shape1','shape2')]
  alpha <- file$alpha
  beta <- file$beta
} 
if(model.type == 1){
  Int_m1 <- file[,c('loc1','loc11')]
  Int_m2 <- file$loc2
  scale <- file[,c('scale1','scale2')]
  shape <- file[,c('shape1','shape2')]
  alpha <- file[,'alpha']
  beta <- file$beta
}  
if(model.type == 2){
  Int_m1 <- file[,c('loc1','loc11', 'loc12')]
  Int_m2 <- file$loc2
  scale <- file[,c('scale1','scale2')]
  shape <- file[,c('shape1','shape2')]
  alpha <- file[,'alpha']
  beta <- file$beta
}
if(model.type == 0){
  boot.est.basin.Int1<-mean(Int_m1)
  boot.se.basin.Int1<- sd(Int_m1)
  
  boot.est.basin.Int2<- mean(Int_m2)
  boot.se.basin.Int2<- sd(Int_m2)
}
  else{
boot.est.basin.Int1<-colMeans(Int_m1)
boot.se.basin.Int1<-apply(Int_m1, 2, sd)

boot.est.basin.Int2<-mean(Int_m2)
boot.se.basin.Int2<- sd(Int_m2)
}
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

#all.boots<-data.frame(Int_m1, Int_m2, scale, shape, alpha, beta)	
#if(model.type==0){
#  names(all.boots)<-c('loc1','loc2','scale1','scale2','shape1','shape2','alpha','beta')
#}else if(model.type==1){
#  names(all.boots)<-c('loc1','loc11','loc2','scale1','scale2','shape1','shape2','alpha','beta')
#} else{
#  names(all.boots)<-c('loc1','loc11','loc12','loc2','scale1','scale2','shape1','shape2','alpha','beta')}


return(list(
  boot.est.basins = boot.est.basins, boot.se.basins = boot.se.basins, se.factor= se.factor))
}
