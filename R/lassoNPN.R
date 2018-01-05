lassoNPN=function(Data_mat){
data_mat=as.matrix(Data_mat)
n=nrow(data_mat)
p=ncol(data_mat)
eta=0.25/n^(0.25)/sqrt(pi*log(n))
WE_mat=matrix(0,n,p)
for(i in 1:p){
mu=mean(data_mat[,i])
sig=sd(data_mat[,i])
temp=rank(data_mat[,i])/(n+1)
temp[temp<eta]=eta
temp[temp>(1-eta)]=1-eta
WE_mat[,i]=mu+sig*qnorm(temp)
}
beta_mat=matrix(0,p-1,p)
for(i in 1:p){
lambda_Liu=2*sqrt(var(WE_mat[,i])*log(p)/n)
flarerun=slim(WE_mat[,-i],WE_mat[,i],lambda=lambda_Liu,nlambda=1, lambda.min.value=NULL,lambda.min.ratio=NULL,rho=1,method="lasso",res.sd=FALSE,prec=1e-5,max.ite=1e5,verbose=FALSE)
temp1=flarerun$beta
zeronorm=sum(temp1!=0)
if(zeronorm>0){
D_i=WE_mat[,-i]
lr=lm(WE_mat[,i]~D_i[,which(temp1!=0)])
beta_mat[which(temp1!=0),i]=lr$coefficients[-1]
}
}
return(beta_mat)
}

