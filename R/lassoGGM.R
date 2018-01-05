lassoGGM=function(Data_mat){	
data_mat=as.matrix(Data_mat)
n=nrow(data_mat)
p=ncol(data_mat)
beta_mat=matrix(0,p-1,p)
for(i in 1:p){
lambda_Liu=2*sqrt(var(data_mat[,i])*log(p)/n)
D_i=data_mat[,-i]
flarerun=slim(D_i,data_mat[,i],rho=1,lambda=lambda_Liu,nlambda=1, lambda.min.value=NULL,lambda.min.ratio=NULL,method="lasso",res.sd=FALSE,prec=1e-5,max.ite=1e5,verbose=FALSE)
temp=flarerun$beta
zeronorm=sum(temp!=0)
if(zeronorm>0){
lr=lm(data_mat[,i]~D_i[,which(temp!=0)])
beta_mat[which(temp!=0),i]=lr$coefficients[-1]
}
}
return(beta_mat)
}


