DNetGGM=function(Data_mat1,Data_mat2,Beta_mat1,Beta_mat2,alpha){
data1=as.matrix(Data_mat1)
data2=as.matrix(Data_mat2)
beta1=as.matrix(Beta_mat1)
beta2=as.matrix(Beta_mat2)

n1=nrow(data1)
n2=nrow(data2)
p=ncol(data1)
alp=alpha

eps1=matrix(0,n1,p)
eps2=matrix(0,n2,p)

cmean1=colMeans(data1)
cmean2=colMeans(data2)

for(j in 1:p){
temp1=data1[,j]-cmean1[j]-(data1[,-j]-matrix(cmean1[-j],n1,p-1,byrow=T))%*%beta1[,j]
temp2=data2[,j]-cmean2[j]-(data2[,-j]-matrix(cmean2[-j],n2,p-1,byrow=T))%*%beta2[,j]
eps1[,j]=temp1
eps2[,j]=temp2
}

T1=matrix(0,p,p)
T2=matrix(0,p,p)
Tstar=matrix(0,p,p)

for(i in 1:(p-1)){
for(j in (i+1):p){
temp1=1/n1*(sum(eps1[,i]*eps1[,j])+sum(eps1[,i]*eps1[,i])*beta1[i,j]+sum(eps1[,j]*eps1[,j])*beta1[j-1,i])
temp2=1/n2*(sum(eps2[,i]*eps2[,j])+sum(eps2[,i]*eps2[,i])*beta2[i,j]+sum(eps2[,j]*eps2[,j])*beta2[j-1,i])
T1[i,j]=1/sqrt(sum(eps1[,i]*eps1[,i])/n1*sum(eps1[,j]*eps1[,j])/n1)*temp1
T2[i,j]=1/sqrt(sum(eps2[,i]*eps2[,i])/n2*sum(eps2[,j]*eps2[,j])/n2)*temp2
}
}

T1=T1+t(T1)
T2=T2+t(T2)

for(i in 1:(p-1)){
for(j in (i+1):p){
rhohat1=ifelse(abs(T1[i,j])>(2*sqrt(log(p)/n1)), T1[i,j], 0)
rhohat2=ifelse(abs(T2[i,j])>(2*sqrt(log(p)/n2)), T2[i,j], 0)
Tstar[i,j]=abs(T1[i,j]-T2[i,j])/sqrt((1-rhohat1^2)^2/n1+(1-rhohat2^2)^2/n2)
}
}

Tstar1=Tstar+t(Tstar)

P0=2*pnorm(1,0,1)-1
P0hat=(sum(Tstar1<1)-p)/(p^2-p)
Q0=sqrt(2)*dnorm(1,0,1)
A=(P0-P0hat)/Q0
t=seq(0,max(Tstar1),by=max(Tstar1)/99)

chooset=numeric(length(t))
for(i in 1:length(t)){
chooset[i]=alp*max(1,sum(Tstar1>t[i])/2)/(p^2-p)/(1-pnorm(t[i])+abs(A*t[i]/sqrt(2))*dnorm(t[i]))
}

t1hat=0
if(sum(chooset>1)>0){
t1hat=t[min(which(chooset>1))]
}
if(sum(chooset>1)==0){
t1hat=2*sqrt(log(p))
}

DNet=ifelse(Tstar1>t1hat,1,0)

return(DNet)
}






