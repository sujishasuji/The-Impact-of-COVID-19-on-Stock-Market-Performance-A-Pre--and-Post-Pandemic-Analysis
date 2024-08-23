#CODE
#load libraries
library(PerformanceAnalytics)
library(quantmod)
library(rugarch)
library(car)
library(FinTS)
library(moments)
library(fGarch)
library(tseries)
library(aTSA)
library(MTS)

A <- read_excel("~/project sujisha/new 15.xlsx")

x=A[,2:16];x
summary(x)
kurtosis(x)
skewness(x)

####Time plot####

for (i in 2:ncol(A)) {
  x1[i]=A[,i]
  xt[i]=ts(x1[i],frequency=52,start=c(2016,1))
  x2=unlist(xt[i])
  p=plot.ts(x2,main=colnames(A[,i]))
  
}

##### Log Return #####
LR=function(x){
  lr={}
  n=length(x)
  for(i in 2:(n-1))
    lr[i-1]=log(x[i])-log(x[i-1])
  return(lr)
}

#### ARCH MODEL ####
arch1=garchFit(~garch(1,0),logreturn);arch1
s1=summary(arch1);s1 
re1=residuals(arch1,standardize=T);re1
lb1=Box.test(re1,lag=round(2*sqrt(n)),type="Ljung-Box");lb1
sw1=shapiro.test(re1);sw1
jb1=jarque.bera.test(re1);jb1
info1=arch1@fit$ics[1];info1
garch11=garchFit(~garch(1,1),logreturn);garch11
s11=summary(garch11);s11
re11=residuals(garch11,standardize=T);re11
lb11=Box.test(re11,lag=round(2*sqrt(n)),type="Ljung-Box");lb11
sw11=shapiro.test(re11);sw11
jb11=jarque.bera.test(re11);jb11
info11=garch11@fit$ics[1];info11
Model=c("arch1","garch11");Model
AIC=c(info1,info11);AIC
LB=c(lb1$p.value,lb11$p.value);LB
SW=c(sw1$p.value,sw11$p.value);SW
JB=c(jb1$p.value,jb11$p.value);JB
Tabel=data.frame(Model,AIC,LB,SW,JB);Tabel
cat("\n The different combinations are given by\n")
print(Tabel)
plot(residuals(garch11),type="l",main="Plot of Garch(1,1) model")
acf(residuals(garch11),main="acf of Garch(1,1) model")
qqnorm(residuals(garch11),main="qqplot of Garch(1,1) model")
qqline(residuals(garch11))

#### Prediction ####
Pred=function(v,f){
  p1=rep(0,52)
  ex=exp(f)
  p1[1]=v*ex[1]
  for(i in 2:52){
    p1[i]=p1[i-1]*ex[i]
  }
  return(p1)
}

### Forecasting ####
P={}
for(i in 2:16){
  x1=A[,i]
  x1=unlist(as.vector(x1))
  v=x1[435]
  x1=LR(x1)
  g= ugarchspec(variance.model = list(garchOrder=c(1,1)),mean.model = list(armaOrder=c(0,1)))
  m =ugarchfit(spec=g, data=x1,out.sample=52)
  m1=ugarchforecast(m, n.ahead = 52, n.roll = 0)
  f=as.numeric(m1@forecast$seriesFor)
  p=Pred(v,f)
  P=cbind(P,p)
}
colnames(P)=colnames(A)[2:16]


### Covid impact###
D1=A$D1
D2=A$D2
D3=A$D3
Dummy=cbind(D2,D3)
g <- ugarchspec(
  variance.model = list(garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), external.regressors = Dummy))
m <- ugarchfit(spec = g, data = A$KARNATAKA, out.sample = 52)
s=sigma(m)