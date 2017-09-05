####using historical data to estimate mu and sigma, then use GBM to simulate price.
library(tseries)
library(lmtest)
library(TSA)
library(fBasics)
library(fGarch)
#library(rugarch)
library(forecast)
library(stats)
library(PerformanceAnalytics)
library(quantmod)
library(car)
library(FinTS)
require(stochvol)
setwd("~/Desktop/Simulation/project")
price=read.table("price.txt",header = TRUE)
price_updated=read.table("price_updated.txt",header = TRUE)
returns=(price_updated[2:233,2]-price_updated[1:232,2])/price_updated[1:232,2]
plot(price)
plot(returns,type='l')
dim(price)
train=price[1:200,2]
train2=price[136:200,2]
train3=price[1:52,2]
xt=log(train[-1])-log(train[-200])

xt2=log(train2[-1])-log(train2[-65])
xt3=log(train3[-1])-log(train3[-52])
mean(xt)
mean(xt3)
plot(price)
points(price_updated$date[225:233],price_updated$close[225:233],col='red',pch='*')
plot(price_updated)
abline(v=2,col='green')
abline(v=52,col='red');abline(v=52,col='green');
xx1=log(price[2:53,2])-log(price[1:52,2])
xx2=log(price[54:109,2])-log(price[53:108,2])
xx3=log(price[110:162,2])-log(price[109:161,2])

ar(xx3)
acf(xx1**2);pacf(xx1**2)
acf(xx3**2);pacf(xx3**2)
plot(price_updated)
abline(v=52,col='green')

abline(v=81,col='green');
abline(v=109,col='green');

abline(v=141,col='red')
abline(v=162,col='red');abline(v=162,col='green')
abline(v=185,col='red')

abline(v=210,col='red')

abline(v=230,col='red')

abline(v=40,col='blue')
abline(v=104,col='blue')
abline(v=169,col='blue')
abline(v=169,col='blue')


acf(xt);acf(xt2);acf(xt3)
pacf(xt);pacf(xt2);pacf(xt3)
acf(xt**2);acf(xt2**2);acf(xt3**2)
pacf(xt**2);pacf(xt2**2);pacf(xt3**2)
Box.test(xt,lag=log(200),type='Ljung')
Box.test(xx3,lag=log(200),type='Ljung')

Box.test(xt2,lag=log(65),type='Ljung')
Box.test(xt3,lag=log(65),type='Ljung')
Box.test(xt3**2,lag=log(65),type='Ljung')

qqnorm(xt)
qqnorm(xt2)
qqnorm(xt3)
qqline(xt)
qqline(xt2)
qqline(xt3)

hist(xt,probability = TRUE)
hist(xt2,probability = TRUE)
densityPlot(as.timeSeries(xt))
densityPlot(as.timeSeries(xt2))


m=sum(xt/199)
m2=sum(xt2/65)
v=sum((xt-m)**2/199)
v2=sum((xt2-m2)**2/65)
sigma2=v*199
sigma22=v2*65
sigma=sqrt(sigma2)
sqrtsigma22=sqrt(sigma22)
r=0.5*sigma2+199*m
r2=0.5*sigma22+65*m2

test=matrix(NA,nrow=10000,ncol=24)
test[,1]=train2[65]*exp((r2-sigma22/2)/65+sqrtsigma22*sqrt(1/65)*rnorm(10000,0,1))

for (i in 2:24){
  test[,i]=mean(test[,i-1])*exp((r2-sigma22/2)/65+sqrtsigma22*sqrt(1/65)*rnorm(10000,0,1))
}




apply(test,2,mean)
plot(price[201:224,2],type='b',pch='o',ylim=c(35,50))
points(apply(test,2,mean),col='red',pch='*')
points(test[1,],col='blue',pch='*')
points(test[2,],col='pink',pch='*')
lines(test[3,],col='red',pch='*')
lines(test[4,],col='red',pch='*')
lines(test[5,],col='red',pch='*')
lines(test[6,],col='red',pch='*')
lines(test[7,],col='red',pch='*')
lines(test[8,],col='red',pch='*')
########################

####time series analysis
setwd("~/Desktop/Simulation/project")
data=read.table("price.txt",header = TRUE)

returns=log(price[2:201,2]/price[1:200,2])

plot(returns,type='l')
abline(h=0,col='red')
adf.test(returns)##no unit root
n=200
acf(returns,lag=100)#ma(27,54)
pacf(returns,lag=85)
eacf(returns)
Box.test(returns,lag=log(n),type='Ljung')
auto.arima(returns)
out=arima(returns,order=c(0,0,0))
coeftest(out)
acf(out$residuals)
pacf(out$residuals)
Box.test(out$residuals,lag=log(n)-4)

acf(out$residuals^2)
pacf(out$residuals^2)
Box.test(out$residuals^2,lag=log(n)-4)


###1/12/2016
xx<-setYuima(data=setData(price$close[1:200]))
ymodel <- setModel(drift=c("mu*x"), diffusion="sigma*x",solve.variable="x")
#grid=as.list(seq(1,200,1))
#ysamp <- setSampling(grid=grid)
yuima <- setYuima(data=setData(price$close[1:200]),model = ymodel)
#set.seed(123)
#yuima <- simulate(yuima, xinit = 1, true.parameter = list(mu = 0.2, sigma = 0.3))

library(yuima)
xx2<-setYuima(data=setData(price$close[150:175]))
ymodel <- setModel(drift=c("6.610696e-05*s"), diffusion="5.126536e-01*s",solve.variable="s")

X=simulate(ymodel)
plot(X)
r2=mean(xt)
sigma22=sd(xt)
prices=c()
prices[1]=price$close[1]



for (i in 1:200){
  prices[i+1]=prices[i]*exp((r2-sigma22/2)+sqrt(sigma22)*rnorm(100000,0,1))
}
test[,1]=price_updated$close[176]*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(100000,0,1)/sqrt(252))
plot(price,ylim=c(10,70))
lines(prices,type='l',col='red')
mod5 <- setModel(drift=c("theta"),diffusion=c("sigma"),solve.variable="x")
set.seed(123)

r2=6.610696e-05
sigma22=5.126536e-01
X <- simulate(mod5, true.p=list(theta=6.610696e-05,sigma=5.126536e-01),sampling=setSampling(n=200,Terminal=200,Initial=0)) 


plot(X)
points(price$close[1:200])


www3=function(test,int){
ymodel <- setModel(drift=c("mu*x"), diffusion="sigma*x",solve.variable="x")
yuima <- setYuima(data=setData(price_updated$close[(test-int)+1:test]),model = ymodel)
param.init <- list(mu=0.5,sigma=0.5)
low.par <-  list(mu=0, sigma=0)
upp.par <-  list(mu=1, sigma=1)
qmle=qmle(yuima, start = param.init,  lower = low.par, upper = upp.par)
sig2=summary(qmle)@coef[1,1]
mu=summary(qmle)@coef[2,1]

testo=array(NA,c(2000,10))

for (k in 1:2000){

    testo[k,1]=price_updated$close[test]*exp((mu-sig2**2/2)/252+sqrt(sig2**2/252)*rnorm(1,0,1))
  
}


for (k in 1:2000){

    for (s in 2:10){
      testo[k,s]=testo[k,s-1]*exp((mu-sig2**2/2)/252+sqrt(sig2**2/252)*rnorm(1,0,1))
    }
  }

#x=(apply(testo,2,f2))>=price_updated$close[(test+1):(test+10)]&(apply(testo,2,f1))<=price_updated$close[(test+1):(test+10)]
#sum(x/10)
rbind(apply(testo,2,mean),apply(testo,2,f1),apply(testo,2,f2))

}


www3(233,60)

www3(100,98)
r2=0.0001066836
sigma22=0.5005550231  
sd(log(price$close[150:175])-log(price$close[149:174]))/sqrt(1/252)









###using the whole dataset.
r2=6.610696e-05
sigma22= 5.048129e-01

r2=6.610696e-05
sigma22=5.126536e-01

###using the last 3 month
r2=6.610696e-05
sigma22=5.015993e-01

##using last 20 days
yuima <- setYuima(data=setData(price_updated$close[213:222]),model = ymodel)
r2= 6.610696e-05 
sigma22=5.126536e-01

r2=0
sigma22=0.5
testss=matrix(NA,nrow=490,ncol=10)
testss[,1]=price_updated$close[233]*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(490,0,1)/sqrt(252))

for (i in 2:10){
  testss[,i]=(testss[,i-1])*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(490,0,1)/sqrt(252))
}


for (i in 2:4){
  test[,i]=mean(test[,i-1])*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(100000,0,1)/sqrt(252))
}
test[,5]=mean(test[,4])*exp((r2-sigma22/2)*3/252+sqrt(sigma22)*sqrt(3)*rnorm(100000,0,1)/sqrt(252))
for (i in 6:9){
  test[,i]=mean(test[,i-1])*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(100000,0,1)/sqrt(252))
}
test[,10]=mean(test[,9])*exp((r2-sigma22/2)*3/252+sqrt(sigma22)*sqrt(3)*rnorm(100000,0,1)/sqrt(252))
for (i in 11:14){
  test[,i]=mean(test[,i-1])*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(100000,0,1)/sqrt(252))
}
test[,15]=mean(test[,14])*exp((r2-sigma22/2)*3/252+sqrt(sigma22)*sqrt(3)*rnorm(100000,0,1)/sqrt(252))
for (i in 16:19){
  test[,i]=mean(test[,i-1])*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(100000,0,1)/sqrt(252))
}
test[,20]=mean(test[,19])*exp((r2-sigma22/2)*3/252+sqrt(sigma22)*sqrt(3)*rnorm(100000,0,1)/sqrt(252))
for (i in 21:24){
  test[,i]=mean(test[,i-1])*exp((r2-sigma22/2)/252+sqrt(sigma22)*rnorm(100000,0,1)/sqrt(252))
}



apply(test,2,mean)
plot(price$close[150:175],ylim=c(40,60),xlim=c(0,40))
points(c(27:37),price_updated[176:186,2],type='b',pch='o',ylim=c(30,70))


f1=function(x){
   quantile(x,0.1)
}

f2=function(x){
  quantile(x,0.9)
}
points(c(1:10),apply(test,2,mean),col='red',pch='*')
points(c(1:10),apply(test,2,f1),col='red',pch='*')
points(c(1:10),apply(test,2,f2),col='red',pch='*')



points(c(1:24),apply(test,2,mean),col='red',pch='*')
points(c(1:10),test[9,],col='blue',pch='*')
points(c(1:24),test[2,],col='blue',pch='*')

points(c(1:24),apply(test,2,f1),col='red',pch='*')
points(c(1:24),apply(test,2,f2),col='red',pch='*')

points(test[1,],col='blue',pch='*')
points(test[2,],col='pink',pch='*')
lines(test[3,],col='red',pch='*')
lines(test[4,],col='red',pch='*')
lines(test[5,],col='red',pch='*')
lines(test[6,],col='red',pch='*')
lines(test[7,],col='red',pch='*')
lines(test[8,],col='red',pch='*')
x=c(1:200)
y=price$close[1:200]

library(KernSmooth)
fit <- locpoly(x, y, bandwidth = 1, gridsize=200, degree=0) 
lines(fit, col = "blue")
newdata=as.data.frame(c(201,202))
predict(fit,newdata=newdata)
fit


###12/2/2016
xt=log(price$close[-1])-log(price$close[-224])
plot(xt,type='l')
ar(xt)
acf(xt)
pacf(xt)

m0=arima(xt,order=c(7,0,0))
acf(m0$residuals)
pacf(m0$residuals)

n=length(xt)##223
p=7#7
d=2
m=n%/%10+p
X=array(0,c(n-p,(p+1)))#row 269,col 12
xd=xt[(p+1-d):(n-d)]#row 10:278
index=order(xd)#返回按顺序排列的xd的index????
xdo=xd[index]#按顺序排列的数据
for (i in 1:(n-p)){#269row
  ind=index[i]+p-d+d
  X[i,]=xt[ind:(ind-p)]
}


res=array(0,c(n-p,2))#269*2  1 column for  residual, 2 for standardized residual
AR=array(0,c(n-p,p+1,3))
for (i in m:(n-p-1)){
  fit=lm(X[1:i,1]~X[1:i,2:(p+1)])
  xnew=as.vector(X[i+1,2:(p+1)])
  xnew=c(1,xnew)
  pred=sum(fit$coef[1:(p+1)]*xnew)#beta*x=prediction
  fit.s=summary(fit)
  res[i+1,1]=X[i+1]-pred
  temp=1+t(xnew)%*%fit.s$cov%*%(xnew)##sigmat^2  REGRESSION LEC4 p13
  res[i+1,2]=res[i+1,1]/sqrt(temp)
  #res[i+1,2]=res[i+1,2]/fit.s$sigma
  AR[i+1,,1]=fit.s$coef[,1]#actual estimation
  AR[i+1,,2]=fit.s$coef[,2]#actual residual,standard error
  AR[i+1,,3]=fit.s$coef[,3]#actual t ratio
}

res.m=lm(res[(m+1):(n-p),2]~X[(m+1):(n-p),2:(p+1)])
res.s=summary(res.m)
summary(res.m)
F=(sum(res[(m+1):(n-p),2]^2)-res.s$sigma^2*(n-p-m-p-1))/(p+1)/res.s$sigma^2####important!
F#adjusted F for alpha0!!!!!!
1-pf(F,p+1,(n-p-m-p-1))

plot(xdo[(m+1):(n-p-m)],res[(m+1):(n-p-m),1],main="Ordinary Predictive Residuals")###xlab is x_pi
plot(xdo[(m+1):(n-p-m)],res[(m+1):(n-p-m),2],main="Standardized Predictive Residuals")###informative plot should be after regime change, the residuals would shift upward or downward.

plot(xdo[(m+1):(n-p-m)],AR[(m+1):(n-p-m),3,3])##t ratio for third column ,informative plots!!
abline(v=c(34.8,70.9))
abline(v=0.01)
#possible threshold!  #try combinations of threshold and fit the 3 models seperately,sum the AIC to see the overall fitness.
## 34.0  34.5  34.8  35.0  35.4  35.6 35.7  36.1 ## first threshold 34.8
## 70.0  70.9  73.0  73.0  73.0  74.0 ## second threshold 70.9


###12/3/2016
xt=log(price_updated$close[214:223])-log(price_updated$close[213:222])
sd(xt)/sqrt(10/252)
