setwd("~/Desktop/565 financial time series/project")
data=read.csv("data.csv")




source("garchM.r")

##USE FIRST 2000 OBSERVATIONS TO FIT THE MODEL
WTI=data[1:2000,4]
exc=data[1:2000,3]
DJIA=data[1:2000,2]
BRE=data[1:2000,5]
n=length(WTI)
return_wti=log(WTI[2:n]/WTI[1:(n-1)])
return_bre=log(BRE[2:n]/BRE[1:(n-1)])
return_djia=log(DJIA[2:n]/DJIA[1:(n-1)])
plot(return_wti,type='l')
abline(h=0,col='red')
adf.test(return_wti)
par(mfrow=c(2,1))
acf(return_wti,lag=60)#ma(3,5,7,12,15,37,55)
pacf(return_wti,lag=60)#ar(3,5,7,12,37,55)
eacf(return_wti)#ma(1),ma(3),arma(1,1)
Box.test(return_wti,lag=log(n),type='Ljung')
auto.arima(return_wti)
out=arima(return_wti,order=c(2,0,2))
coeftest(out)
acf(out$residuals)
pacf(out$residuals)
Box.test(out$residuals,lag=log(n)-4)
polyroot(c(1,-out$coef[1],-out$coef[2]))
polyroot(c(1,out$coef[3],out$coef[4]))
acf(out$residuals^2)
pacf(out$residuals^2)
Box.test(out$residuals^2,lag=log(n)-4)

#fitting best mean component with ARMA(5,0,4) and fixed coefficients
meancom=arima(return_wti,order=c(5,0,4))
meancom=arima(return_wti,order=c(5,0,4),fixed=c(NA,0,NA,NA,NA,NA,0,NA,NA,0),transform.pars=FALSE)
acf(meancom$resi,lag=15)
pacf(meancom$resi,lag=15)
Box.test(meancom$res,lag=log(n)-4)
polyroot(c(1,-meancom$coef[1],-meancom$coef[2],-meancom$coef[3],-meancom$coef[4],-meancom$coef[5]))
polyroot(c(1,meancom$coef[6],meancom$coef[7],meancom$coef[8],meancom$coef[9]))

#######fit ARMA(2,2) GARCH(1,2) 
garchFit(~arma(2,2)+garch(1,2),data=return_wti,trace=F,include.mean = TRUE,cond.dist = "norm")

m=garchFit(~arma(2,2)+garch(1,2),data=return_wti,trace=F,include.mean = TRUE,cond.dist = "std")
m1=garchFit(~arma(2,2)+garch(1,2),data=return_wti,trace=F,include.mean = TRUE,cond.dist = "ged")

sresi=m@residuals/m@sigma.t
acf(sresi)
pacf(sresi)
acf(sresi^2)
pacf(sresi^2)
Box.test(sresi,lag=8,type="Ljung")
Box.test(sresi^2,lag=8,type="Ljung")

plot(m@sigma.t,type='l')
pacf(m@residuals**2)

qqnorm(sresi) 
qqline(sresi)
hist(sresi)
densityPlot(as.timeSeries(sresi))

xx=sqrt(8.76/6.76)*(-50:50)/10  #var(t)=df/df-2
lines(-50:50/10,dt(xx,8.76)*sqrt(8.76/6.76),col="blue")

#################Asymmetric garch
garch.spec = ugarchspec(variance.model = list(garchOrder=c(1,2)),mean.model = list(armaOrder=c(2,2))) 

wti.garch.fit = ugarchfit(spec=garch.spec, data=return_wti,solver.control=list(trace = 1))
wti.garch.fit


# GARCH(1,2)  with Student-t errors
garch.t.spec = ugarchspec(variance.model = list(garchOrder=c(1,2)), 
                            mean.model = list(armaOrder=c(2,2)),
                            distribution.model = "std")
wti.garch.t.fit = ugarchfit(spec=garch.t.spec, data=return_wti)                             
wti.garch.t.fit
plot(wti.garch.t.fit,which='all')
epsi1=wti.garch.t.fit@fit$residuals/wti.garch.t.fit@fit$sigma
acf(epsi1^2)
pacf(epsi1^2)
Box.test(epsi1,lag=log(n))
Box.test(epsi1^2,lag=4)
densityPlot(as.timeSeries(epsi1))

##GARCH(1,2) with skewed student-t errors
garch.sstd.spec = ugarchspec(variance.model = list(garchOrder=c(1,2)), 
                            mean.model = list(armaOrder=c(2,2)),
                            distribution.model = "sstd")
wti.garch.sstd.fit = ugarchfit(spec=garch.sstd.spec, data=return_wti)                             
wti.garch.sstd.fit

# Nelson's egarch model GARCH(1,2) with std
egarch.t.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,2)), mean.model=list(armaOrder=c(2,2)),distribution.model = "std")
wti.egarch.t.fit = ugarchfit(egarch.t.spec, return_wti)
wti.egarch.t.fit
epi=wti.egarch.t.fit@fit$residuals/wti.egarch.t.fit@fit$sigma
acf(epi^2)
# Nelson's egarch model GARCH(1,2) with sstd
egarch.sstd.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,2)), mean.model=list(armaOrder=c(2,2)),distribution.model = "sstd")
wti.egarch.sstd.fit = ugarchfit(egarch.sstd.spec, return_wti)
wti.egarch.sstd.fit
episstd=wti.egarch.sstd.fit@fit$residuals/wti.egarch.sstd.fit@fit$sigma
acf(episstd^2)


#GJR garch(1,2) std
gjrgarch.t.spec = ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,2)),mean.model=list(armaOrder=c(2,2)),distribution.model = "std")
wti.gjrgarch.t.fit = ugarchfit(gjrgarch.t.spec,return_wti)
wti.gjrgarch.t.fit
ugarchforecast(wti.gjrgarch.t.fit,n.ahead = 1,data=return_wti)
#GJR garch(1,2) sstd
gjrgarchsstd.spec = ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,2)),mean.model=list(armaOrder=c(2,2)),distribution.model = "sstd")
wti.gjrgarchsstd.fit = ugarchfit(gjrgarchsstd.spec,return_wti)
wti.gjrgarchsstd.fit
ugarchforecast(wti.gjrgarchsstd.fit,n.ahead = 1,data=return_wti)

# aparch models garch(1,2) std
aparch.t.spec = ugarchspec(variance.model = list(model="apARCH",
                                                     garchOrder=c(1,2)), 
                               mean.model = list(armaOrder=c(2,2)),
                               distribution.model = "std",
                               fixed.pars=list(delta=1))
wti.aparch.t.fit = ugarchfit(spec=aparch.t.spec, data=return_wti)                             
wti.aparch.t.fit
# aparch models garch(1,2) sstd
aparch.sstd.spec = ugarchspec(variance.model = list(model="apARCH",
                                                     garchOrder=c(1,2)), 
                               mean.model = list(armaOrder=c(2,2)),
                               distribution.model = "sstd",
                               fixed.pars=list(delta=1))
wti.aparch.sstd.fit = ugarchfit(spec=aparch.sstd.spec, data=return_wti)                             
wti.aparch.sstd.fit






nic.garch = newsimpact(wti.garch.t.fit)
nic.egarch = newsimpact(wti.egarch.t.fit)
nic.gjrgarch = newsimpact(wti.gjrgarch.t.fit)
nic.aparch = newsimpact(wti.aparch.t.fit)

# compare information criteria
model.list = list(garch = wti.garch.t.fit,
                  egarch = wti.egarch.t.fit,
                  gjrgarch = wti.gjrgarch.t.fit,
                  aparch = wti.aparch.t.fit)
info.mat = sapply(model.list, infocriteria)
rownames(info.mat) = rownames(infocriteria(wti.garch.t.fit))
info.mat

# show news impact curve from estimated standard garch(1,2) and asymmetric garch(1,2)
par(mfrow=c(2,2))
plot(nic.garch$zx, type="l", lwd=2, col="blue", main="GARCH(1,2)", 
     nic.garch$zy, ylab=nic.garch$yexpr, xlab=nic.garch$xexpr)
plot(nic.egarch$zx, type="l", lwd=2, col="blue", main="EGARCH(1,2)", 
     nic.egarch$zy, ylab=nic.egarch$yexpr, xlab=nic.egarch$xexpr)
plot(nic.gjrgarch$zx, type="l", lwd=2, col="blue", main="GJRGARCH(1,2)", 
     nic.gjrgarch$zy, ylab=nic.gjrgarch$yexpr, xlab=nic.gjrgarch$xexpr)
plot(nic.aparch$zx, type="l", lwd=2, col="blue", main="APARCH(1,1,2)", 
     nic.aparch$zy, ylab=nic.aparch$yexpr, xlab=nic.aparch$xexpr)
par(mfrow=c(1,1))


# examine standardized residuals
wti.garch.t.zt = residuals(wti.garch.t.fit)/sigma(wti.garch.t.fit)
qqPlot(coredata(wti.garch.t.zt))
plot(wti.garch.t.fit, which="all")
plot(wti.egarch.t.fit, which="all")
plot(wti.gjrgarch.t.fit, which="all")
plot(wti.aparch.t.fit, which="all")

# compute 100 1-step ahead rolling forecasts
wti=log(data[2:2492,4]/data[1:2491,4])#alldata

wti.garch.t.fit = ugarchfit(garch.t.spec, data=wti[1:2099],out.sample=100) 
wti.garch.t.fcst = ugarchforecast(wti.garch.t.fit, n.roll=100, n.ahead=1)

wti.egarch.t.fit = ugarchfit(egarch.t.spec, data=wti[1:2099],out.sample=100)
wti.egarch.t.fcst = ugarchforecast(wti.egarch.t.fit, n.roll=100, n.ahead=1)

wti.gjrgarch.t.fit = ugarchfit(gjrgarch.t.spec, data=wti[1:2099],out.sample=100)
wti.gjrgarch.t.fcst = ugarchforecast(wti.gjrgarch.t.fit, n.roll=100, n.ahead=1)

wti.aparch.t.fit = ugarchfit(aparch.t.spec, data=wti[1:2099],out.sample=100)
wti.aparch.t.fcst = ugarchforecast(wti.aparch.t.fit, n.roll=100, n.ahead=1)

# compute forecast evaluation statistics
fcst.list = list(garch=wti.garch.t.fcst,
                 egarch=wti.egarch.t.fcst,
                 gjrgarch=wti.gjrgarch.t.fcst,
                 aparch=wti.aparch.t.fcst)
fpm.mat = sapply(fcst.list, fpm)
fpm.mat

########garchM not significant though
source("garchM.r")
garchM(return_wti)
garchM(return_djia)


######stochastic volatility
require(stochvol)
muremove=return_wti-mean(return_wti)
sv1=svsample(muremove)
pre=predict(sv1,1)
htpre=apply(pre,2,median)
v1pre=exp(htpre/2)
apply(sv1$para,2,mean)
sqrt(apply(sv1$para,2,var))
ht=apply(sv1$latent,2,median)
v11=exp(ht/2) ## volatility
ts.plot(v11)####volatility



######rolling forecast   gev block=21
wti=log(data[2:2492,4]/data[1:2491,4])

VaR5E=(1999:2490)*0
for (i in 1999:2490) {
  xx=wti[1:i]
  fit=gev(-xx,block=21)
  par=fit$par.est*c(-1,1,-1)
  VaR5E[i-1998]=par[3]-par[2]/par[1]*(1-(-21*log(1-0.05))**par[1])
}


########rolling forecast   gev blocl=63
wti=log(data[2:2492,4]/data[1:2491,4])

VaR5E2=(1999:2490)*0
for (i in 1999:2490) {
  xx=wti[1:i]
  fit2=gev(-xx,block=63)#length(sp5)=792,n0=21 would be reasonable
  par2=fit2$par.est*c(-1,1,-1)##negative sign!!!!!!
  VaR5E2[i-1998]=par2[3]-par2[2]/par2[1]*(1-(-21*log(1-0.05))**par2[1])
}



##########ARMA(2,2) GARCH(1,2) NORM
var1=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  garch.spec = ugarchspec(variance.model = list(garchOrder=c(1,2)),mean.model = list(armaOrder=c(2,2))) 
  
  wti.garch.t.fit = ugarchfit(spec=garch.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                              
                              fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                              
                              numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                      
                                                      grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                      
                                                      hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  pred1=ugarchforecast(wti.garch.t.fit,n.ahead=1,data=xx)
  var1[i-1998]=as.numeric(pred1@forecast$seriesFor+qnorm(0.05)*pred1@forecast$sigmaFor)
} 



#######ARMA(2,2) GARCH(1,2) STD
var2=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  garch.t.spec = ugarchspec(variance.model = list(garchOrder=c(1,2)), 
                              mean.model = list(armaOrder=c(2,2)),
                              distribution.model = "std")
  wti.garch.t.fit = ugarchfit(spec=garch.t.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                                
                                fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                                
                                numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                        
                                                        grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                        
                                                        hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  pred1=ugarchforecast(wti.garch.t.fit,n.ahead=1,data=xx)
  var2[i-1998]=as.numeric(pred1@forecast$seriesFor+qt(0.05,wti.garch.t.fit@fit$coef[10])*pred1@forecast$sigmaFor)/sqrt(wti.garch.t.fit@fit$coef[10]/(wti.garch.t.fit@fit$coef[10]-2))
} 


#################GARCH(1,2) SSTD
varsstdstan=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  garch.sstd.spec = ugarchspec(variance.model = list(garchOrder=c(1,2)), 
                            mean.model = list(armaOrder=c(2,2)),
                            distribution.model = "sstd")
  wti.garch.sstd.fit = ugarchfit(spec=garch.sstd.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                              
                              fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                              
                              numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                      
                                                      grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                      
                                                      hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  pred1=ugarchforecast(wti.garch.sstd.fit,n.ahead=1,data=xx)
  epsilon=wti.garch.sstd.fit@fit$residuals/wti.garch.sstd.fit@fit$sigma
  mu=mean(epsilon)
  sigma=sqrt(var(epsilon))
  quan=qsstd(0.05, mean = mu, sd = sigma, nu =  wti.garch.sstd.fit@fit$coef[11], xi =  wti.garch.sstd.fit@fit$coef[10])
  quan1=(quan-mu)/var(epsilon)
  varsstdstan[i-1998]=as.numeric(pred1@forecast$seriesFor+quan1*pred1@forecast$sigmaFor)
} 


###########EGARCH STD
predg=ugarchforecast(wti.egarch.t.fit,n.ahead=1,data=wti[1:1999])
varegarch=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  egarch.t.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,2)),
                             mean.model = list(armaOrder=c(2,2)),
                             distribution.model = "std")
  wti.egarch.t.fit= ugarchfit(spec=egarch.t.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                              
                              fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                              
                              numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                      
                                                      grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                      
                                                      hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  
  predg=ugarchforecast( wti.egarch.t.fit,n.ahead=1,data=xx)
  varegarch[i-1998]=as.numeric(predg@forecast$seriesFor+qt(0.05,wti.egarch.t.fit@fit$coef[11])*predg@forecast$sigmaFor)/sqrt(wti.egarch.t.fit@fit$coef[11]/(wti.egarch.t.fit@fit$coef[11]-2))
} 

###########EGARCH SSTD
varegarchsstd=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  egarch.sstd.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,2)),
                             mean.model = list(armaOrder=c(2,2)),
                             distribution.model = "sstd")
  wti.egarch.sstd.fit= ugarchfit(spec=egarch.sstd.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                              
                              fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                              
                              numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                      
                                                      grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                      
                                                      hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  
  predg=ugarchforecast(wti.egarch.sstd.fit,n.ahead=1,data=xx)
  epsilon=wti.egarch.sstd.fit@fit$residuals/wti.egarch.sstd.fit@fit$sigma
  mu=mean(epsilon)
  sigma=sqrt(var(epsilon))
  quan=qsstd(0.05, mean = mu, sd = sigma, nu =  wti.egarch.sstd.fit@fit$coef[12], xi =  wti.egarch.sstd.fit@fit$coef[11])
  quan1=(quan-mu)/var(epsilon)
  varegarchsstd[i-1998]=as.numeric(predg@forecast$seriesFor+quan1*predg@forecast$sigmaFor)
} 

#######GJR STD
vargjr=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  gjrgarch.t.spec = ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,2)),
                               mean.model = list(armaOrder=c(2,2)),
                               distribution.model = "std")
  wti.gjrgarch.t.fit= ugarchfit(spec=  gjrgarch.t.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                                
                                fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                                
                                numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                        
                                                        grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                        
                                                        hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  
  predgjr=ugarchforecast( wti.gjrgarch.t.fit,n.ahead=1,data=xx)
  vargjr[i-1998]=as.numeric(predgjr@forecast$seriesFor+qt(0.05,  wti.gjrgarch.t.fit@fit$coef[11])*predgjr@forecast$sigmaFor)/sqrt(wti.gjrgarch.t.fit@fit$coef[11]/(wti.gjrgarch.t.fit@fit$coef[11]-2))
} 

###############GJR SSTD
vargjrgarchsstd=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  gjrgarch.sstd.spec = ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,2)),
                               mean.model = list(armaOrder=c(2,2)),
                               distribution.model = "sstd")
  wti.gjrgarch.sstd.fit= ugarchfit(spec=  gjrgarch.sstd.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                                
                                fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                                
                                numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                        
                                                        grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                        
                                                        hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  
  predgjr=ugarchforecast( wti.gjrgarch.sstd.fit,n.ahead=1,data=xx)
  epsilon=wti.gjrgarch.sstd.fit@fit$residuals/wti.gjrgarch.sstd.fit@fit$sigma
  mu=mean(epsilon)
  sigma=sqrt(var(epsilon))
  quan=qsstd(0.05, mean = mu, sd = sigma, nu =wti.gjrgarch.sstd.fit@fit$coef[12], xi =  wti.gjrgarch.sstd.fit@fit$coef[11])
  quan1=(quan-mu)/var(epsilon)
  vargjrgarchsstd[i-1998]=as.numeric(predgjr@forecast$seriesFor+quan1*predgjr@forecast$sigmaFor)
} 
write.csv(vargjrgarchsstd,"vargjrgarchsstd.csv")

########APARCH STD
varaparch=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  aparch.t.spec = ugarchspec(variance.model=list(model="apARCH",garchOrder=c(1,2)),
                                 mean.model = list(armaOrder=c(2,2)),
                                 distribution.model = "std")
  wti.aparch.t.fit= ugarchfit(spec=  aparch.t.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                                  
                                  fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                                  
                                  numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                          
                                                          grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                          
                                                          hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  
  predaparch=ugarchforecast(wti.aparch.t.fit,n.ahead=1,data=xx)
  varaparch[i-1998]=as.numeric(predaparch@forecast$seriesFor+qt(0.05,wti.aparch.t.fit@fit$coef[12])*predaparch@forecast$sigmaFor)/sqrt(wti.aparch.t.fit@fit$coef[12]/(wti.aparch.t.fit@fit$coef[12]-2))
} 
##########
########APARCH SSTD
varaparchsstd=(1999:2491)*0
for (i in 1999:2491) {
  xx=wti[1:i]
  aparch.sstd.spec = ugarchspec(variance.model=list(model="apARCH",garchOrder=c(1,2)),
                             mean.model = list(armaOrder=c(2,2)),
                             distribution.model = "sstd")
  wti.aparch.sstd.fit= ugarchfit(spec=  aparch.sstd.spec, data=xx,out.sample = 0, solver = "hybrid", solver.control = list(), 
                              
                              fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                              
                              numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                      
                                                      grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                      
                                                      hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  
  
  predaparch=ugarchforecast(wti.aparch.sstd.fit,n.ahead=1,data=xx)
  epsilon=wti.aparch.sstd.fit@fit$residuals/wti.aparch.sstd.fit@fit$sigma
  mu=mean(epsilon)
  sigma=sqrt(var(epsilon))
  quan=qsstd(0.05, mean = mu, sd = sigma, nu =  wti.aparch.sstd.fit@fit$coef[13], xi =  wti.aparch.sstd.fit@fit$coef[12])
  quan1=(quan-mu)/var(epsilon)
  varaparchsstd[i-1998]=as.numeric(predaparch@forecast$seriesFor+quan1*predaparch@forecast$sigmaFor)
} 
write.csv(varaparchsstd,"varaparchsstd.csv")


###########empirical quantile
VaR5Q=(1999:2491)*0
for (i in 1999:2491){
  xx=wti[1:i]
  VaR5Q[i-1998]=quantile(xx,0.05)
}


###########empirical quantile2 moving window
VaR5Q22=(1999:2491)*0
for (i in 1999:2491){
  xx=wti[i-1998:i]
  VaR5Q22[i-1998]=quantile(xx,0.05)
}



########stochastic volatility

varsv=(1999:2490)*0
for (i in 1999:2490){
  xx=wti[1:i]
  out=arima(xx,order=c(2,0,2))
  muremove=out$residuals
  sv=svsample(muremove)
  pre=predict(sv,1)
  htpre=apply(pre,2,median)
  vlpre=exp(htpre/2)
  arimamean=as.numeric(predict(out,1)$pred)
  varsv[i-1998]=arimamean-1.645*vlpre
}

##### add djia as xreg garch(1,2) std
djia=log(data[2:2492,2]/data[1:2491,2])
varxreg=(1999:2491)*0
for (t in 1998:2491) {
  xx=wti[2:(t+1)]
  yy=djia[1:t]
  out=arima(xx,xreg=yy,order=c(2,0,2))
  garch.t.spec = ugarchspec(variance.model = list(garchOrder=c(1,2)), 
                              mean.model = list(armaOrder=c(2,2)),
                              distribution.model = "std")
  wti.garch.t.fit = ugarchfit(spec=garch.t.spec, data=out$res,out.sample = 0, solver = "hybrid", solver.control = list(), 
                                
                                fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
                                
                                numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, 
                                                        
                                                        grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=.1, 
                                                        
                                                        hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
  pred=predict(out,newxreg = djia[t+1])
  pred$pred
  pred1=ugarchforecast(wti.garch.t.fit,n.ahead=1,data=out$residuals)
  varxreg[t-1997]=as.numeric(pred$pred[1]+qt(0.05,wti.garch.t.fit@fit$coef[10])*pred1@forecast$sigmaFor)/sqrt(wti.garch.t.fit@fit$coef[10]/(wti.garch.t.fit@fit$coef[10]-2))
} 


#####pot 
par(mfrow=c(2,1))
qplot(-return_wti,threshold=0.035)###nonlinear xi unequal to 0
meplot(-return_wti)##to determine threshold
abline(v=0.02)

m1=pot(-return_wti,0.02)
m2=gpd(-return_wti,0.02)

#####pot with threshold 0.02
pot0.02=(1999:2491)*0
for (i in 1999:2491){
  xx=wti[1:i]
  m=pot(-xx,0.02)
  pot0.02[i-1998]=-riskmeasures(m,0.95)[2]
}

#####pot with threshold 0.03
pot0.03=(1999:2491)*0
for (i in 1999:2491){
  xx=wti[1:i]
  m=gpd(-xx,0.03)
  pot0.03[i-1998]=-riskmeasures(m,0.95)[2]
}

#####pot with threshold 0.037
pot0.037=(1999:2491)*0
for (i in 1999:2491){
  xx=wti[1:i]
  m=gpd(-xx,0.037)
  pot0.037[i-1998]=-riskmeasures(m,0.95)[2]
}







plot(wti[2000:2491],type='l',xlab='ahead step',ylab='return',main='one-day 5% rolling forecast VaR')
lines(VaR5E[1:492,2],col='purple')#GEV21
lines(VaR5E2[1:492,2],col='pink')#GEV63
lines(var1[1:492,2],col='green')#GARCH(1,2) NORMAL
lines(var2[1:492,2],col='brown')#GARCH(1,2) STD
lines(vargarchsstd[1:492,2],col=654)#GARCH(1,2) SSTD
lines(varegarch[1:492,2],col='yellow')#EGARCH(1,2) STD
lines(varegarchsstd[1:492,2],col=619)#EGARCH(1,2) SSTD
lines(vargjr[1:492,2],col='522')#GJRGARCH(1,2) STD
lines(vargjrgarchsstd[1:492,2],col=501)#GJRGARCH(1,2) SSTD
lines(varaparch[1:492,3],col=454)#APARCH(1,2) STD
lines(varaparchsstd[1:492,2],col=651)#APARCH(1,2) SSTD
lines(VaR5Q[1:492,2],col='red')#empirical quantile
lines(VaR5Q22[1:492,2],col='blue') #empirical quantile with rolling window
lines(varsv[1:492,2],col=650)#stochastic volatility
lines(varxreg[1:492,2],col=620)#xreg GARCH(1,2) STD
lines(pot0.02[1:492,2],col=610)#POT with threshold=0.02
lines(pot0.03[1:492,2],col=620)#POT with threshold=0.03
lines(pot0.037[1:492,2],col=620)#POT with threshold=0.037











1-(sum(VaR5E[1:492,2]<wti[2000:2491]))/492 #GEV21 0.1382114
1-(sum(VaR5E2[1:492,2]<wti[2000:2491]))/492 #GEV63  0.07520325
1-(sum(var1[1:492,2]<wti[2000:2491]))/492 #GARCH(1,2) NORMAL 0.06504065
1-(sum(var2[1:492,2]<wti[2000:2491]))/492#GARCH(1,2) STD 0.06300813
1-(sum(vargarchsstd[1:492,2]<wti[2000:2491]))/492#GARCH(1,2) SSTD 0.05894309
1-(sum(varegarch[1:492,2]<wti[2000:2491]))/492#EGARCH(1,2) STD 0.05691057
1-(sum(varegarchsstd[1:492,2]<wti[2000:2491]))/492#EGARCH(1,2) SSTD 0.04878049
1-(sum(vargjr[1:492,2]<wti[2000:2491]))/492#GJRGARCH(1,2) STD 0.05894309
1-(sum(vargjrgarchsstd[1:492,2]<wti[2000:2491]))/492#GJRGARCH(1,2) SSTD 0.05894309
1-(sum(varaparch[1:492,3]<wti[2000:2491]))/492#APARCH(1,2) STD 0.05785124
sum(is.na(varaparch[1:492,3]))
varaparch3=varaparch[1:492,3]
varaparch3[is.na(varaparch3)]=10
1-(sum(varaparch3[1:492]<wti[2000:2491]))/484
1-(sum(varaparchsstd[1:492,2]<wti[2000:2491]))/492#APARCH(1,2) SSTD 0.05590062
sum(is.na(varaparchsstd[1:492,2]))
varaparchsstd3=varaparchsstd[1:492,2]
varaparchsstd3[is.na(varaparchsstd3)]=10
1-(sum(varaparchsstd3[1:492]<wti[2000:2491]))/483


1-(sum(VaR5Q[1:492,2]<wti[2000:2491]))/492#empirical quantile 0.09146341
1-(sum(VaR5Q22[1:492,2]<wti[2000:2491]))/492#empirical quantile with rolling window 0.1239837
1-(sum(varsv<wti[2000:2491]))/492#stochastic volatility  0.06910569
1-(sum(varxreg[1:492,2]<wti[2000:2491]))/492##DJIA as xreg GARCH(1,2) STD 0.06097561
1-(sum(pot0.02[1:492,2]<wti[2000:2491]))/492###POT with threshold=0.02 0.09552846
1-(sum(pot0.03[1:492,2]<wti[2000:2491]))/492###POT with threshold=0.03 0.09349593
1-(sum(pot0.037[1:492,2]<wti[2000:2491]))/492###POT with threshold=0.037 0.09146341
