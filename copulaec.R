

install.packages("PerformanceAnalytics")
install.packages("quantmod")
install.packages("car")
install.packages("FinTS")
install.packages("mvtnorm")
install.packages("MASS") # needed for cov.trob
install.packages("mnormt") # needed for dmt
install.packages("sn")
install.packages("copula") 
# load libraries
library(PerformanceAnalytics)
library(quantmod)
library(car)
library(FinTS)
library(mvtnorm)
library(MASS) # needed for cov.trob
library(mnormt) # needed for dmt
library(sn)
library(copula)
library(scatterplot3d)
options(digits=4)


# extract adjusted closing prices
getSymbols('MSFT',src='google',from ="2006-01-03", to = "2016-04-03")
getSymbols('AAPL',src='google',from ="2006-01-03", to = "2016-04-03")
MSFT = MSFT[, "MSFT.Close", drop=F]
AAPL = AAPL[, "AAPL.Close", drop=F]

# calculate log-returns for GARCH analysis
MSFT.ret = CalculateReturns(MSFT, method="log")
AAPL.ret = CalculateReturns(AAPL, method="log")


# remove first NA observation
MSFT.ret = MSFT.ret[-1,]
AAPL.ret = AAPL.ret[-1,]
colnames(MSFT.ret) ="MSFT"
colnames(AAPL.ret) = "AAPL"

# create combined data series
MSFT.AAPL.ret = cbind(MSFT.ret,AAPL.ret)

#
# Bivariate distributions for MSFT and GPSC
#

# plot returns
my.panel <- function(...) {
  lines(...)
  abline(h=0)
}

plot.zoo(MSFT.AAPL.ret, main="Daily Returns", 
         panel=my.panel, col=c("black", "blue"))

# Empirical scatterplots
plot(coredata(MSFT.ret), coredata(AAPL.ret),
     main="Empirical Bivariate Distribution of Returns",
     ylab="AAPL", xlab="MSFT", col="blue")
abline(h=mean(AAPL.ret), v=mean(MSFT.ret))

# normal qq-plots
par(mfrow=c(1,2))
  qqPlot(coredata(MSFT.ret), main="MSFT", ylab="MSFT quantiles")
  qqPlot(coredata(AAPL.ret), main="AAPL", ylab="AAPL quantiles")
par(mfrow=c(1,1))

# simulate multivariate normal returns calibrated to observed data

# estimate mean and covariance
n.obs = nrow(MSFT.AAPL.ret)
mu.hat = apply(MSFT.AAPL.ret, 2, mean)
 = cov(MSFT.AAPL.ret)
Cor.hat = cov2cor()

set.seed(123)
sim.ret = rmvnorm(n.obs, mean=mu.hat, sigma=, method="chol")

# scaterplot of simulated and actual returns

plot(coredata(MSFT.ret), coredata(AAPL.ret),
     main="Empirical vs. Bivariate Normal",
     ylab="AAPL", xlab="MSFT", col="blue")
abline(h=mean(AAPL.ret), v=mean(MSFT.ret))
points(sim.ret[,1], sim.ret[,2], col="red")
legend(x="topleft", legend=c("Empirical", "Normal"), col=c("blue", "red"), pch=1)


# fit multivariate t
# fitting bivariate student's t to data - see SDAFE ch 5

df = seq(2.1,5,.01) # 551 points
n = length(df)
loglik_max = rep(0,n)
for(i in 1:n)
{
  fit = cov.trob(coredata(MSFT.AAPL.ret),nu=df[i])
  loglik_max[i] = sum(log(dmt(coredata(MSFT.AAPL.ret),mean=fit$center,
                              S=fit$cov,df=df[i])))
}

max.lik = max(loglik_max)
v.mle = df[which(loglik_max == max.lik)]
plot(df, loglik_max, type="l", main="Логарифм функции правдоподобия", lwd=2, col="blue")
abline(v=v.mle, lwd=2, col="red")

# extract mle of mu and sigma given v.mle
fit.mle = cov.trob(coredata(MSFT.AAPL.ret),nu=v.mle)
mu.mle.t = fit$center
Sigma.mle.t = fit$cov
Cor.mle.t = cov2cor(Sigma.mle.t)
mu.mle.t
Sigma.mle.t*(v.mle/(v.mle - 2))
Cor.mle.t

# simulate portfolio returns
# use result that Y = mu + sqrt(v/W)*Z is multivariate t

# generate Z ~ N(0, Sigma.mle.t)
set.seed(123)
Z = rmvnorm(n=n.obs, mean=c(0,0), sigma=Sigma.mle.t)
# generate W ~ chi-sq(v.mle)
W = rchisq(n.obs,df=v.mle)
# simulate bivariate t
sim.ret.t = mu.mle.t + sqrt(v.mle/W)*Z
colnames(sim.ret.t) = c("MSFT","AAPL")

# plot simulated data together with actual returns

plot(coredata(MSFT.ret), coredata(AAPL.ret),
     main="Empirical vs. Bivariate t (df=2.61)",
     ylab="AAPL", xlab="MSFT", col="blue")
abline(h=mean(AAPL.ret), v=mean(MSFT.ret))
points(sim.ret.t, col="red")
legend(x="topleft", legend=c("Empirical", "Multivariate t"), col=c("blue", "red"), pch=1)

# fit univariate skew-t distributions

st.msft.fit = st.mle(y=coredata(MSFT.ret))
st.msft.fit$dp
qqPlot(coredata(MSFT.ret), dist="st", location=st.msft.fit$dp[1], scale=st.msft.fit$dp[2], 
       shape=st.msft.fit$dp[3], df=st.msft.fit$dp[4], envelope=FALSE)

st.sp500.fit = st.mle(y=coredata(AAPL.ret))
st.sp500.fit$dp
qqPlot(sp500, dist="st", location=st.sp500.fit$dp[1], scale=st.sp500.fit$dp[2], 
       shape=st.sp500.fit$dp[3], df=st.sp500.fit$dp[4], envelope=FALSE)

# simulate bivariate distribution with zero correlation
st.msft.sim = rst(n.obs, location=st.msft.fit$dp[1], scale=st.msft.fit$dp[2], 
                  shape=st.msft.fit$dp[3], df=st.msft.fit$dp[4])
st.sp500.sim = rst(n.obs, location=st.sp500.fit$dp[1], scale=st.sp500.fit$dp[2], 
                   shape=st.sp500.fit$dp[3], df=st.sp500.fit$dp[4])

plot(st.msft.sim, st.sp500.sim,
     main="Independent skew-t distributions",
     ylab="AAPL", xlab="MSFT", col="blue")
abline(h=mean(AAPL.ret), v=mean(MSFT.ret)) 

########################################################################
# properties of CDFs
########################################################################

x = rnorm(100, mu.hat[1], sqrt()[1])
hist(x)
# should be U[0,1]
u = pnorm(x, mu.hat[1], sqrt()[1])
hist(u)
# should be N(mu, sigma)
x.new = qnorm(u, mu.hat[1], sqrt()[1])
hist(x.new)

########################################################################
# Special Copulas
########################################################################

# independent copula
# copula package indepCoupla() function - S4 
indep.cop = indepCopula(2)
class(indep.cop)
slotNames(indep.cop)
# simulate data from copula
set.seed(123)
u = rCopula(indep.cop, 200)
head(u)

# plots of density and CDF
par(mfrow=c(2,2))
  persp(indep.cop, pcopula, main="CDF", 
        xlab="u", ylab="v", zlab="C(u,v)")
  contour(indep.cop, pcopula, main="CDF", 
          xlab="u", ylab="v")    
  plot(u, main="Simulations", 
       xlab="u", ylab="v")
par(mfrow=c(1,1))


########################################################################
# Dependence measures
########################################################################

# pearson's linear correlation
cor(MSFT.AAPL.ret, method="pearson")[1,2]
# Kendall's tau
cor(MSFT.AAPL.ret, method="kendall")[1,2]
# Spearman's rho
cor(MSFT.AAPL.ret, method="spearman")[1,2]

########################################################################
# Elliptical copulas
########################################################################

# bivariate normal copula with positive dependence
norm.cop.9 = normalCopula(param=0.9, dim=2)
class(norm.cop.9)
slotNames(norm.cop.9)

# plot copula CDF, pdf and contours
par(mfrow=c(2,2))
persp(norm.cop.9, pcopula, main="CDF", 
      xlab="u", ylab="v", zlab="C(u,v)")
persp(norm.cop.9, dcopula, main="pdf", 
      xlab="u", ylab="v", zlab="c(u,v)")
contour(norm.cop.9, pcopula, main="CDF", 
        xlab="u", ylab="v")    
contour(norm.cop.9, dcopula, main="pdf", 
        xlab="u", ylab="v")
par(mfrow=c(1,1))

# compute Kendall's tau, spreaman's rho and tailindex parameters
tau(norm.cop.9)
rho(norm.cop.9)
tailIndex(norm.cop.9)

# bivariate normal copula with negative dependence
norm.cop.m9 = normalCopula(param=-0.9, dim=2)

# plot copula CDF, pdf and contours
par(mfrow=c(2,2))
persp(norm.cop.m9, pcopula, main="CDF", 
      xlab="u", ylab="v", zlab="C(u,v)")
persp(norm.cop.m9, dcopula, main="pdf", 
      xlab="u", ylab="v", zlab="c(u,v)")
contour(norm.cop.m9, pcopula, main="CDF", 
        xlab="u", ylab="v")    
contour(norm.cop.m9, dcopula, main="pdf", 
        xlab="u", ylab="v")
par(mfrow=c(1,1))

tau(norm.cop.m9)
rho(norm.cop.m9)
tailIndex(norm.cop.m9)

# bivariate normal copula with no dependence
norm.cop.0 = normalCopula(param=0, dim=2)

# plot copula CDF, pdf and contours
par(mfrow=c(2,2))
persp(norm.cop.0, pcopula, main="CDF", 
      xlab="u", ylab="v", zlab="C(u,v)")
persp(norm.cop.0, dcopula, main="pdf", 
      xlab="u", ylab="v", zlab="c(u,v)")
contour(norm.cop.0, pcopula, main="CDF", 
        xlab="u", ylab="v")    
contour(norm.cop.0, dcopula, main="pdf", 
        xlab="u", ylab="v")
par(mfrow=c(1,1))

set.seed(123)
u.9 = rCopula(norm.cop.9, 200)
u.m9 = rCopula(norm.cop.m9, 200)
u.0 = rCopula(norm.cop.0, 200)

par(mfrow=c(2,2))
  plot(u.9, main="rho=0.9", 
       xlab="u", ylab="v", col="blue", pch=16)
  plot(u.m9, main="rho=-0.9", 
       xlab="u", ylab="v", col="red", pch=16)
  plot(u.0, main="rho=0", 
       xlab="u", ylab="v", col="black", pch=16)
par(mfrow=c(1,1))

# bivariate student t copula
# bivariate t copulas 
t.cop.9 = tCopula(param=0.9, dim=2, df=4)
t.cop.m9 = tCopula(param=-0.9, dim=2, df=4)
t.cop.0 = tCopula(param=0, dim=2, df=4)

t.cop.5 = tCopula(param=0.5, dim=2, df=4)

# plot copula CDF, pdf and contours
par(mfrow=c(2,2))
persp(t.cop.9, pcopula, main="CDF", 
      xlab="u", ylab="v", zlab="C(u,v)")
persp(t.cop.9, dcopula, main="pdf", 
      xlab="u", ylab="v", zlab="c(u,v)")
contour(t.cop.9, pcopula, main="CDF", 
        xlab="u", ylab="v")    
contour(t.cop.9, dcopula, main="pdf", 
        xlab="u", ylab="v")
par(mfrow=c(1,1))

tau(t.cop.9)
rho(t.cop.9)
tailIndex(t.cop.9)


set.seed(123)
u.9 = rCopula(t.cop.9, 200)
u.m9 = rCopula(t.cop.m9, 200)
u.0 = rCopula(t.cop.0, 200)

par(mfrow=c(2,2))
plot(u.9, main="rho=0.9", 
     xlab="u", ylab="v", col="blue", pch=16)
plot(u.m9, main="rho=-0.9", 
     xlab="u", ylab="v", col="red", pch=16)
plot(u.0, main="rho=0", 
     xlab="u", ylab="v", col="black", pch=16)
par(mfrow=c(1,1))



set.seed(123)
u.2 = rCopula(tCopula(param=0.9, dim=2, df=2), 500)
u.5 = rCopula(tCopula(param=0.9, dim=2, df=5), 500)
u.10 = rCopula(tCopula(param=0.9, dim=2, df=10), 500)
u.30 = rCopula(tCopula(param=0.9, dim=2, df=30), 500)

tailIndex(tCopula(param=0.9, dim=2, df=2))
tailIndex(tCopula(param=0.9, dim=2, df=5))
tailIndex(tCopula(param=0.9, dim=2, df=10))
tailIndex(tCopula(param=0.9, dim=2, df=30))


par(mfrow=c(2,2))
plot(u.2, main="rho=0.9, df=2", 
     xlab="u", ylab="v", col="blue", pch=16)
plot(u.5, main="rho=0.9, df=5", 
     xlab="u", ylab="v", col="red", pch=16)
plot(u.10, main="rho=0.9, df=10", 
     xlab="u", ylab="v", col="black", pch=16)
plot(u.30, main="rho=0.9, df=30", 
     xlab="u", ylab="v", col="green", pch=16)
par(mfrow=c(1,1))

########################################################################
# Archmedian copulas
########################################################################

# Gumbel copula delta = 1
gum.cop.1 = archmCopula(family="gumbel", dim=2, param=1.01)
class(gum.cop.1)
gum.cop.4 = archmCopula(family="gumbel", dim=2, param=4)
gum.cop.10 = archmCopula(family="gumbel", dim=2, param=10)

tau(gum.cop.4)
rho(gum.cop.4)
tailIndex(gum.cop.4)

# plot copula CDF, pdf and contours
par(mfrow=c(2,2))
persp(gum.cop.4, pcopula, main="CDF",
      xlab="u", ylab="v", zlab="C(u,v)")
persp(gum.cop.4, dcopula, main="pdf", 
      xlab="u", ylab="v", zlab="c(u,v)")
contour(gum.cop.4, pcopula, main="CDF", 
        xlab="u", ylab="v")    
contour(gum.cop.4, dcopula, main="pdf", 
        xlab="u", ylab="v")
par(mfrow=c(1,1))

set.seed(123)
u.1 = rCopula(gum.cop.1, 200)
u.4 = rCopula(gum.cop.4, 200)
u.10 = rCopula(gum.cop.10, 200)

# compute tailIndex parameters
tailIndex(gum.cop.1)
tailIndex(gum.cop.4)
tailIndex(gum.cop.10)

par(mfrow=c(2,2))
plot(u.1, main="delta=1", 
     xlab="u", ylab="v", col="blue", pch=16)
plot(u.4, main="delta=4", 
     xlab="u", ylab="v", col="red", pch=16)
plot(u.10, main="delta=10", 
     xlab="u", ylab="v", col="black", pch=16)
par(mfrow=c(1,1))

# Clayton copula delta = 1
clay.cop.0 = archmCopula(family="clayton", dim=2, param=0.01)
class(gum.cop.1)
clay.cop.4 = archmCopula(family="clayton", dim=2, param=4)
clay.cop.10 = archmCopula(family="clayton", dim=2, param=10)

tau(clay.cop.4)
rho(clay.cop.4)
tailIndex(clay.cop.4)

# plot copula CDF, pdf and contours
par(mfrow=c(2,2))
persp(clay.cop.4, pcopula, main="CDF",
      xlab="u", ylab="v", zlab="C(u,v)")
persp(clay.cop.4, dcopula, main="pdf", 
      xlab="u", ylab="v", zlab="c(u,v)")
contour(clay.cop.4, pcopula, main="CDF", 
        xlab="u", ylab="v")    
contour(clay.cop.4, dcopula, main="pdf", 
        xlab="u", ylab="v")
par(mfrow=c(1,1))

set.seed(123)
u.0 = rCopula(clay.cop.0, 200)
u.4 = rCopula(clay.cop.4, 200)
u.10 = rCopula(clay.cop.10, 200)

par(mfrow=c(2,2))
plot(u.0, main="delta=0", 
     xlab="u", ylab="v", col="blue", pch=16)
plot(u.4, main="delta=4", 
     xlab="u", ylab="v", col="red", pch=16)
plot(u.10, main="delta=10", 
     xlab="u", ylab="v", col="black", pch=16)
par(mfrow=c(1,1))

tailIndex(clay.cop.0)
tailIndex(clay.cop.4)
tailIndex(clay.cop.10)

######################################################################
# creating arbitrary distributions
#####################################################################

# use copula function mvdc
args(mvdc)

# bivariate distribution with N(3, 4^2) and t3 margins, and gumbel
# copula with d = 2
my.cop = archmCopula(family="gumbel", dim=2, param=2)
my.margins = c("norm", "t")
my.parms = list(list(mean=3, sd=4), list(df=3))
myBvd = mvdc(copula=my.cop, 
             margins=my.margins, 
             paramMargins=my.parms)
class(myBvd)
slotNames(myBvd)

# plot custom distribution
par(mfrow=c(2,2))
  persp(myBvd, pMvdc, main="CDF", xlim=c(-5,12), ylim=c(-3,3),
        xlab="x", ylab="y", zlab="F(x,y)")
  persp(myBvd, dMvdc, main="pdf",xlim=c(-5,12), ylim=c(-3,3), 
        xlab="x", ylab="y", zlab="c(x,y)")
  contour(myBvd, pMvdc, main="CDF",xlim=c(-5,12), ylim=c(-3,3), 
          xlab="x", ylab="y")    
  contour(myBvd, dMvdc, main="pdf", xlim=c(-5,12), ylim=c(-3,3),
          xlab="x", ylab="y")
par(mfrow=c(1,1))

# simulate observations from distribution
set.seed(123)
myBvd.sim = rMvdc(500, myBvd)
# plot simulations
par(mfrow=c(2,2))
  qqPlot(myBvd.sim[,1], main="X Margin", ylab="X quantiles")
  qqPlot(myBvd.sim[,2], main="Y Margin", ylab="Y quantiles")
  plot(myBvd.sim, main="Simulated X,Y values", 
       xlab="x", ylab="y", pch=16, col="blue")
  abline(h=0, v=3)
par(mfrow=c(1,1))

#
# estimating bivariate distributions defined by margins and copulas
#

# estimate bivariate distn with X ~ N(3,4^2), Y~t(3)
# and gumbel copula with delta = 2 by full MLE
args(fitMvdc)
start.vals = c(1, 2, 4, 5)
names(start.vals) = c("mu", "sigma", "df", "delta")
myBvd.fitMvdc = fitMvdc(myBvd.sim, myBvd, start.vals)
class(myBvd.fitMvdc)
slotNames(myBvd.fitMvdc)
# print results
myBvd.fitMvdc

# simulate values from fitted distribution
param.hat = myBvd.fitMvdc@estimate
param.hat
my.cop.fit = archmCopula(family="gumbel", dim=2, param=param.hat[1])
my.margins.fit = c("norm", "t")
my.parms.fit = list(list(mean=param.hat[2], sd=param.hat[3]), 
                    list(df=param.hat[4]))
myBvd.fit = mvdc(copula=my.cop.fit, 
                 margins=my.margins.fit, 
                 paramMargins=my.parms.fit)
set.seed(123)
myBvd.fit.sim = rMvdc(500, myBvd.fit)

par(mfrow=c(2,2))
  qqPlot(myBvd.fit.sim[,1], main="X Margin", ylab="X quantiles")
  qqPlot(myBvd.fit.sim[,2], main="Y Margin", ylab="Y quantiles")
  plot(myBvd.fit.sim, main="Simulated X,Y values", 
       xlab="x", ylab="y", pch=16, col="blue")
  abline(h=0, v=3)
par(mfrow=c(1,1))


# estimate using two-step IFM procedure
x = myBvd.sim[,1]
y = myBvd.sim[,2]
# step 1: estimate marginal distributions
# X ~ N(mu, sigma^2)
mu.hat = mean(x)
 = sd(x)
mu.hat

# Y ~ t(df=3)
fit.t = fitdistr(y, densfun="t")
df.hat = coef(fit.t)["df"]
df.hat

# transform data to uniform using estimated CDF function
u.hat = pnorm(x, mu.hat, )
v.hat = pt(y, df=df.hat)
# show transformed uniform data
plot(u.hat,v.hat, main="Transformed Data", 
     xlab="u.hat", ylab="v.hat", pch=16, col="blue")

# step 2: estimate copula on transformed uniform observations
args(fitCopula)
fit.ifm = fitCopula(copula=myBvd@copula, 
                    data=cbind(u.hat,v.hat), 
                    method="mpl", start=2)
class(fit.ifm)
slotNames(fit.ifm)
fit.ifm

ifm.parms = c(mu.hat, , df.hat, fit.ifm@estimate)
names(ifm.parms) = c("mu", "sigma", "df", "delta")

# compare MLE and IFM fits
myBvd.fitMvdc@estimate
ifm.parms

# evaluating goodness of fit
gof.test = gofCopula(myBvd@copula, cbind(u.hat,v.hat),
                     estim.method="mpl")
gof.test


############################################################
# Estimate bivariate distn for msft and aapl
############################################################

# fit univariate skew-t distributions
st.msft.fit = st.mle(y=coredata(MSFT.ret))
st.aapl.fit = st.mle(y=coredata(AAPL.ret))
# transform data to uniform
u.hat = pst(coredata(MSFT.ret), dp=st.msft.fit$dp)
v.hat = pst(coredata(AAPL.ret), dp=st.aapl.fit$dp)

plot(u.hat,v.hat, main="Transformed Data", 
     xlab="u.hat", ylab="v.hat", pch=16, col="blue")

# fit by IFM
# set dummy copula objects for fitting. The correlation and df
# parameters will be estimated
n.cop = normalCopula(param=0.5, dim=2)
t.cop = tCopula(param=0.5, dim=2, df=3)

start.vals = 0.5
fit.ncop.ifm = fitCopula(copula=n.cop, 
                    data=cbind(u.hat,v.hat), 
                    start=start.vals)
fit.ncop.ifm

# test goodness of fit of estimated copula
gof.test = gofCopula(n.cop, cbind(u.hat,v.hat),
                     method="SnB",
                     estim.method="mpl")

# Create custom distribution by combining fitted margins
# with fitted copula
MSFT.AAPL.margins = c("st", "st")
MSFT.AAPL.parms = list(list(dp=st.msft.fit$dp), 
                       list(dp=st.aapl.fit$dp))
MSFT.AAPL.n.cop = normalCopula(param=fit.ncop.ifm@estimate, dim=2)
myBvd.MSFT.AAPL.fit = mvdc(copula=MSFT.AAPL.n.cop, 
                           margins=MSFT.AAPL.margins, 
                           paramMargins=MSFT.AAPL.parms)

# simulate from fitted bivariate distn
myBvd.MSFT.AAPL.fit.sim = rMvdc(nrow(MSFT.ret), myBvd.MSFT.AAPL.fit)
par(mfrow=c(2,2))
qqPlot(myBvd.fit.sim[,1], main="MSFT Margin", ylab="MSFT quantiles")
qqPlot(myBvd.fit.sim[,2], main="AAPL Margin", ylab="AAPL quantiles")
plot(myBvd.fit.sim, main="Simulated MSFT,AAPL values", 
     xlab="MSFT", ylab="AAPL", pch=16, col="blue")
abline(h=0,v=0)
par(mfrow=c(1,1))

# scatter of simulated returns from custom distribution

plot(myBvd.MSFT.AAPL.fit.sim, main="Simulated MSFT,AAPL values", 
     xlab="MSFT", ylab="AAPL", pch=16, col="blue")
abline(h=0,v=0)
# overlay actual data
points(coredata(MSFT.AAPL.ret), col="red")
legend(x="topleft", legend=c("Simulated", "Actual"), 
       col=c("blue", "red"), pch=c(16,1))


# now fit with t-copula
start.vals = c(0.5, 3)
names(start.vals) = c("rho.1","df")
fit.tcop.ifm = fitCopula(copula=t.cop, 
                         data=cbind(u.hat,v.hat),
                         method="mpl",
                         start=start.vals, optim.method="L-BFGS-B",
                         lower=c(-0.99, 2),
                         upper=c(0.99, 10))
fit.tcop.ifm
t.cop.fit = tCopula(param=fit.tcop.ifm@estimate[1], dim=2,
                    df=fit.tcop.ifm@estimate[2])
myBvd.tcop.fit = mvdc(copula=t.cop.fit, 
                      margins=MSFT.AAPL.margins, 
                      paramMargins=MSFT.AAPL.parms)

# simulate from fitted bivariate distn
myBvd.fit.sim = rMvdc(nrow(MSFT.ret), myBvd.tcop.fit)
par(mfrow=c(2,2))
qqPlot(myBvd.fit.sim[,1], main="MSFT Margin", ylab="MSFT quantiles")
qqPlot(myBvd.fit.sim[,2], main="AAPL Margin", ylab="AAPL quantiles")
plot(myBvd.fit.sim, main="Simulated MSFT,AAPL values", 
     xlab="MSFT", ylab="AAPL", pch=16, col="blue")
abline(h=0,v=0)
par(mfrow=c(1,1))

plot(myBvd.fit.sim, main="Simulated MSFT,AAPL values", 
     xlab="MSFT", ylab="AAPL", pch=16, col="blue")
abline(h=0,v=0)
# overlay actual data
points(coredata(MSFT.AAPL.ret), col="red")
legend(x="topleft", legend=c("Simulated", "Actual"), 
       col=c("blue", "red"), pch=c(16,1))