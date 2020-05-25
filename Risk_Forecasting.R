# instalacja bibliotek
install.packages(c("tseries", "zoo", "forecast", "FinTS", "rugarch", "ggplot2"))

library("ggplot2")
library("tseries")
library("zoo")
library("forecast")
library("FinTS")
library("rugarch")

#Zgranie do œrodowiska ceny zamkniêcia(Adjusted Close) indeksu S&P 500
#na przedziale czasowym: 2010-04-12 - 2020-04-09
fund = get.hist.quote(instrument="^GSPC", start="2010-04-12",  end="2020-04-09",
                      quote="AdjClose", provider="yahoo", compression="d", retclass="zoo")

#definicja zmiennej r jako lofarytmicznych stóp zwrotu indeksu S&P 500
r <- diff(log(fund))

#Wykres cen indeksu i stóp zwrotu
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(fund, main="Cena za jednostke S&P 500", xlab="", ylab="")
plot(r, main="stopy zwrotu", xlab="", ylab="")


##########################################################
# Ustawienia do analizy                                  #
##########################################################
                       
R  <- coredata(r)
N <- length(r)
p  <- 0.01                              ## poziom tolerancji VaR
H  <- 1                                 ## horyont

mu = sum(R)/N                           ## œrednia ze stóp zwrotu
R0    <- R - mu
M2    <- sum(R0^2)/N                    ## Drugi moment centralny
M3    <- sum(R0^3)/N                    ## Trzeci moment centralny
M4    <- sum(R0^4)/N                    ## Czwarty moment centralny
sig <- sqrt(M2)                         ## zmiennosc/odchylenie standardowe
S   <- M3/(sig^3)                       ## skosnosc
K   <- M4/(sig^4)                       ## kurtoza

S

##############################################################################

##########
# Wykresy#
##########

##############################################################################

#Histogram stóp zwrotu
hist(R, freq=FALSE, col="gray", xlab="stopa zwrotu", ylab = 'czêstoœæ', main="Rozk³ad logarytmicznych stóp zwrotu", breaks = 50)


#T-student
#x <- r
#curve(dt(x, 30), from = -5, to = 5, col = "orange", main= 'T-student plot',
#      xlab = "quantile", ylab = "density", lwd = 2)
#curve(dt(x, 10), from = -5, to = 5, col = "dark green", add = TRUE, lwd = 2)
#curve(dt(x, 5), from = -5, to = 5, col = "sky blue", add = TRUE, lwd = 2)
#curve(dt(x, 1), from = -5, to = 5, col = "grey40", add = TRUE, lwd = 2)
#legend("topleft", legend = paste0("DF = ", c(1, 5, 10, 30)),
#       col = c("grey40", "sky blue", "dark green", "orange"),
#       lty = 1, lwd = 2)


#Kernel Density
#plot(density(r), main="Kernel Density of log returns")
#polygon(density(r), col="red", border="blue")

#QQplot
qqPlot(R, distribution="norm", main = paste('Wykres kwantyl-kwantyl'),
       xlab=paste('Kwantyl Teorytyczny'),
       ylab = paste('Kwantyl rozk³adu'))


#Autokorelacja i autokorelacja cz¹stkowa
par(mfrow=c(2,1), cex = 0.7, bty="l")
acf(R, main = paste('Autokorelacja logarytmicznych stóp zwrotu'),xlab = paste('obserwacja wstecz'),
    ylab = paste('ACF'))
pacf(R, main = paste('Autokorelacja cz¹stkowa logarytmicznych stóp zwrotu'),
     xlab = paste('obserwacja wstecz'),
     ylab = paste('PACF'))

#efekt ARCH, ACF(R^2)
plot(acf(R^2), main = paste('Efekt ARCH logarytmicznych stóp zwrotu'),
     xlab = paste('obserwacja wstecz'),
     ylab = paste('ACF^2'))

#Kod do wygenerowania wykresu przedstawiaj¹cego VaR oraz ES

mu <- 10/1000
sigma <- 30/1000
x   <- seq(mu-3.5*sigma, mu+3.5*sigma, length.out=100)
p_x <- dnorm(x, mean=mu, sd=sigma)
plot(x, p_x, type='l', xlab="Return", ylab="Density of Asset", main="VaR i ES")

VaR <- qnorm(0.05, mean=mu, sd=sigma)

x <-seq(VaR, mu-3.5*sigma, length.out = 50)
xx <- c(VaR, x , mu-3.5*sigma)
px <- c( 0 , dnorm(x, mean=mu, sd=sigma),     0       )
abline(v=mean(x) ,lwd=2, col="black",lty=c(2,2))

polygon(xx,px,col='red')

text(x=VaR+25, y=0.002, labels=paste("VaR = 5%"))
text(x = VaR-20, y = 0.01, labels = paste("ES"))

###############################################################################

##########################################################
# Symulacje Monete Carlo dla rozk³adu normalnego         #
##########################################################

##############################################################################

M      <- 10000  # liczba losowan w symulacjach Monte Carlo
m      <- mean(R)     # Œrednia ze stóp zwrotu
s      <- sd(R)       # odch. st. stóp zwrotu

#Generowanie pseudolosowych liczb dla rozk³.norm. na podstawie
# œredniej i odch. st. logarytmicznych stóp zwrotu
draws0    <- rnorm(M, mean=m, sd=s) 

#Deklaracja funkcji licz¹cej VaR na podstawie symulacji Monte Carlo
my.mc_VaR<-function(x){
#M<-10000
p<-0.01 #kwantyl(0.01), dla którego jest liczona wartoœæ VaR
m<-mean(x)
s<-sd(x)
draws0    <- rnorm(M, mean=m, sd=s)
draws1    <- sort(draws0)                
M0        <- floor(M*p)                               
VAR       <- draws1[M0]
return(VAR)
}

#Aplikacja powy¿szej funkcji do wyliczenia VaR na podstawie 100
# ostatnich obserwacji
w_length = 5000  ## szerokosc okna estymacji
MC_VaR <- rollapply(r, width=w_length, function(x) my.mc_VaR(x),
                by=1, align="right")

#Deklaracja funkcji licz¹cej ES na podstawie symulacji Monte Carlo
my.mc_ES<-function(x){
  #M<-10000
  p<-0.01 #kwantyl(0.01), dla którego jest liczona wartoœæ ES
  m<-mean(x)
  s<-sd(x)
  draws0    <- rnorm(M, mean=m, sd=s)
  draws1    <- sort(draws0)                
  M0        <- floor(M*p)                               
  ES<- mean(draws1[1:M0])
  return(ES)
}

#Aplikacja powy¿szej funkcji do wyliczenia ES na podstawie 100
# ostatnich obserwacji
w_length = 500  ## szerokosc okna estymacji
MC_ES <- rollapply(r, width=w_length, function(x) my.mc_ES(x),
                    by=1, align="right")
###########
#predykcja#
###########

M      <- 10000  # liczba losowan w symulacjach Monte Carlo
m      <- mean(R[2017:2516])     # Œrednia ze stóp zwrotu z ostatnich 500 dni
s      <- sd(R[2017:2516])       # odch. st. stóp zwrotu z ostatnich 500 dni

#Generowanie pseudolosowych liczb dla rozk³.norm. na podstawie
# œredniej i odch. st. logarytmicznych stóp zwrotu z okresu ostatnich 500 dni
draws0    <- rnorm(M, mean=m, sd=s) 

draws1    <- sort(draws0)                 # uporzadkowane obserwacje
M0        <- floor(M*p)                   # obserwacja dla p-tego kwintyla                   
VaR_N0    <- draws1[M0]                   # kwantyl p
# Oczekiwana strata ES
ES_N0     <- mean(draws1[1:M0])           # srednia wartosc w "ogonie"

# Prezentacja wyników na wykresie
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge(r, MC_VaR, MC_ES), plot.type="single",main = paste("Symulacja Monte Carlo VaR i ES 1% dla rozk³.norm."), col=c(1,2,3), ylab="", xlab="" )
legend("topleft", c("VaR", "ES"), lty=1, col=2:3)
plot(density(draws0), main=paste("predykcja dla horyzontu = 1"), xlab="stopa zwrotu", ylab="f. gestosci")
abline(v=c(VaR_N0,ES_N0) ,lwd=2, col=c("red","black"),lty=c(2,2))
legend("right", c("ES","VaR"), lty=c(1,1), col=c("black","red"), bty="n")

################################################################################

#########
# GARCH #
#########

#################################################################################

#First Finding Best Mean Model Using ARIMA(Ewentualnie)
#auto.arima(r, trace = T, 
#           ic = "bic",
#           test = "kpss")

#Specyfikacja modelu GARCH
GARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                        variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                        fixed.pars=list(shape=5), distribution.model = "std")
                       
#Aplikacja modelu GARCH do stóp zwrotu
GARCHfit <- ugarchfit(data = r, spec = GARCHspec, solver="solnp")

#Wspó³czynniki powy¿szego modelu
round(coef(GARCHfit),5)
# VaR i ES (rozklad t-Studenta o 5 stopniach swobody)
#gdzie:
q    <- qdist("std", p=p, shape=5)
qf   <- function(x) qdist("std", p=x, shape=5)

GARCHvar   <- fitted(GARCHfit) + sigma(GARCHfit)*q 
GARCHes    <- fitted(GARCHfit) + sigma(GARCHfit)*(1/p * integrate(qf, 0, p)$value)

#predykcja
GARCHfct   <- ugarchforecast(GARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(GARCHfct))
sT  <- as.numeric(sigma(GARCHfct))
VaR_GARCH  = mT + sT*q
ES_GARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

#Prezentacja wyników
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge( r, GARCHvar, GARCHes), plot.type="single", main=paste(100*p,"% VaR i ES z modelu GARCH(1,1))",  sep=""), col=c(1,2,3), ylab="" , xlab="")
legend("topleft", c("VaR", "ES"), lty=1, col=2:3)
plot(hist(r, breaks = 100), xlim=c(-0.15, 0.10), main=paste("predykcja dla horyzontu = 1"), xlab="stopa zwrotu", ylab="f. gestosci")
abline(v=c(VaR_GARCH, ES_GARCH) ,lwd=2, col=c("red","black"),lty=c(2,2))
legend("right", c("ES","VaR"), lty=c(1,1), col=c("black","red"), bty="n")

################################################################################

##########################################################
# Symulacje Monte Carlo dla modelu GARCH                 #
##########################################################

################################################################################

# liczba losowan w symulacjach Monte Carlo
#M      <- length(r)     

# Specyfikacja modelu: GARCH(1,1) z r. t-studenta
#GARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
#                        variance.model=list(model="sGARCH", garchOrder=c(1,1)),
#                        fixed.pars=list(shape=5), distribution.model = "std")
#estymacja modelu
#garch.fit  <- ugarchfit(data=r, spec=garch.spec)
#M <- length(r)
#H <- 1
#p <- 0.01
#MC_GARCH_VaR <- rollapply(r, width=w_length,
#                          function(w) {
#                            fit <- ugarchfit(data=w, spec=GARCHspec, solver="hybrid")
#                            #frc <- ugarchforecast(fit, n.ahead=1)
#                            garch.sim   <- ugarchsim(fit, n.sim = H, n.start = 0, m.sim = M, startMethod = "sample")
#                            draws <- garch.sim@simulation$seriesSim
#                            draws <- colSums(draws)
#                            draws2 <- sort(draws)                
#                            M0 <- floor(M*p)                          
#                            draws2[M0]
#                          },
#                          by=1, align="right")

#MC_GARCH_ES <- rollapply(r, width=w_length,
#                         function(w) {
#                           fit <- ugarchfit(data=w, spec=GARCHspec, solver="hybrid")
#                           #frc <- ugarchforecast(fit, n.ahead=1)
#                           garch.sim   <- ugarchsim(fit, n.sim = H, n.start = 0, m.sim = M, startMethod = "sample")
#                           draws <- garch.sim@simulation$seriesSim
#                           draws <- colSums(draws)
#                           draws2 <- sort(draws)                
#                           M0 <- floor(M*p)                          
#                           mean(draws2[1:M0])
#                         },
#                         by=1, align="right")
#GARCHvarh   <- fitted(garch.fit) + sigma(garch.fit)*q 
#GARCHesh    <- fitted(garch.fit) + sigma(garch.fit)*(1/p * integrate(qf, 0, p)$value)

#Predykcja
# symulacje MC z modelu GARCH
#garch.sim   <- ugarchsim(garch.fit, n.sim = H, n.start = 0, m.sim = M, startMethod = "sample")
#draws       <- garch.sim@simulation$seriesSim

# Obliczenia VaR
#draws2      <- colSums(draws) 
#draws2      <- sort(draws2)                 # uporzadkowane obserwacje
#M0          <- floor(M*p)                   # obserwacja dla p-tego kwintyla                   
#VaRH_GARCH  <- draws2[M0]                   # porownaj z: quantile(draws0,p)
# Oczekiwana strata ES
#ESH_GARCH   <- mean(draws2[1:M0])           # numeryczne calkowanie dla stop zwrotu




#par(mfrow=c(2,1), cex = 0.7, bty="l")
#plot(merge( r, MC_GARCH_VaR, MC_GARCH_ES), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR i ES z modelu MC GARCH(1,1))",  sep=""), ylab="" )
#legend("topleft", c("VaR", "ES"), lty=1, col=2:3)
#plot(density(draws2), main=paste("predykcja dla horyzontu = 1"), xlab="stopa zwrotu", ylab="f. gestosci")
#abline(v=c(VaRH_GARCH, ESH_GARCH) ,lwd=2, col=c("red","black"),lty=c(2,2))
#legend("right", c("ES","VaR"), lty=c(1,1), col=c("black","red"), bty="n")
#VaRH_GARCH
#ESH_GARCH


#Porównanie 2 metod

plot(hist(r, breaks = 100), xlim=c(-0.15, 0.10),main=paste("Porównanie 2 metod - VaR 1%"), xlab="stopa zwrotu", ylab="f. gestosci")
abline(v=c(VaR_N0, VaR_GARCH) ,lwd=3, col=c("red","black", "blue"),lty=c(2,2))
legend("right", c("MC - rozk³. norm.","GARCH"), lty=c(1,1), col=c("red","black", "blue"), bty="n")
plot(hist(r, breaks = 100), xlim=c(-0.15, 0.10),main=paste("Porównanie 2 metod - ES 1%"), xlab="stopa zwrotu", ylab="f. gestosci")
abline(v=c(ES_N0, ES_GARCH) ,lwd=3, col=c("red","black", "blue"),lty=c(2,2))
legend("right", c("MC - rozk³. norm.","GARCH"), lty=c(1,1), col=c("red","black", "blue"), bty="n")

################################################################################

################
#  Backtesing  #
################

################################################################################

#BACKTEST MC
# empiryczny udzial przekroczen (powinien byc rowny poziomowi tolerancji p)
backN   <- length(MC_VaR)

# przesuwamy indeksy aby wskazywaly na dni prognozy
MCvar <- lag(MC_VaR, -1)
# odpowiadajace prognozom zrealizowane zwroty, rr - realized returns
rr <- r[index(MCvar)]
#Przekroczenia VaR
etaMC  <- tail(rr<=MCvar, backN)
sum(etaMC)
#Wykres przekroczeñ dla VaR
VaRplot(alpha=p, actual=tail(rr,backN), VaR=tail(MCvar,backN))
title(main=paste0("Przekroczenia VaR 1%, Monte Carlo, okno estymacji: ",
                  w_length," obs."), ylab = paste('wartoœæ stopy zwrotu'),
                  xlab = paste('czas'))

#Testy pokrycia(Kupca i Christoffersena)
VaRTest(alpha = 0.01, r[500:2516], MC_VaR)

#BACKTEST GARCH
w_length = 500
GARCH_VaR <- rollapply(r, width=w_length,
                       function(w) {
                         fit <- ugarchfit(data=w, spec=GARCHspec, solver="hybrid")
                         frc <- ugarchforecast(fit, n.ahead=1)
                         quantile(frc, p)
                       },
                       by=1, align="right")
# przesuwamy indeksy aby wskazywaly na dni prognozy
GARCH_VaR <- lag(GARCH_VaR, -1)
# odpowiadajace prognozom zrealizowane zwroty, rr - realized returns
rr <- r[index(GARCH_VaR)]
var <- tail(GARCH_VaR, backN)
act <- tail(rr,backN)
#Przekroczenia VaR
etaGARCH <- coredata(tail(act<var, backN))

#Wykres przekroczeñ VaR
VaRplot(p, act, var)
title(main=paste0("Przekroczenia GARCH, okno estymacji: ",w_length," obs."),
                    xlab = paste('czas'),
                    ylab = paste('wartoœæ stopy zwrotu'))

#Testy pokrycia(Kupca i Christoffersena)
VaRTest(alpha = 0.01, r[501:2516], GARCH_VaR)

# Backtest GARCH MC
# przesuwamy indeksy aby wskazywaly na dni prognozy
#MC_GARCH_VaR <- lag(MC_GARCH_VaR, -1)
# odpowiadajace prognozom zrealizowane zwroty, rr - realized returns
#rr <- r[index(MC_GARCH_VaR)]
#var <- tail(MC_GARCH_VaR, backN)
#act <- tail(rr,backN)
#etaGARCH_MC <- coredata(tail(act<var, backN))
#VaRplot(p, act, var)
#title(main=paste0("GARCH MC, r.t-st. okna estymacji ",w_length," obs."))
#tKupiec(p,etaGARCH_MC)
#tChris1(etaGARCH_MC)
#tChris2(p,etaGARCH_MC)




#ctrl = list(tol = 1e-7, delta = 1e-9)
#res_garch11_roll <- ugarchroll(GARCHspec, r, n.start = 120, refit.every = 1, refit.window = "moving", solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = 0.01, keep.coef = TRUE, solver.control = ctrl, fit.control = list(scale = 1))
#res_garch11_roll <- ugarchroll(spec=GARCHspec, data=R, refit.every=5, forecast.length=2000, refit.window="moving", window.size=300, calculate.VaR=TRUE, VaR.alpha=c(0.01, 0.025, 0.05), solver.control = ctrl, fit.control = list(scale = 1))

#report(res_garch11_roll, type = "VaR", VaR.alpha = 0.01, conf.level = 0.99)
#plot(GARCHfit)




