# Classic SEIR model for transmission dynamics of the novel coronavirus (COVID-19)
# beta is the average number of infected individuals per infectious subject per unit time, 
# alpha is the reciprocal of average latent period, 
# mu is the rate of recovery (reciprocal of duration of the infection).

# SEIR model for Wuhan
# WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
seir.WH = function(N = 1000, SUS = .5*N,  EXP=0, INF = 1, IMM = N-SUS-EXP-INF, days = 200, beta=0.8, alpha=1/6, mu=1/5, 
                TMP = 0.01, CLS=1, D2D=0.5){
  
  M = data.frame(ID = 1:days, SUS = NA, EXP = NA, INF = NA, IMM = NA)
  M[1,] = c(1, SUS, EXP, INF, IMM)
  
  for(i in 2:days){
    
    if (i > 90)  beta <- beta * (1-TMP)^(i-90) # temperature effect Wuhan
    
    if (i < (31+23)) {
      M$SUS[i] = M$SUS[i-1] - beta * M$INF[i-1] * M$SUS[i-1]/N 
      M$EXP[i] = M$EXP[i-1] + beta * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
      M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu 
      M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
    else if (i >= (31+23) & i < (31+31+8)) { # city closure effect
      beta2 = beta * CLS
      M$SUS[i] = M$SUS[i-1] - beta2 * M$INF[i-1] * M$SUS[i-1]/N 
      M$EXP[i] = M$EXP[i-1] + beta2 * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
      M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu 
      M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
    else {
      beta3 = beta * D2D # door to door checking
      M$SUS[i] = M$SUS[i-1] - beta3 * M$INF[i-1] * M$SUS[i-1]/N 
      M$EXP[i] = M$EXP[i-1] + beta3 * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
      M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu 
      M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
  }
  return(invisible(M))
}
# WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW



# SEIR model for Beijing
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
seir.BJ = function(N=1000, SUS=.5*N,  EXP=0, INF=1, IMM=N-SUS-EXP-INF, days=300, beta=0.8, alpha=1/6, mu=1/10, 
                TMP=0.04, IMP=10){
  
  M = data.frame(ID = 1:days,SUS = NA, EXP = NA, INF = NA, IMM = NA)
  M[1,] = c(1, SUS, EXP, INF, IMM)
  
  for(i in 2:days){
    
    if (i > 40)  beta <- beta * (1-TMP)^(i-40) # temperature effect Beijing
    
    if (i < 21){ # stricter control
       M$SUS[i] = M$SUS[i-1] - beta * M$INF[i-1] * M$SUS[i-1]/N 
       M$EXP[i] = M$EXP[i-1] + beta * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
       M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu  + IMP[i] # imported cases every day
       M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
    
    else{
      beta2 = beta*0.85
      M$SUS[i] = M$SUS[i-1] - beta2 * M$INF[i-1] * M$SUS[i-1]/N 
      M$EXP[i] = M$EXP[i-1] + beta2 * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
      M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu  + IMP[i] # imported cases every day
      M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
   }
  }
  return(invisible(M))
}
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB



# Number of confirmed cases in Wuhan
WH = c(41,45,62,121,198,270,375,444,444,549,572,618,698,1590,1905,2261,3215,4109,5142, 6384,8351,10117,11618,
       13603,14982,16902,18454) #Jan15-Feb11

par(mfrow=c(1,1))

# Fig. 3A Match confirmed cases for Wuhan: infectious R0 = beta/mu/(S/N) = 2.3*5/2 = 5.75 
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
DD = seir.WH(N=10000000, beta=2.3, alpha=1/6, mu=1/5, CLS=.5) # Wuhan
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,50), ylim=c(0, 10000),xlab="Days since Dec. 24",ylab="Number of infected cases")
points(c(46:72)-22, WH, pch=16, col="blue")

beta=2.3
for (i in 1:200){
  B = rnorm(1, beta, beta*0.15)
  DD = seir.WH(N=10000000, beta=B, alpha=1/6, mu=1/5) # Wuhan
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.1))
}
DD = seir(N=10000000, beta=2.3, alpha=1/6, mu=1/5) # Wuhan
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(46:72)-22, WH, pch=16, col="blue")
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE





# Fig. 3B Real situation for Wuhan (R0 = 5 in the beginning) - Infectious population
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
DD = seir.WH(N=10000000, beta=10/5, alpha=1/6, mu=1/5, days=200, CLS=.5, D2D=.3) # Wuhan

# plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,120), ylim=c(0, 10000000),xlab="Days since Dec. 1, 2019",ylab="Number of cases")
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(30, 130), ylim=c(0, 700000),xlab="Days since Dec. 1, 2019",ylab="Number of infected cases")
#plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0, 80), ylim=c(0, 1000),xlab="Days since Dec. 1, 2019",ylab="Number of infected cases")
points(c(46:72), WH, pch=16, col="blue", cex=.5)
abline(v=54, lty=2) # city closure
abline(v=70, lty=2) # door to door checking

peak=numeric(200)
beta=10/5

for (i in 1:200){
  B = rnorm(1, beta, beta*0.15)
  DD = seir.WH(N=10000000, beta=B, alpha=1/6, mu=1/5, CLS=.5, D2D=.3) # Wuhan
  # lines(DD$ID, DD$SUS, col=adjustcolor("black", alpha.f=0.2))
  # lines(DD$ID, DD$IMM, col=adjustcolor("green", alpha.f=0.2))
  # lines(DD$ID, DD$EXP, col=adjustcolor("yellow",  alpha.f=0.2))
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.2))
  # peak[i] = DD[DD$INF == max(DD$INF), c('ID')] # peak day
  # peak[i] = DD[DD$INF == max(DD$INF), c('INF')] # peak No. cases
  peak[i] = DD[32, c('INF')]
}

median(peak);mean(peak);sd(peak) # peak of the pandemic wave
DD = seir.WH(N=10000000, beta=beta, alpha=1/6, mu=1/5, CLS=.5, D2D=.3)  # Wuhan
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(46:72)+0, WH, pch=16, col="blue", cex=.8)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


# Fig. 3C Real situation for Wuhan (R0 = 5 in the beginning) - all populations
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
DD = seir.WH(N=10000000, beta=10/5, alpha=1/6, mu=1/5, days=200, CLS=.5, D2D=.3) # Wuhan

plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,140), ylim=c(0, 10000000),xlab="Days since Dec. 1, 2019",ylab="Number of cases")
points(c(46:72), WH, pch=16, col="blue", cex=.5)
abline(v=54, lty=2) # city closure
abline(v=70, lty=2) # door to door checking

peak=numeric(200)
beta=10/5

for (i in 1:200){
  B = rnorm(1, beta, beta*0.15)
  DD = seir.WH(N=10000000, beta=B, alpha=1/6, mu=1/5, CLS=.5, D2D=.3) # Wuhan
  lines(DD$ID, DD$SUS, col=adjustcolor("black", alpha.f=0.2))
  lines(DD$ID, DD$IMM, col=adjustcolor("green", alpha.f=0.2))
  lines(DD$ID, DD$EXP, col=adjustcolor("yellow",  alpha.f=0.2))
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.2))
  # peak[i] = DD[DD$INF == max(DD$INF), c('ID')] # peak day
  # peak[i] = DD[DD$INF == max(DD$INF), c('INF')] # peak No. cases
  peak[i] = DD[32, c('INF')]
}
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


par(mar=c(4.5,4.5,3.5,2))
# Fig. 3D temperature effect for Wuhan
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
DD = seir.WH(N=10000000, beta=10/5, alpha=1/6, mu=1/5, days=200, CLS=.5, TMP=0.01, D2D=.3) # Wuhan
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(30,160), ylim=c(0, 220000),xlab="Days since Dec. 1, 2019",ylab="Number of infected cases")
# points(c(46:67) + 0, WH, pch=16, col="blue", cex=.5)

TMP=0.1
for (i in 1:300){
  TMP = rnorm(1, TMP, TMP*0.3)
  DD = seir.WH(N=10000000, beta=10/5, TMP=TMP, alpha=1/6, mu=1/5, CLS=.5, D2D=.3) # Wuhan
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.1))
}

DD = seir.WH(N=10000000, beta=10/5, alpha=1/6, mu=1/5, CLS=.5, TMP=0.1, D2D=.3) # Wuhan
lines(DD$ID, DD$INF, col="red", lwd=2)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE






BJ = c(1,5,10,14,26,36,51,67,80,96,112,132,156,183,212,228,253,274,297,315,326)
par(mfrow=c(3,1), mar=c(4,4.2,3.2,2))

# Fig. 4A Matching confirmed cases for Beijing
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; imp; mean(imp); sd(imp)
DD = seir.BJ(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=30, IMP=imp) # Beijing

plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,25), ylim=c(0, 400),xlab="Days since Jan. 19, 2020",ylab="Number of infected cases")
points(c(1:21) + 0, BJ, pch=16, col="blue", cex=.8)

for (i in 1:200){
  imp = imp + rnorm(1, 0, 2) # SD 2 was based on the daily report of confirmed cases in Beijing
  DD = seir.BJ(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=30, IMP=imp) # Beijing
  DD = DD[!DD$INF[5]<0,]
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.1))
}

imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; imp; mean(imp); sd(imp)
DD = seir.BJ(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=30, IMP=imp) # Beijing
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(1:21) + 0, BJ, pch=16, col="blue", cex=1)
legend("topleft", "A", bty="n")
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE





# Fig. 4B Simulate R0 for Beijing
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; 
imp = c(imp, rep(0, 170))
DD = seir.BJ(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,150), ylim=c(0, 2000),xlab="Days since Jan. 19, 2020",ylab="Number of infected cases")
points(c(1:length(BJ)) + 0, BJ, pch=16, col="blue", cex=.5)

peak.day = numeric(200)
peak = numeric(200)
beta=2.8/5
TMP=0.001
for (i in 1:200){
  # imp = imp + rnorm(1, 0, 3.85)
  B = rnorm(1, beta, beta*0.15)
  tmp = rnorm(1, TMP, TMP*0.4)
  DD = seir.BJ(N=22000000, beta=B, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=tmp) # Beijing
  peak.day[i] = DD[DD$INF == max(DD$INF), c('ID')]
  peak[i] = DD[DD$INF == max(DD$INF), c('INF')]
  DD$INF[DD$INF<0] <- 0
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.1))
}
median(peak);mean(peak);sd(peak) # peak of the pandemic wave in Beijing
median(peak.day);mean(peak.day);sd(peak.day) 

imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5;imp = c(imp, rep(0, 170))
DD = seir.BJ(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(1:length(BJ)) + 0, BJ, pch=16, col="blue", cex=.8)
legend("topleft", "B", bty="n")
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE




# Fig. 4C Simulate R0 series for Beijing
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; 
imp = c(imp, rep(0, 170))
TMP=0.001
BETA = seq(0.2, 1.6, by=0.2)

DD = seir.BJ(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=TMP) # Beijing
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,130), ylim=c(0, 1500),xlab="Days since Jan. 19, 2020",
     ylab="Number of infected cases", main="")
points(c(1:length(BJ)) + 0, BJ, pch=16, col="blue", cex=.5)

for (i in 1:length(BETA)){
  DD = seir.BJ(N=22000000, beta=BETA[i]*2/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=TMP) # Beijing
  lines(DD$ID, DD$INF, col=rainbow(length(BETA))[i], lwd=2)
  
  legend(121, 300+200*i,  
         expression(paste(R[0])), 
         text.col = rainbow(length(BETA))[i],  
         bty = "n", cex = 0.7)
  legend(124, 300+200*i,  
         paste("= ", BETA[i], sep=""),
         text.col = rainbow(length(BETA))[i],  
         bty = "n", cex = 0.7)
}

legend("topleft", "C", bty="n")
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE




# not run
# no temperature
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; 
imp = c(imp, rep(0, 270))
TMP=0

BETA = seq(0.2, 2, by=0.2)
DD = seir.BJ(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=300, IMP=imp, TMP=TMP) # Beijing
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,300), ylim=c(0, 800000),xlab="自北京发生新型冠状病毒（1月19日）后的天数",
     ylab="新型冠状病毒肺炎发病数", main="病毒活性不受气温升高的影响")
points(c(1:length(BJ)) + 0, BJ, pch=16, col="blue", cex=.5)

for (i in 1:length(BETA)){
  DD = seir.BJ(N=22000000, beta=BETA[i]*2/5, alpha=1/6, mu=1/5, days=300, IMP=imp, TMP=TMP) # Beijing
  lines(DD$ID, DD$INF, col=rainbow(length(BETA))[i], lwd=2)
  
  legend(275, 300000+25000*i,  
         expression(paste(R[0])), 
         text.col = rainbow(length(BETA))[i],  
         bty = "n", cex = 0.7)
  legend(282, 300000+25000*i,  
         paste("= ", BETA[i], sep=""),
         text.col = rainbow(length(BETA))[i],  
         bty = "n", cex = 0.7)
}
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


# Simulate imported cases for Beijing
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; 
imp = c(imp, rep(0, 170))
DD = seir.BJ(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,150), ylim=c(0, 2000),xlab="Days since Jan. 19, 2020",ylab="Number of infected cases")
points(c(1:length(BJ)) + 0, BJ, pch=16, col="blue", cex=.5)

for (i in 1:200){
  imp = imp + rnorm(1, 0, 3.85)
  DD = seir.BJ(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
  DD = DD[!DD$INF[5]<0,]
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.2))
}

imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5;imp = c(imp, rep(0, 170))
DD = seir.BJ(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(1:length(BJ)) + 0, BJ, pch=16, col="blue", cex=.5)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

