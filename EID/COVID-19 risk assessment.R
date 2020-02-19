
# Risk assessment for COVID-19 outbreak in cities in China

# read data
# D = read.csv("D:/coop/2019-nCoV/data/2019_NCoV_city_daily_cases_Feb14.csv", header=T)
D = read.csv("2019_NCoV_city_daily_cases_Feb14.csv", header=T)
head(D)

D = D[!is.na(D$Pop),] # remove cities with incomplete data
# D = D[!D$Diff<0,] # one city withdraw a case


# Checking current risk classification
#########################################################################################
library(randomForest)
D$Risk = as.factor(D$Risk)
RF <- randomForest(Risk ~ Diff + F14 + Pop + GDP + Density, data=D,  ntree=1000, importance=TRUE)
RF
varImpPlot(RF, main='')

imp <- importance(RF); imp = as.data.frame(imp)
par(mfrow=c(1,1))
barplot(imp$MeanDecreaseGini)
#########################################################################################




# Machine learning process
########################################################################################## 
########################################################################################## 

RISK = sample(c(1:4), nrow(D) , rep=T)
DD = data.frame(RISK, Diff = D$Diff, F14 =D$F14, Pop= D$Pop, GDP=D$GDP, Density=D$Density)

# initial values
DD$RISK[D$F14>100 & D$F14<1000] <- 2
DD$RISK[D$Prov=="湖北"] <- 1
DD$RISK[D$F14>10 & D$F14<=100] <- 3
DD$RISK[D$F14<=10] <- 4

# repeat the model
#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
RF <- randomForest(RISK ~ Diff + F14 + Pop + 
                     GDP + Density, data=DD, 
                   ntree=1000, importance=TRUE)
RF

RISK2 = predict(RF, data=DD)
RISK2 = round(RISK2)
DD$RISK = RISK2
#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

# compare 
out = cbind(RISK2, Real=D$Risk)
plot(jitter(out))
plot(jitter(DD$RISK), jitter(as.numeric(D$Risk))) # the same
########################################################################################## 
########################################################################################## 




# partial plot
#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
imp <- importance(RF) 
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)] #sort importance

# Plot partial effects
op <- par(mfrow=c(2, 3),mar=c(4,4,2,2))
for (i in seq_along(impvar)) {
  partialPlot(RF, D, impvar[i], xlab=impvar[i], #Partial effects
              ylab='Class', main='')
}
#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP


# Linear model
#######################################################################
fit = lm(Risk ~ log(Diff+1) + log(F14+1) + Pop + GDP + Density, data=D)
fit = step(fit)
summary(fit)
anova(fit)
#######################################################################
