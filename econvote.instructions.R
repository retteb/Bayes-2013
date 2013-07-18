#########################################################
## Economic voting in Austria: Multinomial logit model ##
#########################################################

## Johannes Karreth
## ICPSR Summer Program 2013

library(foreign)
library(R2jags)

setwd("~/R/Bayes/econvote.mnl")

econ.voting.dat <- read.dta("http://spot.colorado.edu/~joka5204/files/econ.voting.dta")

y <- econ.voting.dat$y
econwor <- econ.voting.dat$econwor
econbet <- econ.voting.dat$econbet
age <- econ.voting.dat$age
income <- econ.voting.dat$income
rural <- econ.voting.dat$rural
gender <- econ.voting.dat$gender
leftrt <- econ.voting.dat$leftrt
N <- length(y)
J <- length(as.numeric(levels(as.factor(y))))  ## number of categories, for our model code

econvote.dat <- list("y", "econwor", "econbet", "age", "income", "rural", "gender", "leftrt", "N", "J")

econvote.model.jags <- function()  {

   for(i in 1:N){
       y[i] ~ dcat(p[i, 1:J])
       
       for (j in 1:J){
           log(q[i,j]) <-  b[1,j] + 
                           b[2,j]*econwor[i] + 
                           b[3,j]*econbet[i] + 
                           b[4,j]*age[i] + 
                           b[5,j]*income[i] + 
                           b[6,j]*rural[i] + 
                           b[7,j]*gender[i] + 
                           b[8,j]*leftrt[i]
           
         p[i,j] <- q[i,j]/sum(q[i,1:J])  ## should be familiar from MLE notes: q is exp(Xb)
           }   # close J loop
          pred[i,1] <- equals(p[i,1], max(p[i,1], p[i,2], p[i,3], p[i,4])) # "1 if p[i,1] = max(p[i,1], p[i,2], p[i,3]), 0 otherwise"
            pred[i,2] <- equals(p[i,2], max(p[i,1], p[i,2], p[i,3], p[i,4]))
            pred[i,3] <- equals(p[i,3], max(p[i,1], p[i,2], p[i,3], p[i,4]))
          pred[i,4] <- equals(p[i,3], max(p[i,1], p[i,2], p[i,3], p[i,4]))
                 
            predcat[i] <- pred[i,1] + 2*pred[i,2] + 3*pred[i,3] + 4*pred[i,4]
       }  # close N loop
           
      for(k in 1:8){
           b[k,1] <- 0          ## MUST set the first set of covariates (for the first outcome category) to 0
               for(j in 2:J){
                    b[k,j] ~ dnorm(0, 0.1)
                   }  # close J loop
               }  # close K loop
       }  # close model loop 

econ.params <- c("b")

econvote.fit <- jags(data=econvote.dat, inits=NULL, econ.params, n.chains=2, n.iter=10, n.burnin=1, model.file=econvote.model.jags)
