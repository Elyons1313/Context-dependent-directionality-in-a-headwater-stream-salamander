
set.seed(1234)
### SIMULATING DATA
library(nimble)
library(ggplot2)
library(tidyverse)
library(basicMCMCplots)
library(coda)
library(magrittr)





### GET THE SAL DATA 
original <- read.csv("2025_Capture_History.csv")

#Only Study Streams
original <- original %>% 
  filter(Stream == "Zigzag")

# Only seen alive (not just pinged by telemetry reader)
original <- original %>%
  filter(POL == "X" | POL == "X?" | POL == "LCA")

original <- original %>%
  filter(Reach != "Inner")

original <- original %>%
  filter(Species != "DF")

# treat metamorphs as larvae
for(i in 1:length(original$Stage)){
  if(is.na(original$Stage[i])) {x = FALSE}
  else if ((original$Stage[i]) == "M") {
    original$Stage[i] <- "L"
  }
}


# convert data to CH

ch <- original[c("FinalID", "CorrLongLoc", "Stage","Primary",
                 "Year", "OldNew", "Date")]
ch$Year <- (ch$Year-2011)

ch$CorrLongLoc <- as.numeric(ch$CorrLongLoc)


ch <- cbind(ch, detect = 1)

Old <- subset(ch, OldNew == "O")
New <- subset(ch, OldNew == "N")


capt.hist <- ch %>%
  spread(Year, detect, fill = 0)
# 
colnames(capt.hist)[7:18] <- paste("Year", seq(1:12), sep="")

capt.hist$CorrLongLoc <- as.numeric(capt.hist$CorrLongLoc)

# capt.hist <- subset(capt.hist, capt.hist$CorrLongLoc > 1000 & capt.hist$CorrLongLoc <=1050 )
# capt.hist$CorrLongLoc <- as.numeric(capt.hist$CorrLongLoc)-1000



capt.hist2 <- capt.hist[,c(1,2,7:18)]

capt.hist2 <- capt.hist2 %>%
  group_by(FinalID) %>%
  dplyr::summarize(Y1 = sum(Year1), Y2 = sum(Year2),
                   Y3 = sum(Year3), Y4 = sum(Year4),
                   Y5 = sum(Year5), Y6 = sum(Year6),
                   Y7 = sum(Year7), Y8 = sum(Year8),
                   Y9 = sum(Year9), Y10 = sum(Year10),
                   Y11 = sum(Year11), Y12 = sum(Year12))

# capt.hist2 <- capt.hist2 %>%
#   group_by(FinalID) %>%
#   dplyr::summarize(CorrLongLoc = mean(CorrLongLoc),
#                    Y1 = sum(Year1), Y2 = sum(Year2), 
#                    Y3 = sum(Year3), Y4 = sum(Year4),
#                    Y5 = sum(Year5), Y6 = sum(Year6),
#                    Y7 = sum(Year7), Y8 = sum(Year8))



capt.hist3 <- subset(capt.hist2,(rowSums(capt.hist2[,2:13])>0))    #### 
# capt.hist_sum3 <- subset(capt.hist2,(rowSums(capt.hist2[,3:14])>0))



capt4 <- as.matrix(capt.hist3[2:13])






T <- K <- 12 # number of years / seasons
lam0 <- 0.3 # Capture Probability
N <- length(capt.hist3$FinalID)
M <- N*3
phi0 <- 0.60   #apparent survival
tau <- 0.5     #App Surv Squared
Mc <- M # upper M for data augmentation



nreps <- K
lam0 <- rep(lam0, T)  ## cap prob in each year
phi<- rep(phi0, T)   ## survival in each year
pmat <- lam <- list()  ## empty list - probability matrix
gamma <- NULL   ## empty vector - fill it 
gamma[1:K] <- N/M 



alAUG <- matrix(0, nrow = M, ncol = T)   ### actual observed dataset
zAUG <- matrix(0, nrow = M, ncol = T)    ### makes 3 copies of empty matrix - this one is the true latent 
rAUG <- matrix(0, nrow = M, ncol = T)   ####  recruitment matrix 




ntot <- length(capt.hist3$FinalID) 





data.in <- as.data.frame(capt.hist3[,2:13])

yAUG2 <- matrix(0, nrow = M, ncol = T)

for(i in 1:nrow(data.in)){
  for(j in 1:ncol(data.in)){
    yAUG2[i,j] <- data.in[i,j]
  }
}




zst <- c(rep(1,ntot), rep(0,(M-ntot)))
zst <- cbind(zst,zst,zst,zst,zst,zst,zst,zst,zst,zst,zst,zst)

constants <- list(M = M, 
                  Time = T,
                  cmr = 8,
                  tel = 10,
                  combo = 12)



inits <- list(psi = runif(1),
              phi = runif(1),
              p.mean = runif(1),
              p.tel = runif(1),
              p.combo = runif(1),
              gamma = runif(12),
              z = zst)



data <- list(Y = yAUG2)



code <- nimbleCode({
  psi ~ dunif(0,1)
  phi  ~ dunif(0,1)
  p.mean ~ dunif(0,1)
  
  for(t in 1:Time){
    N[t] <- sum(z[1:M,t])
    gamma[t] ~ dunif(0,1)
  }
  
  
  
  for(i in 1:M){
    z[i,1] ~ dbern(psi)
    cp[i,1] <- z[i,1] * p.mean    ### just constant for now - fix after adding telemetry
    Y[i,1] ~ dbinom(cp[i,1], 9)    ## Y = Number of encounters
    A[i,1] <- (1-z[i,1])
    # delta[i] <- 1
    
    for(j in 2:cmr){
      a1[i,j] <- sum(z[i, 1:j])
      A[i,j] <- 1- (step(a1[i,j]-1))
      mu[i,j] <- (phi * z[i, j-1]) + (gamma[j] * A[i,j-1])
      z[i,j] ~ dbern(mu[i,j])
      cp[i,j] <- z[i,j] * p.mean
      Y[i,j] ~ dbinom(cp[i,j], 9)
      
    }
    for(k in cmr+1:tel){
      a1[i,k] <- sum(z[i, 1:k])
      A[i,k] <- 1- (step(a1[i,k]-1))
      mu[i,k] <- (phi * z[i, k-1]) + (gamma[k] * A[i,k-1])
      z[i,k] ~ dbern(mu[i,k])
      cp[i,k] <- z[i,k] * p.tel
      Y[i,k] ~ dbinom(cp[i,k], 9)
    }
    
    for(q in tel+1:combo){
      a1[i,q] <- sum(z[i, 1:q])
      A[i,q] <- 1- (step(a1[i,q]-1))
      mu[i,q] <- (phi * z[i, q-1]) + (gamma[q] * A[i,q-1])
      z[i,q] ~ dbern(mu[i,q])
      cp[i,q] <- z[i,q] * p.tel
      Y[i,q] ~ dbinom(cp[i,q], 9)
    }
  }
})




parameters <- c("phi","psi",
                "N","p.mean", "p.tel", "p.combo", 
                "gamma")



reps <- 50000

simtest<-nimbleModel(code= code, name="test", constants = constants, data = data, inits = inits)

confMCMC <- configureMCMC(simtest,
                          monitors = parameters,
                          thin = 1,
                          useConjugacy = TRUE,
                          enableWAIC = FALSE)

nimMCMC <- buildMCMC(confMCMC)
Cnim <- compileNimble(simtest)
CnimMCMC <- compileNimble(nimMCMC, project = simtest)
mcmcout<- runMCMC(CnimMCMC,
                  niter = reps,
                  nburnin = 0.2*reps,
                  nchains = 3,
                  inits = inits,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE,
                  WAIC = FALSE)

mcmcout$summary$all.chains

results <- mcmcout$summary$all.chains
write.csv(results, "Jolly_Seber_Abundance_Zigzag.csv")


#Assess convergence
library(coda)

mcmc.coda <- mcmcout$samples


diag.plots <- function(x,invt){
  xx <- colnames(mcmc.coda$chain1)[
    grep(x,colnames(mcmc.coda$chain1),invert = invt)]
  for(j in xx){
    if(mean(mcmc.coda[,j][[1]])==0){
      next   #skip if parm set to zero
    }
    plot(mcmc.coda[,j], main=j) 
    print(j)
    print(gelman.diag(mcmc.coda[,j]))
  }
}

for(i in c("p", "N")){
  diag.plots(i,invt = FALSE)
}


