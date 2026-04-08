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
  filter(Stream == "Bear")

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
  group_by(FinalID, CorrLongLoc) %>%
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



capt.hist3 <- subset(capt.hist2,(rowSums(capt.hist2[,3:14])>0))    #### 
# capt.hist_sum3 <- subset(capt.hist2,(rowSums(capt.hist2[,3:14])>0))

capt.hist3$CorrLongLoc <- as.numeric(capt.hist3$CorrLongLoc)

capt4 <- as.matrix(capt.hist3[3:14])


e2dist <- function(x, y){ 
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
} 



gridx <- seq(1, 50, 1)
grid <- as.matrix(expand.grid(1,gridx))
J <- dim(grid)[1]


T <- K <- 12 # number of years / seasons
sigma <- 4  #Home Range Size
lam0 <- 0.3 # Capture Probability
N <- length(capt.hist3$FinalID)
M <- N*3
phi0 <- 0.60   #apparent survival
tau <- 0.5     #App Surv Squared
Mc <- M # upper M for data augmentation
ntraps <- J


xl <- 1
yl <- 1
xu <- 1
yu <- 50


ntraps <- dim(grid)[1]
nreps <- K
lam0 <- rep(lam0, T)  ## cap prob in each year
phi<- rep(phi0, T)   ## survival in each year
pmat <- lam <- list()  ## empty list - probability matrix
gamma <- NULL   ## empty vector - fill it 
gamma[1:K] <- N/M 
sx <- rep(1,length(capt.hist3$FinalID)) ### X coords of all indivduals
sy <- capt.hist3$CorrLongLoc ### Y coords of all individuals

pos <- matrix(nrow = length(capt.hist3$FinalID), ncol = ntraps)

for(i in 1:length(capt.hist3$FinalID)){
  for(j in 1:ntraps){
    if(capt.hist3$CorrLongLoc[i] == j)
      pos[i, j] <- 1
    else(
      pos[i,j] <- 0
    )
  }
}

alAUG <- matrix(0, nrow = M, ncol = T)   ### actual observed dataset
zAUG <- matrix(0, nrow = M, ncol = T)    ### makes 3 copies of empty matrix - this one is the true latent 
rAUG <- matrix(0, nrow = M, ncol = T)   ####  recruitment matrix 

S <- cbind(sx, sy)     ### coordinates for all individuals
dmat <- e2dist(S, grid)   ### distance matrix - how far is each salamander from each trap



psi <- exp(-(1 / (2 * sigma * sigma)) * dmat * dmat)



for (t in 1:T){
  lam[[t]] <- lam0[t] * psi             ### prob of being in each location * prob of seeing
  pmat[[t]] <- 1 - exp(-lam[[t]])       ### prob of ***finding*** a given individual at a given trap location 
}


## fill the big array
y <- array(0, dim = c(nrow(capt4), ntraps, K))

for (i in 1:nrow(capt4)){   ### for each individual
  for (j in 1:ntraps){      ### for each trap 
    for (k in 1:K){     ### for each time period
      y[i, j, k] <-  (capt4[i,k] * pos[i,j])   ### individual 1, each trap, each year
      ### random number (n, size , prob)
      ### 49 traps, one value (seen or unseen), prob matrix
    }
  }
}


viz_traps <- grid %>% 
  as_tibble() %>%
  ggplot(aes(x = Var1, y = Var2)) +
  geom_point(pch = 3, size = 20) +
  xlim(xl, xu) +
  ylim(yl, yu)
viz_traps

xlim <- c(1,1)
ylim <- c(1, 50)

viz_traps_ac <- viz_traps + 
  geom_point(data = S, aes(x = sx, y = sy), pch = 16, color = "red")
viz_traps_ac


ntot <- length(y[,1,1]) 


SX <- vector(mode = "numeric", length = M)
SY <- vector(mode = "numeric", length = M)

SX[] <- 1
SY[1:length(capt.hist3$CorrLongLoc)] <- capt.hist3$CorrLongLoc

SY[(ntot+1):M] <- sample.int(50, M-ntot, replace = TRUE)

S <- cbind(SX, SY)     ### coordinates for all individuals
DMAT <- e2dist(S, grid) 




yAUG <- array(0, dim = c(M, J, T))
yAUG[1:ntot, , ] <- y


data.in <- as.data.frame(capt.hist3[,3:14])

yAUG2 <- matrix(0, nrow = M, ncol = T)

for(i in 1:nrow(data.in)){
  for(j in 1:ncol(data.in)){
    yAUG2[i,j] <- data.in[i,j]
  }
}




simdat <- list(y = yAUG, 
               z = zAUG,
               r = rAUG,
               gamma = gamma,
               N = apply(zAUG, 2, sum),
               R = apply(rAUG, 2, sum), 
               SX = SX, 
               SY = SY)

dataugTMJ <- aperm(yAUG, c(3,1,2))



zst <- c(rep(1,ntot), rep(0,M-ntot))
zst <- cbind(zst,zst,zst,zst,zst,zst,zst,zst,zst,zst,zst,zst)

constants <- list(M = M, 
                  Time = T)



inits <- list(psi = runif(1),
              phi = runif(1),
              p.mean = runif(1),
              p.tel = runif(1),
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
    cp[i,1] <- z[i,1] * p.mean
    Y[i,1] ~ dbinom(cp[i,1], 9)    ## Y = Number of encounters
    A[i,1] <- (1-z[i,1])
    # delta[i] <- 1
    
    for(j in 2:8){
      a1[i,j] <- sum(z[i, 1:j])
      A[i,j] <- 1- (step(a1[i,j]-1))
      mu[i,j] <- (phi * z[i, j-1]) + (gamma[j] * A[i,j-1])
      z[i,j] ~ dbern(mu[i,j])
      cp[i,j] <- z[i,j] * p.mean
      Y[i,j] ~ dbinom(cp[i,j], 9)
      # Q[j] <- M/10
    }
    for(j in 9:12){
      a1[i,j] <- sum(z[i, 1:j])
      A[i,j] <- 1- (step(a1[i,j]-1))
      mu[i,j] <- (phi * z[i, j-1]) + (gamma[j] * A[i,j-1])
      z[i,j] ~ dbern(mu[i,j])
      cp[i,j] <- z[i,j] * p.tel
      Y[i,j] ~ dbinom(cp[i,j], 9)
      # Q[j] <- M/10
    }}
  
  }
)




parameters <- c("phi","psi",
                "N","p.mean", 
                "gamma")



reps <- 30000

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
write.csv(results, "Telemetry_Jolly_Seber_Abundance_Bear.csv")


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


