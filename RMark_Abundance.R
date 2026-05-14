set.seed(1234)
### SIMULATING DATA
library(nimble)
library(ggplot2)
library(tidyverse)
library(basicMCMCplots)
library(coda)
library(magrittr)
library(marked)





### GET THE SAL DATA 
original <- read.csv("2025_Capture_History.csv")

# #Only Study Streams
original <- original %>%
  filter(Stream == "Bear"| Stream == "Paradise"| Stream == "Zigzag")


# # # Change this line to look at other streams, remove to look at all three simultaneously
## options: "Bear" , "Paradise", "Zigzag"
original <- original %>%
  filter(Stream == "Paradise")

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
                 "Year", "OldNew", "Date", "Stream")]
ch$Year <- (ch$Year-2011)

ch$CorrLongLoc <- as.numeric(ch$CorrLongLoc)


ch <- cbind(ch, detect = 1)

Old <- subset(ch, OldNew == "O")
New <- subset(ch, OldNew == "N")


capt.hist <- ch %>%
  spread(Year, detect, fill = 0)
# 
colnames(capt.hist)[8:19] <- paste("Year", seq(1:12), sep="")


## corr is the distance from the mainstem (in meters) - used for spatial location
capt.hist$CorrLongLoc <- as.numeric(capt.hist$CorrLongLoc)


## use this line to subset to smaller areas of the stream - "How many salamanders are between 200 and 300m) etc.
## filtering may cause some issues in the following code if certain conditions are not met
## for example - model will not work if there were zero adults in the range selected - set a new range
capt.hist <- subset(capt.hist, capt.hist$CorrLongLoc >= 0 & capt.hist$CorrLongLoc <=1500)

### set to lower bound of filter - first position surveyed is position one 
capt.hist$CorrLongLoc <- as.numeric(capt.hist$CorrLongLoc)-0

## keep relevant columns
capt.hist2 <- capt.hist[,c(1,2,7:15)]


## groups by position as well as ID - useful for OPSCR models or to track individual dispersal
# capt.hist_sum <- capt.hist2 %>%
#   group_by(FinalID, CorrLongLoc) %>%
#   dplyr::summarize(Y1 = sum(Year1), Y2 = sum(Year2), 
#                    Y3 = sum(Year3), Y4 = sum(Year4),
#                    Y5 = sum(Year5), Y6 = sum(Year6),
#                    Y5 = sum(Year7), Y6 = sum(Year8),
#                    Y5 = sum(Year9), Y6 = sum(Year10),
#                    Y5 = sum(Year11), Y6 = sum(Year12))

capt.hist2 <- capt.hist2 %>%
  group_by(FinalID) %>%
  dplyr::summarize(CorrLongLoc = mean(CorrLongLoc),
                   Y1 = sum(Year1), Y2 = sum(Year2), 
                   Y3 = sum(Year3), Y4 = sum(Year4),
                   Y5 = sum(Year5), Y6 = sum(Year6),
                   Y7 = sum(Year7), Y8 = sum(Year8))

capt.hist2$Stage <- NA

### matches stage to individual after converting to capture history
for(i in 1:nrow(capt.hist2)){
  for(j in 1:nrow(capt.hist)){
    if(capt.hist2$FinalID[i]==capt.hist$FinalID[j]){
        if(is.na(capt.hist$Stage[j])) {x = FALSE}
        else(capt.hist2$Stage[i] <- capt.hist$Stage[j])
        }
  }
}


## only keep individuals that were caught at least once (necessary if you change the years)
capt.hist3 <- subset(capt.hist2,(rowSums(capt.hist2[,3:10])>0))    #### 
# capt.hist_sum3 <- subset(capt.hist2,(rowSums(capt.hist2[,3:14])>0))

capt.hist3$CorrLongLoc <- as.numeric(capt.hist3$CorrLongLoc)

capt4 <- as.matrix(capt.hist3[3:10])

##converts to binary - also possible to do in the group_by argument above (change sum() to max())
## used this to make sure it was getting ALL the detections not just the first one
for(i in 1:nrow(capt4)){
  for(j in 1:ncol(capt4)){
    if(capt4[i,j] >= 1)
      capt4[i,j] <- 1
  }
}

### converts everything to a single character of CH, not 12 separate columns - needed for marked
### yes there are way better ways to do this, but this was faster than debugging the apply function I first tried
chs <- vector(length = nrow(capt4))

for(i in 1:length(chs)){
    chs[i] <- paste(as.character(capt4[i,1]), as.character(capt4[i,2]),as.character(capt4[i,3]),as.character(capt4[i,4]),
                    as.character(capt4[i,5]), as.character(capt4[i,6]), as.character(capt4[i,7]), as.character(capt4[i,8]), 
                    sep = "" )
}

capt5 <- cbind(capt4, chs, capt.hist3$Stage)


capt6 <- capt5[,9:10]






# colnames(capt6) <- c("ch")

colnames(capt6) <- c("ch", "stage")

capt6 <- as.data.frame(capt6)
capt6$stage <- as.factor(capt6$stage)


### process the Jolly Seber Model - group by life history stage (can change to sex/stream/etc. etc if desired)
sal.js.proc <- process.data(capt6,
                            model = "JS",
                            groups = "stage")

sal.js.ddl <- make.design.data(sal.js.proc)

### function to fit models (function created by James Patterson - Introduction to CJS in R )
### comments etc in the next section are from his webpage - I don't claim any ownership
### https://jamesepaterson.github.io/jamespatersonblog/2020-04-26_introduction_to_CJS.html
fit.js.sal.models <- function(){
  # Phi formulas
  Phi.dot <- list(formula=~1)
  Phi.time <- list(formula=~time)
  # p formulas
  p.dot <- list(formula=~1)
  # pent formulas. pent estimates MUST SUM to 1 (for each group).
  # This is constained using a Multinomial Logit link
  pent.time <- list(formula=~time)
  pent.sex <- list(formula=~stage)
  pent.dot <- list(formula=~1)
  # Nsuper formulas. Don't confuse "N" from model with predicted population size
  N.sex <- list(formula=~stage)
  N.dot <- list(formula=~1)
  cml <- create.model.list(c("Phi","p", "pent", "N"))
  results <- crm.wrapper(cml, data = sal.js.proc, ddl = sal.js.ddl,
                         external = FALSE, accumulate = FALSE, hessian = TRUE)
  
  return(results)
}


# Run function
sal.js.models <- fit.js.sal.models()


sal.js.models

# Look at estimates of top model (row number on left of model table, or using name)
sal.js.models[[1]]  


# The estimates above are not on probability scale (or in individuals for N)
# (e.g. Phi, p on logit scale, pent on mlogit scale)
# Predict (real) values using top model
sal.js.predicted <- predict(sal.js.models[[1]]) # [[1]] just calls the model row according to the model table.

# Look at predictions of real parameters
sal.js.predicted 


# Abundance (N) is derived from the estimated parameters
# We will estimate population size at each time by making a dataframe of estimates and calculating N
# We will use the predicted estimates from the top-performing model (in this case: "sal.js.predicted")

# NOTE: the below method will have to be adjusted based on your final model and the number of capture events
N.derived.sal <- data.frame(occ = c(1:8), # 7 events
                        Phi = c(rep(sal.js.predicted$Phi$estimate, 7), NA),   # 6 survival estimates all the same
                        Nsuper = rep(sal.js.predicted$N$estimate + nrow(capt6), 8), # Nsuper estimate + number of marked animals
                        pent = c(1-sum(sal.js.predicted$pent$estimate), sal.js.predicted$pent$estimate)) # Sum of all pent must be 1

# Set-up empty vector for calculating N
N.derived.sal$N <- NA

# The inital population size (N[1]) = Nsuper * (1 - sum(all other pent estimates))
# This is because of the link function for estimating pent.
# The sum of all pent parameters MUST equal 1 (therefore, one less must be estimated)
N.derived.sal$N[1] <- (N.derived.sal$Nsuper[1] * N.derived.sal$pent[1])

# Subsequent population sizes are estimated by calculating surviving individuals (N[t-1] * Phi[t]), and
# Adding new births (Nsuper * pent[t])
for(i in 2:nrow(N.derived.sal)){
  N.derived.sal$N[i] <- (N.derived.sal$N[i-1]*N.derived.sal$Phi[i-1]) + (N.derived.sal$Nsuper[i] * N.derived.sal$pent[i])
}

# Look at what we did
N.derived.sal


mean(N.derived.sal$N)


