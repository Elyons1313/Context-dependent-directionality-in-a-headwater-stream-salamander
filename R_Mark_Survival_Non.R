

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

#Only Study Streams
original <- original %>%
  filter(Stream == "Bear"| Stream == "Paradise"| Stream == "Zigzag")
# # 
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

capt.hist$CorrLongLoc <- as.numeric(capt.hist$CorrLongLoc)

capt.hist <- subset(capt.hist, capt.hist$CorrLongLoc >=501 & capt.hist$CorrLongLoc <=1250 )
capt.hist$CorrLongLoc <- as.numeric(capt.hist$CorrLongLoc)-0


capt.hist2 <- capt.hist[,c(1,2,7:15)]

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

for(i in 1:nrow(capt.hist2)){
  for(j in 1:nrow(capt.hist)){
    if(capt.hist2$FinalID[i]==capt.hist$FinalID[j]){
      if(is.na(capt.hist$Stage[j])) {x = FALSE}
      else(capt.hist2$Stage[i] <- capt.hist$Stage[j])
    }
  }
}

capt.hist2$Stream <- NA

for(i in 1:nrow(capt.hist2)){
  for(j in 1:nrow(capt.hist)){
    if(capt.hist2$FinalID[i]==capt.hist$FinalID[j]){
      if(is.na(capt.hist$Stream[j])) {x = FALSE}
      else(capt.hist2$Stream[i] <- capt.hist$Stream[j])
    }
  }
}



capt.hist3 <- subset(capt.hist2,(rowSums(capt.hist2[,3:10])>0))    #### 
# capt.hist_sum3 <- subset(capt.hist2,(rowSums(capt.hist2[,3:14])>0))

capt.hist3$CorrLongLoc <- as.numeric(capt.hist3$CorrLongLoc)


capt.hist3$edges <- 0
for(i in 1:nrow(capt.hist3)){
  if(capt.hist3$Stream[i] == "Paradise" & (capt.hist3$CorrLongLoc[i] <= 100 | capt.hist3$CorrLongLoc[i] >= 1150)){
    capt.hist3$edges[i] <- 1 
  }
  else if (capt.hist3$Stream[i] == "Bear" & (capt.hist3$CorrLongLoc[i] <= 100 | capt.hist3$CorrLongLoc[i] >= 1300)){
    capt.hist3$edges[i] <- 1 
  }
  else if (capt.hist3$Stream[i] == "Zigzag" & (capt.hist3$CorrLongLoc[i] <= 100 | capt.hist3$CorrLongLoc[i] >= 1400)){
    capt.hist3$edges[i] <- 1 
  }
}

edges <- subset(capt.hist3, capt.hist3$edges > 0) 
non <- subset(capt.hist3, capt.hist3$edges == 0) 

capt4 <- as.matrix(non[3:10])

for(i in 1:nrow(capt4)){
  for(j in 1:ncol(capt4)){
    if(capt4[i,j] >= 1)
      capt4[i,j] <- 1
  }
}




chs <- vector(length = nrow(capt4))

for(i in 1:length(chs)){
  chs[i] <- paste(as.character(capt4[i,1]), as.character(capt4[i,2]),as.character(capt4[i,3]),as.character(capt4[i,4]),
                  as.character(capt4[i,5]), as.character(capt4[i,6]), as.character(capt4[i,7]), as.character(capt4[i,8]), 
                  sep = "" )
}

capt5 <- cbind(capt4, chs, non$Stage)


capt6 <- capt5[,9:10]






# colnames(capt6) <- c("ch")

colnames(capt6) <- c("ch", "stage")

capt6 <- as.data.frame(capt6)
capt6$stage <- as.factor(capt6$stage)

# 
# 
# 
# 
# 
# data(dipper, package = "marked")
# 
# # Examine data structure
# head(dipper)
# 
# 
# # Jolly-Seber models (POPAN formulation) are open population models, and 
# # can be used to estimate abundance by including two more parameters than the CJS
# 
# # Additional parameters:
# # Nsuper (or "superpopulation") = total number of individuals available to enter population throughout study
# # pent ("probability of entry") =  the rate at which individuals enter the population from Nsuper (via births and immigration)
# 
# # WARNING: there is no adequate GOF tests for Jolly-Seber models. 
# # One common method: Test equivalent structure of CJS model with R2ucare (previous tutorials).
# 
# # This tests *some* assumptions of Phi and p.
# # Jolly-Seber models have an additional assumption:
# # marked AND unmarked animals have same p (R2ucare doesn't test this)
# # This assumption is required to estimate total abundance (sum of marked and unmarked animals in population)
# 
# # First, process data (Notice model = "JS", previous version = "CJS")
# dipper.js.proc <- process.data(dipper, 
#                                model = "JS", 
#                                groups = "sex")
# 

# 
# # Second, make design data (from processed data)
# dipper.js.ddl <- make.design.data(dipper.js.proc)

# 
# 
# fit.js.dipper.models <- function(){
#   # Phi formulas
#   Phi.dot <- list(formula=~1)
#   Phi.time <- list(formula=~time)
#   # p formulas
#   p.dot <- list(formula=~1)
#   # pent formulas. pent estimates MUST SUM to 1 (for each group).
#   # This is constained using a Multinomial Logit link
#   pent.time <- list(formula=~time)
#   pent.sex <- list(formula=~sex)
#   pent.dot <- list(formula=~1)
#   # Nsuper formulas. Don't confuse "N" from model with predicted population size
#   N.sex <- list(formula=~sex)
#   N.dot <- list(formula=~1)
#   cml <- create.model.list(c("Phi","p", "pent", "N"))
#   results <- crm.wrapper(cml, data = dipper.js.proc, ddl = dipper.js.ddl,
#                          external = FALSE, accumulate = FALSE, hessian = TRUE)
#   
#   return(results)
# }
# 
# 
# # Run function
# dipper.js.models <- fit.js.dipper.models()
# 
# 
# dipper.js.models
# 
# # Look at estimates of top model (row number on left of model table, or using name)
# dipper.js.models[[1]]  # or dipper.js.models[["Phi.dot.p.dot.pent.dot.N.dot"]] or dipper.js.models$Phi.dot.p.dot.pent.dot.N.dot
# 
# 
# # The estimates above are not on probability scale (or in individuals for N)
# # (e.g. Phi, p on logit scale, pent on mlogit scale)
# # Predict (real) values using top model
# dipper.js.predicted <- predict(dipper.js.models[[1]]) # [[1]] just calls the model row according to the model table.
# 
# # Look at predictions of real parameters
# dipper.js.predicted 
# 
# 
# # Abundance (N) is derived from the estimated parameters
# # We will estimate population size at each time by making a dataframe of estimates and calculating N
# # We will use the predicted estimates from the top-performing model (in this case: "dipper.js.predicted")
# 
# # NOTE: the below method will have to be adjusted based on your final model and the number of capture events
# N.derived <- data.frame(occ = c(1:7), # 7 events
#                         Phi = c(rep(dipper.js.predicted$Phi$estimate, 6), NA),   # 6 survival estimates all the same
#                         Nsuper = rep(dipper.js.predicted$N$estimate + nrow(dipper), 7), # Nsuper estimate + number of marked animals
#                         pent = c(1-sum(dipper.js.predicted$pent$estimate), dipper.js.predicted$pent$estimate)) # Sum of all pent must be 1
# 
# # Set-up empty vector for calculating N
# N.derived$N <- NA
# 
# # The inital population size (N[1]) = Nsuper * (1 - sum(all other pent estimates))
# # This is because of the link function for estimating pent.
# # The sum of all pent parameters MUST equal 1 (therefore, one less must be estimated)
# N.derived$N[1] <- (N.derived$Nsuper[1] * N.derived$pent[1])
# 
# # Subsequent population sizes are estimated by calculating surviving individuals (N[t-1] * Phi[t]), and
# # Adding new births (Nsuper * pent[t])
# for(i in 2:nrow(N.derived)){
#   N.derived$N[i] <- (N.derived$N[i-1]*N.derived$Phi[i-1]) + (N.derived$Nsuper[i] * N.derived$pent[i])
# }
# 
# # Look at what we did
# N.derived









###########################################################################################################


sal.js.proc <- process.data(capt6,
                            model = "JS",
                            groups = "stage")

sal.js.ddl <- make.design.data(sal.js.proc)


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
sal.js.models[[1]]  # or dipper.js.models[["Phi.dot.p.dot.pent.dot.N.dot"]] or dipper.js.models$Phi.dot.p.dot.pent.dot.N.dot


# The estimates above are not on probability scale (or in individuals for N)
# (e.g. Phi, p on logit scale, pent on mlogit scale)
# Predict (real) values using top model
sal.js.predicted <- predict(sal.js.models[[1]]) # [[1]] just calls the model row according to the model table.

# Look at predictions of real parameters
sal.js.predicted 


# Abundance (N) is derived from the estimated parameters
# We will estimate population size at each time by making a dataframe of estimates and calculating N
# We will use the predicted estimates from the top-performing model (in this case: "dipper.js.predicted")

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



