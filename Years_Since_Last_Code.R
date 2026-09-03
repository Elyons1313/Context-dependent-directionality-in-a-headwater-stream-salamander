### This script is to answer a question from a reviewer about the distribution of the 
### number of years between detections for a dispersal event - you will need to run the full 
### Publication_Code.Rmd document in order to produce the two necessary dataframes (original and dispersal)
### but it will work on its own after that 

o2 <- original
d2 <- dispersal

o2 <- subset(o2, o2$FinalID %in% dispersal$FinalID)



o3 <- o2 %>%
  group_by(FinalID) %>%
  summarise(works = Year)

o4 <- o3 %>%
  arrange(FinalID, works) %>%
  group_by(FinalID) %>%
  mutate(years_since_last = works - lag(works)) %>%
  ungroup()

o5 <- o4 %>%
  filter(!is.na(years_since_last))

o6 <- o5 %>%
  left_join(dispersal, by = "FinalID")

o7 <- subset(o6, o6$works == o6$Year.x)

o8 <- o7 %>%
  distinct(FinalID, .keep_all = TRUE)

summary(as.factor(o8$years_since_last))

