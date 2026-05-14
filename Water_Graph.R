library(dplyr)
library(ggplot2)

water <- read.csv("Stream_Data_Winsor.csv")

data <- water %>%
  group_by(Stream, Corr_Pos) %>%
  summarize(mean = mean(Wetted_Width, na.rm = TRUE), sd = sd(Wetted_Width, na.rm = TRUE),
            min = min(Wetted_Width,na.rm = TRUE), max = max(Wetted_Width,na.rm = TRUE),
            min2023 = min(Wetted_Width[Year == 2023]),
            min2024 = min(Wetted_Width[Year == 2024]),
            min2025 = min(Wetted_Width[Year == 2025]))
            # min2024 = min(Wetted_Width[Year == 2024 & !is.na(Wetted_Width)]),
            # min2025 = min(Wetted_Width[Year == 2025 & !is.na(Wetted_Width)]))

comp <- water %>%
  group_by(Stream) %>%
  summarize(meanlow = mean(Wetted_Width[Corr_Pos <= 500]),
            sdlow = sd(Wetted_Width[Corr_Pos <= 500]),
            meanup = mean(Wetted_Width[Corr_Pos > 500]),
            sdup = sd(Wetted_Width[Corr_Pos > 500]),
            CVlow = sdlow/meanlow * 100,
            CVup = sdup/meanup * 100)


ggplot(data, aes(x = data$Corr_Pos, y = data$mean))+
  geom_point()
       