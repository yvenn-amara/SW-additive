library(riem)
library(tidyverse)
library(zoo)

weather_data = riem_measures(station = "PAO", date_start = , date_end = )

get_weather(data,station,date_start,date_end){
  weather_data = riem_measures(station = "PAO", date_start = date_start, date_end = date_end)
}

dates = tibble(valid = seq(min(weather_data$valid),max(weather_data$valid),"1 min"))

aa = weather_data %>% 
  dplyr::select(valid,tmpf) %>%
  left_join(dates,.) %>%
  mutate(tmpf = na.spline(tmpf))
