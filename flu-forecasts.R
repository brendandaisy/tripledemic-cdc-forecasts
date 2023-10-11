library(tidyverse)
library(INLA)

# Read in current data to long format (Maya's pipeline here)----------------------
dat0 <- read_csv("data/counts_state_alldisease.csv") |> 
    rename(t=`...1`) |> 
    mutate(date=as.Date("2020-01-05") + lubridate::duration(t-1, "week"), .after=t) # starting from 1st Sunday of 2020 (epiweek 1)

dat <- dat0 |> 
    pivot_longer(ak_covid:last_col(), values_to="count") |> 
    separate_wider_delim(name, "_", names=c("state", "type")) |> 
    mutate(epiweek=as.numeric(str_sub_all(orig_index, 5, 6)), count=round(count))

flu <- filter(dat, type == "flu", date > as.Date("2021-06-01"))
flusight_quantiles <- TODO

###
flu_pred <- prep_model_data(flu)
fit <- fit_current_model(flu_pred)
summary(fit)

summarize_predictions(flu_pred, fit)
plot_predictions(flu_pred, fit, tspan=interval("2022-10-01", "2023-02-01"), file="figs/flu-predictions-10-11.pdf")
