library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)

source("model-fitting.R")

# Read in current data to long format---------------------------------------------
dat0 <- read_csv("data/counts_state_alldisease.csv") |> 
    rename(t=`...1`) |> 
    mutate(date=as.Date("2020-01-05") + lubridate::duration(t-1, "week"), .after=t) # starting from 1st Sunday of 2020 (epiweek 1)

dat <- dat0 |> 
    pivot_longer(ak_covid:last_col(), values_to="count") |> 
    separate_wider_delim(name, "_", names=c("state", "type")) |> 
    mutate(epiweek=as.numeric(str_sub_all(orig_index, 5, 6)), count=round(count))

flu0 <- fetch_flu()

# TODO: setting up `snum` and `us_adj` when there are a variable (< all 52) states 

flu <- flu0 |> 
    rename(state=location_name, count=value) |> 
    group_by(date) |> 
    mutate(t=cur_group_id(), epiweek=epiweek(date), .after=date) |> # add a time counter starting from 1 for earliest week
    ungroup() |> 
    filter(state != "US") |> # make sure US is not in training data
    mutate(snum=as.numeric(fct_inorder(state))) # INLA needs groups as ints starting from 1, so add numeric state code

us0 <- read_sf("data/us-state-boundaries.shp")

us <- us0 |> 
    filter(name %in% unique(flu$state)) |> 
    select(state=name, division, region)

state_order <- fct_inorder(unique(flu$state))

# sort order of states to match their order of appearance in flu data
us <- us |> 
    mutate(state=fct_relevel(state, levels(state_order))) |> 
    arrange(state)

us_adj <- us |> 
    poly2nb() |> 
    nb2mat(style="B", zero.policy=TRUE)

# flu <- filter(dat, type == "flu", date > as.Date("2021-06-01"))
# flusight_quantiles <- TODO

###
flu_pred <- prep_model_data(flu)
fit <- fit_current_model(flu_pred, graph=us_adj)
summary(fit)

summarize_predictions(flu_pred, fit)
plot_seasonal(flu_pred, fit)

plot_predictions(flu_pred, fit, tback=5)
plot_predictions(flu_pred, fit, tspan=interval("2022-10-01", "2023-02-01"))
