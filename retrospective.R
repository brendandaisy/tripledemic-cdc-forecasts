library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)
library(cowplot)

source("model-fitting.R")

add_t_group <- function(df) {
    
}

# Read in current data to long format---------------------------------------------
flu0 <- fetch_flu()

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

###

flu$count_true <- flu$count
tspan1 <- interval(min(flu0$date), "2022-10-01") # the first training window to be plotted

pred_data <- map(duration(0:5, "week"), \(tmax) {
    tspan <- interval(int_start(tspan1), int_end(tspan1) + tmax + duration(5, "week"))
    
    flu_pred <- filter(flu, date %within% tspan) |>
        mutate(count=ifelse(date <= (int_end(tspan1) + tmax), count, NA))
})

fits <- map(pred_data, fit_current_model, graph=us_adj)
ggs <- map2(pred_data, fits, \(df, ft) {
    plot_predictions(df, ft, tback=20, state=c("Vermont")) +
        coord_cartesian(ylim=c(0, NA))
})
plot_grid(plotlist=ggs, nrow=3)



