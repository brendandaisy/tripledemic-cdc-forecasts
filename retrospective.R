library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)
library(cowplot)

# TODO: 1st retrospective question! Does training without the 2021-2022 season improve performance?

source("model-fitting.R")
source("retrospective-helpers.R")

# add_t_group <- function(df) {
#     
# }

# Read in current data to long format---------------------------------------------
flu0 <- fetch_flu()

location_info <- distinct(flu0, location, location_name) # get the location coding for saving final results

flu <- flu0 |> 
    select(-location, location=location_name, count=value) |> 
    filter(location != "US")

# Load the state graph/neighborhood matrix
us <- load_us_graph(flu)
us_adj <- us_adj_mat(us)

###

pred_quants <- summarize_quantiles(pred_samples, nat_samps, forecast_date, quantiles_needed)

tmp <- filter(pred_quants, location == "Ohio", date == "2023-12-30")
weighted_interval_score(450, tmp)

# 4 dimensional path forecast
tmp2 <- filter(pred_samples, location == "Ohio")

# produce n X p (100 X 4) matrix
tmp2 <- do.call(cbind, tmp2$count_samp)
mmed <- hdepthmedian(tmp2)
###

forecast_dates <- rev(seq.Date(ymd("2023-11-04"), ymd("2023-02-01"), by="-4 weeks"))

pred_data <- map(forecast_dates, \(fdate) {
    flu_sub <- filter(flu, date < fdate)
    
    flu_pred <- prep_fit_data(flu_sub, weeks_ahead=4) |> 
        mutate(count_true=filter(flu, date <= (fdate + duration("3 weeks")))$count)
})

fit_results <- map2_dfr(forecast_dates, pred_data, \(fdate, df) {
    fit <- fit_current_model(df, fdate, graph=us_adj)
    pred_samps <- sample_count_predictions(df, fit, fdate, nsamp=600)
    tibble_row(forecast_date=fdate, weeks_ahead=4, fit=list(fit), pred_samps=list(pred_samps))
})



ggs <- map2(pred_data, fits, \(df, ft) {
    plot_predictions(df, ft, tback=20, state=c("Vermont")) +
        coord_cartesian(ylim=c(0, NA))
})
plot_grid(plotlist=ggs, nrow=3)



