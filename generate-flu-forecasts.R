# generate-flu-forecasts.R--------------------------------------------------------
# Main-level script to run and produce weekly forecasts---------------------------
# --------------------------------------------------------------------------------
library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)

source("read-flu-data.R")
source("model-prep-and-fit.R")
source("model-summarize-and-plot.R")

# Data prep and model fitting-----------------------------------------------------
# Read in current data
forecast_date <- ymd("2023-12-02") # Saturday following submission date
flu0 <- fetch_flu()

location_info <- distinct(flu0, location, location_name) # get the location coding for saving final results

flu <- flu0 |> 
    select(-location, location=location_name, count=value)

# Load the state graph/neighborhood matrix
us <- load_us_graph(flu)
us_adj <- us_adj_mat(us)

fit_df <- prep_fit_data(flu, weeks_ahead=4)
fit <- fit_current_model(fit_df, forecast_date, graph=us_adj, joint=TRUE)

# get US predictions
nat_samps <- sample_national(fit_df, fit, forecast_date, nsamp=2000)

# add a column `count_samp` containing posterior predictive samples for each week
# this is the dataframe to be used in all summaries downstream
pred_samples <- sample_count_predictions(fit_df, fit, forecast_date, nsamp=2000)

# Produce the two summary dataframes to be submitted------------------------------
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

cleaned_forecasts_quantiles <- summarize_quantiles(pred_samples, nat_samps, forecast_date, quantiles_needed)
cleaned_forecasts_ratechange <- summarize_rate_change(pred_samples, nat_samps, forecast_date)

# Save the prediction summaries
cleaned_forecasts_quantiles |>
    mutate(output_type_id = as.character(output_type_id)) |>
    bind_rows(cleaned_forecasts_ratechange) |>
    left_join(location_info, by=join_by(location == location_name)) |>  # replace location names with the location numeric ID
    select(-location, location=location.y) |> 
    relocate(location, .after=target_end_date) |> 
    write_csv(paste0("weekly-predictions/", forecast_date,"-UGA_flucast-INfLAenza.csv"))

# Plot the predictions for this week----------------------------------------------
library(cowplot)
theme_set(theme_cowplot())

curr_season_data <- flu |>
    filter(date >= (forecast_date-days(60)))

forecast_df <- cleaned_forecasts_quantiles |> 
    filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |> 
    # left_join(locations, by='location') |> 
    spread(output_type_id, value)

# sample a few forecasts from the joint distribution to show correlations across states and time
ts_ids <- sample(2000, 4)

sampled_timeseries <- pred_samples |> 
    unnest(count_samp) |> 
    group_by(date, location) |> 
    mutate(samp_id=1:n()) |> 
    filter(samp_id %in% ts_ids) |> 
    ungroup()

plots <- unique(flu$location) |> 
    map(plot_state_forecast, 
        curr_season_data = curr_season_data, 
        forecast_df = forecast_df,
        sampled_ts=sampled_timeseries)

plot_grid(plotlist = plots) |> 
    save_plot(filename=paste0("weekly-predictions/prediction-fig-", forecast_date, ".pdf"), base_height=12, base_asp=1.6, bg='white')

cleaned_forecasts_ratechange <- cleaned_forecasts_ratechange |> 
    mutate(output_type_id=fct_relevel(output_type_id, c("large_increase", "increase", "stable", "decrease")))

plots_pmfs <- map(
    unique(flu$location), 
    \(loc) plot_state_pmf(loc, curr_season_data, cleaned_forecasts_ratechange)
)

legend <- pmf_legend(cleaned_forecasts_ratechange)

save_plot(
    filename=paste0("weekly-predictions/ratechange-fig-", forecast_date, ".pdf"), 
    plot=ggdraw(plot_grid(plotlist=plots_pmfs)) + draw_plot(legend, 0.76, -0.4),
    base_height=12, base_asp=1.6, bg='white'
)


# Validate forecast file --------------------------------------------------
library(hubValidations)
sub_validation <- hubValidations::validate_submission(hub_path = '~/projects/FluSight-forecast-hub/',
                                    file_path = 'UGA_flucast-INFLAenza/2023-11-04-UGA_flucast-INFLAenza.csv')
## Want all green checkmarks
sub_validation

## Want to make sure there are no missing required values
sub_validation$req_vals$missing



# plot_predictions(count_pred, state=c("Florida", "Puerto Rico"), tback=100)
# 
# nat_summ <- sample_national(count_pred, quantiles)
# 
# pdf("figs/predictions-11-01.pdf")
# 
# gg <- flu |> 
#     filter(t > 95) |> 
#     group_by(t) |> 
#     summarise(nat_count=sum(count)) |> 
#     ggplot(aes(t, nat_count)) +
#     geom_ribbon(aes(t, ymin=`q25%`, ymax=`q75%`), fill="gray60", alpha=0.6, data=nat_summ, inherit.aes=FALSE) +
#     geom_ribbon(aes(t, ymin=`q2.5%`, ymax=`q97.5%`), fill="gray80", alpha=0.6, data=nat_summ, inherit.aes=FALSE) +
#     geom_point() +
#     geom_line(aes(x=t, y=mean), col="tomato3", data=nat_summ, inherit.aes=FALSE) +
#     labs(x="Weeks", y="Hospitalizations") +
#     theme_bw()
# 
# print(gg)
# 
# pred_dates <- unique(filter(flu_pred, is.na(count))$date)
# year(pred_dates) <- year(pred_dates) - 1
# pred_int <- interval(min(pred_dates), max(pred_dates) + ddays())
# 
# for (s in unique(flu_pred$state)) {
#     ret_pred_dat <- flu_pred |> 
#         filter(date %within% pred_int, state == s) |> 
#         mutate(date=as.Date(date + dyears()))
#     
#     gg <- plot_predictions(count_pred, tback=20, state=s) +
#         geom_point(data=ret_pred_dat, col="maroon1", size=1.05)
#     
#     print(gg)
# }
# 
# dev.off()
# 
# count_pred[which(count_pred$date == "2023-10-28"),]$count <- NA
# nat_summ <- sample_national(count_pred, quantiles)
# 
# pred_summ <- nat_summ |> 
#     select(-mean) |> 
#     pivot_longer(-t, names_to="quantile") |> 
#     mutate(state="National", quantile=0.01*parse_number(quantile))
# 
# pred_summ <- count_pred |> 
#     filter(is.na(count)) |> 
#     select(state, t, contains("%")) |> 
#     pivot_longer(contains("%"), names_to="quantile") |> 
#     mutate(quantile=0.01*parse_number(quantile)) |> 
#     bind_rows(pred_summ)
# 
# write_csv(pred_summ, "predictions-11-01.csv")
# write_csv(flu_pred, "input-pred-data-11-01.csv")
