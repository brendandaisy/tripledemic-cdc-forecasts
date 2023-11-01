library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)

source("model-fitting.R")

# Read in current data to long format---------------------------------------------
flu0 <- fetch_flu()

flu <- flu0 |> 
    rename(state=location_name, count=value) |> 
    group_by(date) |> 
    mutate(t=cur_group_id(), epiweek=epiweek(date), .after=date) |> # add a time counter starting from 1 for earliest week
    ungroup() |> 
    filter(state != "US") |> # make sure US is not in training data
    mutate(snum=as.numeric(fct_inorder(state))) # INLA needs groups as ints starting from 1, so add numeric state code

flu <- flu |> 
    group_by(date) |> 
    mutate(ex_count=(population / sum(population)) * sum(count), .after=count) |> 
    ungroup()

us <- load_us_graph(flu)
us_adj <- us_adj_mat(us)

flu <- flu |> 
    left_join(distinct(us, state, division, region))
    # mutate(count_true=count)

flu_pred <- prep_pred_data(flu)

quantiles <- c(
    0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 
    0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990
)

fit <- fit_current_model(flu_pred, q=quantiles, graph=us_adj)
summary(fit)

plot_predictions(flu_pred, fit, state=c("Florida", "Puerto Rico"), tback=20)

jsamp_fvals <- exp(inla.rjmarginal(5000, fit$selection)$samples)

tslice <- list(1:52, 53:104, 105:156, 157:208, 209:260)
nat_summ <- imap_dfr(tslice, \(idx, t) {
    tsum <- map_dbl(1:5000, \(s) sum(jsamp_fvals[idx, s]))
    qs <- quantile(tsum, quantiles)
    names(qs) <- str_c("q", names(qs))
    tibble_row(t=tlast+t, mean=mean(tsum), !!!qs)
})

pdf("figs/predictions-10-25.pdf")
tlast <- max(filter(flu_pred, !is.na(count))$t)

gg <- flu |> 
    filter(t > 95) |> 
    group_by(t) |> 
    summarise(nat_count=sum(count)) |> 
    ggplot(aes(t, nat_count)) +
    geom_ribbon(aes(t, ymin=`q25%`, ymax=`q75%`), fill="gray60", alpha=0.6, data=nat_summ, inherit.aes=FALSE) +
    geom_ribbon(aes(t, ymin=`q2.5%`, ymax=`q97.5%`), fill="gray80", alpha=0.6, data=nat_summ, inherit.aes=FALSE) +
    geom_point() +
    geom_line(aes(x=t, y=mean), col="tomato3", data=nat_summ, inherit.aes=FALSE) +
    labs(x="Weeks", y="Hospitalizations") +
    theme_bw()

print(gg)

pred_dates <- unique(filter(flu_pred, is.na(count))$date)
year(pred_dates) <- year(pred_dates) - 1
pred_int <- interval(min(pred_dates), max(pred_dates))

for (s in unique(flu_pred$state)) {
    ret_pred_dat <- flu_pred |> 
        filter(date %within% pred_int, state == s) |> 
        mutate(date=as.Date(date + dyears()))
    
    gg <- plot_predictions(flu_pred, fit, tback=20, state=s) +
        geom_point(data=ret_pred_dat, col="maroon1", size=1.05)
    
    print(gg)
}

dev.off()

pred_summ <- nat_summ |> 
    select(-mean) |> 
    pivot_longer(-t, names_to="quantile") |> 
    mutate(state="National", quantile=0.01*parse_number(quantile)) |> 
    bind_rows(summarize_predictions(flu_pred, fit))

write_csv(pred_summ, "predictions-10-25.csv")
write_csv(flu_pred, "input-pred-data-10-25.csv")
