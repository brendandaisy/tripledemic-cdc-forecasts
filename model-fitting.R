library(tidyverse)
library(INLA)

prep_model_data <- function(dat_long) {
    # INLA needs groups as ints starting from 1, so add numeric state code
    dat_long <- mutate(dat_long, snum=as.numeric(factor(state)))
    
    # make a dataframe to hold group info for forecasting
    pred_df <- expand_grid(
        tibble(t=1:5+max(dat_long$t), epiweek=(0:4+last(dat_long$epiweek)) %% 53 + 1),
        distinct(dat_long, state, snum)
    )
    
    bind_rows(dat_long, pred_df) # add to data for counts to be NAs
}

fit_current_model <- function(df, q=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    # the PC priors c(u, a) give the probability that the standard deviation between weeks u is greater than a
    # increasing u increases prior beliefs that there will be large jumps in data between weeks; see `inla.doc("pc.prec")`
    hyper_epwk <- list(theta=list(prior="pc.prec", param=c(0.5, 0.01))) # a priori the seasonal effect should be more smooth
    hyper_wk <- list(theta=list(prior="pc.prec", param=c(1, 0.01))) # more favorable to large jumps
    
    mod <- count ~ 1 + # intercept
        f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk) + # seasonal effect (currently shared over states)
        f(t, model="ar1", group=snum, hyper=hyper_wk, control.group=list(model="exchangeable")) # weekly x state effect
    
    fit <- inla(
        mod, family="poisson", data=df, # poisson regression link
        quantiles=q,
        control.compute=list(dic=FALSE, mlik=FALSE), # don't compute model comparison metrics
        control.predictor=list(link=1, compute=TRUE) # compute quantiles for NA dates (i.e. do the forecasting)
    )
    return(fit)
}

summarize_predictions <- function(df, fit) {
    quantiles <- fit$summary.fitted.values |> 
        as_tibble() |> 
        select(contains("quant"))
    
    df |> 
        bind_cols(quantiles) |> 
        filter(is.na(count)) |> 
        pivot_longer(contains("quant"), names_to="quantile") |> 
        mutate(quantile=parse_number(quantile)) |> 
        select(t, state, quantile, value)
}

plot_predictions <- function(df, fit, tback=10, q=c(0.025, 0.975), tspan=NULL, file=NULL) {
    pred_summ <- fit$summary.fitted.values |> 
        as_tibble() |> 
        select(mean, contains(as.character(q)))
    
    colnames(pred_summ) <- c("mean", "ymin", "ymax")
    tpred <- min(filter(df, is.na(count))$t)
    
    gg <- df |> 
        bind_cols(pred_summ) |> 
        filter(if(is.null(tspan)) t >= tpred - tback else date %within% tspan) |> 
        ggplot(aes(date, count)) +
        geom_ribbon(aes(ymin=ymin, ymax=ymax), col="gray70", alpha=0.6) +
        geom_point(col="steelblue", size=1.05, shape=1) +
        geom_line(aes(y=mean), col="tomato3") +
        facet_wrap(~state, scales="free_y", nrow=3) +
        scale_x_date(date_breaks="4 weeks", date_labels="%b %Y", guide=guide_axis(angle=45)) +
        labs(x="Weeks", y="Hospitalizations") +
        theme_bw() +
        theme(panel.spacing=unit(0, "mm"), axis.title=element_text(size=14))
    
    if (!is.null(file))
        ggsave(file, gg, width=15, height=5)
    gg
}


plot_seasonal <- function(df, fit) {
    ep_wk_fx <- fit$summary.random$epiweek |> 
        as_tibble()
    # mutate(state=rep(unique(dat_pred$state), each=53)) # assumes states are in correct order in dat_pred
    
    ggplot(ep_wk_fx, aes(ID, mean)) +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col="gray70", alpha=0.6) +
        geom_line(col="tomato3") +
        geom_vline(xintercept=last(df$epiweek), col="steelblue", linetype="dashed") +
        # facet_wrap(~state) +
        labs(x="Week of the year", y="Seasonal effect (log scale)") +
        theme_bw() +
        theme(legend.position="none")
}
