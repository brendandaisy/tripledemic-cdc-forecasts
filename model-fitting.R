library(tidyverse)
library(INLA)
library(lubridate)

inla.setOption(inla.mode="classic")

load_us_graph <- function(df, f="data/us-state-boundaries.shp") {
    us0 <- read_sf(f)
    
    us <- us0 |> 
        filter(name %in% unique(df$state)) |> 
        select(state=name, division, region)
    
    state_order <- fct_inorder(unique(df$state))
    
    # sort order of states to match their order of appearance in flu data
    us |> 
        mutate(state=fct_relevel(state, levels(state_order))) |> 
        arrange(state)
}

us_adj_mat <- function(us) {
    us_adj <- us |> 
        poly2nb() |> 
        nb2mat(style="W", zero.policy=TRUE)
    
    us_adj[which(us$state == "Florida"), which(us$state == "Puerto Rico")] <- 1
    us_adj[which(us$state == "Puerto Rico"), which(us$state == "Florida")] <- 1
    # us_adj[which(us$state == "Georgia"), which(us$state == "Puerto Rico")] <- 1
    # us_adj[which(us$state == "Puerto Rico"), which(us$state == "Georgia")] <- 1
    us_adj[which(us$state == "Washington"), which(us$state == "Alaska")] <- 1
    us_adj[which(us$state == "Alaska"), which(us$state == "Washington")] <- 1
    # us_adj[which(us$state == "California"), which(us$state == "Alaska")] <- 1
    # us_adj[which(us$state == "Alaska"), which(us$state == "California")] <- 1
    us_adj[which(us$state == "California"), which(us$state == "Hawaii")] <- 1
    us_adj[which(us$state == "Hawaii"), which(us$state == "California")] <- 1
    return(us_adj)
}

prep_pred_data <- function(df, tpred=5) {
    # make a dataframe to hold group info for forecasting
    pred_df <- expand_grid(
        tibble(
            date=duration(1:tpred, "week") + max(df$date),
            t=1:tpred + max(df$t),
            epiweek=epiweek(date)
        ),
        distinct(df, state, snum, division, region)
    )
    
    bind_rows(df, pred_df) |> # add to data for counts to be NAs
        arrange(t)
}

fit_current_model <- function(df, q=c(0.025, 0.5, 0.975), graph=NULL, joint=FALSE) {
    # the PC priors c(u, a) give the probability that the standard deviation between weeks u is greater than a
    # increasing u increases prior beliefs that there will be large jumps in data between weeks; see `inla.doc("pc.prec")`
    hyper_epwk <- list(prec=list(prior="pc.prec", param=c(0.5, 0.01))) # a priori the seasonal effect should be more smooth
    hyper_wk <- list(theta1=list(prior="pc.prec", param=c(1, 0.01))) # more favorable to large jumps
    
    if (is.null(graph)) {
        sp_control <- list(model="exchangeable")
        mod <- count ~ 1 + # intercept
            f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk) + # seasonal effect (currently shared over states)
            f(t, model="ar1", group=snum, hyper=hyper_wk, control.group=sp_control) # weekly x state effect
    } else {
        mod <- count ~ 1 +
            f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) +
            f(snum, model="besag", graph=graph, hyper=hyper_wk, scale.model=TRUE, group=t, control.group=list(model="ar1"))
    }
    
    pred_idx <- which(is.na(df$count))
    fit <- inla(
        mod, family="poisson", data=df, # poisson regression link
        E=df$ex_count,
        quantiles=q,
        selection=if (!joint) NULL else list(Predictor=pred_idx),
        # lincomb=inla.make.lincombs(Predictor=lc_nat),
        control.compute=list(dic=FALSE, mlik=FALSE, return.marginals.predictor=TRUE),
        control.predictor=list(link=1, compute=TRUE) # compute quantiles for NA dates (i.e. do the forecasting)
    )
    return(fit)
}

summarize_predictions <- function(df, fit) {
    
    map2(fit$marginals.fitted.values, flu_pred$ex_count, \(m, ex) {
        msamp <- inla.rmarginal(1000, m)
        
    })
    
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

plot_predictions <- function(df, fit, state=unique(df$state), tback=10, tspan=NULL, file=NULL) {
    pred_summ <- fit$summary.fitted.values |> 
        as_tibble() |> 
        select(mean, contains("quant"))
    
    tpred <- min(filter(df, is.na(count))$t)
    
    df_plt <- df |> 
        bind_cols(pred_summ) |> 
        filter(
            if(is.null(tspan)) t >= (tpred - tback) else date %within% tspan,
            state %in% !!state
        )
    
    gg <- ggplot(df_plt, aes(date, count)) +
        geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), fill="gray60", alpha=0.6) +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), fill="gray80", alpha=0.6) +
        geom_point(col="steelblue", size=1.05) +
        geom_point(aes(y=count_true), filter(df_plt, is.na(count)), col="maroon1", size=1.05) +
        geom_line(aes(y=mean), col="tomato3") +
        facet_wrap(~state, scales="free_y", nrow=4) +
        scale_x_date(date_breaks="4 weeks", date_labels="%d %b %y", guide=guide_axis(angle=45)) +
        labs(x="Weeks", y="Hospitalizations") +
        theme_bw() +
        theme(panel.spacing=unit(0, "mm"), axis.title=element_text(size=14))
    
    if (!is.null(file))
        ggsave(file, gg, width=14, height=8)
    gg
}

# assumes the most recent week of data is in first row, for current epiweek
plot_seasonal <- function(df, fit) {
    ep_wk_fx <- fit$summary.random$epiweek |> 
        as_tibble()
    # mutate(state=rep(unique(dat_pred$state), each=53)) # assumes states are in correct order in dat_pred
    
    ggplot(ep_wk_fx, aes(ID, mean)) +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col="gray70", alpha=0.6) +
        geom_line(col="tomato3") +
        geom_vline(xintercept=first(df$epiweek), col="steelblue", linetype="dashed") +
        # facet_wrap(~state) +
        labs(x="Week of the year", y="Seasonal effect (log scale)") +
        theme_bw() +
        theme(legend.position="none")
}
