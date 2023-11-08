library(tidyverse)
library(INLA)
library(lubridate)

inla.setOption(inla.mode="classic")

load_us_graph <- function(flu, f="data/us-state-boundaries.shp") {
    us0 <- read_sf(f)
    
    us <- us0 |> 
        filter(name %in% unique(flu$location)) |> 
        select(state=name, division, region)
    
    state_order <- fct_inorder(unique(flu$location))
    
    # sort order of states to match their order of appearance in flu data
    us |> 
        mutate(state=factor(state, levels(state_order))) |> 
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

prep_fit_data <- function(flu, weeks_ahead=4) {
    ret <- flu |> 
        filter(location != "US") |> # make sure US is not in training data
        group_by(date) |> 
        mutate(t=cur_group_id(), epiweek=epiweek(date), .after=date) |> # add a time counter starting from 1 for earliest week
        ungroup() |> 
        mutate(
            snum=as.numeric(fct_inorder(location)), # INLA needs groups as ints starting from 1, so add numeric state code
            ex_lam=population
        )
    
    # flu <- flu |> 
    #     group_by(year, epiweek) |> 
    #     mutate(tot_count=sum(count)) |> 
    #     ungroup(year) |> 
    #     mutate(ex_lam=(population / pop_nat) * mean(tot_count), .after=count) |> 
    #     ungroup() |> 
    #     select(-tot_count)
    
    # make a dataframe to hold group info for forecasting
    pred_df <- expand_grid(
        tibble(
            date=duration(1:weeks_ahead, "week") + max(ret$date),
            t=1:weeks_ahead + max(ret$t),
            epiweek=epiweek(date)
        ),
        distinct(ret, location, snum, population) # makes pairs of new times X each state
    ) |> 
        left_join(distinct(ret, location, epiweek, ex_lam)) # go and find `ex_lam` values for each state and epiweek
    
    bind_rows(ret, pred_df) |> # add to data for counts to be NAs
        arrange(t)
}

fit_current_model <- function(fit_df, forecast_date, q=c(0.025, 0.5, 0.975), graph=NULL, joint=TRUE) {
    # the PC priors c(u, a) give the probability that the standard deviation between weeks u is greater than a
    # increasing u increases prior beliefs that there will be large jumps in data between weeks; see `inla.doc("pc.prec")`
    hyper_epwk <- list(prec=list(prior="pc.prec", param=c(0.5, 0.01))) # a priori the seasonal effect should be more smooth
    hyper_wk <- list(theta1=list(prior="pc.prec", param=c(1, 0.01))) # more favorable to large jumps
    
    if (is.null(graph)) {
        sp_control <- list(model="exchangeable")
        mod <- count ~ 1 + # intercept
            f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) + # seasonal effect (currently shared over states)
            f(t, model="ar1", group=snum, hyper=hyper_wk, control.group=sp_control) # weekly x state effect
    } else {
        mod <- count ~ 1 +
            f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) +
            f(snum, model="besagproper", graph=graph, hyper=hyper_wk, group=t, control.group=list(model="ar1"))
    }
    
    pred_idx <- which(fit_df$date >= forecast_date)
    fit <- inla(
        mod, family="poisson", data=fit_df, # poisson regression link
        E=fit_df$ex_lam,
        quantiles=q,
        selection=if (!joint) NULL else list(Predictor=pred_idx),
        # lincomb=inla.make.lincombs(Predictor=lc_nat),
        control.compute=list(dic=FALSE, mlik=FALSE, return.marginals.predictor=TRUE),
        control.predictor=list(link=1) # compute quantiles for NA dates (i.e. do the forecasting)
    )
    return(fit)
}

# plot_predictions <- function(cound_pred, state=unique(count_pred$state), tback=10, tspan=NULL, file=NULL) {
#     # pred_summ <- fit$summary.fitted.values |> 
#     #     as_tibble() |> 
#     #     select(mean, contains("quant"))
#     
#     tpred <- min(filter(count_pred, is.na(count))$t)
#     
#     df_plt <- filter(
#             count_pred,
#             if(is.null(tspan)) t >= (tpred - tback) else date %within% tspan,
#             state %in% !!state
#         )
#     
#     gg <- ggplot(df_plt, aes(date, count)) +
#         geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="gray60", alpha=0.6) +
#         geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="gray80", alpha=0.6) +
#         geom_point(col="steelblue", size=1.05) +
#         # geom_point(aes(y=count_true), filter(df_plt, is.na(count)), col="maroon1", size=1.05) +
#         geom_line(aes(y=mean), col="tomato3") +
#         facet_wrap(~state, scales="free_y", nrow=4) +
#         scale_x_date(date_breaks="4 weeks", date_labels="%d %b %y", guide=guide_axis(angle=45)) +
#         labs(x="Weeks", y="Hospitalizations") +
#         theme_bw() +
#         theme(panel.spacing=unit(0, "mm"), axis.title=element_text(size=14))
#     
#     if (!is.null(file))
#         ggsave(file, gg, width=14, height=8)
#     gg
# }


