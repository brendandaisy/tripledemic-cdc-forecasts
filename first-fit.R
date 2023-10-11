library(tidyverse)
library(INLA)

dat0 <- read_csv("data/counts_state_alldisease.csv") |> 
    rename(t=`...1`)

dat <- dat0 |> 
    pivot_longer(ak_covid:last_col(), values_to="count") |> 
    separate_wider_delim(name, "_", names=c("state", "type")) |> 
    complete(state, type, t) |> 
    mutate(
        epiweek=as.numeric(str_sub_all(orig_index, 5, 6)),
        count=round(count),
        snum=as.numeric(factor(state)),
        dnum=as.numeric(factor(type))
    )

dat_sub <- filter(
    dat,
    state %in% c("va", "nc", "tn", "sc", "ga", "fl", "al", "la", "tx")
) |> 
    mutate(snum=as.numeric(factor(state)), dnum=as.numeric(factor(type))) # have to redo the indices from 1 for INLA

 ggplot(dat, aes(t, count, col=state, group=state)) +
    geom_line(alpha=0.5) +
    geom_point(alpha=0.7) +
    facet_wrap(~type, scales="free_y") +
    theme_bw()

hyper_epwk <- list(theta=list(prior="loggamma", param=c(1, 0.01))) # more favorable to large jumps
# hyper_st <- list(theta=list(prior="pc.prec", param=c(0.5, 0.01)))
hyper_wk <- list(theta=list(prior="loggamma", param=c(1.5, 0.005))) # precision constrained away from 0; P(prec<1)=0.02%, P(prec<100)=20%. NOTE (10/9): varying this prior has virtually no effect even with a subset of data

# mod2 <- count ~ 1 + f(
#     t, model="rw2", group=tnum, cyclic=FALSE, hyper=hyper1, scale.model=TRUE
#     control.group=list(model="rw1", hyper=hyper2)
# )
# 
# fit1 <- inla(
#     mod1, family="poisson", data=filter(dat, type == "covid"),
#     control.compute=list(dic=TRUE),
#     control.predictor=list(link=1, compute=TRUE)
# )

dat_pred <- expand_grid(
    tibble(t=1:5+max(dat_sub$t), epiweek=1:5+last(dat_sub$epiweek)), 
    tibble(type=unique(dat_sub$type), dnum=unique(dat_sub$dnum)), # allows _predictions_ for states without RSV data
    distinct(dat_sub, state, snum)
)

dat_pred <- dat_sub |> 
    bind_rows(dat_pred) |> 
    mutate(sdnum=as.numeric(factor(str_c(snum, dnum))))
 
mod2 <- count ~ 1 + 
    f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, group=dnum, 
      control.group=list(model="iid", hyper=list(theta=list(initial=1, fixed=TRUE)))) +
    f(t, model="ar1", group=sdnum, hyper=hyper_wk, control.group=list(model="exchangeable"))

fit2 <- inla(
    mod2, family="poisson", data=dat_pred,
    control.compute=list(dic=TRUE),
    control.predictor=list(link=1, compute=TRUE)
)

ep_wk_fx <- fit$summary.random$epiweek |> 
    as_tibble()
    # mutate(state=rep(unique(dat_pred$state), each=53)) # assumes states are in correct order in dat_pred

ggplot(ep_wk_fx, aes(ID, mean)) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col="gray70", alpha=0.6) +
    geom_line() +
    # facet_wrap(~state) +
    labs(x="Week of the year", y="Seasonal effect (log scale)") +
    theme_bw() +
    theme(legend.position="none")

dat_pred |> 
    mutate(
        pred_mean=fit2$summary.fitted.values$mean,
        pred_min=fit2$summary.fitted.values$`0.025quant`,
        pred_max=fit2$summary.fitted.values$`0.975quant`
    ) |> 
    ggplot(aes(t, count)) +
    geom_ribbon(aes(ymin=pred_min, ymax=pred_max), col="gray70", alpha=0.6) +
    geom_point(col="gray30", size=0.95, shape=1, alpha=0.9) +
    geom_line(aes(y=pred_mean, col=type)) +
    facet_grid(type~state, scales="free_y") +
    scale_x_continuous(sec.axis=sec_axis(~.x / 53, "years", breaks=1:3)) +
    theme_bw() +
    theme(panel.spacing=unit(0, "mm"))
    
ggsave("figs/fit-triple-1st-attempt.pdf", width=11.5, height=7.5)
