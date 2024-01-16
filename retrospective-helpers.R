library(tidyverse)
library(INLA)
library(mrfDepth)

pred_cdf <- function(marg, E, max_y=E) {
    marg_lam <- inla.tmarginal(\(x) E*x, marg)
    ret <- double(max_y)
    inla.emarginal(\(x) dpois(2, x), marg_lam)
}

rmse <- function(ytrue, ypred) sqrt(mean((ytrue - ypred)^2))

interval_score <- function(ytrue, lo, hi, alpha) {
    penalty <- case_when(
        ytrue < lo ~ 2/alpha * (lo - ytrue),
        ytrue > hi ~ 2/alpha * (ytrue - hi),
        TRUE ~ 0
    )
    (hi - lo) + penalty
}

weighted_interval_score <- function(ytrue, pred_qs) {
    int_scores <- pred_qs |> 
        filter(quantile != 0.5) |> 
        mutate(alpha=round(pmin(quantile, 1-quantile), 3)) |> 
        group_by(alpha) |> 
        summarize(is=interval_score(ytrue, min(value), max(value), alpha[1]))
    
    K <- nrow(int_scores)
    m <- filter(pred_qs, quantile == 0.5)$value[1]
    w <- int_scores$alpha / 2
    return((0.5 * abs(ytrue - m) / (K+0.5) + sum(w * int_scores$is)))
}

log_score <- function(ytrue, pred_marg, E) {
    marg_lam <- inla.tmarginal(\(x) E*x, pred_marg)
    log(inla.emarginal(\(x) dpois(ytrue, x), marg_lam))
}

# Euclidean norm of a real vector
norm_vec <- function(x) sqrt(sum(x^2))

# assumes rows are samples and cols are dimensions in `pred_samps`
# untested
energy_score <- function(ytrue, pred_samps) {
    M <- nrow(pred_samps)
    norm_acc <- double(M)
    norm_spread <- matrix(nrow=M, ncol=M)
    for (i in seq_len(M)) {
        norm_acc[i] <- norm_vec(pred_samps[i,] - ytrue)
        for (j in seq_len(M)) {
            norm_spread[i, j] <- norm_vec(pred_samps[i,] - pred_samps[j,])
        }
    }
    return(mean(norm_acc), 0.5 * mean(norm_spread))
}

###
