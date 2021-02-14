# PACKAGES <- c("devtools", "dplyr", "ggplot2", "stringr", "lubridate", "lme4")
# install.packages(PACKAGES)
# devtools::install_github("BillPetti/baseballr")

library(dplyr)
library(ggplot2)
library(lme4)
library(splines)
library(baseballr)

set.seed(2020)
PA = 600 # presumptive PA in full season
Q = .975
NSAMP = 10000
BAM_LU_FILEPATH = "~/Downloads/master.csv"

get_exit_data <- function(min_season=2016, max_season=2020){
  #' Wrapper around Bill Petti's package to get tracking data over course of season(s)
  #' If a particular date throws an error within baseballr, the loop simply moves onto the next date
  #' @param min_season : int. Minimum season to scrape
  #' @param max_season : int. Maximum season to scrape
  #' @return data.frame : a dataframe of tracking data
  
  result <- lapply(min_season:max_season, function(s){
    lapply(seq(as.Date(paste0(s, "-03-15")), as.Date(paste0(s, "-10-01")), 1), function(d)
      try(baseballr::scrape_statcast_savant(start_date=d, end_date=d), silent=T)
    ) %>%
      do.call("rbind", .)
  }) %>%
    do.call("rbind", .) %>%
    mutate(season = lubridate::year(game_date))
}

collapse_hr_counts <- function(exit_data_df){
  #' Computes hr, pa, and hr/pa for a player in a season
  #' @param exit_data_df : data.frame. A dataframe output of `get_exit_data()`
  #' @return : data.frame. The collapsed dataframe. 
  exit_data_df %>%
    mutate(season = lubridate::year(game_date)) %>%
    filter(type == 'X') %>%
    group_by(batter, season) %>%
    summarise(pa = n(), hr = sum(events == 'home_run'), hr_pct = mean(events == 'home_run'))%>%
    rename(mlb_id = batter)
}

collapse_exit_data <- function(exit_data_df, n=50){
  #' Extracts each player's max-ev in a given season, and filters players with insufficent tracked batted balls
  #' @param exit_data_df : data.frame. A dataframe of raw tracking data, as retrieved by `get_exit_data()`
  #' @param n : int. Threshold for number of tracked balls
  #' @return : data.frame. A dataframe of the aforementioned max EVs, according to threshold n. 
  
  result <- exit_data_df %>% 
    mutate(season = lubridate::year(game_date)) %>%
    filter(!is.na(launch_speed), 
           launch_speed < 130) %>% # somehow someone's bam_id got listed as an EV???
    group_by(batter, season) %>%
    summarise(ev_max = max(launch_speed, na.rm=T),
              ev_mean = mean(launch_speed, na.rm=T),
              ev_95 = quantile(launch_speed, .95, na.rm=T),
              ev_75 = quantile(launch_speed, .75, na.rm=T),
              ev_50 = quantile(launch_speed, .5, na.rm=T),
              ev_25 = quantile(launch_speed, .25, na.rm=T),
              ev_05 = quantile(launch_speed, .05, na.rm=T),
              n_tracked = length(launch_speed)) %>%
    filter(n_tracked >= n) %>%
    rename(mlb_id = batter)
  
  result
  
}

find_frontier_piecewise <- function(df, yvar='hr'){
  #' Finds a homerun frontier for a given exit velocity via a coarse search
  #' @param df : data.frame. A dataframe of exit velocities and frontier targets. Must have cols `ev_max` and yvar.
  #' @return : data.frame. A dataframe defining the piecewise frontier

  # go through each row:
  # for that row, filter s.t.:
  #     - max exit velo is less than that corresponding to row
  #     - hr% is higher than that corresponding to row
  # If there's nothing there, you're on the "boundary"
  y_boundary_bool = sapply(1:nrow(df), function(i){
    boundary_check_df = df[(df$ev_max < df$ev_max[i]) & (pull(df, eval(yvar)) > as.numeric(df[i, yvar])), ]
    ifelse(nrow(boundary_check_df), F, T)
  })
  
  x_boundary = df$ev_max[y_boundary_bool]
  y_boundary = pull(df[y_boundary_bool, eval(yvar)])
  boundary_df = data.frame(x_boundary, y_boundary) %>% arrange(x_boundary)
  boundary_df = bind_rows(boundary_df, data.frame(x_boundary = max(df$ev_max), y_boundary = max(y_boundary)))
  
  boundary_df
}


plot_frontiers_piecewise <- function(master_df){
  #' Plots piecewise frontiers for a given df that contains: (1) max exit velos, (2) hr season totals, and (3) pa season totals
  #' @param master_df : data.frame. A dataframe with exit velos and hr/pa counts. Cols are `ev_max`, `hr`, `hr_pct`
  #' @return list[gglot, ggplot2], two plots with the piecewise frontiers
  
  # piecewise frontiers
  piecewise_hr_pct = find_frontier_piecewise(master_df, "hr_pct")
  piecewise_hr = find_frontier_piecewise(master_df, "hr")
  
  # smoothed -- no tuning
  loc_lm_hr_pct = loess(y_boundary ~ x_boundary, piecewise_hr_pct, span=.5)
  loc_lm_hr = loess(y_boundary ~ x_boundary, piecewise_hr, span=.5)
  smooth_frontier = data.frame(x_boundary = seq(100, max(piecewise_hr$x_boundary), .01))
  smooth_frontier$yhat_hr_pct = predict(loc_lm_hr_pct, smooth_frontier)
  smooth_frontier$yhat_hr = predict(loc_lm_hr, smooth_frontier)
  
  plt_hr_pct = ggplot(master_df, aes(x = ev_max, y = hr_pct, color = pa)) +
    geom_point() +
    # stat_smooth() +
    geom_path(piecewise_hr_pct, mapping = aes(x = x_boundary, y = y_boundary), color='red', size=1) + 
    geom_path(smooth_frontier, mapping = aes(x = x_boundary, y = yhat_hr_pct), color='purple', size=1)
  
  plt_hr = ggplot(master_df, aes(x = ev_max, y = hr, color = pa)) +
    geom_point() +
    # stat_smooth() +
    geom_path(piecewise_hr, mapping = aes(x = x_boundary, y = y_boundary), color='red', size=1) + 
    geom_path(smooth_frontier, mapping = aes(x = x_boundary, y = yhat_hr), color='purple', size=1)
  
  list(plt_hr_pct, plt_hr)
}


compare_models <- function(df, seas_cut=2019){
  #' Performs a quick holdout CV (no tuning/regularization) due to simplicity of the model, to compare binomial vs. Poisson. 
  #' Function is largely a sanity check. 
  #' @param df : data.frame. A dataframe with train and devset, s.t. one can call glm() on this df
  #' @param seas_cut : int. Season to cut dev set on. Cut is \leq. 
  #' @return : numeric[2]. The -NLL of the binomial and Poisson CV fits, respectively
  #' # TODO: add in bootstrap functionality, since distribution of X is probably whack. 
  holdout_result = list()
  grd = expand.grid(df_ev=1:20, df_seas=1:3)
  for (ii in 1:nrow(grd)){
    df_train = df %>% filter(season <= seas_cut)
    df_test = df %>% filter(season > seas_cut)
    # fit a binomial and a poisson. TODO: neg-binomial???
    glm_binom = glm(cbind(hr, pa-hr) ~ ns(ev_max, df=grd[ii, 1]) + season, df_train, family='binomial')
    glm_pois = glm(hr ~ ns(ev_max, df=grd[ii, 1]) + season + offset(log(pa)), df_train, family='poisson')
    
    pi_hat_binom = predict(glm_binom, newdata=df_test %>% mutate(season=seas_cut), type='response')
    pi_hat_pois = predict(glm_pois, newdata=df_test %>% mutate(season=seas_cut), type='response')/df_test$pa
    
    nll_binom = -sum(df_test$hr * log(pi_hat_binom) - (df_test$pa - df_test$hr) * log(1 - pi_hat_binom))/sum(df_test$pa);
    nll_pois = -sum(df_test$hr * log(pi_hat_pois) - (df_test$pa - df_test$hr) * log(1 - pi_hat_pois))/sum(df_test$pa);
    
    holdout_result[[ii]] = c(nll_binom, nll_pois)
  }
  do.call("rbind", holdout_result) %>%
    cbind(., grd)
}


plot_hr_ceiling <- function(mod, PA=600, season=2019){
  #' Plots mean, as well as 95, and 99 quantiles 
  #' @param mod : glm.lm, a fitted model in the GLM family
  #' @param PA: int. Presumed PA
  #' @param season : int. Season parameter for model, e.g. 2019
  #' @return : ggplot. A plot of the aforementioned
  
  # init features
  xhat = data.frame(ev_max = seq(90, 120, .01), season=season, pa=PA)
  if (mod$call$family == 'poisson'){
    xhat$pi = predict(mod, newdata=xhat, type='response')
    xhat$pi975 = qpois(Q, predict(mod, xhat, type='response'))
  }else if(mod$call$family == 'binomial'){
    xhat$pi = predict(mod, newdata=xhat, type='response')*PA
    xhat$pi975 = qbinom(Q, size=PA, prob=predict(mod, xhat, type='response'))
  }
  
  plt = ggplot(xhat, aes(x = ev_max, y = pi)) + 
    geom_path(size=1) +
    geom_path(mapping=aes(x = ev_max, y = pi975)) + 
    labs(x='Max Exit Velo.', y='HR', title='Max EV vs. Modeled HR: Mean and 97.5 Pctile', subtitle=paste0('Run Env: ', season))
  plt
}

main <- function(){
  
  # Ingest --------------------------------------------------------------------
  # fetch exit data
  exit_data = get_exit_data()
  # compute exit_velo maxes
  ev_maxes = collapse_exit_data(exit_data_df=exit_data, n=100)
  # compute hrs
  hr_counts = collapse_hr_counts(exit_data_df=exit_data)
  # ingest id lookup
  bam_lu_df <- read.csv(BAM_LU_FILEPATH)
  # combine into single dataframe
  master_df <- ev_maxes %>%
    left_join(., bam_lu_df %>% select(mlb_id, mlb_name),
              by=c("mlb_id")) %>% 
    left_join(., hr_counts, by=c("mlb_id", "season")) %>%
    ungroup()

  # Train + Test ------------------------------------------------------------
  # quick comparison of Binomial vs. Poisson: train set is 2016-17; dev set is 2018
  nll_comp = compare_models(master_df, 2018)
  cat('Binomial Best tune:\n\t')
  print(nll_comp[which.min(nll_comp[, 1]), ])
  cat('Poisson Best tune:\n\t')
  print(nll_comp[which.min(nll_comp[, 2]), ])
  
  # train + test: train on 2016-2019, test on 2020
  df_train = master_df %>% filter(season <= 2019); df_test = master_df %>% filter(season == 2020)
  mod_test = glm(hr ~ ns(ev_max, 2) + season + offset(log(pa)), df_train, family='poisson') # see tune from above
  plot_hr_ceiling(mod_test, PA=600, season=2019)
  
  # again -- this is just an approximate prediction interval. NOT RIGOROUS. 
  # more legitimate PI likely requires boostrapping or some Bayesian action, i.e. sample --> fit --> sample theta --> sample/quantile
  # nevertheless, approximate capture rate @ 97.5
  
  # predictions in link space
  yhat_t = as.data.frame(predict(mod_test, newdata=df_test, type='link', se.fit=T))
  # method 1: eval Poison CDF @ .975
  # - 10000 draws of g(X *theta), using asymptotic normal approximation
  # - unlink
  # - 10000 Poisson draws of those
  # - Compute quantile
  yhat_pois_upr = sapply(1:nrow(yhat_t), function(i) 
    rpois(NSAMP, exp(rnorm(NSAMP, yhat_t$fit[i], yhat_t$se.fit[i]))) %>%
      quantile(., Q) %>% 
      as.numeric()
  )
  
  # method 2: Treat expected arrival time in one PA as HR probability, hence a Bernoulli
  # - 10000 draws of g(X *theta), using asymptotic normal approximation
  # - unlink
  # - 10000 *Binomial* draws of those, corresponding to observed PA
  # - Compute quantile
  yhat_binom_upr = sapply(1:nrow(yhat_t), function(i) 
    rbinom(NSAMP, df_test$pa[i], exp(rnorm(NSAMP, yhat_t$fit[i], yhat_t$se.fit[i]))/df_test$pa[i]) %>%
      quantile(., Q) %>%
      as.numeric()
  ) 
  
  # actual (mean) predictions
  yhat = predict(mod_test, newdata=df_test, type='response')
  pi_hat = yhat/df_test$pa
  # store sampled quantiles
  df_test$q_pois = yhat_pois_upr; df_test$q_binom = yhat_binom_upr
  # method 1 capture rate:
  cat('Method 1 (direct Poisson CDF) capture rate: \n\t-', mean(df_test$q_pois > df_test$hr))
  cat('Method 2 (~ binomial) capture rate: \n\t-', mean(df_test$q_binom > df_test$hr))
  
  # Full Model --------------------------------------------------------------
  mod = glm(hr ~ ns(ev_max, 2) + season + offset(log(pa)), master_df, family='poisson')
  plot_hr_ceiling(mod)
  
}
