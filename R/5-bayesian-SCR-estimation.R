#' ---
#' title: "Bayesian SCR inference using random effects"
#' author: ekamau
#' output: html_document
#' ---
#' 

# Notes: Age 1 to 5yo

#+ setup, include=FALSE, echo=FALSE
knitr::opts_chunk$set(collapse = TRUE, echo = FALSE)

# packages
if (!require("pacman")) install.packages("pacman")
p_load(rstan, tidyverse, ggridges, reshape2, viridis, posterior, foreach, bayesplot, 
       patchwork, ggtext, lmtest, sandwich, tidybayes, rstanarm)

# rstan:
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE, threads_per_chain = 1)

# Read data:
here::here()
ind_df <- readRDS(here::here("../../Data/public-devel/", "trachoma_serology_public_data_indiv_devel.rds"))
eus <- unique(ind_df$eu_name)

ind_df2 <- ind_df %>% filter(age_years <= 5) %>% # n = 45073/70681
  rename(cluster_id = cluster_id_public,
         household_id = household_id_public,
         individual_id = individual_id_public) %>% 
  mutate(cluster_id = replace(cluster_id, is.na(cluster_id), 1)) 

eus <- unique(ind_df2$eu_name)

#-------------------------------------------------------------------------------------------------
#' Estimate SCR - fit model
scr_data <- data.frame()
post_scr_df <- list()

scr_1to5 <- foreach(i = eus, .combine = rbind) %do% {
  print(i)
  di <- ind_df2 %>% dplyr::filter(eu_name == i, !is.na(pgp3_pos)) %>%
    group_by(cluster_id, age_years) %>%
    summarize(npos = sum(pgp3_pos), nchild = n(), .groups = "drop",
              prop_seropos = npos / nchild) %>% 
    transform(clusID = as.numeric(factor(cluster_id)))
  
  di <- di[complete.cases(di), ]
  di$age_years <- replace(di$age_years, di$age_years == 0, 1)
  
  stan_data <- list(N = nrow(di), 
                    m = nlevels(as.factor(di$clusID)),
                    clusid = as.numeric(di$clusID),
                    age = di$age_years,
                    z = di$npos, 
                    n = di$nchild)
  
  fit <- rstan::stan(file = here::here("stan/", "constFOI_randomEffects.stan"), data = stan_data, 
                     chains = 4, init = 'random', iter = 3000, warmup = 900, 
                     refresh = 0, seed = 345)
    
  # extract posterior samples
  lambda_draws <- rstan::extract(fit, pars = 'lambda', permuted = TRUE)$lambda
  post_scr_df[[i]] <- lambda_draws
    
  lambda.summary <- as.data.frame(rstan::summary(fit, pars = c('lambda'), 
                                                 probs = c(0.025, 0.50, 0.975))$summary) %>%
    mutate(eu = i)
    
  scr_data <- rbind(scr_data, lambda.summary)
  return(scr_data)
  
}

#' Save posterior samples:
write.csv(scr_data, file = here::here("output/", "bayes_SCR_estimates_summary_1to5yo.csv"), row.names = FALSE)
          
big_scr_df <- dplyr::bind_cols(post_scr_df)
big_scr_df_long <- melt(big_scr_df, variable.name = "eu", value.name = "scr")
write.csv(big_scr_df_long, file = here::here("output/", "posterior_samples_scr-1to5yo.csv"), row.names = FALSE)
saveRDS(big_scr_df_long, file = here::here("output/", "posterior_samples_scr-1to5yo.rds"))

