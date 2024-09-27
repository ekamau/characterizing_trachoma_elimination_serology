#' ---
#' title: "Bayesian seroprevalence inference using random effects"
#' author: ekamau
#' output: html_document
#' ---
#' 

# Notes: Age 1 to 5yo
#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE)

# packages
if (!require("pacman")) install.packages("pacman")
p_load(rstan, tidyverse, ggridges, reshape2, viridis, posterior, foreach, bayesplot, 
       patchwork, ggtext, lmtest, sandwich, tidybayes, rstanarm)

# rstan:
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE, threads_per_chain = 1)

# Read data:
here::here()
ind_df <- readRDS("trachoma_serology_public_data_indiv_devel.rds")
eus <- unique(ind_df$eu_name)

ind_df2 <- ind_df %>% filter(age_years <= 5) %>%
  rename(cluster_id = cluster_id_public,
         household_id = household_id_public,
         individual_id = individual_id_public) %>% 
  mutate(cluster_id = replace(cluster_id, is.na(cluster_id), 1)) 

eus <- unique(ind_df2$eu_name)

#' Calculate seroprevalence
#' Method 1 - using glm
# ref - https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html

seroprev_1to9 <- foreach(eui = eus, .combine = rbind) %do% {
  di <- ind_df2 %>% filter(eu_name == eui, !is.na(pgp3_pos))
  glmfit    <- glm(pgp3_pos ~ 1, data = di, family = binomial(link = "logit"))
  glmfit_rb <- coeftest(glmfit, vcov. = vcovCL(glmfit, cluster = di$cluster_id))
  # transform the log-odds into probability (proportion)
  serop       <- plogis(glmfit_rb[1,1])
  serop_min95 <- plogis(glmfit_rb[1,1] - 1.96*glmfit_rb[1,2])
  serop_max95 <- plogis(glmfit_rb[1,1] + 1.96*glmfit_rb[1,2])
  
  res <- di %>%
    group_by(eu_name) %>%
    summarize(npos = sum(pgp3_pos), nchild = n(), .groups = "drop") %>%
    mutate(ncluster = length(unique(di$cluster_id)),
           logodds = glmfit_rb[1,1],
           logodds_se = glmfit_rb[1,2],
           seroprev = serop,
           seroprev_min95 = serop_min95,
           seroprev_max95 = serop_max95)
  return(res)
  
}

class(seroprev_1to9)
write.csv(as.data.frame(seroprev_1to9), "", row.names = FALSE)

#--------------------------------------------------------------------------------------------
#' Method 2 - fit model using rstanarm
# ref - https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-rstanarm.html
# ref - http://mc-stan.org/rstanarm/articles/binomial.html

#+ warning = FALSE, message = FALSE
post_prev_df <- list()
prev_data <- data.frame()

seroprev_df <- foreach(i = eus, .combine = rbind) %do% {
  print(i)
  di <- ind_df2 %>% dplyr::filter(eu_name == i, !is.na(pgp3_pos)) %>%
    group_by(cluster_id, age_years)
  
  m2 <- rstanarm::stan_glmer(pgp3_pos ~ 1 + (1 | cluster_id), data = di, family = gaussian(link = "identity"),
                             chains = 4, init = 'random', refresh = 0, iter = 3000, warmup = 900)

  # extract all posterior draws
  post_prev_df[[i]] <- as.data.frame(m2 %>% spread_draws(`(Intercept)`))$`(Intercept)`
    
  # get 95% quantile intervals or highest (posterior) density interval
  prev.summary <- (as.data.frame(m2 %>% gather_draws(`(Intercept)`) %>% median_hdi()) %>%
                     mutate(eu = i))[-1]
    
  prev_data <- rbind(prev_data, prev.summary)
  return(prev_data)
    
}

#' Write to files:
colnames(prev_data) <- c("median", "lower", "upper", "width", "point", "interval", "eu")
write.csv(prev_data, file = here::here("output/", "bayesian_seroprev_posterior_summary-1to5yo.csv"), row.names = FALSE)

big_prev_df <- dplyr::bind_cols(post_prev_df)
big_prev_df_long <- melt(big_prev_df, variable.name = "eu", value.name = "prev")

write.csv(big_prev_df_long, file = here::here("output/", "posterior_seroprev_samples_1to5yo.csv"), row.names = FALSE)
saveRDS(big_prev_df_long, file = here::here("output/", "posterior_seroprev_samples_1to5yo.rds"))
