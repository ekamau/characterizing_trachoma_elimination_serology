#' ---
#' title: "Country-level posterior probability LOO analysis - N-1 cross validation"
#' output: html_document
#' author: ekamau
#' ---
#' 
# Notes:
# ref - https://stackoverflow.com/questions/67435996/the-probability-that-x-is-larger-than-or-equal-to-a-given-number

#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE, eval = FALSE)

#+ warning = FALSE, message = FALSE, eval = FALSE
# packages
if (!require("pacman")) install.packages("pacman")
p_load(tidyverse, ggridges, reshape2, viridis, here, patchwork, ggtext, foreach)

# Setup parallel computing:
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(n.cores, type = "FORK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() # check if registered
foreach::getDoParWorkers() # how many workers are available?

#------------------------------------------------------------------------------------------------
# Read Bayesian SCR estimates:
here::here()
df <- read.csv(here("output/", "posterior_samples_scr-1to5yo_AllEUs.csv"))

unpublished_EUs <- c("Niger-Ilela-2022", "Niger-Malbaza-2022", "Niger-Bagaroua-2022", "DRC-Manono-2018", "DRC-Nyemba-2018")

eus_wo_PCR <- c("Sudan-El Seraif-2019", "Sudan-Saraf Omrah-2019", "Sudan-Kotom-2019",
                "Malaysia-Sabah-2015",  "Peru-Amazonia-2020", "Papua New Guinea-Mendi-2015", "Papua New Guinea-Daru-2015", "Papua New Guinea-West New Britain-2015")

df2 <- df[!(df$eu %in% c(unpublished_EUs, eus_wo_PCR)), ] %>% 
  group_by(eu) %>% mutate(median_scr = median(scr, na.rm = TRUE)) %>% 
  ungroup() %>% as.data.frame()

# assign EU categories:
df2 <- df2 %>%
  mutate(trachoma_cat = case_when(startsWith(as.character(eu), "Ghana") ~ "Action not needed",
                                  eu %in% c("Ethiopia-Alefa-2017","Niger-MORDOR/Dosso-2015","Ethiopia-Woreta Town-2017",
                                            "Togo-Anie-2017","Togo-Keran-2017","Morocco-Agdaz-2019","Morocco-Boumalne Dades-2019",
                                            "Gambia-River Regions-2014","Ethiopia-Metema-2021","Ethiopia-Woreta Town-2021",
                                            "Malawi-DHO Nkwazi-2014","Malawi-Kasisi/DHO-2014","Malawi-Luzi Kochilira-2014",
                                            "Malawi-Chapananga-2014") ~ "Action not needed",
                                  eu %in% c("Ethiopia-Andabet-2017","Ethiopia-Goncha-2019","Niger-PRET/Matameye-2013",
                                            "Ethiopia-Ebinat-2019","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016",
                                            "Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016","Tanzania-Kongwa-2018",
                                            "Tanzania-Kongwa-2013","Solomon Islands-Temotu/Rennel/Bellona-2015") ~ "Action needed",
                                  TRUE ~ as.character("Unclassified")),
         trachoma_cat = factor(trachoma_cat, levels = c("Action needed","Action not needed","Unclassified"))
  )

#------------------------------------------------------------------------------------------------
# Calculate posterior probability of each category at each value of the SCR
#' Steps:
# (1) prior probability of 'post elimination / action not needed'
# (2) likelihood - use approx to get the interpolated values of the pdf/density at each SCR estimate 
# (3) calculate marginal probabilities
# (4) use Bayes theorem for posterior inference
# (5) loop through the five different priors

#------------------------------------------------------------------------------------------------
#' (1) Prior probabilities:
prior_elim <- c(1/2, 0.65, 0.7, 0.75, 0.8)
prior_endemic <- c()

for(p in prior_elim){
  print(p)
  prior_endemic <- append(prior_endemic, 1-p, after = length(prior_endemic))
}

# confirm the probabilities add to one:
tibble(prior_elim, prior_endemic) %>% 
  rowwise %>% mutate(total = sum(c(prior_elim, prior_endemic))) %>% 
  as.data.frame()


#------------------------------------------------------------------------------------------------
#' Post elimination group:
# add country variable
df2b <- df2 %>% filter(trachoma_cat != "Unclassified") # n = 34
df3 <- separate(df2b, eu, into = c("country", "Extra"), extra = "merge") %>% select(-Extra)

countries <- unique(df3$country) 
prob_df_elim_LOO <- data.frame()

for(e in 1:(length(countries))){
  print(e)
  print(countries[e])
  # subset dataframe by removing the country - leave one out; N-1 remaining countries
  df4 <- subset(df3, country != countries[e])
  
  # probability of SCR in each category => likelihood:
  #1) group data into categories:
  scr_elim <- df4 %>% filter(trachoma_cat == "Action not needed")
  scr_endemic <- df4 %>% filter(trachoma_cat == "Action needed")
  
  #2) create densities in Fig 2B in manuscript:
  dx_elim <- density(scr_elim$scr, n = 2^12, from = min(scr_elim$scr), to = max(scr_elim$scr))
  dx_endemic <- density(scr_endemic$scr, n = 2^12, from = min(scr_endemic$scr), to = max(scr_endemic$scr))
  
  #3) calculate likelihood and posterior probability:
  ## below computes posterior probability for the post-elimination group:
  A <- foreach(i=prior_elim, k=prior_endemic) %:%
    foreach(c=scr_elim$scr, .combine = 'rbind', .multicombine = TRUE) %dopar% {
      prob_c_elim <- approx(dx_elim$x, dx_elim$y, xout = c)$y
      prob_c_endemic <- approx(dx_endemic$x, dx_endemic$y, xout = c)$y
      
      numerator <- i * prob_c_elim
      denom <- sum((i * prob_c_elim), (k * prob_c_endemic), na.rm = TRUE)
      
      prob_elim_c <- numerator / denom
      data.frame(e, countries[e], i, c, prob_elim_c)
    }
  
  prob_elim_df <- do.call(rbind, A)
  prob_df_elim_LOO <- rbind(prob_df_elim_LOO, prob_elim_df)
  
}

colnames(prob_df_elim_LOO) <- c("LOO_set", "country_left", "prior_set", "c", "prob")

write.csv(prob_df_elim_LOO, "output/probs_scr_elimination-country-N1-2categories.csv", row.names = FALSE)
saveRDS(prob_df_elim_LOO, "output/probs_scr_elimination-country-N1-2categories.rds")

