#' ---
#' title: "Determine thresholds from empirical SCR estimates - 2 categories"
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

# Read Bayesian SCR estimates:
here::here()
df <- read.csv(here("output/", "posterior_samples_scr-1to5yo_AllEUs.csv"))
eus <- unique(df$eu) # check eu labels!

unpublished_EUs <- c("Niger-Ilela-2022", "Niger-Malbaza-2022", "Niger-Bagaroua-2022", "DRC-Manono-2018", "DRC-Nyemba-2018")

left_out_eus <- c("Sudan-El Seraif-2019", "Sudan-Saraf Omrah-2019", "Sudan-Kotom-2019",
                  "Malaysia-Sabah-2015",  "Peru-Amazonia-2020", "Papua New Guinea-Mendi-2015", "Papua New Guinea-Daru-2015", "Papua New Guinea-West New Britain-2015")


df2 <- df[!(df$eu %in% c(unpublished_EUs, left_out_eus)), ] %>% 
  group_by(eu) %>% mutate(median_scr = median(scr, na.rm = TRUE)) %>% 
  ungroup() %>% 
  as.data.frame()
eus <- unique(df2$eu) # check eu labels!

# assign categories - elimination of transmission:
df2 <- df2 %>%
  mutate(trachoma_cat = case_when(startsWith(as.character(eu), "Ghana") ~ "Low",
                                  eu %in% c("Ethiopia-Alefa-2017","Niger-MORDOR/Dosso-2015","Ethiopia-Woreta Town-2017","Togo-Anie-2017","Togo-Keran-2017","Morocco-Agdaz-2019","Morocco-Boumalne Dades-2019","Gambia-River Regions-2014","Ethiopia-Metema-2021","Ethiopia-Woreta Town-2021","Malawi-DHO Nkwazi-2014","Malawi-Kasisi/DHO-2014","Malawi-Luzi Kochilira-2014","Malawi-Chapananga-2014") ~ "Low",
                                  eu %in% c("Ethiopia-Andabet-2017","Ethiopia-Goncha-2019","Niger-PRET/Matameye-2013","Ethiopia-Ebinat-2019","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016","Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016", "Tanzania-Kongwa-2018", "Tanzania-Kongwa-2013", "Solomon Islands-Temotu/Rennel/Bellona-2015") ~ "High",
                                  TRUE ~ as.character("Medium")),
         trachoma_cat = factor(trachoma_cat, levels = c("High","Medium","Low"))
  )

scr_elim <- df2 %>% filter(trachoma_cat == "Low")
scr_interm <- df2 %>% filter(trachoma_cat == "Medium") 
scr_endemic <- df2 %>% filter(trachoma_cat == "High") 


#------------------------------------------------------------------------------------------------
# Posterior probability of each category at each value of the SCR
#-------------------------------------------------------------------------------------------------
#' Steps:
# (1) prior probability of 'post elimination'
# (2) likelihood - use approx to get the interpolated values of the pdf/density at each SCR estimate 
# (3) calculate marginal probabilities
# (4) use Bayes theorem for posterior inference
# (5) loop through the five different priors

#' (1) Prior probabilities:
#' priors:
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

#' likelihood:
# likelihoods - defined as probability densities and not probability:
dx_elim <- density(scr_elim$scr, n = 2^12, from = min(scr_elim$scr), to = max(scr_elim$scr))
dx_endemic <- density(scr_endemic$scr, n = 2^12, from = min(scr_endemic$scr), to = max(scr_endemic$scr))

# post elimination:
A <- foreach(i=prior_elim, k=prior_endemic) %:%
  foreach(c=scr_elim$scr, .combine = 'rbind', .multicombine = TRUE) %dopar% {
    prob_c_elim <- approx(dx_elim$x, dx_elim$y, xout = c)$y
    prob_c_endemic <- approx(dx_endemic$x, dx_endemic$y, xout = c)$y
    
    numerator <- i * prob_c_elim
    denom <- sum((i * prob_c_elim), (k * prob_c_endemic), na.rm = TRUE)
    
    prob_elim_c <- numerator / denom
    data.frame(i, c, prob_elim_c)
  }

class(A)
prob_elim_df <- do.call(rbind, A)
colnames(prob_elim_df) <- c("prior_set", "c", "prob")
write.csv(prob_elim_df, "output/probs_scr_elim-1to5yo-2categories.csv", row.names = FALSE)


# Endemic populations:
C <- foreach(i=prior_elim, k=prior_endemic) %:% # prior
  foreach(c=scr_endemic$scr, .combine = 'rbind', .multicombine = TRUE) %dopar% { 
    prob_c_elim <- approx(dx_elim$x, dx_elim$y, xout = c)$y
    prob_c_endemic <- approx(dx_endemic$x, dx_endemic$y, xout = c)$y
    
    numerator <- k * prob_c_endemic
    denom <- sum((i * prob_c_elim), (k * prob_c_endemic), na.rm = TRUE)
    
    prob_endemic_c <- numerator / denom
    data.frame(k, c, prob_endemic_c)
  }

class(C)
prob_endemic_df <- do.call(rbind, C)
colnames(prob_endemic_df) <- c("prior_set", "c", "prob")
write.csv(prob_endemic_df, "output/probs_scr_endemic-1to5yo-2categories.csv", row.names = FALSE)

#------------------------------------------------------------------------------------------
# Summary table:
#------------------------------------------------------------------------------------------
prob_elim_df$c <- prob_elim_df$c * 100

scr_df1 <- prob_elim_df %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "90%") %>% as.data.frame()

scr_df2 <- prob_elim_df %>% filter(prob >= 0.91 & prob < 0.92) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "91%") %>% as.data.frame() 

scr_df3 <- prob_elim_df %>% filter(prob >= 0.92 & prob < 0.93) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ")%>% mutate(Probability = "92%") %>% as.data.frame()

scr_df4 <- prob_elim_df %>% filter(prob >= 0.93 & prob < 0.94) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "93%") %>% as.data.frame()

scr_df5 <- prob_elim_df %>% filter(prob >= 0.94 & prob < 0.95) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "94%") %>% as.data.frame()

scr_df6 <- prob_elim_df %>% filter(prob >= 0.95 & prob < 0.96) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "95%") %>% as.data.frame()

scr_df7 <- prob_elim_df %>% filter(prob >= 0.96 & prob < 0.97) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "96%") %>% as.data.frame()

scr_df8 <- prob_elim_df %>% filter(prob >= 0.97 & prob < 0.98) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "97%") %>% as.data.frame()

scr_df9 <- prob_elim_df %>% filter(prob >= 0.98 & prob < 0.99) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "98%") %>% as.data.frame()

scr_df10 <- prob_elim_df %>% filter(prob >= 0.99) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = ">99%") %>% as.data.frame()

scr_table_prob <- bind_rows(scr_df1, scr_df2, scr_df3, scr_df4, scr_df5, scr_df6, scr_df7, scr_df8, scr_df9, scr_df10)
head(scr_table_prob)

scr_table_prob2 <- scr_table_prob %>% filter(prior_set == 0.8)

write.csv(scr_table_prob, "output/table_probabilities_for_cutoff-2categories.csv", row.names = FALSE)
write.csv(scr_table_prob2, "output/table_probabilities_for_cutoff-prior0.8-2categories.csv", row.names = FALSE)

#-------------------------------------------------------------------------------------------------
# Summarize / calculate cutoff 'B':
# Near elimination of category posterior probabilities:
#-------------------------------------------------------------------------------------------------
scr_df11 <- prob_endemic_df %>% filter(prob >= 0.5 & prob < 0.51) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "50%",
                                                  trachoma_cat = "Endemic transmission") %>% 
  as.data.frame()

scr_df12 <- prob_endemic_df %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "90%",
                                                  trachoma_cat = "Endemic transmission") %>% 
  as.data.frame()

write.csv(rbind(scr_df11, scr_df12), "output/table_probabilities_2categories-cutoffB.csv", row.names = FALSE)

