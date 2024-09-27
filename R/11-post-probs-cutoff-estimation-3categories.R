#' ---
#' title: "Determine thresholds from empirical SCR estimates - 3 EU categories"
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

unpublished_EUs <- c("Niger-Ilela-2022", "Niger-Malbaza-2022", "Niger-Bagaroua-2022", 
                     "DRC-Manono-2018", "DRC-Nyemba-2018")

eus_wo_PCR <- c("Sudan-El Seraif-2019","Sudan-Saraf Omrah-2019","Sudan-Kotom-2019","Malaysia-Sabah-2015","Peru-Amazonia-2020",
                  "Papua New Guinea-Mendi-2015","Papua New Guinea-Daru-2015","Papua New Guinea-West New Britain-2015","Togo-Anie-2017",
                  "Togo-Keran-2017","Morocco-Agdaz-2019","Morocco-Boumalne Dades-2019","Gambia-River Regions-2014")

df2 <- df[!(df$eu %in% c(unpublished_EUs, eus_wo_PCR)), ] %>% 
  group_by(eu) %>% mutate(median_scr = median(scr, na.rm = TRUE)) %>% 
  ungroup() %>% as.data.frame()
eus <- unique(df2$eu) # check eu labels!

# assign categories - elimination of transmission:
df2 <- df2 %>% mutate(trachoma_cat = case_when(startsWith(as.character(eu), "Ghana") ~ "Elimination",
                                               eu %in% c("Ethiopia-Alefa-2017","Niger-MORDOR/Dosso-2015","Ethiopia-Woreta Town-2017","Ethiopia-Metema-2021","Ethiopia-Woreta Town-2021") ~ "Elimination",
                                               eu %in% c("Ethiopia-Andabet-2017","Niger-PRET/Matameye-2013","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016","Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016",
                                                         "Ethiopia-Ebinat-2019") ~ "Endemic",
                                               TRUE ~ as.character("Near elimination")),
                      trachoma_cat = factor(trachoma_cat, levels = c("Endemic", "Near elimination", "Elimination"))
)

scr_elim <- df2 %>% filter(trachoma_cat == "Elimination") 
scr_interm <- df2 %>% filter(trachoma_cat == "Near elimination") 
scr_endemic <- df2 %>% filter(trachoma_cat == "Endemic") 

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
#' # x = (-b ± √ (b2 - 4ac) )/2a
quadraticRoots <- function(a, b, c) {
  print(paste0(a, "x^2 + ", b, "x + ", c, "."))
  discriminant <- (b^2) - (4*a*c) 
  if(discriminant < 0) {
    return(paste0("This quadratic equation has no real numbered roots."))
  }
  else if(discriminant > 0) {
    x_int_plus <- (-b + sqrt(discriminant)) / (2*a)
    x_int_neg <- (-b - sqrt(discriminant)) / (2*a)
    # return(paste0("The two x-intercepts for the quadratic equation are ",
    #               format(round(x_int_plus, 5), nsmall = 5), " and ",
    #               format(round(x_int_neg, 5), nsmall = 5), "."))
    return(round(x_int_plus, 5))
  }
  else #discriminant = 0  case
    x_int <- (-b) / (2*a)
  #return(paste0("The quadratic equation has only one root. ", x_int))
  return(x_int)
}

quadraticRoots(0.8,0.8,(0.8-1))

#' priors:
prior_elim <- c(1/3, 0.65, 0.7, 0.75, 0.8)
prior_interm <- c()
prior_endemic <- c()

for(p in prior_elim){
  print(p)
  root <- quadraticRoots(a = p, b = p, c = (p-1))
  prior_interm <- append(prior_interm, p * root, after = length(prior_interm))
  prior_endemic <- append(prior_endemic, p * root^2, after = length(prior_endemic))
}

# confirm the probabilities add to one:
tibble(prior_elim, prior_interm, prior_endemic) %>% 
  rowwise %>% mutate(total = sum(c(prior_elim, prior_interm, prior_endemic))) %>% 
  as.data.frame()

#' likelihood:
# likelihoods - defined as probability densities and not probability:
dx_elim <- density(scr_elim$scr, n = 2^12, from = min(scr_elim$scr), to = max(scr_elim$scr))
dx_interm <- density(scr_interm$scr, n = 2^12, from = min(scr_interm$scr), to = max(scr_interm$scr))
dx_endemic <- density(scr_endemic$scr, n = 2^12, from = min(scr_endemic$scr), to = max(scr_endemic$scr))

# post elimination:
A <- foreach(i=prior_elim, j=prior_interm, k=prior_endemic) %:%
  foreach(c=scr_elim$scr, .combine = 'rbind', .multicombine = TRUE) %dopar% {
    prob_c_elim <- approx(dx_elim$x, dx_elim$y, xout = c)$y
    prob_c_interm <- approx(dx_interm$x, dx_interm$y, xout = c)$y
    prob_c_endemic <- approx(dx_endemic$x, dx_endemic$y, xout = c)$y
    
    numerator <- i * prob_c_elim
    denom <- sum((i * prob_c_elim), (j * prob_c_interm), (k * prob_c_endemic), na.rm = TRUE)
    
    prob_elim_c <- numerator / denom
    data.frame(i, c, prob_elim_c)
  }

class(A)
prob_elim_df <- do.call(rbind, A)
colnames(prob_elim_df) <- c("prior_set", "c", "prob")
write.csv(prob_elim_df, "output/probs_scr_elim-1to5yo-3categories.csv", row.names = FALSE)

# Near elimination of group:
B <- foreach(i=prior_elim, j=prior_interm, k=prior_endemic) %:% # prior
  foreach(c=scr_interm$scr, .combine = 'rbind', .multicombine = TRUE) %dopar% { 
    prob_c_elim <- approx(dx_elim$x, dx_elim$y, xout = c)$y
    prob_c_interm <- approx(dx_interm$x, dx_interm$y, xout = c)$y
    prob_c_endemic <- approx(dx_endemic$x, dx_endemic$y, xout = c)$y
    
    numerator <- j * prob_c_interm
    denom <- sum((i * prob_c_elim), (j * prob_c_interm), (k * prob_c_endemic), na.rm = TRUE)
    
    prob_interm_c <- numerator / denom
    data.frame(j, c, prob_interm_c)
  }

class(B)
prob_near_elim_df <- do.call(rbind, B)
colnames(prob_near_elim_df) <- c("prior_set", "c", "prob")
write.csv(prob_near_elim_df, "output/probs_scr_interm-1to5yo-3categories.csv", row.names = FALSE)

# Endemic populations:
C <- foreach(i=prior_elim, j=prior_interm, k=prior_endemic) %:% # prior
  foreach(c=scr_endemic$scr, .combine = 'rbind', .multicombine = TRUE) %dopar% { 
    prob_c_elim <- approx(dx_elim$x, dx_elim$y, xout = c)$y
    prob_c_interm <- approx(dx_interm$x, dx_interm$y, xout = c)$y
    prob_c_endemic <- approx(dx_endemic$x, dx_endemic$y, xout = c)$y
    
    numerator <- k * prob_c_endemic
    denom <- sum((i * prob_c_elim), (j * prob_c_interm), (k * prob_c_endemic), na.rm = TRUE)
    
    prob_endemic_c <- numerator / denom
    data.frame(k, c, prob_endemic_c)
  }

class(C)
prob_endemic_df <- do.call(rbind, C)
colnames(prob_endemic_df) <- c("prior_set", "c", "prob")
write.csv(prob_endemic_df, "output/probs_scr_endemic-1to5yo-3categories.csv", row.names = FALSE)

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

write.csv(scr_table_prob, "output/table_probabilities_3categories.csv", row.names = FALSE)
write.csv(scr_table_prob2, "output/table_probabilities_3categories-prior0.8.csv", row.names = FALSE)

#-------------------------------------------------------------------------------------------------
# Summarize / calculate cutoff 'B':
# Near elimination of category posterior probabilities:
#-------------------------------------------------------------------------------------------------
scr_df11 <- prob_near_elim_df %>% filter(prob >= 0.5 & prob < 0.51) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "50%",
                                                  trachoma_cat = "Near elimination of transmission") %>% 
  as.data.frame()

scr_df12 <- prob_endemic_df %>% filter(prob >= 0.5 & prob < 0.51) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "50%",
                                                  trachoma_cat = "Endemic transmission") %>% 
  as.data.frame()

scr_df13 <- prob_endemic_df %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(prior_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  unite('Range', min:max, sep = " - ") %>% mutate(Probability = "90%",
                                                  trachoma_cat = "Endemic transmission") %>% 
  as.data.frame()

write.csv(rbind(scr_df11, scr_df12, scr_df13), "output/table_probabilities_3categories_cutoffB.csv", row.names = FALSE)

