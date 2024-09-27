#' ---
#' title: "Plot posterior probabilities of SCR estimates - Fig SII"
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

#' (1) Prior probabilities:
#' priors:-- x = (-b ± √ (b2 - 4ac) )/2a
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


#' (2) read data:
prob_elim_df <- read.csv("output/probs_scr_elim-1to5yo-3states.csv")
prob_near_elim_df <- read.csv("output/probs_scr_interm-1to5yo-3states.csv")
prob_endemic_df <- read.csv("output/probs_scr_endemic-1to5yo-3states.csv")

prior_elim <- unique(prob_elim_df$prior_set)
prior_interm <- unique(prob_near_elim_df$prior_set)
prior_endemic <- unique(prob_endemic_df$prior_set)

colPal <- c("#A42820", "#F98400", "cadetblue")
unique(prob_elim_df$prior_set)

dfA <- prob_elim_df %>% filter(prior_set %in% c(prior_elim[5])) %>% 
  mutate(trachoma_cat = "Probable interruption of transmission")
dfB <- prob_near_elim_df %>% filter(prior_set %in% c(prior_interm[5])) %>% 
  mutate(trachoma_cat = "Near interruption of transmission")
dfC <- prob_endemic_df %>% filter(prior_set %in% c(prior_endemic[5])) %>% 
  mutate(trachoma_cat = "Endemic")

# get the coordinates for the intermediate region:
scrA_prob0.9 <- prob_elim_df %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(prior_set) %>% 
  reframe(c = max(c)) %>% # or you can take the min value or mean of min + max
  filter(prior_set == prior_elim[5]) %>% 
  mutate(bound = "lower") %>% as.data.frame()

scrB_prob0.9 <- prob_endemic_df %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(prior_set) %>% 
  reframe(c = max(c)) %>% # or you can take the min value or mean of min + max
  filter(prior_set == prior_endemic[5]) %>% 
  mutate(bound = "upper") %>% as.data.frame()


rbind(dfA, dfB, dfC) %>% 
  ggplot(aes(x = c*100, y = prob)) +
  geom_line(aes(color = trachoma_cat), linewidth = 1) +
  #annotate("text", x = 7, y = 0.5, label = "+", size = 6, hjust = 0.75) +
  #annotate("text", x = 1.04, y = 0.9, label = "+", size = 6) +
  annotate("label", x = 1, y = 1.08, label = "Probable interruption\nof transmission", color = colPal[3], 
           fontface = 2, size = 3.5) +
  annotate("label", x = 4, y = 1.08, label = "Near interruption\nof transmission", color = colPal[2], 
           fontface = 2, size = 3.5) +
  annotate("label", x = 9, y = 1.08, label = "Endemic", color = colPal[1], fontface=2, size = 3.5) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(xlim = c(0, 10), clip = 'off') +
  scale_color_manual(name = "", breaks = c("Probable interruption of transmission", "Near interruption of transmission",
                                           "Endemic"), 
                     values = c("Probable interruption of transmission"="cadetblue", "Near interruption of transmission"=colPal[2],
                                "Endemic"="#A42820")) +
  labs(x = "SCR per 100 person-years", y = "Posterior probability") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        axis.title = element_text(size = 11),
        plot.caption = element_text(hjust = 0.5),
        legend.position = 'none')


ggsave(filename = here("output/", "FigS-II.png"), width = 8, height = 5, dpi = 300)
