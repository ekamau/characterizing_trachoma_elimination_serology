#' title: "Posterior probabilities for held-out EUs - Fig SIII"
#' output: html_document
#' author: ekamau
#' ---
#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE, eval = FALSE)

#+ warning = FALSE, message = FALSE, eval = FALSE
# packages
if (!require("pacman")) install.packages("pacman")
p_load(tidyverse, ggridges, reshape2, viridis, here, patchwork, ggtext, foreach, 
       flextable, cowplot, grid, gridExtra, kableExtra)

# Setup parallel computing:
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(n.cores, type = "FORK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() #check if registered
foreach::getDoParWorkers() #how many workers are available?

#------------------------------------------------------------------------------------------------
#' Plot posterior probabilities per EU:
#------------------------------------------------------------------------------------------------
#' read data:
prob_df1 <- read.csv("output/prob_SCR-heldOutEUs-3states.csv")

# Plot:
df55 <- prob_df1 %>% dplyr::filter(prior_elim == unique(prob_df1$prior_elim)[5] & 
                                                 prior_interm == unique(prob_df1$prior_interm)[5] &
                                                 prior_endemic == unique(prob_df1$prior_endemic)[5])

eu_sorted_SCR <- df55 %>% group_by(eu) %>% summarise(med_scr = median(scr),
                                                     mean_prob_elim = mean(prob_in_post_elim),
                                                     mean_prob_interm = mean(prob_in_interm),
                                                     mean_prob_endemic = mean(prob_in_endemic)) %>% 
  arrange(-med_scr) %>% select(eu) %>% unique() %>% as.data.frame() 

eu_sorted_SCR <- eu_sorted_SCR[,1]

p1 <- df55 %>% group_by(eu) %>% summarise(med_scr = median(scr, na.rm = TRUE),
                                          mean_prob_elim = mean(prob_in_post_elim, na.rm = TRUE),
                                          mean_prob_interm = mean(prob_in_interm, na.rm = TRUE),
                                          mean_prob_endemic = mean(prob_in_endemic, na.rm = TRUE)) %>% 
  rename("Probable interruption of transmission" = "mean_prob_elim",
         "Near interruption of transmission" = "mean_prob_interm",
         "Endemic" = "mean_prob_endemic") %>% 
  reshape2::melt(id.vars = c('eu', 'med_scr'), variable.name = 'mean_prob', 
                 value.name = 'probability') %>%
  arrange(-med_scr) %>% drop_na(probability) %>% as.data.frame() %>% 
  ggplot(aes(x = factor(eu, levels = rev(eu_sorted_SCR)), y = round(probability, digits = 1), fill = forcats::fct_rev(mean_prob))) + 
  geom_bar(stat = "identity", alpha = 0.7,  width = 0.5) +
  scale_fill_manual(name = "", breaks = c("Probable interruption of transmission", 
                                          "Near interruption of transmission", "Endemic"),
                    values = c("Probable interruption of transmission"="cadetblue", 
                               "Near interruption of transmission"="#F98400", "Endemic"="#A42820")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Survey\n", y = "Posterior probability") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 8),
        legend.box.margin = margin(6, 6, 6, 6),
        axis.title = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_flip(ylim = c(0,1))

ggsave(filename = here("output/", "FigS-III.png"), p1, width = 8, height = 5, dpi = 300)
