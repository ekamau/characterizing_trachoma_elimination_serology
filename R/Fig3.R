#' ---
#' title: "Plot posterior probabilities of SCR estimates - Fig 3"
#' output: html_document
#' author: ekamau
#' ---
#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE, eval = FALSE)

#+ warning = FALSE, message = FALSE, eval = FALSE
# packages
if (!require("pacman")) install.packages("pacman")
p_load(tidyverse, ggridges, reshape2, viridis, here, patchwork, ggtext, foreach)

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


#' (2) read data:
prob_elim_df <- read.csv("output/probs_scr_elim-1to5yo-2states.csv")
prob_endemic_df <- read.csv("output/probs_scr_endemic-1to5yo-2states.csv")

prior_elim <- unique(prob_elim_df$prior_set)
prior_endemic <- unique(prob_endemic_df$prior_set)

colPal <- c("#A42820", "cadetblue")
unique(prob_elim_df$prior_set)

dfA <- prob_elim_df %>% filter(prior_set %in% c(prior_elim[5])) %>% 
  mutate(trachoma_cat = "Action not needed")
dfC <- prob_endemic_df %>% filter(prior_set %in% c(prior_endemic[5])) %>% 
  mutate(trachoma_cat = "Action needed")

# get the coordinates for the intermediate region:
scrA_prob0.9 <- prob_elim_df %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(prior_set) %>% 
  reframe(c = max(c)) %>% # or you can take the min value or mean of min + max
  filter(prior_set == prior_elim[5]) %>% 
  mutate(bound = "lower") %>% as.data.frame()

scrB_prob0.9 <- prob_endemic_df %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(prior_set) %>% 
  reframe(c = max(c)) %>% # or you can take the min value or mean of min + max
  filter(prior_set == prior_endemic[5]) %>% 
  mutate(bound = "upper") %>% as.data.frame()


rbind(dfA, dfC) %>% 
  ggplot(aes(x = c*100, y = prob)) +
  geom_line(aes(color = trachoma_cat), linewidth = 1) +
  annotate("rect", fill = "grey", alpha = 0.3, xmin = scrA_prob0.9$c*100, xmax = scrB_prob0.9$c*100, 
           ymin = -Inf, ymax = Inf) +
  geom_point(data = scrA_prob0.9, aes(x = c*100, y = 0.9), shape = 3, size = 3) +
  geom_text(data = scrA_prob0.9, aes(x = c*100, y = 0.9), label="2.2", vjust = -1, hjust = -0.1) +
  geom_point(data = scrB_prob0.9, aes(x = c*100, y = 0.9), shape = 3, size = 3) +
  geom_text(data = scrB_prob0.9, aes(x = c*100, y = 0.9), label="4.5", vjust = -1, hjust = 0.9) +
  annotate("label", x = 1, y = 1.03, label = "Action not needed", color = colPal[2], fontface=2) +
  annotate("label", x = 5.5, y = 1.03, label = "Action needed", color = colPal[1], fontface=2) +
  scale_x_continuous(breaks = seq(0, 12, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(xlim = c(0, 6), clip = 'off') +
  scale_color_manual(name = "", breaks = c("Action not needed", "Action needed"), 
                     values = c("Action not needed"="cadetblue", "Action needed"="#A42820")) +
  labs(x = "SCR per 100 person-years", y = "Posterior probability") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        axis.title = element_text(size = 11),
        plot.caption = element_text(hjust = 0.5),
        legend.position = 'none')


ggsave(filename = here("output/", "Fig3.png"), width = 5, height = 5, dpi = 300)

write.csv(scr_table_prob, "output/output_ages_1to5/table_probabilities_for_cutoff.csv", row.names = FALSE)
write.csv(scr_table_prob2, "output/output_ages_1to5/table_probabilities_for_cutoff-prior0.8.csv", row.names = FALSE)


