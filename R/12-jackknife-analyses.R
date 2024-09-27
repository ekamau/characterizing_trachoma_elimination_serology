#' ---
#' title: "Jack knife analysis for LOO analyses"
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
p_load(tidyverse, ggridges, reshape2, viridis, here, patchwork, ggtext, foreach, ggrepel)

# Setup parallel computing:
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(n.cores, type = "FORK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() # check if registered
foreach::getDoParWorkers() # how many workers are available?

# read data 
scr_table_prob <- read.csv("output/LOO_analyses_part1/table_probabilities_LOO-EU.csv")
scr_table_prob2 <- read.csv("output/LOO_analyses_part2/table_probabilities_LOO-country.csv")

# jackknife analysis:
scr_table_prob_means <- scr_table_prob %>% dplyr::filter(!Data == "All EUs") %>% 
  group_by(Probability) %>% summarise(mean_scr = mean(max)) %>% as.data.frame()

scr_table_prob_means2 <- scr_table_prob2 %>% dplyr::filter(!Data == "All EUs") %>% 
  group_by(Probability) %>% 
  summarise(mean_scr = mean(max),
            mean_scr2 = exp(log(mean(max)))) %>% as.data.frame()

figA <- ggplot() +
  geom_point(data=scr_table_prob, 
             aes(x = max, y = factor(Probability, levels = unique(Probability)), color = as.factor(Data), 
                 alpha = as.factor(Data)), position = position_jitter(height = 0.01), size = 2) +
  geom_point(data=scr_table_prob_means, 
             aes(x = mean_scr, y = factor(Probability, levels = unique(Probability))), color = "black",
             shape = 13, size = 3) +
  #annotate("segment", x = 1, xend = 3, y = 25, yend = 15, colour = "purple", size=3, alpha=0.6)
  scale_color_manual(values = c("All EUs"="#D55E00","N-1"="grey")) +
  scale_alpha_manual(values = c("All EUs"=1,"N-1"=0.4)) +
  scale_x_continuous(breaks = seq(0, 1.6, 0.2)) +
  labs(x = "SCR per 100 PY", y = "Posterior probability", title = "EU level LOO analysis") +
  guides(color = "none", alpha = "none") +
  theme_minimal()

figB <- ggplot() +
  geom_point(data=scr_table_prob2, 
             aes(x = max, y = factor(Probability, levels = unique(Probability)), color = as.factor(Data), 
                 alpha = as.factor(Data)), position = position_jitter(height = 0.01), size = 2) +
  geom_point(data=scr_table_prob_means2, 
             aes(x = mean_scr, y = factor(Probability, levels = unique(Probability))), color = "black",
             shape = 13, size = 3) +
  #annotate("segment", x = 1, xend = 3, y = 25, yend = 15, colour = "purple", size=3, alpha=0.6)
  scale_color_manual(values = c("All EUs"="#D55E00","N-1"="grey")) +
  scale_alpha_manual(values = c("All EUs"=1,"N-1"=0.4)) +
  scale_x_continuous(breaks = seq(0, 1.6, 0.2)) +
  labs(x = "SCR per 100 PY", y = "", title = "Country level LOO analysis") +
  guides(color = "none", alpha = "none") +
  theme_minimal()

figS <- (figA + figB) + plot_annotation(title = 'Post elimination of transmission',
                                        subtitle = 'Prior probability = 0.8')

ggsave(filename = here("output/", "FigS-LOO-analysis-jackknife.png"), 
       figS, width = 8, height = 6, dpi = 640)

