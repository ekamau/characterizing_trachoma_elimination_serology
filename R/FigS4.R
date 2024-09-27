#' ---
#' title: "Sensitivity analysis - prior probability - Fig S4"
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

# Setup parallel computing:
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(n.cores, type = "FORK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() #check if registered
foreach::getDoParWorkers() #how many workers are available?

#------------------------------------------------------------------------------------------------
#' Prior sensitivity: more weight in the action not needed group
#------------------------------------------------------------------------------------------------
prob_elim_df1 <- read.csv("output/elim-prior-sensitivity-2states-v1.csv")
prob_endemic_df1 <- read.csv("output/endemic-prior-sensitivity-2states-v1.csv")
prior_elim1 <- unique(prob_elim_df1$prior_set)
prior_endemic1 <- unique(prob_endemic_df1$prior_set)

#------------------------------------------------------------------------------------------------
# Prior sensitivity - more weight in the action needed group:
#------------------------------------------------------------------------------------------------
prob_elim_df2 <- read.csv("output/elim-prior-sensitivity-2states-v2.csv")
prob_endemic_df2 <- read.csv("output/endemic-prior-sensitivity-2states-v2.csv")
prior_elim2 <- unique(prob_elim_df2$prior_set)
prior_endemic2 <- unique(prob_endemic_df2$prior_set)

#-------------------------------------------------------------------------------------------------
#' Plots:
colPal <- c("black", "#ee7733", "#0077bb", "#BE0032", "#009988")
# legend labels:
prior_labels = c(paste0(round(prior_elim1[1],2)," | ",round(prior_endemic1[1],2)), paste0(round(prior_elim1[2],2)," | ",round(prior_endemic1[2],2)),
                 paste0(round(prior_elim1[3],2)," | ",round(prior_endemic1[3],2)), paste0(round(prior_elim1[4],2)," | ",round(prior_endemic1[4],2)),
                 paste0(round(prior_elim1[5],2)," | ",round(prior_endemic1[5],2)))

# Post elimination:
FigA1 <- ggplot(data = prob_elim_df1, aes(x = c*100, y = prob)) +
  geom_line(aes(color = factor(prior_set, levels = c(0.8,0.7,0.6,0.5,0.4))), linewidth = 0.6) +
  scale_color_viridis(name = "Prior probability\n(Action not needed | Action needed)", discrete = TRUE,
                      labels = prior_labels, breaks = prior_elim1) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "top")) +
  coord_cartesian(xlim = c(0, 6)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "SCR per 100 person-years", y = "Posterior probability", subtitle = "Action not needed") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 11),
        legend.position = "bottom")

# Endemic:
FigA3 <- ggplot(data = prob_endemic_df1, aes(x = c*100, y = prob)) +
  geom_line(aes(color = as.factor(prior_set)), linewidth = 0.6) +
  scale_color_viridis(name = "Prior probability\n(Action not needed | Action needed)", discrete = TRUE,
                      labels = prior_labels, breaks = prior_endemic1) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "top")) +
  coord_cartesian(xlim = c(0, 6)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "SCR per 100 person-years", y = "", subtitle = "Action needed") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 11),
        legend.position = "bottom")


figS6 <- (FigA1 + FigA3) + plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(filename = here("output/", "FigS4.png"), figS6, width = 8, height = 5, dpi = 300)

