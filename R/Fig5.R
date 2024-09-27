#' ---
#' title: "Posterior probabilities for held-out EUs - Fig 5"
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

# left EUs:
eus_wo_PCR <- c("Sudan-El Seraif-2019","Sudan-Saraf Omrah-2019","Sudan-Kotom-2019","Malaysia-Sabah-2015",
                "Peru-Amazonia-2020","Papua New Guinea-Mendi-2015","Papua New Guinea-Daru-2015","Papua New Guinea-West New Britain-2015")

unclassified <- c("Malawi-Ngabu Ngokwe-2014","Malawi-Mkanda Gumba-2014","Ethiopia-Dera-2017","Ethiopia-Debay Tilatgin-2019",
                  "Ethiopia-Machakel-2019", "Vanuatu-Torba/Melampa/Penama/Shefa/Tafea/Sanma-2016")

#------------------------------------------------------------------------------------------------
#' Plot posterior probabilities per EU:
#------------------------------------------------------------------------------------------------
#' read data:
prob_df <- read.csv("output/prob_SCR-heldOutEUs-2states.csv")

# Plot:
df55 <- prob_df %>% dplyr::filter(prior_elim == unique(prob_df$prior_elim)[5] & 
                                    prior_endemic == unique(prob_df$prior_endemic)[5])

eu_sorted_SCR <- df55 %>% group_by(eu) %>% summarise(med_scr = median(scr, na.rm = TRUE),
                                                     mean_prob_elim = mean(prob_in_post_elim, na.rm = TRUE),
                                                     mean_prob_endemic = mean(prob_in_endemic, na.rm = TRUE)) %>% 
  arrange(-med_scr) %>% select(eu) %>% unique() %>% as.data.frame() 

eu_sorted_SCR <- eu_sorted_SCR[,1]

p1 <- df55 %>% group_by(eu) %>% summarise(med_scr = median(scr, na.rm = TRUE),
                                          mean_prob_elim = mean(prob_in_post_elim, na.rm = TRUE),
                                          mean_prob_endemic = mean(prob_in_endemic, na.rm = TRUE)) %>% 
  rename("Action not needed" = "mean_prob_elim",
         "Action needed" = "mean_prob_endemic") %>% 
  reshape2::melt(id.vars = c('eu', 'med_scr'), variable.name = 'mean_prob', 
                 value.name = 'probability') %>%
  arrange(-med_scr) %>%
  as.data.frame() %>% 
  ggplot(aes(x = factor(eu, levels = rev(eu_sorted_SCR)), y = probability, fill = forcats::fct_rev(mean_prob))) + 
  geom_bar(stat = "identity", alpha = 0.7,  width = 0.5) +
  #geom_text(aes(label = round(probability, digits = 2)), size = 2, 
  #          vjust = 0.5, fontface = "bold", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(name = "", breaks = c("Action not needed", "Action needed"),
                    values = c("Action not needed"="cadetblue", "Action needed"="#A42820")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Survey\n", y = "Posterior probability", tag = "A") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 7),
        legend.box.margin = margin(6, 6, 6, 6),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_flip()

xx <- df55 %>% group_by(eu) %>% summarise(med_scr = median(scr)) %>% 
  arrange(-med_scr) %>%
  as.data.frame()

p2 <- gridExtra::tableGrob(round(xx$med_scr*100,1), theme = ttheme_minimal(base_size = 10))
# Set widths/heights to 'fill whatever space I have'
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

# Format table as plot
p3 <- ggplot() + annotation_custom(p2) + theme_bw() + xlab("SCR per\n100 PY")

# Patchwork combine:
figA <- p1 + p3 + plot_layout(ncol = 2) +
  plot_layout(widths = c(4, 1)) & theme(plot.title = element_text(hjust = 0.5))

#------------------------------------------------------------------------------------------------
#' Plot empirical probabilities per EU:
#------------------------------------------------------------------------------------------------
# Read Bayesian SCR estimates:
here::here()
df <- read.csv(here("output/", "posterior_samples_scr-1to5yo.csv"))
eus <- unique(df$eu) # check eu labels!

df0 <- df[(df$eu %in% eu_sorted_SCR), ] %>%
  group_by(eu) %>% mutate(median_scr = median(scr, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame()

# Empirical probabilities:
prob_df2 <- data.frame()
cutoff_elim <- 2.2 # upper limit of 90% posterior probability for 0.8 prior probability
cutoff_endemic <- 4.5 # SCR per 100 PY -> upper limit of 50% posterior probability for 0.33 prior probability 

# get area under the curve or proportion below a certain x value:
for(i in eu_sorted_SCR){
  print(i)
  di <- df %>% dplyr::filter(eu == i)
  prob_elim <- mean(di$scr*100 < cutoff_elim)
  prob_endemic <- mean(di$scr*100 > cutoff_endemic)
  prob_df2 <- rbind(prob_df2, data.frame(i, median(di$scr)*100, prob_elim, prob_endemic))
}

colnames(prob_df2) <- c("eu", "median_scr", "prob_elimination", "prob_endemic")

# Plot:
# (A) SCR density - w/o median line, but showing the threshold as vertical line

figB <- ggplot(df0) +
  geom_density_ridges(aes(x = scr*100, y = fct_reorder(eu, median_scr*100, .desc = FALSE), group = eu),
                      rel_min_height = 0.01, scale = 1) +
  geom_vline(xintercept = cutoff_elim, colour = "cadetblue", linetype = "longdash", linewidth = 0.5) +
  coord_cartesian(xlim = c(0,12)) +
  scale_x_continuous(breaks = seq(0,12,by=2)) +
  labs(x = "SCR per 100\nperson-years", y = "", tag = "B") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# probability of being < cutoff_elimination (right panel)
figC <- prob_df2 %>% select(-prob_endemic) %>% arrange(median_scr) %>% 
  ggplot() +
  geom_bar(aes(x = factor(eu, levels = unique(eu)), y = prob_elimination), fill = "cadetblue", 
           stat = "identity", alpha = 0.7) +
  labs(x = "", y = paste0("Probability SCR ≤", round(cutoff_elim,1)), tag = "C") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.ticks.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank()) +
  coord_flip() 

fig <- p1 + figB + figC
ggsave(filename = "output/Fig5.png", fig, width = 8, height = 6, dpi = 640)
