#' ---
#' title: "Plot thresholds from SCR estimates - LOO analysis (N-1) - Fig 4"
#' output: html_document
#' author: ekamau
#' ---
#' 
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

# read data - all EUs and countries included:
prob_elim_df <- read.csv("output/probs_scr_elim-1to5yo-2states.csv")

#------------------------------------------------------------------------------------------------
# LOO - EUs
#------------------------------------------------------------------------------------------------
#' Plots - post elimination
#' read data - LOO (n-1) analysis:
prob_elim_df2 <- readRDS("output/probs_scr_elimination-eu-N1-2states.rds")
colnames(prob_elim_df2) <- c("LOO_set", "EU_left", "prior_set", "c", "prob")
prior_elim <- unique(prob_elim_df$prior_set)
prior_elim2 <- unique(prob_elim_df2$prior_set)

df1 <- prob_elim_df %>% dplyr::filter(prior_set == unique(prob_elim_df$prior_set)[5]) %>% 
  mutate(data = "All EUs (N=34)", LOO_set = 35, c = c*100)

df2 <- prob_elim_df2 %>% dplyr::filter(prior_set == unique(prob_elim_df2$prior_set)[5]) %>% 
  mutate(data = "N-1", c = c*100) %>% 
  relocate(LOO_set, .after = last_col())

# Fig A:
## Remove the left out endemic EUs
endemic_eus <- c("Ethiopia-Andabet-2017","Ethiopia-Goncha-2019","Niger-PRET/Matameye-2013",
  "Ethiopia-Ebinat-2019","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016",
  "Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016","Tanzania-Kongwa-2018",
  "Tanzania-Kongwa-2013","Solomon Islands-Temotu/Rennel/Bellona-2015")

fig1A <- df2 %>% dplyr::filter(!EU_left %in% endemic_eus) %>% 
  bind_rows(., df1) %>% 
  ggplot(aes(x = c, y = prob, group = LOO_set)) +
  geom_line(aes(color = as.factor(data)), linewidth = 0.65) +
  annotate("text", label = "Malawi-Chapananga-2014", x = 2.4, y = 0.1, size = 3) +
  annotate("text", label = "Morocco-Agdaz-2019", x = 2.5, y = 0.3, size = 3) +
  scale_color_manual(values = c("red", "gray")) +
  coord_cartesian(xlim = c(0, 6)) +
  scale_x_continuous(breaks = seq(0, 6, by=1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "", y = "Posterior probability", tag = "A", color = "Data:") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"))

#------------------------------------------------------------------------------------------------
# to find the discrepant curves/EUs:
xx <- df2 %>% dplyr::filter(!EU_left %in% endemic_eus)
unique(xx$LOO_set)
#[1]  3  5  6  7  8 10 12 13 14 18 19 20 21 22 23 24 25 26 27 28 29 30 31

df2 %>% dplyr::filter(!EU_left %in% endemic_eus) %>% dplyr::filter(!LOO_set == 19) %>% # change here!
  bind_rows(., df1) %>% 
  ggplot(aes(x = c, y = prob, group = LOO_set)) +
  geom_line(aes(color = as.factor(LOO_set)), linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 6)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "", y = "Posterior probability", color = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

df2 %>% dplyr::filter(LOO_set == 5) %>% summarise(eu_name = unique(EU_left))
df2 %>% dplyr::filter(LOO_set == 19) %>% summarise(eu_name = unique(EU_left))

#------------------------------------------------------------------------------------------------
# Summary of posterior probabilities and SCR:
#------------------------------------------------------------------------------------------------
df3 <- df2 %>% dplyr::filter(!EU_left %in% endemic_eus)

scr_df1 <- bind_rows(df1, df3) %>% filter(prob >= 0.5 & prob < 0.51) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.50") %>% as.data.frame()

scr_df2 <- bind_rows(df1, df3) %>% filter(prob >= 0.6 & prob < 0.61) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.60") %>% as.data.frame() 

scr_df3 <- bind_rows(df1, df3) %>% filter(prob >= 0.65 & prob < 0.66) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.65") %>% as.data.frame()

scr_df4 <- bind_rows(df1, df3) %>% filter(prob >= 0.7 & prob < 0.71) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.70") %>% as.data.frame()

scr_df5 <- bind_rows(df1, df3) %>% filter(prob >= 0.75 & prob < 0.76) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.75") %>% as.data.frame()

scr_df6 <- bind_rows(df1, df3) %>% filter(prob >= 0.8 & prob < 0.81) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.80") %>% as.data.frame()

scr_df7 <- bind_rows(df1, df3) %>% filter(prob >= 0.85 & prob < 0.86) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.85") %>% as.data.frame()

scr_df8 <- bind_rows(df1, df3) %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.90") %>% as.data.frame()

scr_df9 <- bind_rows(df1, df3) %>% filter(prob >= 0.95 & prob < 0.96) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.95") %>% as.data.frame()

scr_df10 <- bind_rows(df1, df3) %>% filter(prob >= 0.99) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = ">0.99") %>% as.data.frame()

scr_table_prob <- bind_rows(scr_df1, scr_df2, scr_df3, scr_df4, scr_df5, scr_df6, scr_df7, scr_df8,
                            scr_df9, scr_df10) %>% 
  mutate(Data = case_when(LOO_set == 35 ~ "All EUs (N=34)", TRUE ~ "N-1"))


# jackknife analysis:
scr_table_prob_means <- scr_table_prob %>% dplyr::filter(!Data == "All EUs (N=34)") %>% 
  group_by(Probability) %>% 
  summarise(mean_scr = mean(max)) %>% as.data.frame()

fig1B <- ggplot() +
  geom_point(data=scr_table_prob, 
             aes(x = max, y = factor(Probability, levels = unique(Probability)), color = as.factor(Data), 
                 alpha = as.factor(Data)), position = position_jitter(height = 0.01), size = 2.5) +
  geom_point(data=scr_table_prob_means, 
             aes(x = mean_scr, y = factor(Probability, levels = unique(Probability))), color = "black",
             shape = 13, size = 2) +
  scale_color_manual(values = c("All EUs (N=34)"="red","N-1"="grey")) +
  scale_alpha_manual(values = c("All EUs (N=34)"=1,"N-1"=0.5)) +
  coord_cartesian(xlim = c(0,4)) +
  scale_x_continuous(breaks = seq(0, 4, by=1)) +
  labs(x = "", y = "") +
  guides(color = "none", alpha = "none") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"))

fig1 <- fig1A + fig1B + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

#------------------------------------------------------------------------------------------------
# LOO - countries
#------------------------------------------------------------------------------------------------
prob_elim_df3 <- readRDS("output/probs_scr_elimination-country-N1-2states.rds")
prior_elim3 <- unique(prob_elim_df3$prior_set)

df1b <- prob_elim_df %>% dplyr::filter(prior_set == unique(prob_elim_df$prior_set)[5]) %>% 
  mutate(country_left = "NA", data = "All countries (N=10)", LOO_set = 11, c = c*100) %>% 
  select(country_left, everything())

df3 <- prob_elim_df3 %>% dplyr::filter(prior_set == unique(prob_elim_df3$prior_set)[5]) %>% 
  mutate(data = "N-1", c = c*100) %>% 
  relocate(LOO_set, .after = last_col())

# fig 2A:
fig2A <- bind_rows(df1b, df3) %>% 
  ggplot(aes(x = c, y = prob, group = LOO_set)) +
  geom_line(aes(color = as.factor(data)), linewidth = 0.6) +
  annotate("text", label = "Ethiopia", x = 2.9, y = 0.97, size = 3) +
  annotate("text", label = "Malawi", x = 3.3, y = 0.15, size = 3) +
  scale_color_manual(values = c("red", "grey")) +
  coord_cartesian(xlim = c(0, 6)) +
  scale_x_continuous(breaks = seq(0, 6, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "SCR per 100 person-years", y = "Posterior probability", tag = "B", color = "Data:") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"))

#------------------------------------------------------------------------------------------------
# to find the discrepant curves/countries:
unique(df3$country_left)
df3 %>% dplyr::filter(!country_left == "Malawi") %>% # change here
  bind_rows(., df1b) %>% 
  ggplot(aes(x = c, y = prob, group = LOO_set)) +
  geom_line(aes(color = as.factor(country_left)), linewidth = 0.5) +
  labs(x = "", y = "Posterior probability", color = "") +
  coord_cartesian(xlim = c(0, 6)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "", y = "", color = "Country") +
  theme_minimal() 

#------------------------------------------------------------------------------------------------
# Summary of posterior probabilities and SCR:
#------------------------------------------------------------------------------------------------
scr_df11 <- bind_rows(df1b, df3) %>% filter(prob >= 0.5 & prob < 0.51) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.50") %>% as.data.frame()

scr_df12 <- bind_rows(df1b, df3) %>% filter(prob >= 0.6 & prob < 0.61) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.60") %>% as.data.frame() 

scr_df13 <- bind_rows(df1b, df3) %>% filter(prob >= 0.65 & prob < 0.66) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.65") %>% as.data.frame()

scr_df14 <- bind_rows(df1b, df3) %>% filter(prob >= 0.7 & prob < 0.71) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.70") %>% as.data.frame()

scr_df15 <- bind_rows(df1b, df3) %>% filter(prob >= 0.75 & prob < 0.76) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.75") %>% as.data.frame()

scr_df16 <- bind_rows(df1b, df3) %>% filter(prob >= 0.8 & prob < 0.81) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.80") %>% as.data.frame()

scr_df17 <- bind_rows(df1b, df3) %>% filter(prob >= 0.85 & prob < 0.86) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.85") %>% as.data.frame()

scr_df18 <- bind_rows(df1b, df3) %>% filter(prob >= 0.9 & prob < 0.91) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.90") %>% as.data.frame()

scr_df19 <- bind_rows(df1b, df3) %>% filter(prob >= 0.95 & prob < 0.96) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = "0.95") %>% as.data.frame()

scr_df20 <- bind_rows(df1b, df3) %>% filter(prob >= 0.99) %>% group_by(LOO_set) %>% 
  reframe(min = round(min(c), digits = 2), max = round(max(c), digits = 2)) %>% 
  mutate(Probability = ">0.99") %>% as.data.frame()

scr_table_prob2 <- bind_rows(scr_df11, scr_df12, scr_df13, scr_df14, scr_df15, scr_df16, scr_df17, 
                             scr_df18, scr_df19, scr_df20) %>% 
  mutate(Data = case_when(LOO_set == 11 ~ "All countries (N=10)", TRUE ~ "N-1"))


# jackknife analysis:
scr_table_prob_means2 <- scr_table_prob2 %>% dplyr::filter(!Data == "All countries (N=10)") %>% 
  group_by(Probability) %>% 
  summarise(mean_scr = mean(max),
            mean_scr2 = exp(log(mean(max)))) %>% as.data.frame()

fig2B <- ggplot() +
  geom_point(data=scr_table_prob2, 
             aes(x = max, y = factor(Probability, levels = unique(Probability)), color = as.factor(Data), 
                 alpha = as.factor(Data)), position = position_jitter(height = 0.01), size = 2.5) +
  geom_point(data=scr_table_prob_means2, 
             aes(x = mean_scr, y = factor(Probability, levels = unique(Probability))), color = "black",
             shape = 13, size = 2) +
  scale_color_manual(values = c("All countries (N=10)"="red","N-1"="grey")) +
  scale_alpha_manual(values = c("All countries (N=10)"=1,"N-1"=0.4)) +
  coord_cartesian(xlim = c(0,4)) +
  scale_x_continuous(breaks = seq(0, 4, by=1)) +
  labs(x = "SCR per 100 person-years", y = "") +
  guides(color = "none", alpha = "none") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"))

fig2 <- fig2A + fig2B + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

fig <- (fig1 / fig2) + plot_annotation(theme = theme(plot.title = element_text(size = 12),
                                                     plot.subtitle = element_text(size = 11),
                                                     plot.caption = element_text(hjust = 0.5)))

ggsave(filename = here("output/", "Fig4.png"), fig, width = 8, height = 8, dpi = 300)

