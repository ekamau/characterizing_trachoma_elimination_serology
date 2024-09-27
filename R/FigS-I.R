#' ---
#' title: "Plot Bayesian SCR estimates - Supp Fig I"
#' author: ekamau
#' output: html_document
#' ---
#' 
# Notes: ages 1 to 5 years
#+ setup, include=FALSE, echo=FALSE
knitr::opts_chunk$set(collapse = TRUE, echo = FALSE)

# packages
if (!require("pacman")) install.packages("pacman")
p_load(tidyverse, ggridges, reshape2, viridis, here, patchwork, ggtext)

# Read bayesian SCR estimates:
#+ warning = FALSE, message = FALSE, echo = FALSE, eval = FALSE
df <- read.csv(here("output/", "posterior_samples_scr-1to5yo.csv"))
eus <- unique(df$eu) # check eu labels!

unpublished_EUs <- c("Niger-Ilela-2022", "Niger-Malbaza-2022", "Niger-Bagaroua-2022", "DRC-Manono-2018", "DRC-Nyemba-2018")

eus_wo_PCR <- c("Sudan-El Seraif-2019","Sudan-Saraf Omrah-2019","Sudan-Kotom-2019",
                "Malaysia-Sabah-2015","Peru-Amazonia-2020","Papua New Guinea-Mendi-2015","Papua New Guinea-Daru-2015","Papua New Guinea-West New Britain-2015",
                "Togo-Anie-2017","Togo-Keran-2017","Morocco-Agdaz-2019","Morocco-Boumalne Dades-2019","Gambia-River Regions-2014")

df2 <- df[!(df$eu %in% c(unpublished_EUs)), ] %>% 
  group_by(eu) %>% mutate(median_scr = median(scr, na.rm = TRUE)) %>% 
  ungroup() %>% as.data.frame()
eus <- unique(df2$eu) # check eu labels!

#' Plots:
#+ warning = FALSE, message = FALSE
# assign categories for elimination of transmission:
df2 <- df2 %>% mutate(trachoma_cat = case_when(startsWith(as.character(eu), "Ghana") ~ "Probable interruption of transmission",
                                             eu %in% c("Ethiopia-Alefa-2017","Niger-MORDOR/Dosso-2015","Ethiopia-Woreta Town-2017","Ethiopia-Metema-2021","Ethiopia-Woreta Town-2021") ~ "Probable interruption of transmission",
                                             eu %in% c("Ethiopia-Andabet-2017","Niger-PRET/Matameye-2013","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016","Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016",
                                                            "Ethiopia-Ebinat-2019") ~ "Endemic",
                                             eu %in% eus_wo_PCR ~ "Unclassified",
                                             TRUE ~ as.character("Near interruption of transmission")),
                    trachoma_cat = factor(trachoma_cat, levels = c("Endemic", "Near interruption of transmission", 
                                                                   "Probable interruption of transmission", "Unclassified"))
)

unique(df2$trachoma_cat)
table(df2$trachoma_cat)/8400

colPal <- c("#A42820", "#F98400", "cadetblue") 
  
h <- ggplot(data = df2 %>% filter(trachoma_cat != "Unclassified")) +
  geom_density_ridges(aes(x = scr*100, y = fct_reorder(eu, scr*100, .desc = FALSE), group = eu, 
                          fill = trachoma_cat, color = trachoma_cat), rel_min_height = 0.01,
                      alpha = 0.6, vline_color = "black", quantile_lines = TRUE,
                      quantiles = 2) +
  scale_x_continuous(breaks = seq(0, 22, by = 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 22)) +
  scale_fill_manual(name = "", values = colPal) +
  scale_color_manual(name = "", values = colPal) +
  labs(x = "", y = "Survey") +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 11),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        legend.position = "none",
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill="white"))

h

j <- df2 %>% filter(trachoma_cat != "Unclassified") %>% 
  ggplot(aes(x = scr*100, group = trachoma_cat, fill = trachoma_cat)) +
  geom_density(aes(color = trachoma_cat), alpha = 0.6, trim = TRUE) +
  #annotate("text", label = "density", x = 15, y = 1.2) +
  scale_fill_manual(name = "", breaks = c("Probable interruption of transmission", "Near interruption of transmission", "Endemic"),
                    values = c("Probable interruption of transmission"="cadetblue", "Near interruption of transmission"="#F98400", "Endemic"="#A42820")) +
  scale_color_manual(name = "", breaks = c("Probable interruption of transmission", "Near interruption of transmission", "Endemic"),
                     values = c("Probable interruption of transmission"="cadetblue", "Near interruption of transmission"="#F98400", "Endemic"="#A42820")) +
  scale_x_continuous(breaks = seq(0, 22, by = 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 22)) +
  labs(x = "\nSCR per 100 person-years", y = "Density") +
  guides(fill = guide_legend(nrow = 3, byrow=TRUE)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 8),
        panel.spacing = unit(0.1, "lines"),
        axis.title = element_text(size = 11),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"),
        axis.line.x = element_line(color = "black", linewidth = 0.1),
        axis.text.y = element_blank(),
        legend.position = "bottom") 

j

A <- h + j + plot_layout(ncol = 1, heights = c(4, 1)) + plot_annotation(tag_levels = 'A')
ggsave(filename = here("output/", "FigS-I.png"), A, width = 6, height = 8, dpi = 300)
