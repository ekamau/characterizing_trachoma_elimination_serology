#' ---
#' title: "Plot Bayesian seroprevalence estimates - Fig S3"
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
df <- read.csv(here("output/", "posterior_seroprev_samples_1to5yo.csv"))
eus <- unique(df$eu) # check eu labels!

unpublished_EUs <- c("Niger-Ilela-2022", "Niger-Malbaza-2022", "Niger-Bagaroua-2022", "DRC-Manono-2018", "DRC-Nyemba-2018")

eus_wo_PCR <- c("Sudan-El Seraif-2019", "Sudan-Saraf Omrah-2019", "Sudan-Kotom-2019",
                "Malaysia-Sabah-2015",  "Peru-Amazonia-2020", "Papua New Guinea-Mendi-2015", "Papua New Guinea-Daru-2015", "Papua New Guinea-West New Britain-2015")

df2 <- df[!(df$eu %in% c(unpublished_EUs)), ] %>% 
  group_by(eu) %>% mutate(median_prev = median(prev, na.rm = TRUE)) %>% ungroup() %>% 
  as.data.frame()
eus <- unique(df2$eu) # check eu labels!

#' Plots:
#+ warning = FALSE, message = FALSE
# assign categories for programmatic action:
df2 <- df2 %>%
  mutate(trachoma_cat = case_when(startsWith(as.character(eu), "Ghana") ~ "Action not needed",
                                  eu %in% c("Ethiopia-Alefa-2017","Niger-MORDOR/Dosso-2015","Ethiopia-Woreta Town-2017",
                                            "Togo-Anie-2017","Togo-Keran-2017","Morocco-Agdaz-2019","Morocco-Boumalne Dades-2019",
                                            "Gambia-River Regions-2014","Ethiopia-Metema-2021","Ethiopia-Woreta Town-2021",
                                            "Malawi-DHO Nkwazi-2014","Malawi-Kasisi/DHO-2014","Malawi-Luzi Kochilira-2014",
                                            "Malawi-Chapananga-2014") ~ "Action not needed",
                                  eu %in% c("Ethiopia-Andabet-2017","Ethiopia-Goncha-2019","Niger-PRET/Matameye-2013",
                                            "Ethiopia-Ebinat-2019","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016",
                                            "Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016","Tanzania-Kongwa-2018",
                                            "Tanzania-Kongwa-2013","Solomon Islands-Temotu/Rennel/Bellona-2015") ~ "Action needed",
                                  eu %in% eus_wo_PCR ~ "Unclassified",
                                  TRUE ~ as.character("Unclassified")),
         trachoma_cat = factor(trachoma_cat, levels = c("Action needed","Action not needed","Unclassified"))
  )


# combined plots:
colPal <- c("#A42820", "cadetblue")

h <- df2 %>% filter(trachoma_cat != "Unclassified") %>% 
  ggplot() +
  geom_density_ridges(aes(x = prev*100, y = fct_reorder(eu, prev*100, .desc = FALSE), group = eu, 
                          fill = trachoma_cat, color = trachoma_cat), rel_min_height = 0.01,
                      alpha = 0.6, vline_color = "black", quantile_lines = TRUE, quantiles = 2) +
  scale_x_continuous(breaks = seq(0, 60, by = 10), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 60)) +
  scale_fill_manual(name = "", values = colPal) +
  scale_color_manual(name = "", values = colPal) +
  #ggforce::facet_col(~fct_relevel(trachoma_cat, "Endemic transmission", "Near elimination of transmission", 
                                  #"Post elimination of transmission"), scales = "free_y", space = "free") +
  labs(x = "", y = "Survey") +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 11),
        axis.ticks.x = element_line(linewidth = 0.2, colour = "black"),
        axis.line.x = element_line(color="black", linewidth = 0.1),
        legend.position = "none",
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "white"))

h

j <- df2 %>% filter(trachoma_cat != "Unclassified") %>% 
  ggplot(aes(x = prev*100, group = trachoma_cat, fill = trachoma_cat)) +
  geom_density(aes(color = trachoma_cat), alpha = 0.6, trim = TRUE) +
  #annotate("text", label = "density", x = 40, y = 0.45) +
  scale_fill_manual(name = "", breaks = c("Action not needed","Action needed"), 
                    values = c("Action not needed"="cadetblue", "Action needed"="#A42820")) +
  scale_color_manual(name = "", breaks = c("Action not needed","Action needed"), 
                     values = c("Action not needed"="cadetblue", "Action needed"="#A42820")) +
  scale_x_continuous(breaks = seq(0, 60, by = 10), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "\nSeroprevalence (%)", y = "Density") +
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
B <- h + j + plot_layout(ncol = 1, heights = c(4, 1)) + plot_annotation(tag_levels = 'A')
ggsave(filename = here::here("output/", "FigS3.png"), B, width = 6, height = 8, dpi = 300)
