#-----------------------------------------
# Plot age-seroprevalence curves - Fig 1
# Contributors: Ben Arnold ben.arnold@ucsf.edu
#-----------------------------------------

#-----------------------------------------
# Setup
#-----------------------------------------

#------------------------------
# Source project config file and project functions
#------------------------------
library(here)
source(here("R/0-config.R"))
source(here("R/0-functions.R"))

#-----------------------------------------
# read the individual level dataset
#-----------------------------------------
ind_df <- read_rds("trachoma_serology_public_data_indiv_devel.rds")) 

#-----------------------------------------
# rename id variables to remove "public" for programming ease
# restrict to EUs included in the analysis
#-----------------------------------------
ind_df2 <- ind_df %>%
  filter(!eu_name %in% unpublished_EUs) %>%
  rename(cluster_id = cluster_id_public,
         household_id = household_id_public,
         individual_id = individual_id_public) 

#-----------------------------------------
# create a categorical variable to group EUs: 
# the make_trachoma_cat() function is a shared function in 0-functions.R
#-----------------------------------------
ind_df3 <- make_trachoma_cat(ind_df2) 

#-----------------------------------------
# Age-seroprevalence for 1-9 year olds
# loop over EUs and estimate age Pgp3 seroprev using a spline
#-----------------------------------------
eus <- unique(ind_df3$eu_name)

ageserop_1to9 <- foreach(eui = eus, .combine = rbind) %do% {
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3_pos))
  gamfit <- mgcv::gam(pgp3_pos ~ s(age_years,bs="cr",k=3), data=di,family=binomial(link="logit"))
  # get predicted values
  dpred <- di %>%
    select(study_id,country,district,eu,eu_desc,eu_name,trachoma_cat, location_name,age_years) %>%
    distinct()
  dpred$logodds <- predict(gamfit, newdata = dpred, type = "link", se.fit = FALSE)
  
  # transform the log-odds into probability (proportion)
  dpred$serop <- plogis(dpred$logodds)
  
  return(dpred)
  
}

ageserop_1to9 <- ageserop_1to9 %>%
  arrange(study_id, eu_name, age_years)

#-----------------------------------------
# format a new trachoma category variable that includes a carriage return for a
# condensed plot, only first 3 cats!
#-----------------------------------------
ageserop_1to9 <- ageserop_1to9 %>%
  mutate(trachoma_cat_plot = factor(trachoma_cat, levels = levels(ind_df3$trachoma_cat), 
                                    labels = c("Action needed","Action not\nneeded","Unclassified"))
         )


#-----------------------------------------
# plot age-seroprevalence, faceted by trachoma category
#-----------------------------------------
tcats <- c("Action needed","Action not\nneeded","Unclassified")
pcols <- c("#A42820", "cadetblue",  "black") 
d_shade_box <- data.frame(minX = rep(1,3), 
                          maxX = rep(5,3),
                          trachoma_cat_plot = factor(tcats,levels=tcats))

plot_ageseroprev <- ggplot(data = ageserop_1to9, aes(x = age_years, y = serop, group = eu_name, color = trachoma_cat_plot)) +
  geom_rect(data=d_shade_box, inherit.aes = FALSE, aes(xmin = minX, xmax = maxX, ymin = 0, ymax = 1), alpha=0.1) +
  geom_line(alpha=0.6) +
  scale_x_continuous(breaks = seq(1,9,by=1), labels=seq(1,9,by=1)) +
  scale_y_continuous(breaks = seq(0,0.9,by=0.1), labels=seq(0,90,by=10)) +
  coord_cartesian(ylim = c(0,0.8)) +
  scale_color_manual(values = pcols) +
  labs(x="Age in years", y = "Pgp3 IgG seroprevalence (%)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor.x = element_blank(), 
    strip.text = element_text(size=12)) +
  facet_grid(.~trachoma_cat_plot)
  

plot_ageseroprev

# count number of eus per category
ageserop_1to9 %>% group_by(trachoma_cat) %>% summarize(count = n_distinct(eu_name))
ageserop_1to9 %>% group_by(trachoma_cat_plot) %>% summarize(count = n_distinct(eu_name))

ggsave(file = here("output/","Fig1.png"), plot_ageseroprev, device = "png",
       width=180, height = 100, units = "mm")


#-----------------------------------------
# Session Info
#-----------------------------------------
sessionInfo()

