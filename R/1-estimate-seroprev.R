#-----------------------------------------
# 1-estimate-seroprev.R

# Estimate the seroprevalence at the EU level for 1-9y olds, 1-5y olds, 1-3y olds
# Use robust SEs for 95% CIs
#
# Contributors: Ben Arnold ben.arnold@ucsf.edu
#-----------------------------------------

#-----------------------------------------
# Source project config file and project functions
#-----------------------------------------
library(here)
source(here("R/0-config.R"))
source(here("R/0-functions.R"))

#-----------------------------------------
# read the individual level dataset
#-----------------------------------------
ind_df <- read_rds("/trachoma_serology_public_data_indiv_devel.rds")
ind_df %>% dplyr::filter(!eu_name %in% unpublished_EUs) %>% nrow() # n = 63911 (1-9yo)

ind_df %>% dplyr::filter(age_years <= 5, !eu_name %in% unpublished_EUs) %>% nrow() # n = 41168 (1-5yo)

ind_df %>% dplyr::filter(age_years <= 3, !eu_name %in% unpublished_EUs) %>% nrow() # n = 24353 (1-3yo)

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
# create a categorical variable to group EUs 
# the remaining EUs that cannot be classified are set to "Unclassified"
# 
# the make_trachoma_cat() function is a shared function in 0-functions.R
#-----------------------------------------
ind_df3 <- make_trachoma_cat(ind_df2) 

#-----------------------------------------
# Estimate seroprevalence 1-9y
#-----------------------------------------
#---------------------------------
# exclude TAITU, PRET, MORDOR because they only included children <5
#---------------------------------
ind_df1to9 <- ind_df3 %>% filter(!study_id %in% c("TAITU2018","PRET2013","MORDOR2015"))

#-----------------------------------------
# loop over EUs and estimate mean seroprev using a GLM with robust SEs to account for clustering
#-----------------------------------------
eus <- unique(ind_df1to9$eu_name)

seroprev_ests_1to9 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df1to9 %>% filter(eu_name == eui, !is.na(pgp3_pos))
  seroprev_est <- estimate_seroprev_glm(serop=di$pgp3_pos, agey=di$age_years, id=di$cluster_id,
                                        variance="robust")
  
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           logodds = seroprev_est$logodds,
           logodds_se=seroprev_est$logodds_se,
           seroprev=seroprev_est$seroprev,
           seroprev_min95=seroprev_est$seroprev_min95,seroprev_max95=seroprev_est$seroprev_max95,
           seroprev_ci_method=seroprev_est$seroprev_ci_method)
  return(res)
  
}

#-----------------------------------------
# Plot estimates
#-----------------------------------------

# arrange estimates by seroprev
seroprev_ests_1to9_2 <- seroprev_ests_1to9 %>%
  mutate(eu_name = factor(eu_name,levels=seroprev_ests_1to9$eu_name[order(seroprev_ests_1to9$seroprev)]))

plot_seroprev_ests_1to9 <- ggplot(data=seroprev_ests_1to9_2,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=seroprev_min95*100, ymax=seroprev_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=seroprev*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,70,by=10)) +
  labs(x="",y="Seroprevalence (%)") +
  coord_flip(ylim=c(0,70)) +
  theme_minimal()

#-----------------------------------------
# Estimate seroprevalence 1-5y
#-----------------------------------------
#-----------------------------------------
# loop over EUs and estimate mean seroprev using a GLM with robust SEs to account for clustering
#-----------------------------------------
eus <- unique(ind_df3$eu_name)

seroprev_ests_1to5 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3_pos), age_years <= 5 )
  seroprev_est <- estimate_seroprev_glm(serop=di$pgp3_pos, agey=di$age_years, id=di$cluster_id,
                                        variance="robust")
  
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           logodds = seroprev_est$logodds,
           logodds_se=seroprev_est$logodds_se,
           seroprev=seroprev_est$seroprev,
           seroprev_min95=seroprev_est$seroprev_min95,seroprev_max95=seroprev_est$seroprev_max95,
           seroprev_ci_method=seroprev_est$seroprev_ci_method)
  return(res)
  
}

#-----------------------------------------
# Plot estimates
#-----------------------------------------
# arrange estimates by seroprev
seroprev_ests_1to5_2 <- seroprev_ests_1to5 %>%
  mutate(eu_name = factor(eu_name,levels=seroprev_ests_1to5$eu_name[order(seroprev_ests_1to5$seroprev)]))

plot_seroprev_ests_1to5 <- ggplot(data=seroprev_ests_1to5_2,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=seroprev_min95*100, ymax=seroprev_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=seroprev*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,60,by=10)) +
  labs(x="",y="Seroprevalence (%)") +
  coord_flip(ylim=c(0,60)) +
  theme_minimal()


#-----------------------------------------
# Estimate seroprevalence 1-3y
#-----------------------------------------
#-----------------------------------------
# loop over EUs and estimate mean seroprev using a GLM with robust SEs to account for clustering
#-----------------------------------------
eus <- unique(ind_df3$eu_name)

seroprev_ests_1to3 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n")
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3_pos), age_years <= 3 )
  seroprev_est <- estimate_seroprev_glm(serop=di$pgp3_pos, agey=di$age_years, id=di$cluster_id,
                                        variance="robust")

  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           logodds = seroprev_est$logodds,
           logodds_se=seroprev_est$logodds_se,
           seroprev=seroprev_est$seroprev,
           seroprev_min95=seroprev_est$seroprev_min95,seroprev_max95=seroprev_est$seroprev_max95,
           seroprev_ci_method=seroprev_est$seroprev_ci_method)
}

#-----------------------------------------
# Plot estimates
#-----------------------------------------
# arrange estimates by seroprev
seroprev_ests_1to3_2 <- seroprev_ests_1to3 %>%
  mutate(eu_name = factor(eu_name,levels=seroprev_ests_1to5$eu_name[order(seroprev_ests_1to5$seroprev)]))

plot_seroprev_ests_1to3 <- ggplot(data=seroprev_ests_1to3_2,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=seroprev_min95*100, ymax=seroprev_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=seroprev*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,50,by=10)) +
  labs(x="",y="Seroprevalence (%)") +
  coord_flip(ylim=c(0,50)) +
  theme_minimal()


#-----------------------------------------
# save estimates
#-----------------------------------------

write_rds(seroprev_ests_1to9, file=here("output","trachoma-sero-thresholds-seroprev-estimates-1to9.rds"))
write_csv(seroprev_ests_1to9, file=here("output","trachoma-sero-thresholds-seroprev-estimates-1to9.csv"))

write_rds(seroprev_ests_1to5, file=here("output","trachoma-sero-thresholds-seroprev-estimates-1to5.rds"))
write_csv(seroprev_ests_1to5, file=here("output","trachoma-sero-thresholds-serorpev-estimates-1to5.csv"))

write_rds(seroprev_ests_1to3, file=here("output","trachoma-sero-thresholds-seroprev-estimates-1to3.rds"))
write_csv(seroprev_ests_1to3, file=here("output","trachoma-sero-thresholds-serorpev-estimates-1to3.csv"))


#-----------------------------------------
# Session Info
#-----------------------------------------
sessionInfo()

