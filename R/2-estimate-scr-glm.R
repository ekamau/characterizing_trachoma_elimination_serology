#-----------------------------------------
# 2-estimate-scr-glm.R
#
# Estimate the SCR using a catalytic model with no sero-reversion.  
# Estimate SCR using GLM with a complementary log-log link and robust SEs
#
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
ind_df <- read_rds("/trachoma_serology_public_data_indiv_devel.rds")

#-----------------------------------------
# rename id variables to remove "public" for programming ease
# restrict to EUs included in the analysis
#-----------------------------------------
ind_df2 <- ind_df %>%
  filter(!eu_name %in% unpublished_EUs) %>%
  rename(cluster_id = cluster_id_public,
         household_id = household_id_public,
         individual_id = individual_id_public) 

dim(ind_df2) # n=63911

#-----------------------------------------
# create a categorical variable to group EUs into three levels of transmission:
# the remaining EUs that cannot be classified are set to "Unclassified"
# 
# the make_trachoma_cat() function is a shared function in 0-trachoma-sero-thresholds-functions.R
#-----------------------------------------
ind_df3 <- make_trachoma_cat(ind_df2) 

#-----------------------------------------
# Estimate SCRs with GLM model 1 - 9 year olds
#-----------------------------------------
#---------------------------------
# exclude TAITU, PRET, MORDOR because they only included children <5
#---------------------------------
ind_df1to9 <- ind_df3 %>% filter(!study_id %in% c("TAITU2018","PRET2013","MORDOR2015"))

#---------------------------------
# loop over EUs and estimate the SCR using a GLM with complementary log-log link
# and robust SEs the estimate_scr_glm() function is in the 0-functions.R script
#---------------------------------
eus <- unique(ind_df1to9$eu_name)

scr_ests_1to9 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df1to9 %>% filter(eu_name == eui, !is.na(pgp3_pos))
  scr_est <- estimate_scr_glm(serop=di$pgp3_pos, agey=di$age_years, id=di$cluster_id,
                              variance="robust")
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           scr=scr_est$scr,
           logscr_se=scr_est$logscr_se,
           scr_min95=scr_est$scr_min95,scr_max95=scr_est$scr_max95,
           scr_ci_method=scr_est$scr_ci_method)
  return(res)
  
}

#-----------------------------------------
# Plot estimates
#-----------------------------------------
#---------------------------------
# arrange estimates by SCR
#---------------------------------
scr_ests_1to9 <- scr_ests_1to9 %>%
  mutate(eu_name = factor(eu_name,levels=scr_ests_1to9$eu_name[order(scr_ests_1to9$scr)]))

plot_scr_ests_1to9 <- ggplot(data=scr_ests_1to9,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=scr_min95*100, ymax=scr_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=scr*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,22,by=2)) +
  labs(x="",y="Seroconversion Rate per 100 child-years") +
  coord_flip(ylim=c(0,22)) +
  theme_minimal()


#-----------------------------------------
# Estimate SCRs with SIR model - 1 - 5 year olds
#-----------------------------------------
#---------------------------------
# loop over EUs and estimate the SCR using a GLM with complementary log-log link
# and robust SEs the estimate_scr_glm() function is in the 0-functions.R script
#---------------------------------
eus <- unique(ind_df3$eu_name)
scr_ests_1to5 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3_pos), age_years <= 5)
  scr_est <- estimate_scr_glm(serop=di$pgp3_pos, agey=di$age_years, id=di$cluster_id,
                              variance="robust")
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           scr=scr_est$scr,
           logscr_se=scr_est$logscr_se,
           scr_min95=scr_est$scr_min95,scr_max95=scr_est$scr_max95,
           scr_ci_method=scr_est$scr_ci_method)
  return(res)
  
}

#-----------------------------------------
# Plot estimates
#-----------------------------------------
#---------------------------------
# arrange estimates by SCR
#---------------------------------
scr_ests_1to5 <- scr_ests_1to5 %>%
  mutate(eu_name = factor(eu_name,levels=scr_ests_1to5$eu_name[order(scr_ests_1to5$scr)]))

plot_scr_ests_1to5 <- ggplot(data=scr_ests_1to5,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=scr_min95*100, ymax=scr_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=scr*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,24,by=4)) +
  labs(x="",y="Seroconversion Rate per 100 child-years") +
  coord_flip(ylim=c(0,24)) +
  theme_minimal()


#-----------------------------------------
# Estimate SCRs with SIR model - 1 - 3 year olds
#-----------------------------------------
#---------------------------------
# loop over EUs and estimate the SCR using a GLM with complementary log-log link
# and robust SEs the estimate_scr_glm() function is in the 0-functions.R script
#---------------------------------
eus <- unique(ind_df3$eu_name)
scr_ests_1to3 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3_pos), age_years <= 3)
  scr_est <- estimate_scr_glm(serop=di$pgp3_pos,
                              agey=di$age_years,
                              id=di$cluster_id,
                              variance="robust"
  )
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           scr=scr_est$scr,
           logscr_se=scr_est$logscr_se,
           scr_min95=scr_est$scr_min95,scr_max95=scr_est$scr_max95,
           scr_ci_method=scr_est$scr_ci_method)
  return(res)
  
}

#-----------------------------------------
# Plot estimates
#-----------------------------------------
#---------------------------------
# arrange estimates by SCR
#---------------------------------
scr_ests_1to3 <- scr_ests_1to3 %>%
  mutate(eu_name = factor(eu_name,levels=scr_ests_1to3$eu_name[order(scr_ests_1to3$scr)]))

plot_scr_ests_1to3 <- ggplot(data=scr_ests_1to3,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=scr_min95*100, ymax=scr_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=scr*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,30,by=2)) +
  labs(x="",y="Seroconversion Rate per 100 child-years") +
  coord_flip(ylim=c(0,30)) +
  theme_minimal()


#-----------------------------------------
# save estimates for later use
#-----------------------------------------
write_rds(scr_ests_1to9, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-1to9.rds"))
write_csv(scr_ests_1to9, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-1to9.csv"))

write_rds(scr_ests_1to5, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-1to5.rds"))
write_csv(scr_ests_1to5, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-1to5.csv"))

write_rds(scr_ests_1to3, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-1to3.rds"))
write_csv(scr_ests_1to3, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-1to3.csv"))


#-----------------------------------------
# Session Info
#-----------------------------------------
sessionInfo()

