#-----------------------------------------
# 4-estimate-seroprev-scr-pgp3ct694.R
#
# Estimate the seroprevalence and SCR at the EU level for Pgp3 + CT694
# Use robust SEs for 95% CIs
#
# Contributors: Ben Arnold ben.arnold@ucsf.edu
#-----------------------------------------

#-----------------------------------------
# Source project config file
# and project functions
#-----------------------------------------
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
ind_df <- ind_df %>%
  filter(!eu_name %in% unpublished_EUs) %>%
  rename(cluster_id = cluster_id_public,
         household_id = household_id_public,
         individual_id = individual_id_public) 

#-----------------------------------------
# limit to studies and samples that measured both Pgp3 and CT694:
# This omits Gambia2014, Malawi2014, Kongwa2018, Solomon2015, Vanuatu2016
# and some samples from the other studies
#-----------------------------------------
table(ind_df$study_id, ind_df$pgp3ct694_pos, useNA = "ifany")
ind_df2 <- ind_df %>% filter(!is.na(pgp3ct694_pos))

#-----------------------------------------
# create a categorical variable to group EUs:
# make_trachoma_cat() function is a shared function 0-functions.R
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
# loop over EUs and estimate mean seroprev using a GLM with
# robust SEs to account for clustering
#-----------------------------------------
eus <- unique(ind_df1to9$eu_name)

seroprev_ests_1to9 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df1to9 %>% filter(eu_name == eui, !is.na(pgp3ct694_pos))
  seroprev_est <- estimate_seroprev_glm(serop=di$pgp3ct694_pos, agey=di$age_years,
                                        id=di$cluster_id, variance="robust")
  
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3ct694_pos), nchild=n(), .groups="drop") %>%
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
  scale_y_continuous(breaks=seq(0,50,by=10)) +
  labs(x="",y="Seroprevalence (%)") +
  coord_flip(ylim=c(0,50)) +
  theme_minimal()

#-----------------------------------------
# Estimate seroprevalence 1-5y
#-----------------------------------------
#-----------------------------------------
# loop over EUs and estimate mean seroprev using a GLM with
# robust SEs to account for clustering
#-----------------------------------------
eus <- unique(ind_df3$eu_name)

seroprev_ests_1to5 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3ct694_pos), age_years <= 5 )
  seroprev_est <- estimate_seroprev_glm(serop=di$pgp3ct694_pos, agey=di$age_years,
                                        id=di$cluster_id, variance="robust")
  
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3ct694_pos), nchild=n(), .groups="drop") %>%
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
  scale_y_continuous(breaks=seq(0,50,by=10)) +
  labs(x="",y="Seroprevalence (%)") +
  coord_flip(ylim=c(0,50)) +
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
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3ct694_pos), age_years <= 3 )
  seroprev_est <- estimate_seroprev_glm(serop=di$pgp3ct694_pos, agey=di$age_years,
                                        id=di$cluster_id, variance="robust")
  
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3ct694_pos), nchild=n(), .groups="drop") %>%
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
seroprev_ests_1to3_2 <- seroprev_ests_1to3 %>%
  mutate(eu_name = factor(eu_name,levels=seroprev_ests_1to3$eu_name[order(seroprev_ests_1to3$seroprev)]))

plot_seroprev_ests_1to3 <- ggplot(data=seroprev_ests_1to3_2,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=seroprev_min95*100, ymax=seroprev_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=seroprev*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,50,by=10)) +
  labs(x="",y="Seroprevalence (%)") +
  coord_flip(ylim=c(0,50)) +
  theme_minimal()


#-----------------------------------------
# for EU estimates with 0 samples that were seropositive by both Pgp3 and CT694
# estimate the upper 95% CI using an exact binomial distribution
#-----------------------------------------
seroprev_ests_1to9_b <- seroprev_ests_1to9 %>%
  mutate(seroprev = ifelse(npos==0, 0, seroprev),
         seroprev_min95 = ifelse(npos==0,0,seroprev_min95)) %>%
  rowwise() %>%
  mutate(seroprev_max95 = ifelse(npos==0, binom.test(x=npos,n=nchild)$conf.int[2],seroprev_max95),
         seroprev_ci_method = ifelse(npos==0, "Exact binomial", seroprev_ci_method))

seroprev_ests_1to5_b <- seroprev_ests_1to5 %>%
  mutate(seroprev = ifelse(npos==0, 0, seroprev),
         seroprev_min95 = ifelse(npos==0,0,seroprev_min95)) %>%
  rowwise() %>%
  mutate(seroprev_max95 = ifelse(npos==0, binom.test(x=npos,n=nchild)$conf.int[2],seroprev_max95),
         seroprev_ci_method = ifelse(npos==0, "Exact binomial", seroprev_ci_method))

seroprev_ests_1to3_b <- seroprev_ests_1to3 %>%
  mutate(seroprev = ifelse(npos==0, 0, seroprev),
         seroprev_min95 = ifelse(npos==0,0,seroprev_min95)) %>%
  rowwise() %>%
  mutate(seroprev_max95 = ifelse(npos==0, binom.test(x=npos,n=nchild)$conf.int[2],seroprev_max95),
         seroprev_ci_method = ifelse(npos==0, "Exact binomial", seroprev_ci_method))


#-----------------------------------------
# save estimates for later use
#-----------------------------------------

write_rds(seroprev_ests_1to9_b, file=here("output","trachoma-sero-thresholds-seroprev-estimates-pgp3ct694-1to9.rds"))
write_csv(seroprev_ests_1to9_b, file=here("output","trachoma-sero-thresholds-seroprev-estimates-pgp3ct694-1to9.csv"))

write_rds(seroprev_ests_1to5_b, file=here("output","trachoma-sero-thresholds-seroprev-estimates-pgp3ct694-1to5.rds"))
write_csv(seroprev_ests_1to5_b, file=here("output","trachoma-sero-thresholds-serorpev-estimates-pgp3ct694-1to5.csv"))

write_rds(seroprev_ests_1to3_b, file=here("output","trachoma-sero-thresholds-seroprev-estimates-pgp3ct694-1to3.rds"))
write_csv(seroprev_ests_1to3_b, file=here("output","trachoma-sero-thresholds-serorpev-estimates-pgp3ct694-1to3.csv"))


#-----------------------------------------
# Estimate SCRs with GLM model 1 - 9 year olds
#-----------------------------------------
#---------------------------------
# exclude TAITU, PRET, MORDOR because they only included children <5
#---------------------------------
ind_df1to9 <- ind_df3 %>% filter(!study_id %in% c("TAITU2018","PRET2013","MORDOR2015"))

#---------------------------------
# loop over EUs and estimate the SCR using a GLM with complementary log-log link and robust SEs
# estimate_scr_glm() function is in the 0-functions.R
#---------------------------------
eus <- unique(ind_df1to9$eu_name)

scr_ests_1to9 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df1to9 %>% filter(eu_name == eui, !is.na(pgp3ct694_pos))
  scr_est <- estimate_scr_glm(serop=di$pgp3ct694_pos, agey=di$age_years,
                              id=di$cluster_id, variance="robust"
  )
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3ct694_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           scr=scr_est$scr,
           logscr_se=scr_est$logscr_se,
           scr_min95=scr_est$scr_min95,
           scr_max95=scr_est$scr_max95,
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
  scale_y_continuous(breaks=seq(0,20,by=2)) +
  labs(x="",y="Seroconversion Rate per 100 child-years") +
  coord_flip(ylim=c(0,16)) +
  theme_minimal()

#-----------------------------------------
# Estimate SCRs with GLM model 1 - 5 year olds
#-----------------------------------------
#---------------------------------
# loop over EUs and estimate the SCR using a GLM with complementary log-log link and robust SEs
# estimate_scr_glm() function is in the 0-functions.R script
#---------------------------------
eus <- unique(ind_df3$eu_name)
scr_ests_1to5 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3ct694_pos), age_years <= 5)
  scr_est <- estimate_scr_glm(serop=di$pgp3ct694_pos, agey=di$age_years,
                              id=di$cluster_id, variance="robust")
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3ct694_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           scr=scr_est$scr,
           logscr_se=scr_est$logscr_se,
           scr_min95=scr_est$scr_min95,
           scr_max95=scr_est$scr_max95,
           scr_ci_method=scr_est$scr_ci_method)
  return(res)
  
}

#---------------------------------
# for EUs with zero positives
# set the SCR and its 95% CI to zero
#---------------------------------
table(scr_ests_1to5$eu_name, scr_ests_1to5$npos==0)
scr_ests_1to5 <- scr_ests_1to5 %>%
  mutate(scr = ifelse(npos==0,0,scr),
         logscr_se = ifelse(npos==0,0,logscr_se),
         scr_min95 = ifelse(npos==0,0,scr_min95),
         scr_max95 = ifelse(npos==0,0,scr_max95))

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
  scale_y_continuous(breaks=seq(0,20,by=2)) +
  labs(x="",y="Seroconversion Rate per 100 child-years") +
  coord_flip(ylim=c(0,16)) +
  theme_minimal()


#-----------------------------------------
# Estimate SCRs with GLM model 1 - 3 year olds
#-----------------------------------------
#---------------------------------
# loop over EUs and estimate the SCR using a GLM with complementary log-log link and robust SEs
# estimate_scr_glm() function is in the 0-functions.R script
#---------------------------------
eus <- unique(ind_df3$eu_name)
scr_ests_1to3 <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3ct694_pos), age_years <= 3)
  scr_est <- estimate_scr_glm(serop=di$pgp3ct694_pos, agey=di$age_years, id=di$cluster_id,
                              variance="robust")
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3ct694_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           scr=scr_est$scr,
           logscr_se=scr_est$logscr_se,
           scr_min95=scr_est$scr_min95,
           scr_max95=scr_est$scr_max95,
           scr_ci_method=scr_est$scr_ci_method)
  return(res)
  
}

#---------------------------------
# for EUs with zero positives set the SCR and its 95% CI to zero
#---------------------------------
table(scr_ests_1to3$eu_name, scr_ests_1to5$npos==0)
scr_ests_1to3 <- scr_ests_1to3 %>%
  mutate(scr = ifelse(npos==0,0,scr),
         logscr_se = ifelse(npos==0,0,logscr_se),
         scr_min95 = ifelse(npos==0,0,scr_min95),
         scr_max95 = ifelse(npos==0,0,scr_max95))


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
  scale_y_continuous(breaks=seq(0,20,by=2)) +
  labs(x="",y="Seroconversion Rate per 100 child-years") +
  coord_flip(ylim=c(0,20)) +
  theme_minimal()

#-----------------------------------------
# save estimates for later use
#-----------------------------------------
write_rds(scr_ests_1to9, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-pgp3ct694-1to9.rds"))
write_csv(scr_ests_1to9, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-pgp3ct694-1to9.csv"))

write_rds(scr_ests_1to5, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-pgp3ct694-1to5.rds"))
write_csv(scr_ests_1to5, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-pgp3ct694-1to5.csv"))

write_rds(scr_ests_1to3, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-pgp3ct694-1to3.rds"))
write_csv(scr_ests_1to3, file=here("output","trachoma-sero-thresholds-scr-estimates-glm-pgp3ct694-1to3.csv"))

#-----------------------------------------
# Session Info
#-----------------------------------------
sessionInfo()

