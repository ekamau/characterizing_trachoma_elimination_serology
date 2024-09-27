#-----------------------------------------
# 3-estimate-scr-sis.R
# Estimate the SCR using a reversible catalytic model (SIS)
# assume a fixed seroreversion rate of 6 per 100 child-years
# based on estiates from the Kongwa, Tanzania study
#
# Contributors: Ben Arnold ben.arnold@ucsf.edu
#-----------------------------------------

#-----------------------------------------
# Preamble
#-----------------------------------------
library(here)
source(here("R/0-config.R"))
source(here("R/0-functions.R"))

#-----------------------------------------
# read the individual level dataset
#-----------------------------------------
ind_df <- read_rds(here(box_data_path,"/public-devel/trachoma_serology_public_data_indiv_devel.rds")) 

#-----------------------------------------
# rename id variables to remove "public" for programming ease
# restrict to EUs included in the analysis
#-----------------------------------------
ind_df2 <- ind_df %>%
  filter(!eu_name %in% unpublished_EUs) %>%
  rename(cluster_id = cluster_id_public,
         household_id = household_id_public,
         individual_id = individual_id_public) 

dim(ind_df2)

#-----------------------------------------
# create a categorical variable to group EUs into three levels of transmission:
# the remaining EUs that cannot be classified are set to "Unclassified"
# 
# make_trachoma_cat() function is a shared function in 0-functions.R
#-----------------------------------------
ind_df3 <- make_trachoma_cat(ind_df2) 

# limit to children 1-5 years
ind_df3 <- ind_df3 %>% filter(age_years <=5)

#-----------------------------------------
# Estimate SCRs with SIS model limit to 1-5 years old
#-----------------------------------------
#---------------------------------
# loop over EUs and estimate the SCR using a reversible catalytic model. 
# Use a bootstrap, resampling clusters, to estimate the 95% CI the estimate_scr_sis() function
# is in the 0-functions.R script
#---------------------------------
eus <- unique(ind_df3$eu_name)

set.seed(147123)
sis_ests <- foreach(eui = eus, .combine = rbind) %do% {
  cat("\n\n----------------------------------------\n",
      "Estimate for:", eui,
      "\n----------------------------------------\n") 
  di <- ind_df3 %>% filter(eu_name == eui, !is.na(pgp3_pos))
  scr_est <- estimate_scr_sis(serop=di$pgp3_pos, agey=di$age_years, id=di$cluster_id,
                              srr=0.06, variance="bootstrap", nboots = 100)
  res <- di %>%
    group_by(study_id,country,eu_name,trachoma_cat) %>%
    summarize(npos=sum(pgp3_pos), nchild=n(), .groups="drop") %>%
    mutate(ncluster=length(unique(di$cluster_id)),
           scr=scr_est$scr,
           scr_min95=scr_est$scr_min95,scr_max95=scr_est$scr_max95,
           scr_ci_method=scr_est$scr_ci_method,
           sis_srr = scr_est$srr,
           sis_err_warn_msg=scr_est$sis_err_warn_msg)
  return(res)
  
}

#-----------------------------------------
# Plot estimates
#-----------------------------------------
#---------------------------------
# arrange estimates by SCR
#---------------------------------
sis_ests2 <- sis_ests %>%
  mutate(eu_name = factor(eu_name,levels=sis_ests$eu_name[order(sis_ests$scr)]))

#---------------------------------
# plot estimates
#---------------------------------
plot_sis_ests <- ggplot(data=sis_ests2,aes(x=eu_name)) +
  geom_errorbar(aes(ymin=scr_min95*100, ymax=scr_max95*100),width=0.5,color=cbpal[4]) +
  geom_point(aes(y=scr*100),pch=19,color=cbpal[4]) +
  scale_y_continuous(breaks=seq(0,20,by=2)) +
  labs(x="",y="Seroconversion Rate per 100 child-years") +
  coord_flip(ylim=c(0,18)) +
  theme_minimal()

#-----------------------------------------
# save estimates for later use
#-----------------------------------------

write_rds(sis_ests, file=here("output","trachoma-sero-thresholds-scr-estimates-sis-1to5.rds"))
write_csv(sis_ests, file=here("output","trachoma-sero-thresholds-scr-estimates-sis-1to5.csv"))

#-----------------------------------------
# Session Info
#-----------------------------------------
sessionInfo()

