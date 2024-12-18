#----------------------------------
# Characterizing trachoma elimination using serology
#
# Download .rds versions of public datasets from osf.io 
# and save them in the local working repository:
# ~/trachoma-serology/data
#
# data are available here: https://osf.io/ykjc4/
#----------------------------------


#----------------------------------
# preamble
#----------------------------------
library(here)

# source configuration file
source(here("R","0-config.R"))

# https://osf.io/9p3jn
df_v4 <- osf_retrieve_file("9p3jn") %>%
  osf_download(path=here("data"), conflicts = "overwrite", progress = TRUE)

#----------------------------------
# Session info
#----------------------------------
sessionInfo()
