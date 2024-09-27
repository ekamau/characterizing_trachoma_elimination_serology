#-------------------------------
# 0-trachoma-sero-meta-RUNALL.sh
#
# run all of the analysis scripts
# in batch mode to save log files
#-------------------------------

R CMD BATCH 1-trachoma-sero-thresholds-estimate-seroprev.R 
R CMD BATCH 2-trachoma-sero-thresholds-estimate-scr-sir.R
R CMD BATCH 3-trachoma-sero-thresholds-estimate-scr-sis.R
R CMD BATCH 4-trachoma-sero-thresholds-pool-seroprev.R
R CMD BATCH 5-trachoma-sero-thresholds-pool-scr.R
Rscript -e "rmarkdown::render('6-trachoma-sero-thresholds-scr-agerange-comparison.Rmd')"
Rscript -e "rmarkdown::render('7-trachoma-sero-thresholds-scr-method-comparison.Rmd')"
