# Last updated: Feb 25, 2023

library(tidyverse)
library(lubridate)
library(haven)
library(stargazer)
library(xtable)
library(lmtest)
library(sandwich)
library(broom)
library(whitestrap)
library(fixest)
library(ivreg)
library(car)


# ddir_json <- 
ddir_matt <- './data/'
ddir <- ddir_matt

df <- read.csv(file.path(ddir, 'final_data.csv'))[,-1] %>%
  mutate(date = ymd(date))


## variable column extraction ----

### time fe + instrument int ----
indiv_mth_fe <- colnames(df)[str_detect(colnames(df),'mthfe_')]
mth_fe <- paste(indiv_mth_fe,
                collapse = ' + ')
mth.aftexpl <- paste(colnames(df)[str_detect(colnames(df),'mthfe.aftexpl_')],
                collapse = ' + ')
mth.aftexpl.dist_to_ref <- 
  paste(colnames(df)[str_detect(colnames(df),'mthfe.aftexpl.dist_to_ref_')],
        collapse = ' + ')

### county fe + instrument int ----
indiv_cnty_fe <- colnames(df)[str_detect(colnames(df), "cntyfe_")]
cnty_fe <- paste(indiv_cnty_fe,
                collapse = ' + ')
cnty.aftexpl <- paste(colnames(df)[str_detect(colnames(df),'cntyfe.aftexpl_')],
                     collapse = ' + ')

## functions ----



## preliminary analysis ----


### regressions ----

# see if distance to refinery has relevance on its own
ols_lprice_dist_to_ref_pre_expl <-
  lm(lprice ~ dist_to_refinery,
     data = df %>% filter(Date > my('Feb 2015')))

stargazer(ols_lprice_dist_to_ref_pre_expl,
          type='latex',digits=3, column.sep.width = "5pt",
          dep.var.labels.include = F,
          title='Log Price vs Distance to Refinery')

ols_base <- feols(lsales ~ lprice, data = df)

# baseline regressions
# base_formula <- paste(
#   'lsales ~ ',
#   paste(mth_fe,
#         mth.aftexpl_fe,
#         mth.aftexpl.dist_to_ref_fe,
#         sep = ' + '),
#   '| lprice ~',
#   'aftexpl + aftexpl.dist_to_ref'
#   )

tsls_base <- feols(as.formula(base_formula),
                   fixef = indiv_mth_fe,
                   data = df)
etable(ols_base,tsls_base,
       tex = T,
       title = 'Baseline Regressions')


# testing if county has effect
## saturated first stage
fm <- paste(
  'lprice ~ ',
  paste(
    mth.aftexpl,
    mth.aftexpl.dist_to_ref,
    'aftexpl + aftexpl.dist_to_ref',
    sep = '+')
)
lm_county_test <- feols(
  as.formula(fm),
  fixef = c(indiv_mth_fe,indiv_cnty_fe),
  data = df
)

etable(lm_county_test,
       tex = T,
       title = 'Test for County Effect')

stargazer(ols_base,tsls_base,
          type='latex',digits=3, column.sep.width = "5pt",
          dep.var.labels.include = F,
          model.names = T,
          title='Baseline Regressions')


### graphs ----

# Dingel Plot (looks bad)
ggplot(data = df, 
       mapping = aes(
         x = aftexpl + aftexpl.dist_to_ref,
         y = lprice),
       ) +
  geom_point()
  


