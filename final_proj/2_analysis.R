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
  mutate(Date = ymd(Date))


## variable column extraction ----

mth_fe <- paste(colnames(df)[str_detect(colnames(df),'fe_')],
                collapse = ' + ')
mth.aftexpl_fe <- paste(colnames(df)[str_detect(colnames(df),'fe.aftexpl_')],
                collapse = ' + ')
mth.aftexpl.dist_to_ref_fe <- 
  paste(colnames(df)[str_detect(colnames(df),'fe.aftexpl.dist_to_ref_')],
        collapse = ' + ')


## functions ----



## preliminary analysis ----


### regressions ----

ols_lprice_dist_to_ref_pre_expl <-
  lm(lprice ~ dist_to_refinery,
     data = df %>% filter(Date > my('Feb 2015')))

stargazer(ols_lprice_dist_to_ref_pre_expl,
          type='latex',digits=3, column.sep.width = "5pt",
          dep.var.labels.include = F,
          title='Log Price vs Distance to Refinery')

ols_base <- lm(lsales ~ lprice, data = df)

base_formula <- paste(
  'lsales',
  X,
  '| lprice ~',
  'aftexpl + aftexpl.dist_to_ref'
  )

tsls_base <- feols(as.formula(base_formula),
                   data = df)
make_table(ols_lprice_dist_to_ref_pre_expl,
           'Base')


##
