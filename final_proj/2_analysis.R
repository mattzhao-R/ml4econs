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



## functions ----

make_table <- function(reg,title){
  stargazer(reg,
            type='latex',digits=3, column.sep.width = "5pt",
            dep.var.labels.include = F,
            title=title)
}

## preliminary analysis ----


### regressions ----

ols_lprice_dist_to_ref_pre_expl <-
  lm(lprice ~ dist_to_refinery,
     data = df %>% filter(Date > my('Feb 2015')))

make_table(ols_lprice_dist_to_ref_pre_expl,
           'Log Price vs Distance to Refinery')


##
