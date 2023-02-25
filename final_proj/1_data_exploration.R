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

# data preparation ----

## prep functions ----
make_longdf <- function(filename){
  df <- read.csv(file.path(ddir, filename)) %>%
    mutate(Date = mdy(Date),
           month = month(Date),
           year = year(Date)) %>%
    pivot_longer(cols=!c(Date,month,year),
                 names_to = 'county',
                 values_to = 'value')
    df
}

make_covariates <- function(df){
  df <- df %>%
    mutate(lprice = log(price),
           lsales = log(sales),
           quarter = quarter(Date),
           month = month(Date),
           year = year(Date),
           aftexpl = ifelse((Date > my('Feb 2015')),1,0),
           t_expl = as.numeric(Date - my('Feb 2015')),
           t_expl_sq = t_expl ** 2)
  df
}


## Data Processing - making dataframe for analysis ----

m_sales_df <- make_longdf('m_county_sales_2012_2022.csv') %>%
  rename(sales = value)
w_price_df <- make_longdf('w_county_price_idxsf_2000_2023.csv') %>%
  rename(price = value)
m_price_df <- w_price_df %>%
  group_by(year,month,county) %>%
  summarise(price = mean(price), .groups = 'keep')
dist_to_torr <- read.csv(file.path(ddir, 'dist_to_refinery.csv')) %>%
  rowwise() %>%
  mutate(county = str_replace_all(county,' ','.'))

temp <- left_join(m_price_df,m_sales_df,
                    by = c('year','month','county')) %>%
  left_join(dist_to_torr, by = 'county') %>%
  drop_na()
m_merged_df <- make_covariates(temp)

### adding indicators ----

pre_interact_df <- m_merged_df %>%
  mutate(aftexpl.dist_to_ref = dist_to_refinery * aftexpl)

county <- pre_interact_df$county
mth <- paste0('mth',pre_interact_df$month)
yr <- paste0('yr',pre_interact_df$year)
aftexpl <- pre_interact_df$aftexpl
aftexpl.dist_to_ref <- pre_interact_df$aftexpl.dist_to_ref
dist_to_torr <- pre_interact_df$dist_to_refinery

analysis_df <- 
  bind_cols(pre_interact_df,
  data.frame(i(county, ref = T)),
  data.frame(i(mth, ref = T)),
  data.frame(i(yr, ref = T)),
  data.frame(i(mth, i.yr, ref = T)),
  data.frame(i(mth, aftexpl, ref = T)),
  data.frame(i(mth, aftexpl.dist_to_ref, ref = T))
  )

write.csv(analysis_df,file.path(ddir,'final_data.csv'))


