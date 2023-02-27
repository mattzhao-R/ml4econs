# Last updated: Feb 25, 2023

library(tidyverse)
library(lubridate)
library(haven)
library(stargazer)
library(xtable)
library(broom)
library(fixest)
library(whitestrap)
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
  drop_na() %>%
  filter(year <= 2018)
m_merged_df <- make_covariates(temp)

### adding FEs and interactions ----

pre_interact_df <- m_merged_df %>%
  mutate(aftexpl.dist_to_ref = dist_to_refinery * aftexpl,
         time = as.character(year * 100 + month)) %>%
  rename(date = Date) %>%
  ungroup() %>%
  dplyr::select(county,time,everything())
time_map <- setNames(seq(1,length(unique(pre_interact_df$time))), 
                     unique(pre_interact_df$time))
pre_interact_df$time <- time_map[pre_interact_df$time]


county <- pre_interact_df$county
time <- paste0('t',pre_interact_df$time)
# mth <- paste0('mth',pre_interact_df$month)
# yr <- paste0('yr',pre_interact_df$year)
aftexpl <- pre_interact_df$aftexpl
aftexpl.dist_to_ref <- pre_interact_df$aftexpl.dist_to_ref
aftexpl.dist_to_ref_i <- paste0('aftexpl.dist_to_ref',
                              pre_interact_df$aftexpl.dist_to_ref)

analysis_df <- 
  bind_cols(pre_interact_df,
  data.frame(i(county, ref = T)),
  data.frame(i(time, ref = T)),
  data.frame(i(time, aftexpl, ref = T)),
  data.frame(i(time, aftexpl.dist_to_ref, ref = T)),
  data.frame(i(aftexpl.dist_to_ref_i, ref = T)),
  # data.frame(i(mth, i.yr, ref = T)),
  data.frame(i(county, aftexpl, ref = T))
  )

### renaming columns ----

cn <- colnames(analysis_df)

#### month-year fixed effects
# mthyrfe <- cn[str_detect(cn, "mth")]
# 
# mthyrfe_pos <- match(mthyrfe, colnames(analysis_df))
# mthfe_names <- paste0('mthyrfe_',mthyrfe)
# colnames(analysis_df)[mthyrfe_pos] <- mthfe_names

# mthfe.aftexpl_pos <- match(mthfe.aftexpl, colnames(analysis_df))
# mthfe.aftexpl_names <- paste0('mthfe.aftexpl_',
#                               str_extract(mthfe.aftexpl,
#                                           "[a-z]{3}[0-9]{1,2}"))
# colnames(analysis_df)[mthfe.aftexpl_pos] <- mthfe.aftexpl_names

#### time fe + interaction w instrument

timefe_allint <- cn[str_detect(cn, "^t\\d+")]
timefe <- timefe_allint[1:83]
timefe.aftexpl <- timefe_allint[84:166]
timefe.aftexpl.dist_to_ref <- timefe_allint[167:249]


timefe_pos <- match(timefe, cn)
timefe_names <- paste0('timefe_',str_extract(timefe,"^t\\d+"))
colnames(analysis_df)[timefe_pos] <- timefe_names

timefe.aftexpl_pos <- match(timefe.aftexpl, cn)
timefe.aftexpl_names <- paste0('timefe.aftexpl_',str_extract(timefe.aftexpl,"^t\\d+"))
colnames(analysis_df)[timefe.aftexpl_pos] <- timefe.aftexpl_names

timefe.aftexpl.dist_to_ref_pos <- match(timefe.aftexpl.dist_to_ref, cn)
timefe.aftexpl.dist_to_ref_names <- paste0('timefe.aftexpl.dist_to_ref_',
                                           str_extract(timefe.aftexpl.dist_to_ref,
                                                       "^t\\d+"))
colnames(analysis_df)[timefe.aftexpl.dist_to_ref_pos] <- timefe.aftexpl.dist_to_ref_names



#### county fe + interaction w instrument

cntyfe_allint <- colnames(analysis_df)[str_detect(colnames(analysis_df), "^[A-Z]")]
cntyfe <- cntyfe_allint[1:55]
cntyfe.aftexpl <- cntyfe_allint[56:110]

cntyfe_pos <- match(cntyfe, colnames(analysis_df))
cntyfe_names <- paste0('cntyfe_',
                       str_extract_all(cntyfe,
                                       "[[:alpha:]]+(?:\\.[[:alpha:]]+)?+(?:\\.[[:alpha:]]+)?"))
colnames(analysis_df)[cntyfe_pos] <- cntyfe_names

cntyfe.aftexpl_pos <- match(cntyfe.aftexpl, colnames(analysis_df))
cntyfe.aftexpl_names <- paste0('cntyfe.aftexpl_',
                       str_extract_all(cntyfe.aftexpl,
                                       "[[:alpha:]]+(?:\\.[[:alpha:]]+)?+(?:\\.[[:alpha:]]+)?"))
colnames(analysis_df)[cntyfe.aftexpl_pos] <- cntyfe.aftexpl_names

# #### creating month-year.instrument interactions
# 
# mthyrfe_colnames <- colnames(analysis_df)[str_detect(colnames(analysis_df), "mthyrfe")]
# 
# mthyr.aftexpl_df <- analysis_df %>% 
#   select(year,month,county,aftexpl,aftexpl.dist_to_ref,
#          all_of(contains('mthyrfe_'))
#          )
# 
# mthyr.aftexpl_df[mthyrfe_colnames] <- 
#   mthyr.aftexpl_df[mthyrfe_colnames] * 
#   mthyr.aftexpl_df$aftexpl
# mthyr.aftexpl_df <- mthyr.aftexpl_df %>%
#   select(all_of(contains('mthyrfe_')))
# colnames(mthyr.aftexpl_df) <- 
#   str_replace(colnames(mthyr.aftexpl_df),
#               'mthyrfe_',
#               'mthyr.aftexpl')
# 
# mthyr.aftexpl.dist_to_ref_df <- analysis_df %>% 
#   select(year,month,county,aftexpl,aftexpl.dist_to_ref,
#          all_of(contains('mthyrfe_'))
#   )
# mthyr.aftexpl.dist_to_ref_df[mthyrfe_colnames] <- 
#   mthyr.aftexpl.dist_to_ref_df[mthyrfe_colnames] * 
#   mthyr.aftexpl.dist_to_ref_df$aftexpl.dist_to_ref
# mthyr.aftexpl.dist_to_ref_df <- mthyr.aftexpl.dist_to_ref_df %>%
#   select(all_of(contains('mthyrfe_')))
# colnames(mthyr.aftexpl.dist_to_ref_df) <- 
#   str_replace(colnames(mthyr.aftexpl.dist_to_ref_df),
#               'mthyrfe_',
#               'mthyr.aftexpl.dist_to_ref')

# final_data_df <- 
#   bind_cols(
#     analysis_df,
#     mthyr.aftexpl_df,
#     mthyr.aftexpl.dist_to_ref_df
#     )

# export for analysis
write.csv(analysis_df,file.path(ddir,'final_data.csv'))


