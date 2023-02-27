# Last updated: Feb 27, 2023

library(MASS)
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
library(plm)
source('3_estimators.R')

# ddir_json <- 
ddir_matt <- './data/'
ddir <- ddir_matt

df <- read.csv(file.path(ddir, 'final_data.csv'))[,-1] %>%
  mutate(date = ymd(date))


## variable column extraction ----

### time fe + instrument int ----
indiv_mthyr_fe <- colnames(df)[str_detect(colnames(df),'mthyrfe_')]
mthyr_fe <- paste(indiv_mthyr_fe,
                collapse = ' + ')
indiv_mthyr.aftexpl <- colnames(df)[str_detect(colnames(df),'mthyr.aftexplmth')]
mthyr.aftexpl <- paste(indiv_mthyr.aftexpl,
                collapse = ' + ')
indiv_mthyr.aftexpl.dist_to_ref <- colnames(df)[str_detect(colnames(df),'mthyr.aftexpl.dist_to_ref')]
mthyr.aftexpl.dist_to_ref <- 
  paste(indiv_mthyr.aftexpl.dist_to_ref,
        collapse = ' + ')

### county fe + instrument int ----
indiv_cnty_fe <- colnames(df)[str_detect(colnames(df), "cntyfe_")]
cnty_fe <- paste(indiv_cnty_fe,
                collapse = ' + ')
indiv_cnty.aftexpl <- colnames(df)[str_detect(colnames(df),'cntyfe.aftexpl_')]
cnty.aftexpl <- paste(indiv_cnty.aftexpl,
                     collapse = ' + ')

## functions ----



## preliminary analysis ----



### regressions ----

#### see if distance to refinery has relevance on its own ----
ols_lprice_dist_to_ref_pre_expl <-
  lm(lprice ~ dist_to_refinery,
     data = df %>% filter(date > my('Feb 2015')))

stargazer(ols_lprice_dist_to_ref_pre_expl,
          type='latex',digits=3, column.sep.width = "5pt",
          dep.var.labels.include = F,
          title='Log Price vs Distance to Refinery')

# ols_base <- feols(lsales ~ lprice, 
#                   fixef = c(indiv_mthyr_fe,indiv_cnty_fe),
#                   data = df)
# 
# etable(ols_base,
#        tex = T,
#        title = 'Baseline Regressions')


#### saturated first stage ----
fm <- paste(
  'lprice ~ ',
  paste(
    mthyr.aftexpl,
    mthyr.aftexpl.dist_to_ref,
    'aftexpl + aftexpl.dist_to_ref',
    sep = '+')
)
lm_old_fstg <- feols(
  as.formula(fm),
  fixef = c(indiv_mthyr_fe,indiv_cnty_fe),
  data = df
)
fm2 <- 'lprice ~ aftexpl + aftexpl.dist_to_ref'
lm_new_fstg <- feols(
  as.formula(fm2),
  fixef = c(indiv_mthyr_fe,indiv_cnty_fe),
  data = df
)
fm3 <- paste(
  'lprice ~ ',
  paste(
    mthyr.aftexpl,
    'aftexpl + aftexpl.dist_to_ref',
    sep = '+')
)
lm_new_fstg_wmthyr.aftexpl <- feols(
  as.formula(fm3),
  fixef = c(indiv_mthyr_fe,indiv_cnty_fe),
  data = df
)

etable(lm_old_fstg,lm_new_fstg,lm_new_fstg_wmthyr.aftexpl,
       tex = T,
       title = 'First Stages',
       fitstat = 'f',
       style.tex = style.tex('aer'),
       keep = "!mthyr")


# checking for interaction significance in first stage (ovb)
## regress aftexpl.county on U (fstg residuals) 

betapvals_aftexpl.cnty_on_U <- c()
U <- lm_new_fstg$residuals
# beta <- c()

for (i in 1:length(indiv_cnty.aftexpl)) {
  y <- indiv_cnty.aftexpl[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  # temp_beta <- temp_lm$coefficients['U']
  # beta <- c(beta,
  #           temp_beta)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_aftexpl.cnty_on_U <- c(betapvals_aftexpl.cnty_on_U,
                                   beta_pval)
}

sum(betapvals_aftexpl.cnty_on_U < 0.01) / length(betapvals_aftexpl.cnty_on_U)


## regress mthyr.aftexpl, mthyr.aftexpl.dist_to_ref on U

betapvals_mthyr.aftexpl_on_U <- c()
U <- lm_new_fstg$residuals

for (i in 1:length(indiv_mthyr.aftexpl)) {
  y <- indiv_mthyr.aftexpl[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_mthyr.aftexpl_on_U <- c(betapvals_mthyr.aftexpl_on_U,
                                   beta_pval)
}

sum(betapvals_mthyr.aftexpl_on_U < 0.01) / length(betapvals_mthyr.aftexpl_on_U)


betapvals_mthyr.aftexpl.dist_to_expl_on_U <- c()
U <- lm_new_fstg$residuals

for (i in 1:length(indiv_mthyr.aftexpl.dist_to_ref)) {
  y <- indiv_mthyr.aftexpl.dist_to_ref[i]
  temp_fm <- paste(y,'~ U')
  temp_lm <- lm(as.formula(temp_fm),data=df)
  
  beta_pval <- summary(temp_lm)$coefficients[2,4]
  betapvals_mthyr.aftexpl.dist_to_expl_on_U <- c(betapvals_mthyr.aftexpl.dist_to_expl_on_U,
                                    beta_pval)
}

sum(betapvals_mthyr.aftexpl.dist_to_expl_on_U < 0.01) / length(betapvals_mthyr.aftexpl.dist_to_expl_on_U)


### Estimators (OLS, TSLS, JIVE, RJIVE, Post-Lasso) ----

y <- data.matrix(df$lsales)
n <- length(y)
X <- cbind(rep(1,n),data.matrix(df %>%
                   select(lprice,
                          contains('mthyrfe_'),contains('cntyfe_')))
           )
Z <- cbind(rep(1,n),data.matrix(df %>%
                   select(aftexpl,aftexpl.dist_to_ref,
                          contains('mthyrfe_'),contains('cntyfe_')))
           )

#### OLS ----
fm_ols <- 'lsales ~ lprice'
ols_est <- feols(
  as.formula(fm_ols),
  fixef = c(indiv_mthyr_fe,indiv_cnty_fe),
  data = df
)

#### TSLS ----
# X_hat <- lm_new_fstg$fitted.values
# fm_tsls <- paste(
#   'lsales ~ ',
#   'lprice ~',
#   paste(
#     mthyr_fe,
#     cnty_fe,
#     'aftexpl + aftexpl.dist_to_ref',
#     sep = '+')
# )
fm_tsls <- paste(
  'lsales ~',
  mthyr_fe,
  cnty_fe,
  '| lprice ~', 
  'aftexpl + aftexpl.dist_to_ref')
tsls_est <- feols(
  as.formula(fm_tsls),
  # fixef = c(indiv_mthyr_fe,indiv_cnty_fe),
  data = df
)
# tsls.est(y,X,Z,SE=T)



#### jive ----
# from estimators script
jive.est(y,X,Z,SE=T)


#### rjive ----



#### Post Lasso ----



#### Table ----

etable(ols_est,tsls_est,
       tex = T,
       title = 'Main Results',
       fitstat = 'f',
       style.tex = style.tex('aer')
       )

### graphs ----

# Dingel Plot (looks bad)
ggplot(data = df, 
       mapping = aes(
         x = aftexpl + aftexpl.dist_to_ref,
         y = lprice),
       ) +
  geom_point()
  


