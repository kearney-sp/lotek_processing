#### Prep ####
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(qpcR)
library(parameters)

df <- read.csv('C:/Users/sean.kearney/OneDrive - USDA/Documents/Projects/GPS_v_hetgen/data/full_grazing_time_gridded.csv')
df <- df %>%
  filter((Biomass <= 4.0 & Biomass >= -4.0) &
           (CP <= 4.0 & CP >= -4.0) &
           (PC_div <= 4.0 & PC_div >= -4.0)) %>%
  mutate(season_str = factor(season_str, levels=c('early', 'mid', 'late')),
         PC_dmt = relevel(factor(PC_dmt), ref='C3_C4_mix')) %>%
  mutate(pred_zprob = NA,
         pred_response = NA) %>%
  mutate(grazing_secs = case_when(grazing_secs < 30 ~ as.integer(0),
                                  TRUE ~ grazing_secs))

df$grz_lag_early <- NA
df$grz_lag_mid <- NA
df[df$season_str == 'late', 
   'grz_lag_mid'] <- left_join(df[df$season_str == 'late',], 
                               df[df$season_str == 'mid', c('UTM_X',
                                                            'UTM_Y', 
                                                            'Year',
                                                            'grazing_rel_freq')] %>%
                                 group_by(Year, UTM_X, UTM_Y) %>%
                                 summarise(grz_lag_mid=mean(grazing_rel_freq)),
                               by=c('UTM_X', 'UTM_Y', 'Year')) %>% 
  pull(grz_lag_mid.y)
df[df$season_str %in% c('mid', 'late'), 
   'grz_lag_early'] <- left_join(df[df$season_str %in% c('mid', 'late'),], 
                               df[df$season_str == 'early', c('UTM_X',
                                                            'UTM_Y', 
                                                            'Year',
                                                            'grazing_rel_freq')] %>%
                                 group_by(Year, UTM_X, UTM_Y) %>%
                                 summarise(grz_lag_early=mean(grazing_rel_freq)),
                               by=c('UTM_X', 'UTM_Y', 'Year')) %>%
  pull(grz_lag_early.y)

df_aic <- data.frame()
df_sub_cor <- data.frame()

avg_grazing_secs <- df %>% group_by(Year, week, Pasture) %>% 
  summarise(n = n(), grazing_secs = mean(grazing_wkly_sum)) %>% 
  mutate(expected_secs = grazing_secs/n) %>% 
  ungroup() %>%
  dplyr::select(grazing_secs) %>% 
  colMeans() %>% unname()

avg_expected_secs <- df %>% group_by(Year, week, Pasture) %>% 
  summarise(n = n(), grazing_secs = mean(grazing_wkly_sum)) %>% 
  mutate(expected_secs = grazing_secs/n) %>% 
  ungroup() %>%
  dplyr::select(expected_secs) %>% 
  colMeans() %>% unname()

#### 2016 Early ####
df_sub <- df %>% filter(Year == 2016 & season_str == 'early') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response",
                          "grz_lag_early", "grz_lag_mid"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div + 
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))

akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div + 
                                TPC_c +
                                offset(log(grazing_wkly_sum)) +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2016,
                           season='early',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2016 & 
          df_aic$season == 'early', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2016 & 
                 df_aic$season == 'early',])
}

modFNL_2016_early <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2016
  df_tmp_cor['season'] = 'early'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2016 Mid ####
df_sub <- df %>% filter(Year == 2016 & season_str == 'mid') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response",
                          "grz_lag_mid"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                     grz_lag_avg + grz_lag_early +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                       grz_lag_avg + grz_lag_early +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                         grz_lag_avg + grz_lag_early +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                        grz_lag_avg + grz_lag_early +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                          grz_lag_avg + grz_lag_early +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div +
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))



akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                                grz_lag_avg + grz_lag_early +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div +
                                TPC_c +
                                offset(log(grazing_wkly_sum)) +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2016,
                           season='mid',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2016 & 
          df_aic$season == 'mid', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2016 & 
                 df_aic$season == 'mid',])
}

modFNL_2016_mid <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2016
  df_tmp_cor['season'] = 'mid'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2016 Late ####
df_sub <- df %>% filter(Year == 2016 & season_str == 'late') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                     grz_lag_avg + grz_lag_early + grz_lag_mid +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                       grz_lag_avg + grz_lag_early + grz_lag_mid +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                         grz_lag_avg + grz_lag_early + grz_lag_mid +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                        grz_lag_avg + grz_lag_early + grz_lag_mid +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                          grz_lag_avg + grz_lag_early + grz_lag_mid +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div +
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))


akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                                grz_lag_avg + grz_lag_early + grz_lag_mid +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div +
                              TPC_c +
                                offset(log(grazing_wkly_sum)) +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2016,
                           season='late',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2016 & 
          df_aic$season == 'late', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2016 & 
                 df_aic$season == 'late',])
}

modFNL_2016_late <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2016
  df_tmp_cor['season'] = 'late'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2017 Early ####
df_sub <- df %>% filter(Year == 2017 & season_str == 'early') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response",
                          "grz_lag_early", "grz_lag_mid"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~.,
                   family=truncated_nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div + 
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))

akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div + 
                            TPC_c +
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2017,
                           season='early',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2017 & 
          df_aic$season == 'early', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2017 & 
                 df_aic$season == 'early',])
}

modFNL_2017_early <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2017
  df_tmp_cor['season'] = 'early'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2017 Mid ####
df_sub <- df %>% filter(Year == 2017 & season_str == 'mid') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response",
                          "grz_lag_mid"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                     grz_lag_avg + grz_lag_early + 
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                       grz_lag_avg + grz_lag_early +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                         grz_lag_avg + grz_lag_early +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                        grz_lag_avg + grz_lag_early +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                          grz_lag_avg + grz_lag_early +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div +
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))



akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                                grz_lag_avg + grz_lag_early +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div +
                                TPC_c +
                                offset(log(grazing_wkly_sum)) +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2017,
                           season='mid',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2017 & 
          df_aic$season == 'mid', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2017 & 
                 df_aic$season == 'mid',])
}

modFNL_2017_mid <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2017
  df_tmp_cor['season'] = 'mid'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2017 Late ####
df_sub <- df %>% filter(Year == 2017 & season_str == 'late') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                     grz_lag_avg + grz_lag_early + grz_lag_mid +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                       grz_lag_avg + grz_lag_early + grz_lag_mid +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                         grz_lag_avg + grz_lag_early + grz_lag_mid +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                        grz_lag_avg + grz_lag_early + grz_lag_mid +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                          grz_lag_avg + grz_lag_early + grz_lag_mid +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div +
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))


akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                                grz_lag_avg + grz_lag_early + grz_lag_mid +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div +
                                TPC_c +
                                offset(log(grazing_wkly_sum)) +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2017,
                           season='late',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2017 & 
          df_aic$season == 'late', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2017 & 
                 df_aic$season == 'late',])
}

modFNL_2017_late <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2017
  df_tmp_cor['season'] = 'late'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2018 Early ####
df_sub <- df %>% filter(Year == 2018 & season_str == 'early') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response",
                          "grz_lag_early", "grz_lag_mid"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div + 
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))

akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))


modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div + 
                                offset(log(grazing_wkly_sum)) +
                              TPC_c +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2018,
                           season='early',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2018 & 
          df_aic$season == 'early', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2018 & 
                 df_aic$season == 'early',])
}

modFNL_2018_early <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2018
  df_tmp_cor['season'] = 'early'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2018 Mid ####
df_sub <- df %>% filter(Year == 2018 & season_str == 'mid') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response",
                          "grz_lag_mid"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                     grz_lag_avg + grz_lag_early +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                       grz_lag_avg + grz_lag_early +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                         grz_lag_avg + grz_lag_early +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                        grz_lag_avg + grz_lag_early +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                          grz_lag_avg + grz_lag_early +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div +
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))

akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                                grz_lag_avg + grz_lag_early +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div +
                              TPC_c +
                                offset(log(grazing_wkly_sum)) +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2018,
                           season='mid',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))


if(df_aic[df_aic$Year == 2018 & 
         df_aic$season == 'mid', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2018 & 
                 df_aic$season == 'mid',])
}

modFNL_2018_mid <- modFULL_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                                  pred_secs = predict(mod_test, ., 
                                                      type='response', 
                                                      allow.new.levels=TRUE),
                                  pred_rel = pred_secs / avg_expected_secs,
                                  pred_rel = case_when(pred_rel < 0 ~ 0,
                                                       pred_rel > 40 ~ 40,
                                                       TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2018
  df_tmp_cor['season'] = 'mid'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}

#### 2018 Late ####
df_sub <- df %>% filter(Year == 2018 & season_str == 'late') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(c("pred_zprob", "pred_response"))) %>%
  na.omit(.)

mod0_nb <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                     grz_lag_avg + grz_lag_early + grz_lag_mid +
                     offset(log(grazing_wkly_sum)) +
                     (1|Pasture) + (1|week),
                   zi=~0,
                   family=nbinom2,
                   data=df_sub,
                   control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nb)
plot(simulationOutput)
testZeroInflation(simulationOutput)

mod0_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                       grz_lag_avg + grz_lag_early + grz_lag_mid +
                       offset(log(grazing_wkly_sum)) +
                       (1|Pasture) + (1|week),
                     zi=~.,
                     family=truncated_nbinom2,
                     data=df_sub,
                     control = glmmTMBControl(parallel = 10))
simulationOutput <- simulateResiduals(fittedModel = mod0_nbzi)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(mod0_nb, mod0_nbzi)

modVEG_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                         grz_lag_avg + grz_lag_early + grz_lag_mid +
                         poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                         offset(log(grazing_wkly_sum)) +
                         (1|Pasture) + (1|week),
                       zi=~.,
                       family=truncated_nbinom2,
                       data=df_sub,
                       control = glmmTMBControl(parallel = 10))

modPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                        grz_lag_avg + grz_lag_early + grz_lag_mid +
                        PC_dmt * PC_div +
                        offset(log(grazing_wkly_sum)) +
                        (1|Pasture) + (1|week),
                      zi=~.,
                      family=truncated_nbinom2,
                      data=df_sub,
                      control = glmmTMBControl(parallel = 10))

modFULL_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                          grz_lag_avg + grz_lag_early + grz_lag_mid +
                          poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                          PC_dmt * PC_div +
                          offset(log(grazing_wkly_sum)) +
                          (1|Pasture) + (1|week),
                        zi=~.,
                        family=truncated_nbinom2,
                        data=df_sub,
                        control = glmmTMBControl(parallel = 10))


akaike.weights(c(AIC(mod0_nbzi), 
                 AIC(modVEG_nbzi),
                 AIC(modPC_nbzi),
                 AIC(modFULL_nbzi)))

modFULL_TPC_nbzi <- glmmTMB(formula = grazing_secs ~ dFence + dTank + 
                                grz_lag_avg + grz_lag_early + grz_lag_mid +
                                poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                                PC_dmt * PC_div +
                              TPC_c + 
                                offset(log(grazing_wkly_sum)) +
                                (1|Pasture) + (1|week),
                              zi=~.,
                              family=truncated_nbinom2,
                              data=df_sub,
                              control = glmmTMBControl(parallel = 10))

anova(modFULL_nbzi, modFULL_TPC_nbzi)

all_aic <- c(AIC(mod0_nbzi), 
             AIC(modVEG_nbzi),
             AIC(modPC_nbzi),
             AIC(modFULL_nbzi),
             AIC(modFULL_TPC_nbzi))

df_aic <- rbind(df_aic, 
                data.frame(Year=2018,
                           season='late',
                           models=c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC'),
                           AIC=all_aic,
                           AICd=akaike.weights(all_aic)$deltaAIC,
                           AICw=akaike.weights(all_aic)$weights))

if(df_aic[df_aic$Year == 2018 & 
          df_aic$season == 'late', 'AICd'][5] > 0.0){
  print('NOTE: Full model + TPC not best.')
  print(df_aic[df_aic$Year == 2018 & 
                 df_aic$season == 'late',])
}

modFNL_2018_late <- modFULL_TPC_nbzi

mod_names <- c('Null', 'Veg', 'PC', 'PC+Veg', 'PC+Veg+TPC')
mods <- list(mod0_nbzi, modVEG_nbzi, modPC_nbzi, modFULL_nbzi, modFULL_TPC_nbzi)
for(i in seq(1:5)){
  mod_name <- mod_names[i]
  mod_test <- mods[[i]]
  df_tmp_cor <- df_sub %>% mutate(grazing_wkly_sum = avg_grazing_secs,
                    pred_secs = predict(mod_test, ., 
                                        type='response', 
                                        allow.new.levels=TRUE),
                    pred_rel = pred_secs / avg_expected_secs,
                    pred_rel = case_when(pred_rel < 0 ~ 0,
                                         pred_rel > 40 ~ 40,
                                         TRUE ~ pred_rel)) %>%
    group_by(week, Pasture) %>%
    mutate(pred_rel_bins = cut(pred_rel,
                               breaks=quantile(pred_rel,
                                               probs=c(0, 0.05, 0.10, 0.20,
                                                       0.40, 0.60, 0.80,
                                                       0.90, 0.95, 1.0),
                                               type=6),
                               include.lowest = TRUE,
                               labels=seq(1, 9))) %>%
    ungroup() %>%
    summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
              pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
  df_tmp_cor['models'] = mod_name
  df_tmp_cor['Year'] = 2018
  df_tmp_cor['season'] = 'late'
  df_sub_cor <- rbind(df_sub_cor, df_tmp_cor)
}


#### Plot model fits ####
avg_grazing_secs <- df %>% group_by(Year, week, Pasture) %>% 
  summarise(n = n(), grazing_secs = mean(grazing_wkly_sum)) %>% 
  mutate(expected_secs = grazing_secs/n) %>% 
  ungroup() %>%
  dplyr::select(grazing_secs) %>% 
  colMeans() %>% unname()

avg_expected_secs <- df %>% group_by(Year, week, Pasture) %>% 
  summarise(n = n(), grazing_secs = mean(grazing_wkly_sum)) %>% 
  mutate(expected_secs = grazing_secs/n) %>% 
  ungroup() %>%
  dplyr::select(expected_secs) %>% 
  colMeans() %>% unname()

pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
df_pred <- data.frame()
for(yr in unique(df$Year)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_test <- get(paste('modFNL', yr, seas, sep='_'))
    if(seas == 'late'){
      col_sub = c("pred_zprob", "pred_response")
    } else if(seas == 'mid'){
      col_sub = c("pred_zprob", "pred_response", 'grz_lag_mid')
    } else {
      col_sub = c("pred_zprob", "pred_response", 'grz_lag_early', 'grz_lag_mid')
    }
    df_sub <- df %>% filter(Year == yr & season_str == seas) %>%
      mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
      dplyr::select(-any_of(col_sub)) %>%
      na.omit(.) %>%
      mutate(grazing_wkly_sum = avg_grazing_secs,
             pred_z = predict(mod_test, ., 
                                 type='zprob', 
                                 allow.new.levels=TRUE),
             pred_secs = predict(mod_test, ., 
                                 type='response', 
                                 allow.new.levels=TRUE)) %>%
      
      group_by(week, Pasture) %>%
      mutate(pred_rel = case_when(pred_z > 0.75 ~ 0,
                                  pred_secs < 30 ~ 0,
                                  TRUE ~ pred_secs / mean(pred_secs)),
             pred_rel2 = pred_secs / mean(pred_secs),
             pred_rel = case_when(pred_rel < 0 ~ 0,
                                  pred_rel > 40 ~ 40,
                                  TRUE ~ pred_rel),
             pred_rel2 = case_when(pred_rel2 < 0 ~ 0,
                                   pred_rel2 > 40 ~ 40,
                                   TRUE ~ pred_rel2)) %>%
      mutate(pred_rel_bins = cut(pred_rel2,
                                 breaks=quantile(pred_rel2,
                                                 probs=c(0, 0.05, 0.10, 0.20,
                                                         0.40, 0.60, 0.80,
                                                         0.90, 0.95, 1.0),
                                                 type=6),
                                 include.lowest = TRUE,
                                 labels=seq(1, 9)))
    df_pred <- rbind(df_pred, df_sub %>%
                       dplyr::select(-any_of(c('grz_lag_early', 'grz_lag_mid'))))
    stepi = stepi + 1
  }
}
close(pb)

df_pred %>%
  group_by(Year, season_str, pred_rel_bins) %>%
  summarise(mean = mean(grazing_rel_freq, na.rm=TRUE),
            sd = sd(grazing_rel_freq, na.rm=TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se) %>%
  ggplot(aes(x=pred_rel_bins, y=mean)) +
  facet_wrap(Year~season_str) +
  geom_point() +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=0.2) +
  geom_smooth(aes(x=as.numeric(pred_rel_bins), y=mean), method=lm) +
  geom_hline(yintercept=1.0, linetype='dashed', color='black')


df_pred_cor <- df_pred %>%
  group_by(Year, season_str) %>%
  summarise(cor = cor(as.numeric(pred_rel_bins), grazing_rel_freq),
            pval = cor.test(as.numeric(pred_rel_bins), grazing_rel_freq)[[3]])
df_pred_cor

df_pred %>%
  group_by(Year, season_str) %>%
  summarise(cor = cor(as.numeric(pred_rel), grazing_rel_freq),
            pval = cor.test(as.numeric(pred_rel), grazing_rel_freq)[[3]])

#### Maps ####
write.csv(df_pred, 
          'C:/Users/Sean.Kearney/OneDrive - USDA/Documents/Projects/GPS_v_hetgen/data/pred_all_gridded.csv')
yr = 2017
seas = 'early'
past = '15E'

df_pred_sub <- df_pred %>%
  filter(Year == yr, Pasture == past) %>%
  group_by(season_str, UTM_X, UTM_Y) %>%
  filter(week == max(week)) %>%


library(sp)
dat = data.frame(res = residuals(simGrp), x = df_sub_grp$UTM_X, 
                 y = df_sub_grp$UTM_Y)
coordinates(dat) = ~ x + y


#### Plot grazing lag ####
pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
lag_zi_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'zero_inflated') %>%
      filter(str_detect(Parameter, 'grz_lag')) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('avg', seas, sub('grz_lag_', '', Parameter))) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    lag_zi_coef_all <- rbind(lag_zi_coef_all, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)
lag_zi_coef_all <- lag_zi_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')),
         Class = factor(Class, levels=c('early', 'mid', 'late')))

pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
lag_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'conditional') %>%
      filter(str_detect(Parameter, 'grz_lag')) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('avg', seas, sub('grz_lag_', '', Parameter))) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    lag_coef_all <- rbind(lag_coef_all, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)
lag_coef_all <- lag_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')),
         Class = factor(Class, levels=c('early', 'mid', 'late')))

ggplot(data=lag_zi_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Lag of grazing (Certain zero)') +
  xlim(0.5, 1.5)

ggplot(data=lag_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Lag of grazing (Conditional)') +
  xlim(0.75, 1.25)

#### Plot dFence and dTank ####
pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
dist_zi_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'zero_inflated') %>%
      filter(str_detect(Parameter, 'dFence') | 
               str_detect(Parameter, 'dTank') ) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('d', '', Parameter)) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    dist_zi_coef_all <- rbind(dist_zi_coef_all, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)
dist_zi_coef_all <- dist_zi_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')),
         Class = factor(Class, levels=c('Tank', 'Fence')))

pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
dist_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'conditional') %>%
      filter(str_detect(Parameter, 'dFence') | 
               str_detect(Parameter, 'dTank') ) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('d', '', Parameter)) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    dist_coef_all <- rbind(dist_coef_all, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)
dist_coef_all <- dist_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')),
         Class = factor(Class, levels=c('Tank', 'Fence')))

ggplot(data=dist_zi_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Distance to feature (Certain zero)')

ggplot(data=dist_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Distance to feature (Conditional)')

#### Plot Biomass and CP ####
avg_grazing_secs <- df %>% group_by(Year, week, Pasture) %>% 
  summarise(n = n(), grazing_secs = mean(grazing_wkly_sum)) %>% 
  mutate(expected_secs = grazing_secs/n) %>% 
  ungroup() %>%
  dplyr::select(grazing_secs) %>% 
  colMeans() %>% unname()

avg_expected_secs <- df %>% group_by(Year, week, Pasture) %>% 
  summarise(n = n(), grazing_secs = mean(grazing_wkly_sum)) %>% 
  mutate(expected_secs = grazing_secs/n) %>% 
  ungroup() %>%
  dplyr::select(expected_secs) %>% 
  colMeans() %>% unname()


cp_v <- c()
bm_v <- c()
for(cp in seq(-2.0, 2.0, 0.1)){
  for(bm in seq(-2.0, 2.0, 0.1)){
    cp_v <-c(cp_v, cp)
    bm_v <- c(bm_v, bm)
  }
}

pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
df_plot2 <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    df_sub <- df %>% filter(Year == yr & season_str == seas)
    df_tmp <- data.frame(Year = yr,
                         season_str = seas,
                         Pasture = NA,
                         week = NA,
                         grazing_wkly_sum = avg_grazing_secs,
                         TPC_c = 'Flat Plains',
                         PC_dmt = 'C3_C4_mix',
                         PC_div = mean(df_sub$PC_div),
                         dFence = mean(df_sub$dFence),
                         dTank = mean(df_sub$dTank),
                         grz_lag_avg = 1.0,
                         grz_lag_early = 1.0,
                         grz_lag_mid = 1.0,
                         CP = cp_v,
                         Biomass = bm_v)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp$pred_zprob <- predict(mod_tmp, newdata=df_tmp, type='zprob')
    df_tmp$pred_response <- predict(mod_tmp, newdata=df_tmp, type='response')
    df_plot2 <- rbind(df_plot2, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)

df_plot2$season_str <- factor(df_plot2$season_str, 
                              levels=c('early', 'mid', 'late'))

df_plot2 <- df_plot2 %>% 
  mutate(CP_grp = case_when(CP < -1.5 ~ 'V. Low',
                            CP < -1.0 & CP>=-1.5 ~ 'Low',
                            CP < -0.5 & CP>=-1.0 ~ 'Mod low',
                            CP < 0.5 & CP>=-0.5 ~ 'Mod',
                            CP < 1.0 & CP>=0.5 ~ 'Mod hi',
                            CP < 1.5 & CP>=1.0 ~ 'Hi',
                            CP >= 1.5 ~ 'V. Hi')) %>%
  mutate(CP_grp = factor(CP_grp, levels=c('V. Low', 
                                          'Low', 
                                          'Mod low',
                                          'Mod', 
                                          'Mod hi',
                                          'Hi',
                                          'V. Hi')))

ggplot(data=df_plot2, mapping=aes(x=Biomass, y=pred_zprob, color=CP_grp)) +
  geom_smooth(mapping=aes(color=CP_grp),
              method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_grid(Year~season_str) +
  geom_hline(yintercept = 0.0, linetype='dashed', color='white') +
  scale_color_viridis_d(direction=1) +
  ggtitle('Certain zero') + 
  theme_dark()

ggplot(data=df_plot2, mapping=aes(x=Biomass, y=pred_response, color=CP_grp)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_grid(Year~season_str) +
  geom_hline(yintercept = avg_expected_secs, linetype='dashed', color='white') +
  scale_color_viridis_d(direction=1) +
  ggtitle('Conditional') + 
  theme_dark()

#### Plot PC_dmt ####
pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
pc_zi_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'zero_inflated') %>%
      filter(str_detect(Parameter, 'PC_dmt') & 
               str_detect(Parameter, ':', negate=TRUE)) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('PC_dmt', '', Parameter)) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    pc_zi_coef_all <- rbind(pc_zi_coef_all, df_tmp)
    stepi = stepi + 1
    }
}
close(pb)
pc_zi_coef_all <- pc_zi_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')))

pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
pc_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'conditional') %>%
      filter(str_detect(Parameter, 'PC_dmt') & 
               str_detect(Parameter, ':', negate=TRUE)) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('PC_dmt', '', Parameter)) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    pc_coef_all <- rbind(pc_coef_all, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)
pc_coef_all <- pc_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')))

ggplot(data=pc_zi_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Plant Community (Certain zero)')

ggplot(data=pc_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Plant Community (Conditional)')


#### Plot PC_div ####
div_v <- c()
for(div in seq(-2.0, 2.0, 0.1)){
  div_v <-c(div_v, div)
  }

pb = txtProgressBar(min = 0, max = 9*length(unique(df$PC_dmt)), initial = 0) 
stepi=0
df_plot3 <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    for(pcc in unique(df$PC_dmt)){
      setTxtProgressBar(pb,stepi)
      df_sub <- df %>% filter(Year == yr & season_str == seas)
      df_tmp <- data.frame(Year = yr,
                           season_str = seas,
                           Pasture = NA,
                           week = NA,
                           grazing_wkly_sum = avg_grazing_secs,
                           TPC_c = 'Flat Plains',
                           PC_dmt = pcc,
                           PC_div = div_v,
                           dFence = mean(df_sub$dFence),
                           dTank = mean(df_sub$dTank),
                           grz_lag_avg = 1.0,
                           grz_lag_early = 1.0,
                           grz_lag_mid = 1.0,
                           CP = mean(df_sub$CP),
                           Biomass = mean(df_sub$Biomass))
      mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
      df_tmp$pred_zprob <- predict(mod_tmp, newdata=df_tmp, type='zprob')
      df_tmp$pred_response <- predict(mod_tmp, newdata=df_tmp, type='response')
      df_plot3 <- rbind(df_plot3, df_tmp)
      stepi = stepi + 1
    }
    
  }
}
close(pb)

df_plot3$season_str <- factor(df_plot3$season_str, 
                              levels=c('early', 'mid', 'late'))

df_plot3 <- df_plot3 %>% 
  mutate(CP_grp = case_when(CP < -1.5 ~ 'V. Low',
                            CP < -1.0 & CP>=-1.5 ~ 'Low',
                            CP < -0.5 & CP>=-1.0 ~ 'Mod low',
                            CP < 0.5 & CP>=-0.5 ~ 'Mod',
                            CP < 1.0 & CP>=0.5 ~ 'Mod hi',
                            CP < 1.5 & CP>=1.0 ~ 'Hi',
                            CP >= 1.5 ~ 'V. Hi')) %>%
  mutate(CP_grp = factor(CP_grp, levels=c('V. Low', 
                                          'Low', 
                                          'Mod low',
                                          'Mod', 
                                          'Mod hi',
                                          'Hi',
                                          'V. Hi')))

ggplot(data=df_plot3, mapping=aes(x=PC_div, y=pred_zprob, color=PC_dmt)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_grid(Year~season_str) +
  geom_hline(yintercept = 0.0, linetype='dashed', color='white') +
  ggtitle('Certain zero') + 
  theme_dark()

ggplot(data=df_plot3, mapping=aes(x=PC_div, y=pred_response, color=PC_dmt)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_grid(Year~season_str) +
  geom_hline(yintercept = avg_expected_secs, linetype='dashed', color='white') +
  ggtitle('Conditional') + 
  theme_dark()

#### Plot TPC ####
pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
tpc_zi_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'zero_inflated') %>%
      filter(str_detect(Parameter, 'TPC_c')) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('TPC_c', '', Parameter)) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    tpc_zi_coef_all <- rbind(tpc_zi_coef_all, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)
tpc_zi_coef_all <- tpc_zi_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')),
         Class = factor(Class, levels=c('Open Slopes', 'Lowlands',
                                        'Highlands', 'Other')))

pb = txtProgressBar(min = 0, max = 9, initial = 0) 
stepi=0
tpc_coef_all <- data.frame()
for(yr in seq(2016, 2018, 1)){
  for(seas in unique(df$season_str)){
    setTxtProgressBar(pb,stepi)
    mod_tmp <- get(paste('modFNL', yr, seas, sep='_'))
    df_tmp <- model_parameters(mod_tmp, exponentiate=TRUE) %>%
      filter(Effects == 'fixed' & Component == 'conditional') %>%
      filter(str_detect(Parameter, 'TPC_c')) %>% 
      mutate(Year = yr,
             season = seas,
             Class = sub('TPC_c', '', Parameter)) %>%
      dplyr::select(Year, season, Class, Coefficient, CI_low, CI_high)
    tpc_coef_all <- rbind(tpc_coef_all, df_tmp)
    stepi = stepi + 1
  }
}
close(pb)
tpc_coef_all <- tpc_coef_all %>%
  mutate(season = factor(season, levels=c('early', 'mid', 'late')),
         Class = factor(Class, levels=c('Open Slopes', 'Lowlands',
                                        'Highlands', 'Other')))

ggplot(data=tpc_zi_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Topographic Position (Certain zero)')

ggplot(data=tpc_coef_all, aes(y=Class, x=Coefficient)) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.2) +
  facet_grid(Year~season) +
  geom_vline(xintercept = 1.0, linetype="dashed") +
  theme_bw() +
  ggtitle('Topographic Position (Conditional)')

#### Summary and spatial autocor ####
yr=2016
seas='early'

mod_test <- get(paste('modFNL', yr, seas, sep='_'))
if(seas == 'late'){
  col_sub = c("pred_zprob", "pred_response")
} else if(seas == 'mid'){
  col_sub = c("pred_zprob", "pred_response", 'grz_lag_mid')
} else {
  col_sub = c("pred_zprob", "pred_response", 'grz_lag_early', 'grz_lag_mid')
}
df_sub <- df %>% filter(Year == yr & season_str == seas) %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  dplyr::select(-any_of(col_sub)) %>%
  na.omit(.)

summary(mod_test)
model_parameters(mod_test, exponentiate = FALSE)

df_sub_grp <- df_sub %>%
  group_by(ID) %>%
  summarise(week=mean(week)) %>%
  ungroup()
df_sub_grp$UTM_X = sapply(strsplit(df_sub_grp$ID, '_'), function(x){
  as.integer(x[1])})
df_sub_grp$UTM_Y = sapply(strsplit(df_sub_grp$ID, '_'), function(x){
  as.integer(x[2])})

simOut <- simulateResiduals(fittedModel = mod_test)
plot(simOut)
simGrp <- recalculateResiduals(simOut, group=df_sub$ID)

# add a spatial variogram
library(gstat)
library(sp)
dat = data.frame(res = residuals(simGrp), x = df_sub_grp$UTM_X, 
                 y = df_sub_grp$UTM_Y)
coordinates(dat) = ~ x + y
vario = variogram(res~1, data = dat, alpha=c(0,45,90,135), 
                  width=10, cutoff=500)
plot(vario, ylim = c(-1,1))
