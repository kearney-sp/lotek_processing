df <- read.csv('C:/Users/sean.kearney/OneDrive - USDA/Documents/Projects/GPS_v_hetgen/data/training_grazing_time_gridded.csv')
df$Biomass_sq = df$Biomass^2
df$grz_lag_avg_sq = df$grz_lag_avg^2
df$grz_lag_avg_cu = df$grz_lag_avg^3

library(piecewiseSEM)
library(tidyverse)
library(nlme)
library(ctsem)

df_sub <- df %>% filter(Year == 2017 & season_str == 'early') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  na.omit(.)
#  group_by(week) %>% slice_sample(n=500)

df_wide <- ctLongToWide(datalong=df_sub, 
                        id="ID",
                        time="week", 
                        manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                        TIpredNames = c("dFence", "dTank"))

df_wide <- as.data.frame(df_wide)
df_wide$ID <- row.names(df_wide)

df_wide <- df_wide %>% left_join(df_sub[,c('ID', 'Pasture')], by='ID')
form <- as.formula(paste0(
  'Biomass_T13 ~ ', 
  paste0('grazing_rel_freq_T', seq(0, 13), collapse = ' + ')))

test <- lme(form, random=~1|Pasture, data=df_wide)
summary(test)

veg_sem <- psem(
  lme(grz_lag_avg ~ dFence + dTank,
      random=~1|week/Pasture,
      data = df_sub),
  lme(Biomass ~ grz_lag_avg + dFence + dTank,
      random=~1|week/Pasture,
      data = df_sub),
  lme(CP ~ grz_lag_avg + dFence + dTank,
      random=~1|week/Pasture,
      data = df_sub),
  lme(grazing_rel_freq ~ dFence + dTank + grz_lag_avg + 
        Biomass + Biomass_sq + CP + Biomass:CP + Biomass_sq:CP,
      random=~1|week/Pasture,
      data=df_sub),
  Biomass %~~% CP,
  data = df_sub
)

PC_veg_sem <- psem(
  lme(grz_lag_avg ~ dFence + dTank + PC_dmt + PC_div,
      random=~1|Pasture/week,
      data = df_sub),
  lme(Biomass ~ grz_lag_avg + dFence + dTank + PC_dmt + PC_div,
      random=~1|Pasture/week,
      data = df_sub),
  lme(CP ~ grz_lag_avg + dFence + dTank + PC_dmt + PC_div,
      random=~1|Pasture/week,
      data = df_sub),
  lme(grazing_rel_freq ~ dFence + dTank + grz_lag_avg + 
        Biomass + Biomass_sq + CP + Biomass:CP + Biomass_sq:CP +
        PC_dmt + PC_div,
      random=~1|Pasture/week,
      data=df_sub),
  Biomass %~~% CP,
  data = df_sub
)

summary(veg_sem, .progressBar = TRUE)
summary(PC_veg_sem, .progressBar = TRUE)


library(lme4)

test <- glmer(grazing_rel_freq ~ dFence + dTank + grz_lag_avg + 
        Biomass + CP + Biomass:CP +
        PC_dmt + PC_div + (1|Pasture) + (1|week),
      family=poisson,
      data=df_sub,
      control = glmerControl(calc.derivs = FALSE),
      nAGQ = 0)


library(glmmTMB)

df_sub <- df %>% filter(Year == 2017 & season_str == 'early') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  na.omit(.)

test <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg + 
                  poly(Biomass, 2) + CP + Biomass:CP +
                  PC_dmt + PC_div + 
                  offset(log(grazing_wkly_sum)) +
                  (1|Pasture) + (1|week),
                zi=~.,
                family=nbinom2,
                data=df_sub,
                control = glmmTMBControl(parallel = 10))
summary(test)
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = test)
plot(simulationOutput)
testZeroInflation(simulationOutput)

anova(veg_sem, PC_veg_sem)

df_sub$pred_zprob <- predict(test, type='zprob')
df_sub$pred_response <- predict(test, type='response') 
df_sub <- df_sub %>%
  group_by(Pasture) %>%
  mutate(pred_rel = pred_response * grazing_wkly_sum / grazing_wkly_sum / n()) %>%
  ungroup()

library(tidyverse)
library(ggplot2)
df_sub <- df_sub %>% 
  mutate(CP_grp = case_when(CP_orig < 6 ~ 'V. Low',
                            CP_orig < 7 & CP_orig>=6 ~ 'Low',
                            CP_orig < 8 & CP_orig>=7 ~ 'Mod low',
                            CP_orig < 9 & CP_orig>=8 ~ 'Mod',
                            CP_orig < 10 & CP_orig>=9 ~ 'Mod hi',
                            CP_orig >= 10 ~ 'Hi')) %>%
  mutate(lag_grp = case_when(grz_lag_avg < 0.90 ~ 'Low',
                             grz_lag_avg < 1.1 & grz_lag_avg >= 0.90 ~ 'Avg',
                             grz_lag_avg >= 1.1 ~ 'Hi'))
ggplot(data=df_sub, mapping=aes(x=Biomass_orig, y=pred_rel, color=CP_grp)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_wrap(~Year)

ggplot(data=df_sub, mapping=aes(x=Biomass_orig, y=pred_zprob, color=CP_grp)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_wrap(~Year)





df_sub <- df %>% filter(Year == 2017 & season_str == 'mid') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  na.omit(.)

mod_Veg_mid <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg + 
                  poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                  offset(log(grazing_wkly_sum)) +
                  (1|Pasture) + (1|week),
                zi=~.,
                family=nbinom2,
                data=df_sub,
                control = glmmTMBControl(parallel = 10))

mod_VegPC_mid <- glmmTMB(formula = grazing_secs ~ dFence + dTank + grz_lag_avg + 
                           poly(Biomass, 2) + CP + poly(Biomass, 2):CP +
                           PC_dmt + PC_div + 
                           offset(log(grazing_wkly_sum)) +
                           (1|Pasture) + (1|week),
                         zi=~.,
                         family=nbinom2,
                         data=df_sub,
                         control = glmmTMBControl(parallel = 10))
anova(mod_Veg_mid, mod_VegPC_mid)
simulationOutput <- simulateResiduals(fittedModel = test2)
plot(simulationOutput)
testZeroInflation(simulationOutput)

df_sub$pred_zprob <- predict(mod_VegPC_mid, type='zprob')
df_sub$pred_response <- predict(mod_VegPC_mid, type='response') 
df_sub <- df_sub %>%
  group_by(Pasture) %>%
  mutate(pred_rel = pred_response * grazing_wkly_sum / grazing_wkly_sum / n()) %>%
  ungroup()

df_sub <- df_sub %>% 
  mutate(CP_grp = case_when(CP_orig < 6 ~ 'V. Low',
                            CP_orig < 7 & CP_orig>=6 ~ 'Low',
                            CP_orig < 8 & CP_orig>=7 ~ 'Mod low',
                            CP_orig < 9 & CP_orig>=8 ~ 'Mod',
                            CP_orig < 10 & CP_orig>=9 ~ 'Mod hi',
                            CP_orig >= 10 ~ 'Hi')) %>%
  mutate(lag_grp = case_when(grz_lag_avg < 0.90 ~ 'Low',
                             grz_lag_avg < 1.1 & grz_lag_avg >= 0.90 ~ 'Avg',
                             grz_lag_avg >= 1.1 ~ 'Hi'))
ggplot(data=df_sub, mapping=aes(x=Biomass_orig, y=pred_rel, color=CP_grp)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_wrap(~Year)

ggplot(data=df_sub, mapping=aes(x=Biomass_orig, y=pred_zprob, color=CP_grp)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_wrap(~Year)

ggplot(data=df_sub, mapping=aes(x=Biomass_orig, y=pred_rel, color=PC_dmt)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_wrap(~Year)

ggplot(data=df_sub, mapping=aes(x=Biomass_orig, y=pred_zprob, color=PC_dmt)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_wrap(~Year)















grz_sem_grp <- multigroup(grz_sem, 'Year')

grz_sem_grp

bm_func <- function(x, b1, b2){b1*x + b2*(x^2)}
b1 <- -0.700823
b2 <- -0.151348
x <- seq(-3.0, 3.0, 0.1)
y <- bm_func(x, b1, b2)
plot(x, y)

sem.predict(grz_sem)




df_sub <- df %>% filter(Year == 2017 & Pasture == '15E') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  na.omit(.)
test <- lme(grazing_rel_freq ~ Biomass + Biomass_sq + CP + dFence + dTank + grz_lag_avg,
            random=~1|week/Pasture/ID,
            data=df_sub)
summary(test)


library(mgcv)
test2 <- gamm(grazing_rel_freq ~ Biomass + Biomass_sq + CP + dFence + dTank + grz_lag_avg,
              random=list(week=~1),
              correlation=corExp(form = ~ UTM_X + UTM_Y),
              data=df_sub)
summary(test2$gam)
plot(test2)

library(sp)
library(gstat) 
df_spat <- df_sub
# convert simple data frame into a spatial data frame object
coordinates(df_spat)= ~ UTM_X + UTM_Y
# create a bubble plot with the random values
bubble(df_spat, zcol='grazing_rel_freq', fill=TRUE, do.sqrt=FALSE, maxsize=3)
TheVariogram=variogram(dFence~1, data=df_spat)
plot(TheVariogram)
TheVariogramModel <- vgm(psill=2, model="Gau", nugget=1, range=100)
FittedModel <- fit.variogram(TheVariogram, model=TheVariogramModel)    
plot(TheVariogram, model=FittedModel)
FittedModel
