df <- read.csv('C:/Users/sean.kearney/OneDrive - USDA/Documents/Projects/GPS_v_hetgen/data/training_grazing_time_gridded.csv')
df$Biomass_sq = df$Biomass^2
df$grz_lag_avg_sq = df$grz_lag_avg^2
df$grz_lag_avg_cu = df$grz_lag_avg^3

library(piecewiseSEM)
library(tidyverse)
library(lme4)
library(nlme)
library(MASS)
library(splines)

df_sub <- df %>% filter(Year == 2017 & Pasture %in% c('15E')) %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  na.omit(.)

groups <- df_sub %>%
  group_by(ID) %>%
  arrange(Pasture) %>%
  filter(row_number()==1) %>%
  pull(Pasture)

groups <- paste0("g", groups)

ids <- df_sub %>%
  group_by(ID) %>%
  arrange(ID) %>%
  filter(row_number()==1) %>%
  pull(ID)
ids_samp <- sample(ids, size=50)
df_sub <- df_sub %>%
  filter(ID %in% ids_samp)

library(ctsemOMX)
df_wide <- ctLongToWide(datalong=df_sub, 
                        id="ID",
                        time="week", 
                        manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                        TIpredNames = c("dFence", "dTank"))

wide <- ctIntervalise(datawide = df_wide,
                      Tpoints = 14,
                      n.manifest = 3,
                      n.TDpred = 0, 
                      n.TIpred = 2,
                      manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                      TIpredNames = c("dFence", "dTank"))

example1model <- ctModel(n.latent = 3, n.manifest = 3, Tpoints = 14,
                         manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                         latentNames = c("Biomass", "CP", "grazing_rel_freq"),
                         TIpredNames = c("dFence", "dTank"),
                         LAMBDA = diag(3),
                         TRAITVAR = "auto",
                         MANIFESTVAR = 'free')
#example1model$DRIFT[1, 2] = "0" #remove effect of CP on Biomass
#example1model$DRIFT[2, 1] = "0" #remove effect of Biomass on CP
example1fit <- ctFit(dat = wide, ctmodelobj = example1model)
groupfit <- ctMultigroupFit(dat=wide, ctmodelobj = example1model, 
                            groupings=groups)

plot(groupfit)


print(example1model)
plot(example1fit,
     standardiseCR = TRUE,
     randomImpulse= TRUE,
     experimentalImpulse = TRUE)
summary(example1fit)
ctPlot(example1fit, 'CR', xlim=c(0, 15))

library(ctsem)
df_sub <- df %>% filter(Year == 2017 & Pasture == '15E') %>%
  mutate(ID = paste(UTM_X, UTM_Y, sep="_")) %>%
  na.omit(.)
#  group_by(week) %>% slice_sample(n=500)

ids <- df_sub %>%
  group_by(ID) %>%
  arrange(ID) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  group_by(Pasture) %>%
  sample_n(500) %>%
  pull(ID)

df_sub <- df_sub %>%
  filter(ID %in% ids)

df_wide <- ctLongToWide(datalong=df_sub, 
                        id="ID",
                        time="week", 
                        manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                        TIpredNames = c("dFence", "dTank"))

wide <- ctIntervalise(datawide = df_wide,
                      Tpoints = 14,
                      n.manifest = 3,
                      n.TDpred = 0, 
                      n.TIpred = 2,
                      manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                      TIpredNames = c("dFence", "dTank"))

df_long <- ctWideToLong(datawide = wide,
                     Tpoints = 14,
                     n.manifest = 3,
                     n.TDpred = 0,
                     n.TIpred = 2,
                     manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                     TIpredNames = c("dFence", "dTank"))

long <- ctDeintervalise(datalong = df_long, id='id', dT='dT')

example1model <- ctModel(type='stanct',
                         n.latent = 3, 
                         n.manifest = 3,
                         Tpoints = 14,
                         manifestNames = c("Biomass", "CP", "grazing_rel_freq"),
                         latentNames = c("Biomass", "CP", "grazing_rel_freq"),
                         TIpredNames = c("dFence", "dTank"),
                         tipredDefault = TRUE,
                         LAMBDA = diag(3), 
                         TRAITVAR = NULL,
                         MANIFESTMEANS = matrix(c(0, 0, 0)),
                         #MANIFESTVAR = matrix(c(0, 0, 0), nrow=3, ncol = 3)
                         )
example1model$pars[example1model$pars$matrix == 'DRIFT', 
                   c('dFence_effect', 'dTank_effect')] = TRUE
example1model$pars
print(example1model)

example1fit <- ctStanFit(dat = long, cores=10, 
                         optimize = TRUE, nopriors = TRUE,
                         ctstanmodel = example1model)
summary(example1fit)
plot(example1fit)

ctModelLatex(example1model)

grz_sem <- psem(
  lme(grz_lag_avg ~ dFence + dTank,
      random=~1|Year/week/Pasture,
      data = df_sub),
  lme(Biomass ~ grz_lag_avg + grz_lag_avg_sq + grz_lag_avg_cu + dFence + dTank,
      random=~1|Year/week/Pasture,
      data = df_sub),
  lme(grazing_rel_freq ~ Biomass + Biomass_sq + CP + dFence + dTank + grz_lag_avg,
      random=~1|Year/week/Pasture,
      data=df_sub),
  data = df_sub
)

summary(grz_sem, .progressBar = TRUE)

grz_sem_grp <- multigroup(grz_sem, 'Year')

grz_sem_grp

bm_func <- function(x, b1, b2){b1*x + b2*(x^2)}
b1 <- -0.15361430
b2 <- -0.06883238
x <- seq(-3.0, 3.0, 0.1)
y <- bm_func(x, b1, b2)
plot(x, y)

sem.predict(grz_sem)


test <- lm(grazing_rel_freq ~ dFence + dTank + grz_lag_avg,
            data=df_sub)
summary(test)

df_sub$resid <- test$residuals
df_sub$diff <- df_sub$grazing_rel_freq - df_sub$grz_lag_avg
df_sub <- df_sub %>%
  mutate(binary = case_when(grazing_rel_freq > 1.0 ~ 1,
                            TRUE ~ 0))

for (yr in c(2016, 2017, 2018)){
  
  test2 <- lme(grazing_rel_freq ~ Biomass + Biomass_sq + CP + Biomass:CP + Biomass_sq:CP + 
                    dFence + dTank + grz_lag_avg + grz_lag_avg:Biomass + grz_lag_avg:Biomass_sq,
                  random=~1|Pasture/week,
                  data=df_sub[df_sub$Year == yr, ])
  df_sub[df_sub$Year == yr, "pred"] <- predict(test2)
}


library(ggplot2)
df_sub <- df_sub %>% 
  mutate(CP_grp = case_when(CP < 6 ~ 'V. Low',
                            CP < 7 & CP>=6 ~ 'Low',
                            CP < 8 & CP>=7 ~ 'Mod low',
                            CP < 9 & CP>=8 ~ 'Mod',
                            CP < 10 & CP>=9 ~ 'Mod hi',
                            CP >= 10 ~ 'Hi')) %>%
  mutate(lag_grp = case_when(grz_lag_avg < 0.90 ~ 'Low',
                             grz_lag_avg < 1.1 & grz_lag_avg >= 0.90 ~ 'Avg',
                             grz_lag_avg >= 1.1 ~ 'Hi'))
ggplot(data=df_sub, mapping=aes(x=Biomass_orig, y=pred, color=CP_grp)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs="cs", k=3)) +
  facet_wrap(~Year)

test3 <- glmer(binary ~ Biomass + Biomass_sq + CP + Biomass:CP + dFence + dTank + (1|Pasture),
               family='binomial',
              data=df_sub)
summary(test3)
