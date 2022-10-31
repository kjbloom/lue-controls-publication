#########################################################
###    Inter-model comparisons: P-model versus NACP   ###
### Producing Figure 2 in Bloomfield et al. 2022, GCB ###
#########################################################


## What happens when we apply our statistical model to other TBM predictions relating to NACP sites included in main dataset?

library(tidyverse)
library(lme4)
library(visreg)
library(zoo)
library(gridExtra)
library(grid)


## This is at the daily time-step. 

nacp <- read_csv(file = "input_data/df_nacp_gpp.csv") %>% 
  rename(CAN.IBIS = 'CAN-IBIS', CLM.CN = 'CLM-CN') # dashes in the naming convention cause some problems later


# we capture the unique model names for later steps:
tbm_list <- names(nacp)[c(7:29)]


## The NACP download has been updated to include LAI estimates.
# so we need to differentiate GPP and LAI columns.  We append a suffix to the relevant column headings:

names(nacp)[7:31] <- str_c(names(nacp)[7:31], "_gpp")


## NEXT we read in the LAI estimates:

nacp2 <- read_csv(file = "input_data/df_nacp_LAI.csv") %>% 
  # dashes in the naming convention cause some problems later
  rename(CAN.IBIS = 'CAN-IBIS', CLM.CN = 'CLM-CN') %>%
  # remove ambiguous column (-9999 passim)
  select(- c(PCT_FILLED, X31))  


names(nacp2)[6:30] <- str_c(names(nacp2)[6:30], "_lai")


## now we'll join these NACP estimates together; the two files are exactly the same length

nacp <- nacp %>% 
  inner_join(nacp2, by = c("sitename", "YEAR", "FDOY"))


rm(nacp2) # clear out the separate LAI file to avoid any confusion later on


## and to accompany that we have the cleaned up input files used for the main analysis; 
# here we are interested in ddf_grow which is daily data restricted to a site-specific growing season

source("load_GPP-dependencies.R")

# clear out objects we don't need here:
rm(ddf, ddf_2, ddf_3)


# So we want to merge these two dataframes - nacp and ddf_grow - using unique combination of Site_Year_Day

# To make that exercise easier, let's standardise the names of the key variables.  'sitename' is good

ddf_grow <- ddf_grow %>%
  rename(YEAR = year, FDOY = day) %>%
  # and we retain only the variables we need here - specifically we revert to Iabs derived from FPAR
  select(sitename, YEAR, FDOY,
         gpp_obs, gpp_mod_FULL, 
         ppfd_fluxnet2015, ppfd_pot, temp_day_fluxnet2015, vpd_day_kPa, 
         fapar_spl, sm_lim) 


# and then we try an inner-join which retains only those rows that match: 

trevor <- nacp %>%
  inner_join(ddf_grow, by = c("sitename", "YEAR", "FDOY"))



## Many of the TBMs may not have generated daily estimates for these sites, can we weed these out?

summary(trevor)

# there are FOUR that return '-9999' passim for GPP: AGROIBIS, DNDC, EPIC, SIBCROP

# and SEVEN such for LAI: AGROIBIS, DNDC, ECLUEEDCM, ED2, EPIC, SIBCROP, TRIPLEX



######
### 'How does our preferred statistical model perform when applied to these alternative TBM predictions?'
######

# There are multiple TBMs included here, so we need to loop through those.

# A bit of tidying up to do: organise on the 15-day timestep & calculate LUE for the TBM predictions

source("functions/averaging_timesteps.R")

# the 15-day averaging function (using rollapply from zoo) does not recognise breaks in Site or Year, so embed that within a looped function:


######
### Here is the fully nested function: TBM/Site/Year .. ###
######


## The visreg command needs access to the model object AND the data used to generate the model.  

## We use a loop to create a unique dataframe for each TBM and then create the accompanying model tied to that specific tibble.  Each model with it's unique dataset, is then available for the visreg() function.


# first we prune our model list to exclude those TBMs that generated NULL values throughout:
tbm_list2 <- tbm_list[-c(1, 7, 8, 10, 11, 19, 22)]


for (k in 1 : length(tbm_list2)) {
  
  tbm.k <- tbm_list2[k] 
  
  tbm.gpp <- str_c(tbm.k, "_gpp")
  tbm.lai <- str_c(tbm.k, "_lai")

  df_name <- paste("trev", tbm.k, sep = "_")
  
  trev2 <- trevor %>%
    # let's keep only the forcing variables and the tbm of interest.  For validation we also need inferred GPP.
    select(sitename, YEAR, gpp_obs,
           ppfd_fluxnet2015, ppfd_pot, temp_day_fluxnet2015, vpd_day_kPa, sm_lim,
           all_of(tbm.gpp), all_of(tbm.lai)) %>%
    
    # let's filter out prediction gaps and we need both GPP and LAI:
    filter(get(tbm.gpp) != -9999 & get(tbm.lai) != -9999 )
  
  
  AllSites <- unique(trev2$sitename)
  
  SitesYears <- list()
  
  for (i in 1:length(AllSites)) {
    Site.i <- AllSites[i]
    ddf.comb.i <- subset(trev2, sitename == Site.i)
    AllYears <- unique(ddf.comb.i$YEAR)
    
    site.yrs <- list() # a nested list to capture all years for a given site
    
    for (j in 1:length(AllYears)) {
      Year.j <- AllYears[j]
      
      site_yr.j <- subset(ddf.comb.i, YEAR == Year.j)
      
      steps <- site_yr.j %>%
        
        # and we calculate fAPAR based on LAI using an extinction coefficient of 0.5
        mutate(fAPAR = 1 - exp(-0.5 * get(tbm.lai)),
               I_abs = fAPAR * ppfd_fluxnet2015) %>% 
        
        # and exclude any very low values
        filter(I_abs > 0.1) %>% 
        
        
        mutate(tDay.ma = mafun(temp_day_fluxnet2015),
               vpdDay.ma = mafun(vpd_day_kPa), 
               splash.ma = mafun(sm_lim), 
               Iabs_sum = sum_fun(I_abs),
               gpp_obs_sum = sum_fun(gpp_obs),
               tbm_sum = sum_fun(get(tbm.gpp)),
               # and our workings for a cloudiness index (here I opt for 15-day totals rather than daily averages)
               Rg_sum = sum_fun(ppfd_fluxnet2015),
               Rp_sum = sum_fun(ppfd_pot),
               Cloud_index = (1 - (Rg_sum / Rp_sum))
               )
      
      
      site.yrs[[j]] <- steps
      
      rm(steps)
    }
    
    SitesYears[[i]] <- bind_rows(site.yrs)
    
    rm(site.yrs)  # clear this out before looping back to the next site
  }
  
  
  trev3 <- bind_rows(SitesYears) 

  
  ## Now we tidy up this df to retain only the averaged/summed rows:
  
  trev4 <- trev3 %>%
    
    # drop the raw daily columns
    select(-c(3:11)) %>%
    
    # omit rows with missing values (e.g. days 1:14):
    na.omit() %>%
    
    # next our revised method for calculating the LUE index; and a step to overwrite zero soilm estimates:
    mutate(lue_TBM = tbm_sum / (mass_C * Iabs_sum),
           lue_obs = gpp_obs_sum / (mass_C * Iabs_sum),
           splash.ma = ifelse(splash.ma < 0.001, 0.001, splash.ma)) %>% 
    
    # finally restrict ourselves to reasonable LUE estimates:
    filter(lue_obs < 0.12 & lue_TBM > 0.001)
  
  
  assign(df_name, trev4) # and a unique name to this TBM specific dataset
  
  
  # and clean up the environment before looping back through
  rm(AllSites, AllYears, ddf.comb.i, df_name, i, j, k, Site.i, site_yr.j, SitesYears, 
     tbm.k, tbm.gpp, tbm.lai, 
     trev2, trev3, trev4, Year.j)
  
  }



## That generates 16 tibbles of differing length.
## For comparative purposes we also need to provide the tailored fit for our preferred empirical model applied to the GPP observed in this NACP overlap dataset, also the P-model.

# so here equivalent averaging and calculations for the 'observations' and P-model.  And we proceed from the tibble that merges Beni's workings with the NACP simulations:

trev2 <- trevor %>%
  # let's keep only the forcing variables 
  select(sitename, YEAR, 
         gpp_obs, gpp_mod_FULL,
         ppfd_fluxnet2015, ppfd_pot, temp_day_fluxnet2015, vpd_day_kPa, fapar_spl, sm_lim) 


AllSites <- unique(trev2$sitename)
  
SitesYears <- list()
  
for (i in 1:length(AllSites)) {
    Site.i <- AllSites[i]
    ddf.comb.i <- subset(trev2, sitename == Site.i)
    AllYears <- unique(ddf.comb.i$YEAR)
    
    site.yrs <- list() # a nested list to capture all years for a given site
    
    for (j in 1:length(AllYears)) {
      Year.j <- AllYears[j]
      
      site_yr.j <- subset(ddf.comb.i, YEAR == Year.j)
      
      steps <- site_yr.j %>%
        
        # let's filter out any negative GPP estimates
        filter(gpp_obs >= 0) %>%
        
        # we calculate I_abs and exclude very low values; here I follow Beni's advice and prefer fapar_splined:
        mutate(I_abs = ppfd_fluxnet2015 * fapar_spl) %>%
        filter(I_abs > 0.1) %>%
        
        mutate(tDay.ma = mafun(temp_day_fluxnet2015),
               vpdDay.ma = mafun(vpd_day_kPa), 
               splash.ma = mafun(sm_lim), 
               Iabs_sum = sum_fun(I_abs), 
               gpp_obs_sum = sum_fun(gpp_obs),
               gpp_full_sum = sum_fun(gpp_mod_FULL),
               Rg_sum = sum_fun(ppfd_fluxnet2015),
               Rp_sum = sum_fun(ppfd_pot),
               Cloud_index = (1 - (Rg_sum / Rp_sum))) %>% 
        
        # next our revised method for calculating the LUE index (sum BEFORE averaging):
        mutate(lue_obs = gpp_obs_sum / (mass_C * Iabs_sum), 
               lue_mod_FULL = gpp_full_sum / (mass_C * Iabs_sum))
        
        # finally we filter out unreasonable LUE estimates:
        #filter(lue_obs < 0.12 & lue_mod_FULL > 0.001)
      
      
      site.yrs[[j]] <- steps
      
      rm(steps)
    }
    
    SitesYears[[i]] <- bind_rows(site.yrs)
    
    rm(site.yrs)  # clear this out before looping back to the next site
  }
  
  
trev3 <- bind_rows(SitesYears)
  
## Now we tidy up this df to retain only the averaged rows:
  
trev_emp <- trev3 %>%
  
  # drop the raw daily columns
  select(-c(3:10)) %>%
    
  # omit rows with missing values (e.g. days 1:14):
  na.omit() %>%
    
  # next a step to overwrite zero soilm estimates:
  mutate(splash.ma = ifelse(splash.ma < 0.001, 0.001, splash.ma))
  
  
# and clean up the environment 
rm(AllSites, AllYears, ddf.comb.i, i, j, Site.i, site_yr.j, SitesYears, Year.j, trev2, trev3)
  



### ... next we fit our statistical model (main script) to these tbm-specific datasets ###

Form <- formula(lue_TBM ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                  (1 | sitename / YEAR))

# the BEPS tibble only provides data for a single site and so we need to modify the random term there:
Form_yr <- formula(lue_TBM ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                     (1 | YEAR))

FM_beps <- glmer(Form_yr, family = Gamma(link = "log"), 
                data = trev_BEPS,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_biomebgc <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_BIOMEBGC,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_ibis <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_CAN.IBIS, 
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 


FM_clm <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_CLM.CN,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_class <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_CNCLASS,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

# convergence warning; try restarting from the previous fit 
ss.class <- getME(FM_class, c("theta", "fixef")) # extract components from an earlier attempt

FM_class2 <- update(FM_class, start = ss.class, control = glmerControl(optCtrl = list(maxfun = 2e5)))

# warning message persists, proceed with caution

FM_dlem <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_DLEM,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))



FM_ecosys <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_ECOSYS, 
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_isam <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_ISAM,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_isolsm <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_ISOLSM,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_lotec <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_LOTEC,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
# warning of boundary singular fit??


FM_lpj <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_LPJ_wsl,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_orchid <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_ORCHIDEE,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_sib <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_SIB,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_sibcasa <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_SIBCASA,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_ssib2 <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_SSIB2,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


FM_teco <- glmer(Form, family = Gamma(link = "log"), 
                 data = trev_TECO,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))




#  and here is the same model applied to the FLUXNET observations:

FM_empirical <- glmer(lue_obs ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                        (1 | sitename / YEAR), 
                      family = Gamma(link = "log"), 
                      data = trev_emp,
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# and the P-model:

FM_Pmodel <- glmer(lue_mod_FULL ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                     (1 | sitename / YEAR), 
                   family = Gamma(link = "log"), 
                   data = trev_emp, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))



###########################################
### AND NOW for a huge multi-TBM figure ###
###########################################

### To aid comparision, we hold the axes' ranges constant for all TBMs:

lue_range <- c(0, 0.08) # this will truncate a few plots (e.g. LPJ), but the fitted lines are unaffected
T_range <- c(0, 30)


# we create separate visreg objects to generate the data but suppress the plot (plot = FALSE)

# and now a looped function to run through those iterations; here the nine models included in the Figure 3 of the main text:

fig3mods <- c("FM_empirical", "FM_Pmodel", "FM_isam", "FM_lpj", "FM_class", 
              "FM_ibis", "FM_lotec", "FM_orchid", "FM_ssib2")


fig3labs <- c("Empirical", "P-model", "ISAM", "LPJ", "CN-CLASS", "Can-IBIS", "LoTEC", "ORCHIDEE", "SSiB2")


T3mods <- list()

for(i in 1:length(fig3mods)) {
  
  vmod <- fig3mods[i]
  vlabel <- fig3labs[i]
  vT_mod <- visreg(get(vmod), "tDay.ma", type = "conditional", trans = exp, plot = FALSE)
  
  # create a tibble of those predictions and add a column identifying the relevant TBM
  T3mods[[i]] <- vT_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vT_mod)
  
  }


# and we bind those listed tibbles into a single dataframe:
Tfits3 <- bind_rows(T3mods)


# a quick clean out of intermediate objects from the looped function
rm(i, T3mods)



### creating grobs for each separate TBM_effect ###

# We tailor some unique formatting for the first (Empirical) and last (Ensemble) panels, but can loop through the other TBMs:

# we start with the temperature effects:

Tfigs <- list()

# first the Empirical model:

Tfigs[[1]] <- visreg(get(fig3mods[1]), "tDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                     gg = TRUE,
                     points = list(pch = 1),
                     line = list(lty = 1, col = "red")
                     ) + 
  labs(title = fig3labs[1], x = "Temp.", y = " ") +
  scale_x_continuous(limits = T_range) +
  scale_y_continuous(limits = lue_range) +
  theme(panel.background = element_rect(fill = rgb(0.1, 0.5, 0.1, alpha=0.1), color = 'green'),
        # here we tailor the plot margins; sequence is: top, right, bottom, left
        plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(fig3mods)) {
  
  Tfigs[[i]] <- visreg(get(fig3mods[i]), "tDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                       gg = TRUE,
                       points = list(pch = 1),
                       line = list(lty = 1, col = "red"),
                       yaxt = 'n') + 
    labs(title = fig3labs[i], x = " ", y = " ") +
    scale_x_continuous(limits = T_range) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))
  
  }


# and then we create a grob for the final ensemble plot:

Tfigs[[10]] <- ggplot(data = Tfits3, aes(tDay.ma, visregFit)) +
  geom_line(data = subset(Tfits3, nacpMod != "Empirical"), aes(group = nacpMod), col = "grey") +
  geom_line(data = subset(Tfits3, nacpMod == "Empirical"), col = "red") +
  labs(title = "Ensemble", x = "Temp.", y = " ") +
  scale_x_continuous(limits = T_range) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.5, 0.2, 0.1, 0.1), "cm"))




### now the equivalent plot for VPD effects ## 

D3mods <- list()

for(i in 1:length(fig3mods)) {
  
  vmod <- fig3mods[i]
  vlabel <- fig3labs[i]
  vD_mod <- visreg(get(vmod), "vpdDay.ma", type = "conditional", trans = exp, plot = FALSE)
  
  D3mods[[i]] <- vD_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vD_mod)
  
}


Dfits3 <- bind_rows(D3mods)


## and now the VPD grobs; we drop the model labels here:


Dfigs <- list()

# first the Empirical model:

Dfigs[[1]] <- visreg(get(fig3mods[1]), "vpdDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                     gg = TRUE,
                     points = list(pch = 1),
                     line = list(lty = 1, col = "red")
                     ) + 
  labs(x = "VPD", y = " ") +
  scale_x_continuous(limits = c(0, 3.5)) +
  scale_y_continuous(limits = lue_range) +
  theme(panel.background = element_rect(fill = rgb(0.1, 0.5, 0.1, alpha=0.1), color = 'green'),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(fig3mods)) {
  
  Dfigs[[i]] <- visreg(get(fig3mods[i]), "vpdDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                       gg = TRUE,
                       points = list(pch = 1),
                       line = list(lty = 1, col = "red"),
                       yaxt = 'n') + 
    labs(x = " ", y = " ") +
    scale_x_continuous(limits = c(0, 3.5)) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  }


# and then we create a grob for the final ensemble plot:

Dfigs[[10]] <- ggplot(data = Dfits3, aes(vpdDay.ma, visregFit)) +
  geom_line(data = subset(Dfits3, nacpMod != "Empirical"), aes(group = nacpMod), col = "grey") +
  geom_line(data = subset(Dfits3, nacpMod == "Empirical"), col = "red") +
  labs(x = "VPD", y = " ") +
  scale_x_continuous(limits = c(0, 3.5)) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"))



### and next for Soil moisture term ##

Sm3mods <- list()

for(i in 1:length(fig3mods)) {
  
  vmod <- fig3mods[i]
  vlabel <- fig3labs[i]
  vSm_mod <- visreg(get(vmod), "splash.ma", type = "conditional", trans = exp, plot = FALSE)
  
  Sm3mods[[i]] <- vSm_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vSm_mod)
  
}


Smfits3 <- bind_rows(Sm3mods)


## and now the Soilm grobs:


Smfigs <- list()

# first the Empirical model:

Smfigs[[1]] <- visreg(get(fig3mods[1]), "splash.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                     gg = TRUE,
                     points = list(pch = 1),
                     line = list(lty = 1, col = "red")
                     ) + 
  labs(x = "Soilm", y = " ") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c("0", " ", "0.4", " ", "0.8", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(panel.background = element_rect(fill = rgb(0.1, 0.5, 0.1, alpha=0.1), color = 'green'),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(fig3mods)) {
  
  Smfigs[[i]] <- visreg(get(fig3mods[i]), "splash.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                       gg = TRUE,
                       points = list(pch = 1),
                       line = list(lty = 1, col = "red"),
                       yaxt = 'n') + 
    labs(x = " ", y = " ") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c("0", " ", "0.4", " ", "0.8", " ")) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  
}


# and then we create a grob for the final ensemble plot:

Smfigs[[10]] <- ggplot(data = Smfits3, aes(splash.ma, visregFit)) +
  geom_line(data = subset(Smfits3, nacpMod != "Empirical"), aes(group = nacpMod), col = "grey") +
  geom_line(data = subset(Smfits3, nacpMod == "Empirical"), col = "red") +
  labs(x = "Soilm", y = " ") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c("0", " ", "0.4", " ", "0.8", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"))


### And finally include conditional plots for the effect of Cloudiness_Index

Cloud3mods <- list()

for(i in 1:length(fig3mods)) {
  
  vmod <- fig3mods[i]
  vlabel <- fig3labs[i]
  vCloud_mod <- visreg(get(vmod), "Cloud_index", type = "conditional", trans = exp, plot = FALSE)
  
  Cloud3mods[[i]] <- vCloud_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vCloud_mod)
  
}


Cloudfits3 <- bind_rows(Cloud3mods)


## and now the cloudiness grobs:


Cloudfigs <- list()

# first the Empirical model:

Cloudfigs[[1]] <- visreg(get(fig3mods[1]), "Cloud_index", type = "conditional", trans = exp, 
                         rug = F, partial = T, 
                      gg = TRUE,
                      points = list(pch = 1),
                      line = list(lty = 1, col = "red")
                      ) + 
  labs(x = "Cloudiness Index", y = " ") +
  scale_x_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2), labels = c(" ", "0.4", "0.6 ", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(panel.background = element_rect(fill = rgb(0.1, 0.5, 0.1, alpha=0.1), color = 'green'),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(fig3mods)) {
  
  Cloudfigs[[i]] <- visreg(get(fig3mods[i]), "Cloud_index", type = "conditional", trans = exp, 
                           rug = F, partial = T, 
                        gg = TRUE,
                        points = list(pch = 1),
                        line = list(lty = 1, col = "red"),
                        yaxt = 'n') + 
    labs(x = " ", y = " ") +
    scale_x_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2), labels = c(" ", "0.4", "0.6 ", " ")) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  }


# and then we create a grob for the final ensemble plot:

Cloudfigs[[10]] <- ggplot(data = Cloudfits3, aes(Cloud_index, visregFit)) +
  geom_line(data = subset(Cloudfits3, nacpMod != "Empirical"), aes(group = nacpMod), col = "grey") +
  geom_line(data = subset(Cloudfits3, nacpMod == "Empirical"), col = "red") +
  labs(x = "Cloudiness Index", y = " ") +
  scale_x_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2), labels = c(" ", "0.4", "0.6 ", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"))




# and we arrange those all together:


fig4 <- grid.arrange(
         arrangeGrob(grobs = c(Tfigs, Dfigs, Smfigs, Cloudfigs), 
                     ncol = 10, nrow = 4, 
                     # need to give a little more space to the first column (yaxis labels) and first row (titles)
                     widths = c(1, rep(0.85, 9)), 
                     heights = c(1, 0.8, 0.8, 0.8), 
                     left = textGrob(label = expression("LUE "~ (molC~mol^{-1}~photons)), rot = 90, 
                                     gp = gpar(fontface = "bold", fontsize = 15)),
                     # reduce the margin around the left axis label
                     padding = unit(0.1, "line")
                     )
         )



######
### And now we need the companion plot for Supplementary ###
######

# first looped function to run through the nine models relegated to SI:

figSImods <- c("FM_beps", "FM_biomebgc", "FM_clm", "FM_dlem", "FM_ecosys", 
              "FM_isolsm", "FM_sib", "FM_sibcasa", "FM_teco")


figSIlabs <- c("BEPS", "Biome-BGC", "CLM-CN", "DLEM", "Ecosys", "ISOLSM", "SiB", "SiBCASA", "TECO")


TSImods <- list()

for(i in 1:length(figSImods)) {
  
  vmod <- figSImods[i]
  vlabel <- figSIlabs[i]
  vT_mod <- visreg(get(vmod), "tDay.ma", type = "conditional", trans = exp, plot = FALSE)
  
  # create a tibble of those predictions and add a column identifying the relevant TBM
  TSImods[[i]] <- vT_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vT_mod)
  
}


# and we bind those listed tibbles into a single dataframe:
TfitsSI <- bind_rows(TSImods)


### creating grobs for each separate TBM_effect ###


Tfigs_SI <- list()

# We still need some special axis features for the first column:

Tfigs_SI[[1]] <- visreg(get(figSImods[1]), "tDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                     gg = TRUE,
                     points = list(pch = 1),
                     line = list(lty = 1, col = "red")
                     ) + 
  labs(title = figSIlabs[1], x = "Temp.", y = " ") +
  scale_x_continuous(limits = T_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(figSImods)) {
  
  Tfigs_SI[[i]] <- visreg(get(figSImods[i]), "tDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                       gg = TRUE,
                       points = list(pch = 1),
                       line = list(lty = 1, col = "red"),
                       yaxt = 'n') + 
    labs(title = figSIlabs[i], x = " ", y = " ") +
    scale_x_continuous(limits = T_range) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))
  
}


# and then we create a grob for the final ensemble plot:

Tfigs_SI[[10]] <- ggplot(data = TfitsSI, aes(tDay.ma, visregFit)) +
  geom_line(data = TfitsSI, aes(group = nacpMod), col = "grey") +
  #geom_line(data = subset(TfitsSI, nacpMod == "Empirical"), col = "red") +
  labs(title = "Ensemble", x = "Temp.", y = " ") +
  scale_x_continuous(limits = T_range) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.5, 0.2, 0.1, 0.1), "cm"))




### now the equivalent plot for VPD effects ## 

Dmods_SI <- list()

for(i in 1:length(figSImods)) {
  
  vmod <- figSImods[i]
  vlabel <- figSIlabs[i]
  vD_mod <- visreg(get(vmod), "vpdDay.ma", type = "conditional", trans = exp, plot = FALSE)
  
  Dmods_SI[[i]] <- vD_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vD_mod)
  
}


Dfits_SI <- bind_rows(Dmods_SI)


## and now the VPD grobs; we drop the model labels here:


Dfigs_SI <- list()

# first the Empirical model:

Dfigs_SI[[1]] <- visreg(get(figSImods[1]), "vpdDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                     gg = TRUE,
                     points = list(pch = 1),
                     line = list(lty = 1, col = "red")
                     ) + 
  labs(x = "VPD", y = " ") +
  scale_x_continuous(limits = c(0, 3.5)) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(figSImods)) {
  
  Dfigs_SI[[i]] <- visreg(get(figSImods[i]), "vpdDay.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                       gg = TRUE,
                       points = list(pch = 1),
                       line = list(lty = 1, col = "red"),
                       yaxt = 'n') + 
    labs(x = " ", y = " ") +
    scale_x_continuous(limits = c(0, 3.5)) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}


# and then we create a grob for the final ensemble plot:

Dfigs_SI[[10]] <- ggplot(data = Dfits_SI, aes(vpdDay.ma, visregFit)) +
  geom_line(data = Dfits_SI, aes(group = nacpMod), col = "grey") +
  #geom_line(data = subset(Dfits_SI, nacpMod == "Empirical"), col = "red") +
  labs(x = "VPD", y = " ") +
  scale_x_continuous(limits = c(0, 3.5)) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"))



### and now for Soil moisture term ##

Smmods_SI <- list()

for(i in 1:length(figSImods)) {
  
  vmod <- figSImods[i]
  vlabel <- figSIlabs[i]
  vSm_mod <- visreg(get(vmod), "splash.ma", type = "conditional", trans = exp, plot = FALSE)
  
  Smmods_SI[[i]] <- vSm_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vSm_mod)
  
}


Smfits_SI <- bind_rows(Smmods_SI)


## and now the Soilm grobs:


Smfigs_SI <- list()


Smfigs_SI[[1]] <- visreg(get(figSImods[1]), "splash.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                      gg = TRUE,
                      points = list(pch = 1),
                      line = list(lty = 1, col = "red")
                      ) + 
  labs(x = "Soilm", y = " ") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c("0", " ", "0.4", " ", "0.8", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(figSImods)) {
  
  Smfigs_SI[[i]] <- visreg(get(figSImods[i]), "splash.ma", type = "conditional", trans = exp, rug = F, partial = T, 
                        gg = TRUE,
                        points = list(pch = 1),
                        line = list(lty = 1, col = "red"),
                        yaxt = 'n') + 
    labs(x = " ", y = " ") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c("0", " ", "0.4", " ", "0.8", " ")) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  
}


# and then we create a grob for the final ensemble plot:

Smfigs_SI[[10]] <- ggplot(data = Smfits_SI, aes(splash.ma, visregFit)) +
  geom_line(data = Smfits_SI, aes(group = nacpMod), col = "grey") +
  #geom_line(data = subset(Smfits3, nacpMod == "Empirical"), col = "red") +
  labs(x = "Soilm", y = " ") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c("0", " ", "0.4", " ", "0.8", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"))




### and here for our Cloudiness Index:

Cloudmods_SI <- list()

for(i in 1:length(figSImods)) {
  
  vmod <- figSImods[i]
  vlabel <- figSIlabs[i]
  vCloud_mod <- visreg(get(vmod), "Cloud_index", type = "conditional", trans = exp, plot = FALSE)
  
  Cloudmods_SI[[i]] <- vCloud_mod$fit %>% mutate(nacpMod = vlabel)
  
  rm(vmod, vlabel, vCloud_mod)
  }


Cloudfits_SI <- bind_rows(Cloudmods_SI)


## and now the Cloudiness Index grobs:


Cloudfigs_SI <- list()


Cloudfigs_SI[[1]] <- visreg(get(figSImods[1]), "Cloud_index", type = "conditional", 
                            trans = exp, rug = F, partial = T, 
                         gg = TRUE,
                         points = list(pch = 1),
                         line = list(lty = 1, col = "red")
                         ) + 
  labs(x = "Cloudiness Index", y = " ") +
  scale_x_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2), labels = c(" ", "0.4", "0.6 ", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


for(i in 2:length(figSImods)) {
  
  Cloudfigs_SI[[i]] <- visreg(get(figSImods[i]), "Cloud_index", type = "conditional", 
                              trans = exp, rug = F, partial = T, 
                           gg = TRUE,
                           points = list(pch = 1),
                           line = list(lty = 1, col = "red"),
                           yaxt = 'n') + 
    labs(x = " ", y = " ") +
    scale_x_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2), labels = c(" ", "0.4", "0.6 ", " ")) +
    scale_y_continuous(limits = lue_range) +
    theme(axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  
}


# and then we create a grob for the final ensemble plot:

Cloudfigs_SI[[10]] <- ggplot(data = Cloudfits_SI, aes(Cloud_index, visregFit)) +
  geom_line(data = Cloudfits_SI, aes(group = nacpMod), col = "grey") +
  labs(x = "Cloudiness Index", y = " ") +
  scale_x_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2), labels = c(" ", "0.4", "0.6 ", " ")) +
  scale_y_continuous(limits = lue_range) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', color = 'blue'),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"))




# and we arrange those all together:


figS7 <- grid.arrange(
         arrangeGrob(grobs = c(Tfigs_SI, Dfigs_SI, Smfigs_SI, Cloudfigs_SI), 
                     ncol = 10, nrow = 4, 
                     # we give a little more space to the first column and first row 
                     widths = c(1, rep(0.85, 9)), 
                     heights = c(1, 0.8, 0.8, 0.8), 
                     left = textGrob(label = expression("LUE "~ (molC~mol^{-1}~photons)), rot = 90, 
                                     gp = gpar(fontface = "bold", fontsize = 15)),
                     # reduce the margin around the left axis label
                     padding = unit(0.1, "line")
         )
       )




######
### Table S5: performance criteria for the various TBMs. 
## to include the Nash-Sutcliffe model efficiency coefficient
######


# here we call in a tailored function
source("functions/nash-sutcliff-ec.R")

source("functions/analyse_modobs.R")


## Prepare a table providing performance metrics for the TBMs

# We develop a list of the 16 individual tibbles:


df_list <- list()

for (i in 1 : length(tbm_list2)) {
  
  df.i <- tbm_list2[i]
  
  df_list[i] <- paste("trev", df.i, sep = "_")
  
  }

df_list <- str_c(df_list)


# and now loop through each tibble generating the performance metrics

nacp_metrics <- list()

for (i in 1 : length(tbm_list2)) {
  
  df_name <- df_list[i]
  
  nacp_metrics[[i]] <- with(get(df_name), metRics(lue_obs, lue_TBM))
  
  }


df_metrics <- bind_rows(nacp_metrics) %>% 
  mutate(TBM = tbm_list2)

# clean up our environment
rm(df_list, df_name, df.i, i, nacp_metrics)


## Who are the winners and losers?

tbm_list2[which.max(df_metrics$Rsquared)] # so on Rsquared, the best performing TBM is ORCHIDEE

tbm_list2[which.max(df_metrics$nsEC)] # whilst on Nash Sutcliffe EC it is SIBCASA


# We can visualise that with Beni's function


with(trev_SIBCASA, analyse_modobs(lue_TBM, lue_obs, 
                               heat = T, 
                               plot.title = "SiBCASA",
                               xlab = "Predicted LUE (TBM)",
                               ylab = expression("Inferred LUE"~ (molC~mol^{-1}~photons))))




### and here is the table; notice that none of the models achieve an EC of > 0.5 

view(df_metrics)
    

### END ###

  
  
