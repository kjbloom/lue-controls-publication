#############################################################
### Reading in GPP estimates and environmental variables  ###
#############################################################


## load the packages we shall need

library(tidyverse)
library(lubridate)
library(zoo)   # infrastructure for time series (Zeileis's Ordered Observations)


mass_C <- 12.0107  # molecular mass of carbon
BandG <- 0.889 # Bristow et al (1985), their B coefficient (Eqn 6) as calculated from the available data here


# read in the daily eddy-covariance data and P-model simulations provided by Beni Stocker (Stocker et al. 2020) together with site details and forcing variables such as temperature, VPD etc.

# see README file for column descriptions:

ddf <- read_csv("input_data/df_fluxnet-pmodel-forcings_publication.csv") 

summary(ddf)

# lots of missing values there (e.g. GPP-eddycovariance and diffuse radiation). 
# start with a clean up exercise:


ddf_2 <- ddf %>% 
  
  # we exclude any sites explicitly identified as C4 vegetation together with two Tropical_Grassland sites:
  filter(c4 == "FALSE" | is.na(c4), sitename != "AU-DaP", sitename != "AU-Stp") %>%
  
  # we drop the negative gpp measurements:
  filter(gpp_obs >= 0) %>%
  
  # we calculate a daily value for absorbed solar radiation and exclude very low values:
  mutate(I_abs = ppfd_fluxnet2015 * fapar_spl) %>%
  filter(I_abs > 0.1) %>%
  
  # and we derive an alternative Iabs based on EVI rather than FPAR:
  mutate(Iabs_evi = ppfd_fluxnet2015 * evi) %>%
  
  # we gap-fill the diffuse radiation fraction using Bristow et al. (1985) and the B coefficient calculated using the available data:
  mutate(Ttotal = ppfd_fluxnet2015 / ppfd_pot,
         diffuse_fraction = PPFD_DIF / ppfd_pot,
         # knock off crazy values,
         diffuse_fraction = ifelse(diffuse_fraction <= 0, NA, diffuse_fraction),
         diffuse_fraction = ifelse(diffuse_fraction >= 1, NA, diffuse_fraction),
         # then gap fill
         diffuse_fraction = ifelse(is.na(diffuse_fraction), 
                                   Ttotal * (1 - exp(0.6 * (1 - BandG / Ttotal) / (BandG - 0.4))),
                                   diffuse_fraction)) %>% 
  
  # and again set fraction limits between zero and unity:
  mutate(diffuse_fraction = ifelse(diffuse_fraction <= 0, NA, diffuse_fraction),
         diffuse_fraction = ifelse(diffuse_fraction >= 1, NA, diffuse_fraction)) %>% 
  
  # and let's convert the vpd values to kPa as more manageable:
  mutate(vpd_day_kPa = vpd_day_fluxnet2015 / 1000) %>%
  
  # next we can create new columns in the dataframe that use the 'date' to assign Day, Week, Month:
  mutate(day = yday(date), 
         week = week(date), 
         month = month(date), 
         year = year(date)) %>%
  
  # and drop some redundant columns, including the c4 identifier:
  select(-c(vpd_day_fluxnet2015, c4, Ttotal)) 
  


# coerce vegetation class and climate zone as factors (we leave sitename as a character variable for now):
ddf_2$veg <- as.factor(ddf_2$veg)
ddf_2$kgclim <- as.factor(ddf_2$kgclim)


## after all that we have lost about 2/5 of the observation rows. 

summary(ddf_2)

# So some missing values related to SW_IN_POT, SPLASH simulations and P-model runs.


######
### We restrict our data to the growing season BEFORE we calculate time-step averages
######

##  A modified approach for each site_year based on a threshold criteria suggested by Lasslop et al. 2012.  Here we simply exclude those rows that fall below the Lasslop threshold.

## Setting the growing season needs to be tailored for each Site_Year combination and for that we can use a nested loop function with Site at the highest level. 

# First we need to define a variable for the sites:

AllSites <- unique(ddf_2$sitename)

SitesYears <- list()

for (i in 1:length(AllSites)) {
  Site.i <- AllSites[i]
  ddf.comb.i <- subset(ddf_2, sitename == Site.i)
  AllYears <- unique(ddf.comb.i$year)
  
  site.yrs <- list() # a nested list to capture all years for a given site
  
  for (j in 1:length(AllYears)) {
    Year.j <- AllYears[j]
    
    site_yr.j <- subset(ddf.comb.i, year == Year.j)
    
    top <- quantile(site_yr.j$gpp_obs, 0.95, na.rm = T)
    toe <- quantile(site_yr.j$gpp_obs, 0.05, na.rm = T) # we then set this as our zero for the threshold
    
    grow <- subset(site_yr.j, gpp_obs > toe + (0.2 * (top - toe)))
    
    site.yrs[[j]] <- grow
    
    rm(grow)
    
  }
  
  SitesYears[[i]] <- bind_rows(site.yrs)
  
  rm(site.yrs)  # clear this out before looping back to the next site
}


ddf_grow <- bind_rows(SitesYears)


# this "growing-season" tibble is approaching two thirds the size of the earlier, pruned dataset. 


######
### Next we time-average the variables and employ a 15-day window e.g. Reichstein et al. (2005), Keenan et al. (2019).
######

source("functions/averaging_timesteps.R")


## Here we embed that averaging function within a second looped function.  We do this in two stages because later functions averaging the underlying data at a range of time-steps require access to the ddf_grow tibble.

# To avoid confusion, we remove the named (and now redundant) objects created in the earlier loop:
rm(SitesYears, Site.i, ddf.comb.i, AllYears, Year.j, site_yr.j, i, j, top, toe)


SitesYears <- list()

for (i in 1:length(AllSites)) {
  Site.i <- AllSites[i]
  ddf.comb.i <- subset(ddf_grow, sitename == Site.i)
  AllYears <- unique(ddf.comb.i$year)
  
  site.yrs <- list() # a nested list to capture all years for a given site
  
  for (j in 1:length(AllYears)) {
    Year.j <- AllYears[j]
    
    site_yr.j <- subset(ddf.comb.i, year == Year.j)
    
    steps <- site_yr.j %>% 
      
      # now for our moving averages tied to each Site_Year:
      mutate(tDay.ma = mafun(temp_day_fluxnet2015),
             prec.ma = mafun(prec_fluxnet2015),
             alpha.ma = mafun(meanalpha),
             vpdDay.ma = mafun(vpd_day_kPa),
             patm.ma = mafun(patm_fluxnet2015),
             ppfd.ma = mafun(ppfd_fluxnet2015),
             fpar.ma = mafun(fapar_spl),
             evi.ma = mafun(evi),
             splash.ma = mafun(sm_lim),
             co2.ma = mafun(co2_site_filled),
             co2_ML.ma = mafun(co2_mauna),
             Rdiff.ma = mafun(diffuse_fraction),
             gpp_obs_sum = sum_fun(gpp_obs),
             Iabs_sum = sum_fun(I_abs),
             Iabs2_sum = sum_fun(Iabs_evi),
             gpp_full_sum = sum_fun(gpp_mod_FULL),
             prec_sum = sum_fun(prec_fluxnet2015),
             # and our workings for a cloudiness index (here I opt for 15-day totals rather than daily averages)
             Rg_sum = sum_fun(ppfd_fluxnet2015),
             Rp_sum = sum_fun(ppfd_pot),
             Cloud_index = (1 - (Rg_sum / Rp_sum))) %>%
      
      
      # next our method for calculating the LUE index (sum BEFORE averaging):
      mutate(lue2_obs = gpp_obs_sum / (mass_C * Iabs2_sum), # EVI basis
             lue3_obs = gpp_obs_sum / (mass_C * Iabs_sum), # fPAR basis
             lue2_mod_FULL = gpp_full_sum / (mass_C * Iabs_sum)) %>%
      
      # finally we filter out unreasonable LUE estimates:
      filter(lue2_obs < 0.12, lue2_mod_FULL > 0.001)
    
    site.yrs[[j]] <- steps
    
    rm(steps)
  }
  
  SitesYears[[i]] <- bind_rows(site.yrs)
  
  rm(site.yrs)  # clear this out before looping back to the next site
}


ddf_3 <- bind_rows(SitesYears) %>% 
  
  # prune out redundant columns.  What we want to retain are the time-averaged variables:
  select(sitename, date, elv, veg, day, week, month, year,
       tDay.ma, prec.ma, alpha.ma, vpdDay.ma, patm.ma, ppfd.ma, fpar.ma, evi.ma, splash.ma, co2.ma, co2_ML.ma, 
       Rdiff.ma, Cloud_index,
       lue2_obs, lue3_obs, lue2_mod_FULL)



### Before proceeding, let's clean up the work-space and remove some of the intermediate objects we created:

rm(AllSites, AllYears, BandG, ddf.comb.i, i, j, mafun, Site.i, site_yr.j, SitesYears, sum_fun, Year.j)


###### End



