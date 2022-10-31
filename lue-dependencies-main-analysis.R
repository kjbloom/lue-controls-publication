################################################################################################
###                 Environmental controls on Light Use Efficiency                           ###
## script supporting the analysis presented by Bloomfield et al. (2022) Global Change Biology ##
################################################################################################

### import the cleaned up file of inferred and simulated LUE values accompanied by the forcing variables; all averaged at 15-day timestep

source("load_GPP-dependencies.R")

## load additional packages to be employed:

library(nlme)
library(visreg)
library(grid)
library(gridExtra)


## See text file "input_data/README_sessionInfo.txt" for details of package versions etc.


summary(ddf_3)

## Some cleaning up of the imported tibble, we need complete cases for the statistical models to run:

ben_cc <- ddf_3 %>% 
  
  # remove some variables not needed here:
  select(-c(alpha.ma, patm.ma, co2_ML.ma)) %>% 
  
  # omit rows with missing values:
  na.omit() %>%
  
  # drop redundant factor levels (e.g. Sites)
  droplevels()


## replace zero soil moisture estimates with a nominal value to aid later computation (e.g. log transformation)

ben_cc$splash.ma <- with(ben_cc, ifelse(splash.ma < 0.001, 0.001, splash.ma))


# reorder the vegetation classifications - broadly following the distribution in the dataset:
fct_count(ben_cc$veg, sort = TRUE)

ben_cc$veg <- factor(ben_cc$veg, levels = c("ENF", "EBF", "DBF", "MF", "GRA", "SAV", "WSA", "CSH", "OSH"))


# after all of that we have a data-frame of about 8040 complete cases.


###########################################
### A PREDICTIVE MODEL: model selection ###
###########################################

## Notice here that our response variable is 'lue3_obs' driven by FPAR

## to help numerical estimation we scale some of the candidate variables:

ben_cc <- ben_cc %>%
  mutate(temp_scale = scale(tDay.ma, center = T, scale = T),
         elv_scale = scale(elv, center = T, scale = T),
         ppfd_scale = scale(ppfd.ma, center = T, scale = T),
         co2_scale = scale(co2.ma, center = T, scale = T))


##Step1#  
# We start with a 'beyond optimal' model including as many explanatory variables as possible and then try to whittle that down.  Here we include a quadratic term for daytime_temperature, daytime_vpd, soil moisture, elevation, CO2 and diffuse radiation.

# Notice that within the poly() function we can specify raw = TRUE to avoid orthogonal polynomials which are difficult to interpret, but that choice makes no difference to selection steps (AIC scores etc are identical).  The raw option appears to create some difficulties with specifying formulae on subsequent steps, so we shall stick with the default (orthogonal) setting for now and ask for raw only at the final interpretation stage.

# we begin with a simple linear model:

lm1 <- lm(lue3_obs ~ poly(temp_scale, degree = 2, raw = F) + log(vpdDay.ma) + splash.ma + 
            elv_scale + co2_scale + Cloud_index, 
          data = ben_cc)


anova(lm1)

# so the output suggests that the elevation term might be dropped:

lm2 <- update(lm1, .~. - elv_scale)

# how do the models compare?
anova(lm1, lm2)

# so we proceed with the simpler model

anova(lm2)

# and all of those terms are significant and the model explains ca. 45% of the variation in LUE.  From the F_values it looks like vpd is easily the strongest determinant.

# We have longitudinal data, so we need to reflect that in our model structure - repeat observations for a given site.  Something like Year nested in Site.  


##Step2# 
# Restate our model with gls to allow direct comparisons, default method here is REML:

Form <- formula(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + splash.ma + co2_scale + Cloud_index)

M02 <- gls(Form, method = "REML", data = ben_cc)


# We have a continuous response variable, so the distribution options here are normal (i.e. Gaussian) or Gamma distributions.  Both of these error distributions accept "log" as the link function.
# The default method for these models is 'glm.fit' based on reweighted least squares; but that method is NOT directly comparable (anova()) with REML.  We rely then on AIC for our selection.  Let's try some of these non-normal families and links and see if they offer any improvement on the standard approach:


M03 <- glm(Form, family = gaussian(link = "log"), data = ben_cc)

plot(M03, which = 1) # still wedged 

M04 <- glm(Form, family = Gamma(link = "log"), data = ben_cc)

plot(M04, which = 1) # that offers some improvement


# How does model performance compare?  

AIC(M02, M03, M04)

# The Akaike scores indicate that the gamma model performs better here.


##Step3# 
# Choice of variance structure: 
# we consider a range of Random terms taking us to a generalised linear MIXED model (glmm).
# A number of packages are available, but we'll use lme4:

detach(package:nlme)
library(lme4)

# First a random term providing random intercepts for Year.  The control argument is designed to help the model reach convergence, in part by increasing the number of iterations.

M05_yr <- glmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + splash.ma + co2_scale + Cloud_index + 
                  (1 | year), 
                family = Gamma(link = "log"), 
                data = ben_cc,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# For comparison, a term providing random intercepts for individual Sites

M05_st <- glmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + splash.ma + co2_scale + Cloud_index + 
                  (1 | sitename), 
                family = Gamma(link = "log"), 
                data = ben_cc,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# What about Year nested in Site?

M05_nest <- glmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + splash.ma + co2_scale + Cloud_index + 
                    (1 | sitename / year), 
                  family = Gamma(link = "log"), 
                  data = ben_cc,
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))



##Step4# 
# How do these glmm variants compare against a simpler glm with the same fixed formula?:

M05_glm <- glm(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + splash.ma + co2_scale + Cloud_index, 
               family = Gamma(link = "log"), 
               data = ben_cc)


AIC(M05_glm, M05_yr, M05_st, M05_nest)

# all right, it looks clear that the nested glmm performs best and offers a big improvement on the glm.  


##Step5# 
# With the variance structure settled and the preferred random term in place, we move to selecting the optimal fixed structure.

# We revert to a maximal model (beyond optimal) - and include a term for the interaction between diffuse fraction and canopy density:

M06 <- glmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + log(splash.ma) + 
               elv_scale + co2_scale + Cloud_index:fpar.ma +
               (1 | sitename / year), 
             family = Gamma(link = "log"), 
             data = ben_cc,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

# that fails to converge

ss1 <- getME(M06, c("theta", "fixef")) # extract components from the earlier attempt

M06b <- update(M06, 
                  start = ss1,
                  control = glmerControl(optCtrl = list(maxfun = 2e4)))

# convergence warnings persist; proceed with caution

anova(M06b) 
# the output suggests that we might drop CO2.  The diffuse radiation term looks important.


## Model selection steps.  Can we safely drop any terms without affecting performance unduly?
# The weakest term appears to be CO2:

M07 <- update(M06b, .~. - co2_scale)

# how does the simpler model compare?

anova(M06b, M07)

# So that confirms that we are justified in dropping the CO2 term.

anova(M07)

# The weakest candidate now looks to be elevation.  Can that be dropped?

M08 <- update(M07, .~. - elv_scale)

# That runs to convergence.  How do the nested variants compare?

anova(M07, M08)

# So that's a marginal call, but on BIC and in the interests of parsimony we can justify dropping.  Can we go further?

anova(M08)

# It looks like the diffuse radiation term is now the weakest candidate.  Can that be dropped?

M09 <- update(M08, .~. - Cloud_index:fpar.ma)

anova(M08, M09)

# All right, so that CLEARLY indicates that the diffuse radiation term should be retained.

# What about substituting Rdiffuse for Cloud Index?

M10 <- glmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + log(splash.ma) + Rdiff.ma:fpar.ma +
               (1 | sitename / year), 
             family = Gamma(link = "log"), 
             data = ben_cc,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


anova(M08, M10)
# and that confirms that the Cloud_Index variant performs much better

# or how about dropping the FPAR interaction?

M11 <- glmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index +
               (1 | sitename / year), 
             family = Gamma(link = "log"), 
             data = ben_cc,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

# no warnings there

anova(M08, M11)

# so that indicates that the interaction with FPAR is NOT required

anova(M11)


# can we dispense with the soil moisture term?

M12 <- update(M11, .~. - log(splash.ma))

anova(M11, M12) # confirms that we need to retain the soilm term


## So to this point we have T + D + Sm + CI.  But do we need the transformations?

M13 <- glmer(lue3_obs ~ temp_scale + vpdDay.ma + splash.ma + Cloud_index + 
               (1 | sitename / year), 
             family = Gamma(link = "log"), 
             data = ben_cc)


anova(M11, M13)

# and that confirms that the transformed explanatory variables give a much better fit.

# As a final step we could verify that the generalised design is really superior to a simpler mixed-effects model.
# effectively this operates as family = gaussian(link = "identity"), but such usage is deprecated and we use lmer directly:

M14 <- lmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index +
              (1 | sitename / year),
            data = ben_cc,
            # to be consistent with preceding steps we fit with maximum likelihood
            REML = FALSE)


# that runs with no difficulties

anova(M11, M14)

# and that leaves no doubt that the generalised model is preferred. 


###
## So we finish the selection protocol by specifying our preferred model:

Form <- formula(lue3_obs ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                  (1 | sitename / year))

M_fin <- glmer(Form, 
               family = Gamma(link = "log"), 
               data = ben_cc,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# to aid interpretation we have dropped the scaling terms, but we satisfy ourselves that this is the same model as M11 in the stepped selection above:

anova(M11, M_fin)

# the outputs are effectively identical including the coefficients.  We proceed.

# Remember that the default poly() function computes ORTHOGONAL polynomials



####################################
### RESULTS - figures and tables ###
####################################

# a standard argument for formatting axis labels:
kj.lab_size <- theme(axis.title.x=element_text(size=15), 
                     axis.title.y=element_text(size=15), 
                     axis.text.x=element_text(size=12), 
                     axis.text.y=element_text(size=12))


# We start by printing a summary of the final empirical model:

summary(M_fin)


# to ease plotting we add the fitted values to our tibble.  And we do this in two forms - with and without the random effects:

ben_cc <- ben_cc %>%
  mutate(M_fit = as.numeric(fitted(M_fin)),
         M_fit.FE = predict(M_fin, type = "response", re.form = NA)) # ignoring random effects



## Some site information from the FLUXNET2015 files:

MySites <- unique(ben_cc$sitename)

load("input_data/metainfo_Tier1_sites_kgclimate_fluxnet2015.RData")

sites_loc <- metainfo_Tier1_sites_kgclimate_fluxnet2015 %>% 
  select(sitename, lon, lat, year_start, year_end, classid, koeppen_code) %>% 
  
  mutate(period = str_c(year_start, year_end, sep = "-")) %>% 
  
  # restrict this to the sites retained for our final analysis
  subset(sitename %in% MySites) %>%
  
  # for some reason there are a handful of duplicated site references in the Tier1 file
  distinct() 




### Model evaluation ###

### Figure 1 ## is a 'Goodness of fit' style plot, and we adopt a heat-map style here:

source("functions/analyse_modobs.R")

op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

with(ben_cc, analyse_modobs(M_fit, lue3_obs, heat = T, 
                            plot.title = "Empirical model - Goodness of fit",
                            xlab = "Predicted LUE (full model)",
                            ylab = expression("Inferred LUE"~ (molC~mol^{-1}~photons)),
                            xlim = c(0, 0.09)))
mtext("(a)", side = 3, line = 0.5, adj = -0.15, cex = 1.4)


with(ben_cc, analyse_modobs(M_fit.FE, lue3_obs, heat = T, 
                            plot.title = "Ignoring random effects of Site and Year",
                            xlab = "Predicted LUE (fixed effects only)",
                            ylab = " ",
                            xlim = c(0, 0.09)))
mtext("(b)", side = 3, line = 0.5, adj = -0.15, cex = 1.4)

par(op)



### Figure 02 ## model experiments.  See separate script: "lue-dependencies_model-experiments.R"##


### Figure 03 ## diagnostic plots testing the empirical model's assumption of equality of variance:

E01 <- resid(M_fin, type = "deviance") 


fig3a <- ggplot(data = ben_cc) +
  geom_boxplot(aes(x = veg, y = E01), varwidth = F) + 
  kj.lab_size + 
  labs(x = "", y = "Model residuals (log scale)") +
  geom_text(data=NULL, x=Inf, y=Inf, label="(a)", size=5, colour="black", hjust=1.5, vjust=1.5)



## And as a companion plot we investigate the site random intercepts 

rr_fin <- ranef(M_fin, condVar = T)

rrfin_site <- tibble::rownames_to_column(rr_fin$sitename, "sitename")

# Let's tie these intercepts back to the site information:

site_int <- merge(sites_loc, rrfin_site, by = "sitename")

names(site_int)[9] <- "Intercept"

# and reorder the PFT classes to match with the main dataset:
site_int$classid <- factor(site_int$classid, 
                           levels = c("ENF", "EBF", "DBF", "MF", "GRA", "SAV", "WSA", "CSH", "OSH"))


## And now for a boxplot organised by vegetation class:

fig3b <- ggplot(data = site_int) + 
  geom_boxplot(aes(x = classid, y = Intercept), varwidth = F) +  
  kj.lab_size + 
  scale_y_continuous(limits=c(-1, 1), expand=c(0,0)) +
  labs(x = "Vegetation class", 
       y = "Random site intercepts (log scale)",
       title = "greenness index is FPAR") +
  geom_hline(yintercept = 0, colour = "blue", linetype=2, size=1.2) +
  geom_text(data=NULL, x=Inf, y=Inf, label="(b)", size=5, colour="black", hjust=1.5, vjust=1.5)


# So a somewhat altered pattern there, but consistently shrublands sitting below the general trend.  

## How to interpret these site intercept terms given our log transformation?
# The population intercept is -4.61 (log scale), which translates to an LUE around 0.010 mol/mol

exp(summary(M_fin)$coefficients[1,1])

which.max(site_int$Intercept)

# the top site from our model's random term is AR-SLu with an intercept of + 0.87

# that implies (all other model terms at zero) that the LUE at AR-SLu would be ~ 0.024 mol/mol

exp(summary(M_fin)$coefficients[1,1] + rr_fin$sitename[1,1])



## And what happens to those patterns if we adopt MODIS_EVI instead of FPAR ? ##

# So we replicate our final model structure, but with the revised LUE ratio based on EVI (lue2_obs):

M_EVI <- glmer(lue2_obs ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                 (1 | sitename / year), 
               family = Gamma(link = "log"), 
               data = ben_cc)

# we repeat the boxplot of site intercept terms organised by PFT:


rr_EVI <- ranef(M_EVI, condVar = T)

rrEVI_site <- tibble::rownames_to_column(rr_EVI$sitename, "sitename")

# Let's merge these into our existing tibble:

site_int <- site_int %>% 
  left_join(rrEVI_site, by = "sitename") %>% 
  rename(EVI_intercept = '(Intercept)')


## And the figure:

fig3c <- ggplot(data = site_int) + 
  geom_boxplot(aes(x = classid, y = EVI_intercept), varwidth = F) +  
  kj.lab_size + 
  scale_y_continuous(limits=c(-1, 1), expand=c(0,0)) +
  labs(x = "Vegetation class", 
       y = " ",
       title = "greenness index is EVI") +
  geom_hline(yintercept = 0, colour = "blue", linetype=2, size=1.2) +
  geom_text(data=NULL, x=Inf, y=Inf, label="(c)", size=5, colour="black", hjust=1.5, vjust=1.5)




### Figure 04 ## NACP comparisons.  See separate script: "lue-dependencies_NACP-model-comparisons.R" ##



#################################
### SUPPLEMENTARY INFORMATION ###
#################################


### Figure S3:  sample gap-filled curves for FPAR over time 

# The most important vegetation classes in the dataset ( cumulatively 0.8) are: ENF >> GRA > DBF > EBF
# We provide an example for each, here we use the daily dataset before filters, growing season or time-averaging:

grid.arrange(
  #ENF - DAVOS has a long history: 
  ddf %>% 
    filter(sitename == "CH-Dav") %>% 
    ggplot(mapping = aes(date, fapar_spl)) + 
    geom_point() + 
    labs(title = "CH-Dav, ENF", 
         x = "Date", 
         y = NULL) +
    scale_y_continuous(limits = c(0, 1)) + 
    kj.lab_size,
  
  
  #GRA - AT-Neu:
  ddf %>% 
    filter(sitename == "AT-Neu") %>% 
    ggplot(mapping = aes(date, fapar_spl)) + 
    geom_point() + 
    labs(title = "AT-Neu, GRA",
         x = "Date", 
         y = NULL) + 
    scale_y_continuous(limits = c(0, 1)) + 
    kj.lab_size,
  
  
  #DBF - US-MMS:
  ddf %>% 
    filter(sitename == "US-MMS") %>% 
    ggplot(mapping = aes(date, fapar_spl)) + 
    geom_point() + 
    labs(title = "US-MMS, DBF",
         x = "Date", 
         y = NULL) + 
    scale_y_continuous(limits = c(0, 1)) + 
    kj.lab_size,
  
  
  #EBF - AU-Tum 
  ddf %>% 
    filter(sitename == "AU-Tum") %>% 
    ggplot(mapping = aes(date, fapar_spl)) +
    geom_point() + 
    labs(title = "AU-Tum, EBF", 
         x = "Date", 
         y = NULL) + 
    scale_y_continuous(limits = c(0, 1)) + 
    kj.lab_size,
  
  
  ncol = 2, nrow = 2, 
  left = textGrob(label = "MODIS_FPAR", rot = 90, 
                  gp = gpar(fontface = "bold", fontsize = 15))
  )



### Figure S4: comparing fPAR variants ###

# We are interested to see how substituting EVI for FPAR might affect the different biomes.  

ggplot(data = ben_cc, aes(fpar.ma, evi.ma)) +
  geom_point(alpha = 1/10) +
  labs(y = "Enhanced vegetation index (MODIS)",
       x = "Fraction of photosynthetically active radiation (MODIS)") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_abline(intercept = 0, slope = 1.0, colour = "red", linetype = 2, size = 1.2) +
  facet_wrap(~ veg, ncol = 3) +
  kj.lab_size



### Figure S5: Choice of time-step ###

## Notice that we get convergence warnings for the model fit to the Daily dataframe ##

## WEEKLY

wkdf <- ddf_grow %>%
  # let's simplify some of the parameter names to start
  rename(T_day = temp_day_fluxnet2015, D_day = vpd_day_kPa) %>%
  
  group_by(sitename, year, week) %>% 
  
  # for our explanatory variables we want a mix of averages and sums
  summarise_at(vars(T_day, D_day, sm_lim, gpp_obs, I_abs, ppfd_fluxnet2015, ppfd_pot), 
               list(sum = sum, mean = mean), na.rm = T) %>%
  
  # here is our new method for calculating LUE using the weekly sums
  mutate(lue_obs = gpp_obs_sum / (mass_C * I_abs_sum),
         Cloud_index = (1 - (ppfd_fluxnet2015_sum / ppfd_pot_sum))) %>%
  
  # and we'll filter out any unreasonable ratios and sub-zero temperatures
  filter(lue_obs < 0.12, T_day_mean > 0) %>%
  
  # despite the NA argument above, we still generate a handful of missing values - we exclude those
  na.omit()


summary(wkdf)


# and we need to sort out any near-zero forcing values; replace with a negligible value.  Over-writing these very low VPD and soilm averages helps the model to converge:

wkdf$D_day_mean <- with(wkdf, ifelse(D_day_mean < 0.001, 0.001, D_day_mean))
wkdf$sm_lim_mean <- with(wkdf, ifelse(sm_lim_mean < 0.001, 0.001, sm_lim_mean))


# now for our model formula - tweaked here for the slightly altered syntax:

FM_wk <- glmer(lue_obs ~ poly(T_day_mean, 2) + log(D_day_mean) + log(sm_lim_mean) + Cloud_index + 
                 (1 | sitename / year), 
               family = Gamma(link = "log"), 
               data = wkdf,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# and that runs without generating warnings

summary(FM_wk)



## DAILY ##

day_df <- ddf_grow %>%
  
  # select the columns we need
  select(sitename, year, gpp_obs, 
         ppfd_fluxnet2015, temp_day_fluxnet2015, sm_lim, ppfd_pot, I_abs, vpd_day_kPa) %>% 
  
  # let's simplify some of the parameter names to start
  rename(T_day = temp_day_fluxnet2015, D_day = vpd_day_kPa) %>%
  
  # here we revert to fPAR:
  mutate(lue_obs = gpp_obs / (mass_C * I_abs),
         Cloud_index = (1 - (ppfd_fluxnet2015 / ppfd_pot))) %>%
  
  # and we'll filter out any unreasonable ratios
  filter(lue_obs < 0.12, T_day > 0) %>%
  
  # despite the NA argument above, we still generate a handful of missing values - we exclude those
  na.omit()



summary(day_df)


# and we need to sort out any zero VPD values; replace with a negligible value:

day_df$D_day <- with(day_df, ifelse(D_day < 0.001, 0.001, D_day))
day_df$sm_lim <- with(day_df, ifelse(sm_lim < 0.001, 0.001, sm_lim))


# now for our model formula and this takes some minutes to run:

FM_day <- glmer(lue_obs ~ poly(T_day, 2) + log(D_day) + log(sm_lim)  + Cloud_index + (1 | sitename / year), 
                family = Gamma(link = "log"), 
                data = day_df,
                control = glmerControl(optCtrl = list(maxfun = 2e4)))


# this generates convergence warnings and we should proceed with caution.


### We produce a multi-panel figure of these conditional plots comparing each effect over the time-steps ###  

# setting uniform axes ranges:

lue_range <- c(0, 0.10) 
T_range <- c(0, 40)
D_range <- c(0, 7)
Sm_range <- c(0, 1)
CI_range <- c(0, 1)


## We'll create the individual grobs before trying to arrange them:

Tfig1 <- visreg(FM_day, "T_day", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")
                ) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = T_range) +
  scale_y_continuous(limits = lue_range) +
  geom_text(data=NULL, x=Inf, y=Inf, label="Daily", size=4, colour="black",
              hjust=1, vjust=1) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


Tfig2 <- visreg(FM_wk, "T_day_mean", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")
                ) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = T_range) +
  scale_y_continuous(limits = lue_range) +
  geom_text(data=NULL, x=Inf, y=Inf, label="Weekly averages", size=4, colour="black",
              hjust=1, vjust=1) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


Tfig3 <- visreg(M_fin, "tDay.ma", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")
                ) + 
  labs(x = "Temperature", y = NULL) +
  scale_x_continuous(limits = T_range) +
  scale_y_continuous(limits = lue_range) +
  geom_text(data=NULL, x=Inf, y=Inf, label="15-day averages", size=4, colour="black",
              hjust = 1, vjust = 1) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


## now for VPD

Dfig1 <- visreg(FM_day, "D_day", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = D_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


Dfig2 <- visreg(FM_wk, "D_day_mean", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = D_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


Dfig3 <- visreg(M_fin, "vpdDay.ma", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = "VPD", y = NULL) +
  scale_x_continuous(limits = D_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


## now for soil moisture:

Smfig1 <- visreg(FM_day, "sm_lim", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = Sm_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


Smfig2 <- visreg(FM_wk, "sm_lim_mean", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = Sm_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


Smfig3 <- visreg(M_fin, "splash.ma", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = "Soil moisture", y = NULL) +
  scale_x_continuous(limits = Sm_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


## and finally, Cloudiness Index:

CIfig1 <- visreg(FM_day, "Cloud_index", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = CI_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


CIfig2 <- visreg(FM_wk, "Cloud_index", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = CI_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))


CIfig3 <- visreg(M_fin, "Cloud_index", type = "conditional", trans = exp, rug = F, partial = T,
                gg = TRUE,
                points = list(pch = 1),
                line = list(lty = 1, col = "red")) + 
  labs(x = "Cloudiness Index", y = NULL) +
  scale_x_continuous(limits = CI_range) +
  scale_y_continuous(limits = lue_range) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))



Fig_S5 <- grid.arrange(
  arrangeGrob(
    grobs = list(Tfig1, Tfig2, Tfig3, 
                 Dfig1, Dfig2, Dfig3, 
                 Smfig1, Smfig2, Smfig3, 
                 CIfig1, CIfig2, CIfig3), 
    as.table = FALSE, 
    ncol = 4, nrow = 3, 
    left = textGrob(label = expression("LUE "~ (molC~mol^{-1}~photons)), rot = 90, 
                    gp = gpar(fontface = "bold", fontsize = 15)), 
    # reduce the margin around the left axis label 
    padding = unit(0.1, "line")
    )
  )



### Figure-S6, cross-validation exercise: see separate script "lue-dependencies_cross-validation.R" ###


### Figure-S7, NACP companion plot: see separate script "lue-dependencies_NACP-model-comparisons.R" ###


### Figure-S8, looking at seasonal patterns to observed and modelled LUE 
# and we restrict this to temperate/boreal sites.  

ben_cc %>%
  left_join(sites_loc) %>% 
  filter(lat > 23) %>% # this puts us above the Tropic of Cancer
  mutate(bias = M_fit - lue3_obs) %>% 
  ggplot(aes(day, bias)) + 
  geom_point(na.rm = T, alpha = 1/10) + 
  geom_smooth(se = F) +
  labs(title = "Sites above the Tropic of Cancer", 
       caption = "15-day averages",
       x = "Day of Year",
       y = expression("LUE bias"~ (molC~mol^{-1}~photons))) +
  kj.lab_size



### Supplementary Figure 09 ## 
# In our workings, we assume a linear dependence on light, but is that borne out?

# Here we extend the preferred model (M11 above) to include an additional explanatory term of PPFD.


M15 <- glmer(lue3_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + ppfd_scale + 
               (1 | sitename / year), 
             family = Gamma(link = "log"), 
             data = ben_cc,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# How does that compare with our preferred empirical model above?
anova(M_fin, M15)

# all right, so that indicates the PPFD term is valuable and it's inclusion improves model performance.  But it is the weakest of the explanatory terms:

anova(M15)


# What does the conditional plot look like?:

visreg(M15, "ppfd_scale", type = "conditional", trans = exp, rug = F, partial = T, gg = T,
       points = list(pch = 1),
       line = list(lty = 1, col = "red"),
       xlab = expression("PPFD (scaled, underlying range is 1 to 70"~mol~m^{-2}~day^{-1}~")"),
       ylab = expression("LUE "~(molC~mol^{-1}~photons))) +
  ylim(0, 0.1) +
  labs(title = "Conditional response to PPFD") +
  kj.lab_size




### Table S4: Here we demonstrate that the choice of FPAR versus EVI in calculating LUE makes little (or no) difference to the design of the empirical model.  So the key step here is to substitute the response (dependent) variable and here we want lue2_obs: 


## We repeat the selection protocol with the new response variable:

M01_evi <- glmer(lue2_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + log(splash.ma) + 
                   elv_scale + co2_scale + Cloud_index:fpar.ma +
                   (1 | sitename / year), 
                 family = Gamma(link = "log"), 
                 data = ben_cc,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


# so that reaches convergence

ss.evi <- getME(M01_evi, c("theta", "fixef")) # extract components from the earlier attempt

M01_evi.b <- update(M01_evi, 
                    start = ss.evi,
                    control = glmerControl(optCtrl = list(maxfun = 2e4)))


#  so the warning persists, proceed with caution:


summary(M01_evi.b) 
# the output suggests that we can drop CO2, elevation looks important perhaps.


## Model selection steps.  Can we safely drop any terms without affecting performance unduly?
# The weakest term appears to be CO2:

M02_evi <- update(M01_evi.b, .~. - co2_scale)

# how does the simpler model compare?

anova(M01_evi.b, M02_evi)

# So that confirms that we are justified in dropping the CO2 term.

anova(M02_evi)

# The weakest candidate now looks to be the elevation term.  Can that be dropped?

M03_evi <- update(M02_evi, .~. - elv_scale)

# how do the nested variants compare?

anova(M02_evi, M03_evi)

# So that's a marginal call, but on BIC and in the interests of parsimony we could justify dropping.  Can we go further?

summary(M03_evi)

# Can we dispense with the FPAR interaction term?

M04_evi <- glmer(lue2_obs ~ poly(temp_scale, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index +
                   (1 | sitename / year), 
                 family = Gamma(link = "log"), 
                 data = ben_cc,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


anova(M03_evi, M04_evi)

# so that indicates that the interaction with FPAR can be safely dropped

anova(M04_evi)


# what about soil moisture?

M05_evi <- update(M04_evi, .~. - log(splash.ma))

anova(M04_evi, M05_evi)


# And that confirms that we should retain the soil term.
# So to this point we have T + D + Sm + CI.  But do we need the transformations?

M06_evi <- glmer(lue2_obs ~ temp_scale + vpdDay.ma + splash.ma + Cloud_index + (1 | sitename / year), 
                 family = Gamma(link = "log"), 
                 data = ben_cc)


anova(M04_evi, M06_evi)

# and that confirms that the transformed responses give a much better fit.




####################################
### Bits and pieces for the text ###
####################################


## What proportion of our soilm stress estimates fall in the transition zone > 0 but < 1?  Remember that we jittered zero splash scores to 0.001 in our initial steps above (to accommodate log transformations):

ggplot(data = ben_cc, aes(splash.ma)) +
  geom_histogram(binwidth = 0.05) 

# so the distribution for this averaged, growing season dataset is heavily left-skewed.

ben_cc %>% 
  mutate(stress = if_else(splash.ma == 0.001, "full", if_else(splash.ma == 1, "absent", "transit"))) %>% 
  group_by(stress) %>% 
  summarise(count = n())



###########
### END ###
###########



