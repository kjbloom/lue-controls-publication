##########################################################################################
### Model experiments testing the functional form of environmental dependencies of LUE ###
###               Producing Figure 2 in Bloomfield et al. 2022, GCB                    ###
##########################################################################################


## Calculate the functional forms for each of the empirical model's explanatory variables - temperature, vpd, soil moisture and Cloudiness_index. This is done using conditional plots available within the visreg package.

## Two versions are produced: 
#(i) adopting eddy-covariance inferred values of LUE as response 
#(ii) adopting P-model simulations of LUE as response variable

## The dataset is also used to conduct P-model experiments: successive iterations look at individual effects on simulated LUE whilst holding All other model inputs at median value (this mimics the working of the conditional plots).  There is no step here for Cloudiness_index because the P-model currently has no formulation (or input) for diffuse radiation.


library(lme4)
library(visreg)
library(gridExtra)
library(grid)
library(rpmodel) # this is v1-0-4 


# Read in the cleaned-up data-set used in the main analysis - confined to growing season and time-averaged:

source("load_GPP-dependencies.R")


# Calculate median values for the explanatory variables necessary to the P-model.  

T_med <- signif(median(ddf_3$tDay.ma), digits = 4)
D_med <- signif(median(ddf_3$vpdDay.ma), digits = 5) * 1000  # we convert these back to Pa
fapar_med <- signif(median(ddf_3$fpar.ma), digits = 3)
patm_med <- signif(median(ddf_3$patm.ma), digits = 5)
soilm_med <- signif(median(ddf_3$splash.ma), digits = 3)
Ca_med <- signif(median(ddf_3$co2_ML.ma), digits = 4) # annual global averages, per Mauna Loa
ai_med <- signif(median(ddf_3$alpha.ma), digits = 3) # ratio of Actual to Potential Evapotranspiration



## and we want regular sequences for the explanatory variables retained in our final empirical model:

range(ddf_3$tDay.ma)

t_seq <- seq(from = -20, to = 36, length.out = 100)

range(ddf_3$vpdDay.ma)

vpd_seq <- seq(from = 20, to = 5000, length.out = 100)

soilm_seq <- seq(from = 0, to = 1, length.out = 100)


## and now we construct a dummy data-frame:

modex <- tibble(
  tc = t_seq,
  vpd = vpd_seq,
  soilm = soilm_seq,
  co2 = Ca_med,
  fapar = fapar_med,
  patm = patm_med,
  meanalpha = ai_med 
  )



## And now the three experimental iterations; 
# this is the 'FULL' implementation of the P-model (Stocker et al. 2020):

# (i) Allowing only T to vary:

mod_T <- modex %>% 
  mutate(vpd = D_med, soilm = soilm_med) %>% 
  mutate(out_pmodel = purrr::pmap(., rpmodel,
                                  
                                  ppfd = NA, # optional if we want GPP estimates - ignore here
                                  kphio = 0.0870, # tailored for the FULL implementation 
                                  beta = 146,
                                  apar_soilm = 0,
                                  bpar_soilm = 0.685, # Table 1, Stocker et al. 2020
                                  c4 = FALSE,
                                  method_optci = "prentice14",
                                  method_jmaxlim = "wang17",
                                  do_ftemp_kphio = TRUE, 
                                  do_soilmstress = TRUE,
                                  returnvar = "lue"
                                  ) 
         )





# (ii) Allowing only VPD to vary:

mod_D <- modex %>% 
  mutate(tc = T_med, soilm = soilm_med) %>% 
  mutate(out_pmodel = purrr::pmap(., rpmodel,
                                  
                                  ppfd = NA,
                                  kphio = 0.0870, # tailored for the FULL implementation 
                                  beta = 146,
                                  apar_soilm = 0,
                                  bpar_soilm = 0.685, # Table 1, Stocker et al. 2020
                                  c4 = FALSE,
                                  method_optci = "prentice14",
                                  method_jmaxlim = "wang17",
                                  do_ftemp_kphio = TRUE, 
                                  do_soilmstress = TRUE,
                                  returnvar = "lue" 
                                  ) 
         )




# (iii) Allowing only soil moisture stress to vary:

mod_Sm <- modex %>% 
  mutate(tc = T_med, vpd = D_med) %>% 
  mutate(out_pmodel = purrr::pmap(., rpmodel,
                                  
                                  ppfd = NA,
                                  kphio = 0.0870, # tailored for the FULL implementation 
                                  beta = 146,
                                  apar_soilm = 0,
                                  bpar_soilm = 0.685, # Table 1, Stocker et al. 2020
                                  c4 = FALSE,
                                  method_optci = "prentice14",
                                  method_jmaxlim = "wang17",
                                  do_ftemp_kphio = TRUE, 
                                  do_soilmstress = TRUE,
                                  returnvar = "lue" 
                                  ) 
         )


## Now combine all of those runs into a single tibble.  As a first step un-nest the output so that the model estimates appear in separate (and conventional) columns:
  
  
mod_T <- mod_T %>% 
  mutate(out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>% 
  unnest(out_pmodel) %>% 
  mutate(lue_Tvar = lue / mass_C)


mod_D <- mod_D %>% 
  mutate(out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>% 
  unnest(out_pmodel) %>% 
  mutate(lue_Dvar = lue / mass_C)


mod_Sm <- mod_Sm %>% 
  mutate(out_pmodel = purrr::map(out_pmodel, ~as_tibble(.))) %>% 
  unnest(out_pmodel) %>% 
  mutate(lue_Smvary = lue / mass_C)



# and now combine these into a single object - we simply bind the columns:

modex_02 <- cbind(modex, mod_T[ ,9], mod_D[ ,9], mod_Sm[ ,9])



######
### Now we compare that with our empirical model: and we recreate the final statistical model explaining variation in inferred GPP; see model selection steps in companion script "lue-dependencies-main-analysis.R'


# and I need a complete case dataset for these runs:

ben_cc <- ddf_3 %>%
  
  # we prune out redundant columns (remember to retain 'year' here for site_year constructs)
  select(sitename, year,
         lue3_obs, lue2_mod_FULL,
         tDay.ma, vpdDay.ma, splash.ma, Cloud_index) %>%
  
  # replace zero soilm estimates with a nominal value to allow log transformation:
  mutate(splash.ma = ifelse(splash.ma < 0.001, 0.001, splash.ma),
         temp_scale = scale(tDay.ma, center = T, scale = T)) %>% 
  
  # omit rows with missing values:
  na.omit() %>%
  
  # drop redundant factor levels (e.g. Sites)
  droplevels()


# here's the final model construct.  

M_fin <- glmer(lue3_obs ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                 (1 | sitename / year), 
               family = Gamma(link = "log"), 
               data = ben_cc,
               control = glmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))



## Next rather than 'observed' GPP, we are interested in the functional forms when our response variable is generated by the P-model and here the response variable again relies on FPAR ##


ff_full <- glmer(lue2_mod_FULL ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index + 
                   (1 | sitename / year), 
                 family = Gamma(link = "log"), 
                 data = ben_cc,
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5)))


## Notice how much weaker the vpd effect is here - compared with our empirical model (M_fin above). 
## AND the proxy for diffuse radiation (Cloud_index) has NO role ##


summary(ff_full)



### Figure 2 ###

### We can visualise those effects with plots of the models' partial residuals.  And for presentation purposes, we want a 4 x 2 plot.  To combine those plots effectively, we include a ggplot argument in the visreg function:


kj.lab_size <- theme(axis.title.x=element_text(size=15), 
                     axis.title.y=element_text(size=15), 
                     axis.text.x=element_text(size=12), 
                     axis.text.y=element_text(size=12))


## the top panel shows the effects as they relate to the final empirical model:

plot_a <- visreg(M_fin, "tDay.ma", type = "conditional", trans = exp, rug = F, partial = T, gg = T, 
                 points = list(pch = 1, cex = 1.0),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", 
                 xlab = " ") + 
  kj.lab_size + 
  xlim(0, 36) +
  ylim(0, 0.10) +
  geom_text(data=NULL, x=0, y=Inf, label="FLUXNET 2015: Inferred", size=5, colour="black", hjust=0, vjust=1.5)



plot_b <- visreg(M_fin, "vpdDay.ma", type = "conditional", trans = exp, rug = F, partial = T, gg = T,
                 points = list(pch = 1),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", xlab = " ") + 
  kj.lab_size + 
  xlim(0, 5) +
  ylim(0, 0.10) 


plot_c <- visreg(M_fin, "splash.ma", type = "conditional", trans = exp, rug = F, partial = T, gg = T, 
                 points = list(pch = 1),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", xlab = " ") + 
  kj.lab_size + 
  xlim(0, 1) +
  ylim(0, 0.10)



plot_d <- visreg(M_fin, "Cloud_index", type = "conditional", trans = exp, rug = F, partial = T, gg = T, 
                 points = list(pch = 1),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", xlab = " ") + 
  kj.lab_size + 
  xlim(0, 1) +
  ylim(0, 0.10)




## the bottom panel shows equivalent plots with P-model estimates as response variable:

plot_e <- visreg(ff_full, "tDay.ma", type = "conditional", trans = exp, rug = F, partial = T, gg = T, 
                 points = list(pch = 1),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", 
                 xlab = "Daytime temperature (°C)") + 
  kj.lab_size + 
  xlim(0, 36) +
  ylim(0, 0.05) +
  geom_text(data=NULL, x=0, y=Inf, label="P-model ('FULL')", size=5, colour="black", hjust=0, vjust=1.5) +
  # and we add a line to show the response according to our model experiment:
  geom_smooth(data = modex_02, aes(tc, lue_Tvar), 
              method = "gam", formula = y ~ s(x, bs = "cs"), 
              colour = "blue")



plot_f <- visreg(ff_full, "vpdDay.ma", type = "conditional", trans = exp, rug = F, partial = T, gg = T, 
                 points = list(pch = 1),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", 
                 xlab = "Daytime VPD (kPa)") + 
  kj.lab_size + 
  xlim(0, 5) +
  ylim(0, 0.05) +
  # and we superimpose the line for the modelling experiment
  geom_smooth(data = modex_02, aes(vpd / 1000, lue_Dvar), 
              method = "gam", formula = y ~ s(x, bs = "cs"), 
              colour = "blue")
  


plot_g <- visreg(ff_full, "splash.ma", type = "conditional", trans = exp, rug = F, partial = T, gg = T, 
                 points = list(pch = 1),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", 
                 xlab = "Soil moisture stress (SPLASH)") + 
  kj.lab_size + 
  xlim(0, 1) +
  ylim(0, 0.05) +
  geom_smooth(data = modex_02, aes(soilm, lue_Smvary), 
              method = "gam", formula = y ~ s(x, bs = "cs"), 
              colour = "blue")


# the P-model has no explicit diffuse radiation term and so confine ourselves here to the conditonal plot
plot_h <- visreg(ff_full, "Cloud_index", type = "conditional", trans = exp, rug = F, partial = T, gg = T, 
                 points = list(pch = 1),
                 line = list(lty = 1, col = "red"),
                 ylab = " ", 
                 xlab = "Cloudiness Index") + 
  kj.lab_size + 
  xlim(0, 1) +
  ylim(0, 0.05)




       
fig2 <- grid.arrange(
         arrangeGrob(grobs = list(plot_a, plot_b, plot_c, plot_d, plot_e, plot_f, plot_g, plot_h), 
                     ncol = 4, nrow = 2, 
                     left = textGrob(label = expression("LUE "~ (molC~mol^{-1}~photons)), rot = 90, 
                                     gp = gpar(fontface = "bold", fontsize = 15)
                                     ),
                     # reduce the margin around the left axis label
                     padding = unit(0.1, "line")
                     )
         )



### END ###


