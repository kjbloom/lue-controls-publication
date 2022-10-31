###########################################################
###   Cross-validation of the final empirical model     ###
### Supplementary figure 6, Bloomfield et al. 2020, GCB ###
###########################################################


## load the necessary packages:

library(lme4)
library(reshape2)


# Read in the cleaned-up data-set used in the main analysis - confined to growing season and time-averaged:

source("load_GPP-dependencies.R")


# remove some items not needed here
rm(ddf, ddf_2, ddf_grow, mass_C)


## to evaluate the predicted values against the observations we employ Beni Stocker's function:

source("functions/analyse_modobs.R")


## We have missing values in the core file - especially related to flux site measures of soil moisture; so a first step is to prune the dataframe to the variables required for the final model.

summary(ddf_3)


ben_cc <- ddf_3 %>%
  
  # we retain only the variables of interest:
  select(sitename, year, tDay.ma, vpdDay.ma, splash.ma, Cloud_index, lue3_obs) %>%
  
  # omit rows with missing values:
  na.omit() %>%
  
  # drop redundant factor levels (e.g. Sites)
  droplevels()
  
 

#I'm going to replace zero soilm estimates with a nominal value to aid later computation (e.g. log transformation)
ben_cc$splash.ma <- with(ben_cc, ifelse(splash.ma < 0.001, 0.001, splash.ma))


# We have a final agreed model structure (M_fin, main analysis).  The cross-validation exercise, however, in which we sequentially drop one site from the data used to train the model, is unable to generate site-level random effects and so here we adopt a simpler GLM that excludes any mixed effects:


form <- formula(lue3_obs ~ poly(tDay.ma, 2) + log(vpdDay.ma) + log(splash.ma) + Cloud_index)


## Now the idea is to sequentially drop one site from the dataset.  The remaining df becomes 'training' and the single excluded site is 'evaluation'.

# We run our model for the training df and then apply those coefficients to create a predicted LUE for the single site under evaluation.

# We repeat those steps, generating predictions for each evaluation site in turn.  And combine all the evaluation predictions into a single dataframe.
  

all_sites <- unique(ben_cc$sitename)

oob_preds <- list()


for (i in seq_along(all_sites)) {
  
  eval_site <- all_sites[i]
  
  train_df <- ben_cc %>%
    subset(sitename != eval_site)
  
  eval_df <- ben_cc %>%
    subset(sitename == eval_site)
  
  train_mod <- glm(form, 
                   family = Gamma(link = "log"), 
                   data = train_df)
  
  oob_preds[[i]] <- predict(train_mod, newdata = eval_df, re.form = NULL, type = "response") 
  
  }



## We use another Hadley Wickham function to convert this list into a dataframe:
pred_lue <- melt(oob_preds)


ben_cc <- ben_cc %>%  
  mutate(pred_cv = pred_lue$value)


# Now for the figure:

with(ben_cc, analyse_modobs(pred_cv, lue3_obs, heat = T, 
                            plot.title = "Statistical model cross-validation (no random term)",
                            xlab = "Predicted LUE (glm)",
                            ylab = expression("Measured LUE"~ (molC~mol^{-1}~photons)),
                            xlim = c(0, 0.10)))  


# clean-out
rm(oob_preds, pred_lue, eval_df, train_df, train_mod, eval_site, form, i)


### END ###



   










