### There are different ways to average these variables (e.g. weekly, monthly); many references seem to employ a 15 day window e.g. Reichstein et al. (2005), Keenan et al. (2019).  

## Colin prefers that we do not employ a moving average with overlaps and so width and by arguments now set to the same 15 day span:

mafun <- function (x) {
  rollapply(x, width = 15, by = 15, mean, na.rm = T, fill = NA, partial = F, align = "right")
}


# and a version for sums, instead of means (this for our revised LUE calculation):

sum_fun <- function (x) {
  rollapply(x, width = 15, by = 15, sum, na.rm = T, fill = NA, partial = F, align = "right")
}


