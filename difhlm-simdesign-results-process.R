#08mar2018
#process results.rds object for later plotting and analysis
#creates 'res' object as well as some constants

# library(dplyr)
# library(reshape2)

# read and subset data
readRDS(file = "difhlm-simdesign-results.rds") %>%
  as.data.frame %>%
  select(-contains("WARNING")) %>%
  rename(rc_fail = `ERROR: .Error in getModels(d = dat, g5 = 0, ts_type = "ts", difitem = it) : \n  rc failed to converge it25\n`,
         ri_fail = `ERROR: .Error in getModels(d = dat, g5 = 0, ts_type = "ts", difitem = it) : \n  ri failed to converge it25\n`,
         rc_ri_fail = `ERROR: .Error in getModels(d = dat, g5 = 0, ts_type = "ts", difitem = it) : \n  rc and ri failed to converge it25\n`) %>% 
  mutate(rc_fail = ifelse(is.na(rc_fail), 0, rc_fail),
         ri_fail = ifelse(is.na(ri_fail), 0, ri_fail),
         rc_ri_fail = ifelse(is.na(rc_ri_fail), 0, rc_ri_fail),
         rc_fail_tot = rc_fail+rc_ri_fail,
         ri_fail_tot = ri_fail+rc_ri_fail,
         fail_tot = rc_fail+ri_fail+rc_ri_fail,
         tot_reps=fail_tot+200) -> res


# some constants

t1_limit <- 0.08 # very close to the SE(p) = 0.05+2*sqrt((0.05*0.95)/200)

condition_vars <- c('dif_lab', 'sample_size', 'bval_lab', 'difvar_lab', 'icc_lab',
                    'impact_lab')

condition_vars_nodif <- c('sample_size', 'bval_lab', 'difvar_lab', 'icc_lab',
                          'impact_lab')


# label variables

res$sample_size <- factor(paste0("J=", res$J, ", N=", res$N),
                          levels = c("J=25, N=10", "J=100, N=10",
                                     "J=25, N=40", "J=100, N=40"))
res$dif_lab    <- paste0("DIF=", res$dif)
res$bval       <- factor(res$bval, levels = c("low", "med", "high"))
res$bval_lab   <- factor(paste0("b=", res$bval),
                         levels = c("b=low", "b=med", "b=high"))
res$difvar_lab <- paste0("Var(DIF)=", res$difvar)
res$icc_lab    <- paste0("ICC=", res$icc)
res$impact_lab <- paste0("Impact=", res$impact)
res$truedif    <- ifelse(res$dif=="with", -0.6, 0)
res$truedifvar <- ifelse(res$difvar=="with", 0.80, 0)

res$`Sample Size`  <- res$sample_size
res$Difficulty     <- res$bval
res$`DIF Variance` <- res$difvar
res$Impact         <- res$impact
res$ICC            <- res$icc
res$DIF            <- res$dif


################################################################################
# create data frame with DIF type i error rates and power
res %>% filter(dif=="no") %>%
  select(condition_vars_nodif,
         mn.lr_dif_sig, mn.ri_x_g_dif_sig, mn.rc_x_g_dif_sig) %>%
  rename(LR=mn.lr_dif_sig,
         RI=mn.ri_x_g_dif_sig,
         RC=mn.rc_x_g_dif_sig) %>%
  melt(., id.vars = condition_vars_nodif,
       variable.name = 'model', value.name = 't1') -> dif_t1

# power rates
res %>% filter(dif=="with") %>%
  select(condition_vars_nodif,
         mn.lr_dif_sig, mn.ri_x_g_dif_sig, mn.rc_x_g_dif_sig) %>%
  rename(LR=mn.lr_dif_sig,
         RI=mn.ri_x_g_dif_sig,
         RC=mn.rc_x_g_dif_sig) %>%
  melt(., id.vars = condition_vars_nodif,
       variable.name = 'model', value.name = 'pwr') -> dif_pwr

# merge together
dif_t1 %>%
  merge(x = ., y = dif_pwr, all = TRUE, by = c(condition_vars_nodif, 'model')) %>%
  mutate(inflated = ifelse(t1>t1_limit, 1, 0),
         inflated_star = ifelse(t1>t1_limit, "*", "")) -> difsig


rm(dif_pwr, dif_t1)
