### generate "fake" data
library(tidyverse)

## make event times file 

make_event_times = function(n){
  ids = seq(1:n)
  anes_start = as.POSIXct("2019-01-01 00:00", tz = "UTC") + 
    runif(n, min = 0, max = 24 * 60)
  cpb_start = anes_start + runif(n, min = 2 * 60 * 60, max = 4 * 60 * 60)
  axc_start = cpb_start + runif(n, min = 0.25 * 60 * 60, max = 0.5 * 60 * 60)
  axc_end = axc_start + runif(n, min = 2 * 60 * 60, max = 4 * 60 * 60)
  cpb_end = axc_end + runif(n, min = 0.25 * 60 * 60, max = 0.5 * 60 * 60)
  anes_end = cpb_end + runif(n, min = 2 * 60 * 60, max = 4 * 60 * 60)

  eventtimes_wide = 
    tibble(id = ids,
           'start_ADT_0-Anesthesia_0' = anes_start,
           'start_ADT_0-Anesthesia_0-CPB_0' = cpb_start,
           'start_ADT_0-Anesthesia_0-CPB_0-AXC_0' = axc_start,
           'end_ADT_0-Anesthesia_0' = anes_end,
           'end_ADT_0-Anesthesia_0-CPB_0' = cpb_end,
           'end_ADT_0-Anesthesia_0-CPB_0-AXC_0' = axc_end)
           
}

eventtime_wide = make_event_times(n = 100)

if(!dir.exists(here::here("data", "processed"))){
  dir.create(here::here("data", "processed"))
}

write_csv(eventtime_wide, here("data/processed/event_start_stop_times_cohort.csv"))


for(i in 1:nrow(eventtime_wide)){
  temp = eventtime_wide[i,]
  id = temp$id
  start = temp$`start_ADT_0-Anesthesia_0` - as.period(runif(1, max = 10), units = "minutes")
  end = temp$`end_ADT_0-Anesthesia_0` +  as.period(runif(1, max = 10), units = "minutes")
  
  ts = seq.POSIXt(start, end, by = "min")
  
  # > mean(hemo_data$val_MAP, na.rm = TRUE)
  # [1] 74.21569
  # > sd(hemo_data$val_MAP, na.rm = TRUE)
  # [1] 15.91696
  # >   mean(hemo_data$val_CVP, na.rm = TRUE)
  # [1] 8.083316
  # >   sd(hemo_data$val_CVP, na.rm = TRUE)
  # [1] 4.994828
  df = 
    tibble(timestamp = ts,
           ADT_0.Anesthesia_0 = NA_real_,
           ADT_0.Anesthesia_0.CPB_0 = NA_real_,
           ADT_0.Anesthesia_0.CPB_0.AXC_0 = NA_real_,
           CPBFlow = NA_real_,
           CVP = rnorm(n = length(ts), mean = 8.1, sd = 5),
           DBP = NA_real_,
           HR = NA_real_,
           MAP = rnorm(n = length(ts), mean = 74, sd = 15.9),
           SBP = NA_real_,
           MPP = NA_real_,
           id = i)
  
  # inject missingness 
  
  prop_miss = runif(2, 0, .1)
  
  
  na_inds_map = sample(length(ts), length(ts) * prop_miss[1])
  na_inds_cvp = sample(length(ts), length(ts) * prop_miss[2])
  
  df[na_inds_map,which(colnames(df) == "MAP")] <- NA_real_
  df[na_inds_cvp,which(colnames(df) == "CVP")] <- NA_real_
  
  fname = paste0(i, ".csv")
  readr::write_csv(df, here::here("data", "raw", "hemodynamics", "individual_timeseries", fname))
  
}
    
n = nrow(eventtime_wide)
covariates_analysis =
  tibble(id = seq(1:n),
         cat_surgcat16 = NA_character_,
         cat_surgcat4 = NA_character_,
         cat_surgcat4b = NA_character_,
         cat_surgcat7 = NA_character_,
         bin_preopinotrop = rbernoulli(n, p = 0.2),
         bin_afib = rbernoulli(n, p = 0.1) * 1,
         val_egfr = rnorm(n, 80, 30),
         val_age = rnorm(n, 60, 13),
         cat_gender = c(rep("Female", n * .5), rep("Male", n * .5)),
         val_bmi = rnorm(n, 10, 30),
         bin_hypertn = rbernoulli(n, p = 0.25) * 1,
         bin_dm = rbernoulli(n, p = 0.3) * 1,
         bin_copd = rbernoulli(n, p = 0.15) * 1,
         bin_chf = rbernoulli(n, p = 0.1) * 1,
         bin_priormi = rbernoulli(n, p = 0.35) * 1,
         bin_cva = rbernoulli(n, p = 0.2)* 1,
         bin_emergent = rbernoulli(n, p = 0.08)* 1,
         bin_redo = rbernoulli(n, p = 0.2)* 1,
         bin_pvd = rbernoulli(n, p = 0.18)* 1,
         val_hct = rbernoulli(n, p = 0.2)* 1,
         val_hdef = rnorm(n, 39, 6),
         bin_statin = rbernoulli(n, p = 0.5)* 1,
         bin_acearb = rbernoulli(n, p = 0.4)* 1,
         bin_betablocker = rbernoulli(n, p = 0.6)* 1,
         bin_iabp =rbernoulli(n, p = 0.08)* 1,
         val_perfustm = rnorm(n, 110, 50), 
         val_visscore = rnorm(n, 0.4, 0.4), 
         val_creatlst = rnorm(n, 1.1, 0.4), 
         cat_transbc = c(rep("1", 30), rep("0", 30), rep(">1", 40)),
         val_transplasma = NA_character_,
         val_transplatelet = NA_character_,
         cat_disloctn = NA_character_,
         val_predmort = rnorm(n, 0.02, 0.03), 
         val_predrenf = rnorm(n, 0.02, 0.03), 
         cat_mt30stat = NA_character_,
         val_crystalloid = rnorm(n, 2000, 2000), 
         cat_aki7d = NA_character_,
         bin_aki48h = rbernoulli(n, p = 0.15)* 1,
         cat_race = c(rep("Asian", 10), rep("Black", 30), rep("Caucasian", 50), rep("Other", 10)),
         bin_readmit = rbernoulli(n, p = 0.15)* 1,
         dtm_anesstart = eventtime_wide$`start_ADT_0-Anesthesia_0`,
         dtm_anesend = eventtime_wide$`end_ADT_0-Anesthesia_0`,
         dtm_hospadm = NA_real_,
         dtm_hospdis = NA_real_,
         val_los = NA_real_
  )

covariates_analysis <-
  covariates_analysis %>% 
  mutate(val_proctime = as.numeric(difftime(dtm_anesend, dtm_anesstart, units = "mins")))

write_csv(covariates_analysis, here("data/processed/covariates_cohort.csv"))




    