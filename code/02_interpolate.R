library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(purrr)
library(stringr)
library(haven)
library(here)
rm(list = ls())
rle2 <- function (x)  {
  stopifnot("'x' must be a vector of an atomic type" = is.atomic(x))
  
  n <- length(x)
  if (n == 0L) {
    return(structure(list(
      lengths = integer(), values = x)
    ), class = 'rle')
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Where does next value not equal current value?
  # i.e. y will be TRUE at the last index before a change
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y <- (x[-1L] != x[-n])
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Since NAs are never equal to anything, NAs in 'x' will lead to NAs in 'y'.
  # These current NAs in 'y' tell use nothing - Set all NAs in y to FALSE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y[is.na(y)] <- FALSE
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # When a value in x is NA and its successor is not (or vice versa) then that
  # should also count as a value change and location in 'y' should be set to TRUE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y <- y | xor(is.na(x[-1L]), is.na(x[-n]))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Any TRUE locations in 'y' are the end of a run, where the next value
  # changes to something different
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  i <- c(which(y), n)
  
  structure(list(
    lengths = diff(c(0L, i)),
    values  = x[i]
  ), class = 'rle')
}
`%notin%` <- Negate(`%in%`)
options(dplyr.summarise.inform=FALSE)

# now do interpolation 
all_hemo = readRDS(here("data", "processed", "hemo_timeseries_cohort.rds"))

# CVP interpolation 
cvp_approx <-
  all_hemo %>%
  group_by(id) %>%
  mutate(first_meas = min(time[!is.na(CVP)]),
         last_meas = max(time[!is.na(CVP)])) %>%
  filter(time >= first_meas & time <= last_meas) %>% # filter to measurements after first valid one 
  mutate(
    time_num = seq(1, n()), # x's for interpolation 
    rollCVP = case_when( # conditional rolling mean - previous 3 values if before NA, next 3 if after NA, NA if NA 
      is.na(CVP) ~ NA,
      is.na(lead(CVP, 1)) &
        !is.na(lag(CVP, 1)) ~ zoo::rollmean(
          CVP,
          k = 3,
          fill = NA,
          align = "right"
        ),
      !is.na(lead(CVP, 1)) &
        is.na(lag(CVP, 1)) ~ zoo::rollmean(
          CVP,
          k = 3,
          fill = NA,
          align = "left"
        )
    ),
    CVP_approx = approx( # do the interpolation 
      x = time_num,
      y = rollCVP,
      xout = time_num,
      method = "linear"
    )$y
  )

# run length encoding to find length of NA stretch 
cvp_approx$rle <- rep(rle2(cvp_approx$CVP)$lengths, times = rle2(cvp_approx$CVP)$lengths)

# get rid of values if the length of missingness is > 15 min 
cvp_approx=cvp_approx  %>%
  mutate(
    CVP_final = ifelse(is.na(CVP) & rle <= 15, CVP_approx, CVP))

# repeat the process for MAP 
MAP_approx <-
  all_hemo %>%
  group_by(id) %>%
  mutate(first_meas = min(time[!is.na(MAP)]),
         last_meas = max(time[!is.na(MAP)])) %>%
  filter(time >= first_meas & time <= last_meas) %>% # filter to measurements after first valid one 
  mutate(
    time_num = seq(1, n()), # x's for interpolation 
    rollMAP = case_when( # conditional rolling mean - previous 3 values if before NA, next 3 if after NA, NA if NA 
      is.na(MAP) ~ NA,
      is.na(lead(MAP, 1)) &
        !is.na(lag(MAP, 1)) ~ zoo::rollmean(
          MAP,
          k = 3,
          fill = NA,
          align = "right"
        ),
      !is.na(lead(MAP, 1)) &
        is.na(lag(MAP, 1)) ~ zoo::rollmean(
          MAP,
          k = 3,
          fill = NA,
          align = "left"
        )
    ),
    MAP_approx = approx( # do the interpolation 
      x = time_num,
      y = rollMAP,
      xout = time_num,
      method = "linear"
    )$y
  )

# run length encoding to find length of NA stretch 
MAP_approx$rle <- rep(rle2(MAP_approx$MAP)$lengths, times = rle2(MAP_approx$MAP)$lengths)

# get rid of values if the length of missingness is > 15 min 
MAP_approx=MAP_approx %>%
  mutate(
    MAP_final = ifelse(is.na(MAP) & rle <= 15, MAP_approx, MAP))

# join back with original df 
hemo_interp <-
  all_hemo %>%
  left_join(., MAP_approx %>% dplyr::select(id, time, MAP_approx, MAP_final, rle) %>% rename(rle_map = rle), by = c("id" = "id", "time" = "time")) %>%
  left_join(., cvp_approx %>% dplyr::select(id, time, CVP_approx, CVP_final, rle) %>% rename(rle_cvp = rle), by = c("id" = "id", "time" = "time"))

hemo_interp <-
  hemo_interp %>% 
  select(id, time, cat_cpb = CPB, cat_anes = anes, 
         val_CPBFlow = CPBFlow, val_DBP = DBP, val_HR = HR, 
         val_MAP = MAP_final, val_CVP = CVP_final, 
         val_SBP = SBP, val_MPP = MPP)

saveRDS(hemo_interp, here("data", "processed", "hemo_timeseries_interp_cohort.rds"))

