library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(purrr)
library(stringr)
library(haven)
library(here)
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

# get filenames of individual time series 
file_names <-
  list.files(here::here("data/raw/hemodynamics/individual_timeseries"), full.names = T)

# make sure to only get files that are csvs 
csvs <- file_names[grep(".csv", file_names)]

# function to read in files and make an id column
read.files = function(x) {
  read.csv(x) %>% 
    mutate(id = sub(".*individual_timeseries/", "", sub(".csv.*", "", x)))
}

# read in all files and bind rows of resulting list
# make time POSIXct 

all_hemo_timeseries <-
  map(.x = csvs,
      .f = read.files)

all_hemo_timeseries <-
  map(.x = csvs,
      .f = read.files) %>%
  list_rbind() %>%
  mutate(
    time = ymd_hms(timestamp)
  ) %>% 
  select(-timestamp)

if(!dir.exists(here::here("data", "processed"))){
  dir.create(here::here("data", "processed"))
}

# save one dataframe that is just all hemodynamic time series row-bound 
saveRDS(all_hemo_timeseries, here("data/processed/hemo_timeseries_all.rds"))

# create another dataframe with individuals eligible for the study 


# create a wide dataset, so that there is one row for each person with the start and end of each event
eventtimes_wide = read_csv(here("data/processed/event_start_stop_times_cohort.csv"))

cohort_ids = unique(eventtimes_wide$id)
# remove people with no anesthesia because it will throw error for next part 
no_anesthesia <-
  eventtimes_wide %>%
  filter(is.na(`start_ADT_0-Anesthesia_0`) |
           is.na(`end_ADT_0-Anesthesia_0`)) %>%
  pull(id)


cohort_ids_anes <- cohort_ids[cohort_ids %notin% no_anesthesia]

# function to create time series with rows for every minute between start and end of anesthesia,
# since Zach's original processing allows for rows to be left out if they don't have a recorded measurement

create_full_ts <-
  function(subject) {
    # get row for the subject
    times <- eventtimes_wide %>%
      filter(id == subject)
    # time sequence from start to end of anesthesia by minute
    
    # we don't exclude for multiple bypass, so we just find the first bypass time here
    start_cpb <-
      min(
        times$`start_ADT_0-Anesthesia_0-CPB_0`,
        na.rm = TRUE
      )
    # and the last bypass time
    end_cpb <-
      max(
        times$`end_ADT_0-Anesthesia_0-CPB_0`,
        na.rm = TRUE
      )
    # create data frame with row for every minute between start and end of anesthesia
    # create indicators for whether each timestamp is pre, intra, or post bypass
    # by definition all times are intra-anesthesia
    # if there's no bypass, these people will be excluded later but for now we just make CPB NA
    
    all_hemo_timeseries %>% 
      filter(id == subject) %>% 
      mutate(
        CPB = case_when(
          is.infinite(start_cpb) ~ NA,
          time < start_cpb ~ "pre",
          time > end_cpb ~ "post",
          time >= start_cpb &
            time <= end_cpb ~ "intra"
        ),
        anes = "intra"
      )

  }

# create the "full" time series data set 
full_timeseries <-
  map_dfr(.x = cohort_ids_anes,
      .f = create_full_ts) 

# join with the hemodynamics time series 
saveRDS(full_timeseries, file = here("data/processed/hemo_timeseries_cohort.rds"))


# now do interpolation 
all_hemo <- readRDS(here("data/processed/hemo_timeseries_cohort.rds"))

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

saveRDS(hemo_interp, here("data/processed/hemo_timeseries_interp_cohort.rds"))

