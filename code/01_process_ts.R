library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(purrr)
library(stringr)
library(haven)
library(here)
rm(list = ls())
`%notin%` <- Negate(`%in%`)
options(dplyr.summarise.inform=FALSE)

# get filenames of individual time series 
file_names =
  list.files(here::here("data/raw/hemodynamics/individual_timeseries"), full.names = T)

# make sure to only get files that are csvs 
csvs <- file_names[grep(".csv", file_names)]

# function to read in files and make an id column
read.files = function(x) {
  read_csv(x) %>% 
    mutate(id = as.character(sub(".*individual_timeseries/", "", sub(".csv.*", "", x))))
}

# read in all files and bind rows of resulting list
# make time POSIXct 

all_hemo_timeseries =
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
saveRDS(all_hemo_timeseries, here("data", "processed", "hemo_timeseries_all.rds"))

# create another dataframe with individuals eligible for the study 


# create a wide dataset, so that there is one row for each person with the start and end of each event
eventtimes_wide = read_csv(here("data", "raw", "event_start_stop_times.csv"))

cohort_ids = unique(eventtimes_wide$id)
# remove people with no anesthesia because it will throw error for next part (we don't have this in the simulated data)
no_anesthesia =
  eventtimes_wide %>%
  filter(is.na(`start_ADT_0-Anesthesia_0`) |
           is.na(`end_ADT_0-Anesthesia_0`)) %>%
  pull(id)


cohort_ids_anes = cohort_ids[cohort_ids %notin% no_anesthesia]

# function to create time series with rows for every minute between start and end of anesthesia,

create_full_ts =
  function(subject) {
    # get row for the subject
    times = eventtimes_wide %>%
      filter(id == subject)
    # time sequence from start to end of anesthesia by minute
    
    # we don't exclude for multiple bypass, so we just find the first bypass time here
    start_cpb =
      min(
        times$`start_ADT_0-Anesthesia_0-CPB_0`,
        na.rm = TRUE
      )
    # and the last bypass time
    end_cpb =
      max(
        times$`end_ADT_0-Anesthesia_0-CPB_0`,
        na.rm = TRUE
      )
    
    timeseq =
      seq.POSIXt(times$`start_ADT_0-Anesthesia_0`,
                 times$`end_ADT_0-Anesthesia_0`,
                 by = "mins")
    # create data frame with row for every minute between start and end of anesthesia
    # create indicators for whether each timestamp is pre, intra, or post bypass
    # by definition all times are intra-anesthesia
    # if there's no bypass, these people will be excluded later but for now we just make CPB NA
    df = tibble(id = rep(subject, length(timeseq)),
                     time = timeseq) %>%
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
    df %>% 
      mutate(id = as.character(id))
    

}

# create the "full" time series data set 
full_timeseries =
  map_dfr(.x = cohort_ids_anes,
      .f = create_full_ts) 

full_ts = 
  full_timeseries %>% 
  left_join(all_hemo_timeseries, 
            by = c("id", "time"))

# join with the hemodynamics time series 
saveRDS(full_ts, file = here("data", "processed", "hemo_timeseries_cohort.rds"))

