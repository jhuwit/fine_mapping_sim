library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(purrr)
library(stringr)
library(haven)
library(readxl)
library(here)
`%notin%` <- Negate(`%in%`)
options(dplyr.summarise.inform=FALSE)

# read in hemodynamics data 
all_hemo <- readRDS(here("data/processed/hemo_timeseries_interp_cohort.rds"))

starting_cohort <-
  all_hemo %>% 
  select(id) %>% 
  distinct() %>% 
  pull(id)



# read in covariates data
covariates <- read_csv(here("data/processed/covariates_cohort.csv")) 

cohort_event_start_stop_times = read_csv(here("data/processed/event_start_stop_times_cohort.csv")) %>% 
  filter(id %in% starting_cohort)



renal_failure <-
  covariates %>%
  filter(val_creatlst > 4) %>%
  pull(id)




get_missing_anesthesia <- function(event_df){
  event_df %>%
    filter(is.na(`start_ADT_0-Anesthesia_0`) |
             is.na(`end_ADT_0-Anesthesia_0`)) %>%
    pull(id)
}

missing_anesthesia = get_missing_anesthesia(cohort_event_start_stop_times)

# no bypass 
get_missing_cpb <- function(hemodf){
  hemodf %>% 
    group_by(id) %>% 
    summarize(cpb = sum(cat_cpb == "intra", na.rm = TRUE)) %>% 
    filter(cpb == 0) %>% 
    pull(id) %>% 
    unique()
}

missing_cpb = get_missing_cpb(all_hemo)



# no non-missimg MAP during anesthesia 

missing_map_all <-
  all_hemo %>% 
  filter(cat_anes == "intra") %>% 
  group_by(id) %>% 
  summarize(non_na_MAP = sum(!is.na(val_MAP))) %>% 
  filter(non_na_MAP == 0) %>% 
  pull(id)

# no non-missing CVP during anesthesia 
missing_cvp_all <-
  all_hemo %>% 
  filter(cat_anes == "intra") %>% 
  group_by(id) %>% 
  summarize(non_na_CVP = sum(!is.na(val_CVP))) %>% 
  filter(non_na_CVP == 0) %>% 
  pull(id)


# missing MAP or CVP during pre, intra, or post bypass 
get_missing_period <- function(variable, period, hemodf){
  hemodf %>% 
    rename(var = all_of(variable)) %>% 
    filter(cat_anes == "intra", cat_cpb == period) %>% 
    group_by(id) %>% 
    summarize(non_na = sum(!is.na(var))) %>% 
    filter(non_na == 0) %>% 
    pull(id) 
}

missing_map_precpb <- get_missing_period("val_MAP", "pre", all_hemo)
missing_map_intracpb <- get_missing_period("val_MAP", "intra", all_hemo)
missing_map_postcpb <- get_missing_period("val_MAP", "post", all_hemo)

missing_cvp_precpb <- get_missing_period("val_CVP", "pre", all_hemo)
missing_cvp_intracpb <- get_missing_period("val_CVP", "intra", all_hemo)
missing_cvp_postcpb <- get_missing_period("val_CVP", "post", all_hemo)


# last one: both >50% missing and <20 minutes of MAP or CVP in one period
# denominator is minutes between first and last arterial BP measurement 
# get df that's only between first and last arterial BP measurement 

between_abp <- 
  all_hemo %>%
  group_by(id) %>%
  mutate(start_MAP = min(time[!is.na(val_MAP)]),
         end_MAP = max(time[!is.na(val_MAP)])) %>%
  filter(time <= end_MAP & time >= start_MAP)

# summarize the percent missing and number of minutes present for each id and period
missing_summary <- 
  between_abp %>%
  group_by(id, cat_cpb) %>%
  summarize(
    elapsed = as.numeric(difftime(max(time), min(time), units = "mins"))+1,
    na_map = sum(is.na(val_MAP)),
    na_cvp = sum(is.na(val_CVP)),
    present_map = sum(!is.na(val_MAP)),
    present_cvp = sum(!is.na(val_CVP))
  ) %>%
  ungroup() %>% 
  mutate(across(starts_with("na"), ~.x/elapsed*100, .names = "{.col}_pct")) %>% 
  rowwise() %>% 
  mutate(missing_map_50_20 = ifelse(na_map_pct > 50 & present_map < 20, 1, 0),
         missing_cvp_50_20 = ifelse(na_cvp_pct > 50 & present_cvp < 20, 1, 0))


missing_map_50_20_precpb <- 
  missing_summary %>% 
  filter(cat_cpb == "pre" & missing_map_50_20 == 1) %>% 
  pull(id)

missing_map_50_20_intracpb <- 
  missing_summary %>% 
  filter(cat_cpb == "intra" & missing_map_50_20 == 1) %>% 
  pull(id)

missing_map_50_20_postcpb <- 
  missing_summary %>% 
  filter(cat_cpb == "post" & missing_map_50_20 == 1) %>% 
  pull(id)

missing_cvp_50_20_precpb <- 
  missing_summary %>% 
  filter(cat_cpb == "pre" & missing_cvp_50_20 == 1) %>% 
  pull(id)

missing_cvp_50_20_intracpb <- 
  missing_summary %>% 
  filter(cat_cpb == "intra" & missing_cvp_50_20 == 1) %>% 
  pull(id)

missing_cvp_50_20_postcpb <- 
  missing_summary %>% 
  filter(cat_cpb == "post" & missing_cvp_50_20 == 1) %>% 
  pull(id)


# summary of missingness 
missing_summary %>% 
  filter(!is.na(cat_cpb)) %>% 
  group_by(cat_cpb) %>% 
  summarize(mean_elapsed = mean(elapsed),
            sd = sd(elapsed), 
            mean_miss_map = mean(na_map),
            sd_miss_map = sd(na_map),
            mean_miss_cvp = mean(na_cvp),
            sd_miss_cvp = sd(na_cvp),
            mean_pct_map = mean(na_map_pct),
            sd_pct_map = sd(na_map_pct),
            mean_pct_cvp = mean(na_cvp_pct),
            sd_pct_cvp = sd(na_cvp_pct)) %>% 
  kableExtra::kable(align = "llll", digits = 2) %>% 
  kableExtra::kable_styling()


# put all criteria in one dataset 
exclusion_summary <-
  data.frame(
    id = starting_cohort
  ) %>% 
  mutate(
    is_missing_anesthesia = 0,
    has_multiple_anesthesia = 0,
    is_missing_cpb = ifelse(id %in% missing_cpb, 1, 0),
    has_multiple_axc = 0,
    has_ecmo = 0,
    is_missing_map_all = ifelse(id %in% missing_map_all, 1, 0),
    is_missing_cvp_all = ifelse(id %in% missing_cvp_all, 1, 0),
    is_missing_map_precpb = ifelse(id %in% missing_map_precpb, 1, 0),
    is_missing_map_intracpb = ifelse(id %in% missing_map_intracpb, 1, 0),
    is_missing_map_postcpb = ifelse(id %in% missing_map_postcpb, 1, 0),
    is_missing_cvp_precpb = ifelse(id %in% missing_cvp_precpb, 1, 0),
    is_missing_cvp_intracpb = ifelse(id %in% missing_cvp_intracpb, 1, 0),
    is_missing_cvp_postcpb = ifelse(id %in% missing_cvp_postcpb, 1, 0),
    is_missing_map_50_20_precpb = ifelse(id %in% missing_map_50_20_precpb, 1, 0),
    is_missing_map_50_20_intracpb = ifelse(id %in% missing_map_50_20_intracpb, 1, 0),
    is_missing_map_50_20_postcpb = ifelse(id %in% missing_map_50_20_postcpb, 1, 0),
    is_missing_cvp_50_20_precpb = ifelse(id %in% missing_cvp_50_20_precpb, 1, 0),
    is_missing_cvp_50_20_intracpb = ifelse(id %in% missing_cvp_50_20_intracpb, 1, 0),
    is_missing_cvp_50_20_postcpb = ifelse(id %in% missing_cvp_50_20_postcpb, 1, 0),
    has_renal_failure = ifelse(id %in% renal_failure, 1, 0),
    has_dialysis = 0
  ) %>%
  mutate(excluded = case_when(if_any(contains("missing"), ~ . == 1) ~ 1,
                              if_any(contains("has"), ~ . == 1) ~ 1,
                              TRUE ~ 0))


write_csv(exclusion_summary, here("data/processed/exclusion_summary_cohort.csv"))
write_csv(exclusion_summary, here("data/analytic/exclusion_summary_cohort.csv"))

excluded <- exclusion_summary %>% 
  filter(excluded == 1) %>% 
  pull(id)

covariates_filtered <-
  covariates %>% 
  filter(id %notin% excluded) 

write_csv(covariates_filtered, here("data/processed/covariates_post_exclusion.csv"))
write_csv(covariates_filtered, here("data/analytic/covariates_post_exclusion.csv"))


hemo_interp_filtered <-
  all_hemo %>% 
  filter(id %notin% excluded) 

saveRDS(hemo_interp_filtered, 
        here("data/processed/hemo_timeseries_interp_post_exclusion.rds"))
saveRDS(hemo_interp_filtered, 
        here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))


write.csv(hemo_interp_filtered, 
        here("data/processed/hemo_timeseries_interp_post_exclusion.csv"))
