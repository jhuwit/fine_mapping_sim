library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
set.seed(123)
options(dplyr.summarise.inform=F)


# cohort = "full" 
cohort = "cabg"
# read in covariate data 
covar_data = read_csv(here::here("data/analytic/covariates_post_exclusion.csv"))
# read in hemo data 
hemo_data = readRDS(here::here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))

# make categorical variables factors for logistic reg 
covar_data = covar_data %>% 
  mutate(across(starts_with("cat"), as.factor))

# filter to just CABG ids 
if(cohort == "cabg"){
  covar_data = covar_data %>% filter(cat_surgcat4 == "1_CABG")
  hemo_data = hemo_data %>% filter(id %in% covar_data$id)
}

# get covariates for analysis 
primary_and_mediators = colnames(covar_data)[c(6:27, 32, 33, 49)]
# source regression functions 
source(here::here("code/CABG/analysis/utilities.R"))

map_data = readRDS(here::here("results/data_cabg_manuscript/map_data.rds"))
cvp_data = readRDS(here::here("results/data_cabg_manuscript/cvp_data.rds"))
bivar_bricks = readRDS(here::here("results/data_cabg_manuscript/bricks_data.rds"))


set.seed(2123)
# regression for 48h outcome 
# univariate 
map_48h = 
  run_regression_full(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 exp = TRUE,
                 hemo_type = "MAP") 

map_48h %>% 
  rename(regression_type = variable,
         predictor = term) %>% 
  select(regression_type, predictor, estimate, p.value, conf.low, conf.high, AIC) %>% 
  gt::gt() %>% 
  fmt_number(columns = -c(regression_type, predictor), decimals = 3) %>% 
  fmt_number(columns = AIC, decimals = 0) %>% 
  gt::gtsave(here::here("results", "tables_cabg_manuscript", "all_map_covariate_estimates.docx"))

# cvp univariate 48h 
cvp_48h = 
  run_regression_full(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 exp = TRUE,
                 hemo_type = "CVP")

cvp_48h %>% 
  rename(regression_type = variable,
         predictor = term) %>% 
  select(regression_type, predictor, estimate, p.value, conf.low, conf.high, AIC) %>% 
  gt::gt() %>% 
  fmt_number(columns = -c(regression_type, predictor), decimals = 3) %>% 
  fmt_number(columns = AIC, decimals = 0) %>% 
  gt::gtsave(here::here("results", "tables_cabg_manuscript", "all_cvp_covariate_estimates.docx"))


## bricks 
shell_df = readRDS(here::here("results/data_cabg_manuscript/shell_data.rds"))
shell_df_3 = readRDS(here::here("results/data_cabg_manuscript/shell_data_threezones.rds"))


regression_df = 
  shell_df %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>% 
  mutate(across(starts_with("group"), ~.x / 5))


run_multi = function(var, df){
  tempdf = df %>% select(id, all_of(var)) %>% 
    left_join(covar_data %>% select(id, bin_aki48h, all_of(primary_and_mediators)), by = "id") %>% 
    select(-id) %>% 
    mutate(across(bin_aki48h, as.factor))
  m = 
    logistic_reg() %>%
    set_engine("glm") %>%
    fit(bin_aki48h ~ .,
        data = tempdf) 
  
  m %>% 
    tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
    mutate(group = var) %>% 
    mutate(AIC = m$fit$aic)
}


multi_results = map(.x = regression_df %>% select(starts_with("group")) %>% colnames,
                    .f = run_multi,
                    df = regression_df) %>% 
  bind_rows()

regression_df_3 = 
  shell_df_3 %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>% 
  mutate(across(starts_with("zone"), ~.x / 5))


multi_results_3 = map(.x = regression_df_3 %>% select(starts_with("zone")) %>% colnames,
                    .f = run_multi,
                    df = regression_df_3) %>% 
  bind_rows()

multi_results_3 %>% 
  rename(regression_type = group,
       predictor = term) %>% 
  select(regression_type, predictor, estimate, p.value, conf.low, conf.high, AIC) %>% 
  gt::gt() %>% 
  fmt_number(columns = -c(regression_type, predictor), decimals = 3) %>% 
  fmt_number(columns = AIC, decimals = 0) %>% 
  gt::gtsave(here::here("results", "tables_cabg_manuscript", "zones3_covariate_estimates.docx"))

  
  
multi_results = 
  multi_results %>% 
  mutate(group = paste0("zone_", sub(".*group\\_", "", group)))

multi_results %>% 
  rename(regression_type = group,
         predictor = term) %>% 
  select(regression_type, predictor, estimate, p.value, conf.low, conf.high, AIC) %>% 
  gt::gt() %>% 
  fmt_number(columns = -c(regression_type, predictor), decimals = 3) %>% 
  fmt_number(columns = AIC, decimals = 0) %>% 
  gt::gtsave(here::here("results", "tables_cabg_manuscript", "zones5_covariate_estimates.docx"))



