library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
set.seed(123)
options(dplyr.summarise.inform=F)

if(!dir.exists(here::here("results", "data_cabg_manuscript"))){
  dir.create(here::here("results", "data_cabg_manuscript"),
             recursive = TRUE)
}
# cohort = "full" 
cohort = "cabg"

# read in covariate data 
covar_data = read_csv(here::here("data/analytic/covariates_post_exclusion.csv"),
                      col_types = list(id = col_character()))
# read in hemo data 
hemo_data = readRDS(here::here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))

# make categorical variables factors for logistic reg 
covar_data = covar_data %>% 
  mutate(across(starts_with("cat"), as.factor))

# get covariates for analysis 
# primary_and_mediators = colnames(covar_data)[c(6:27, 32, 33, 49)]
primary_and_mediators = colnames(covar_data)[c(6:27)]

# source regression functions 
source(here::here("code/analysis/utilities.R"))

map_thresholds = seq(45, 115, 5)
map_thresholds_bivar = seq(45, 115, 10)
cvp_thresholds = seq(0, 20, 2)

# create time in range regression predictors 
map_data = get_ranges(hemo_data = hemo_data,
                      thresholds = map_thresholds,
                      hemo_variable = "MAP") %>% 
  select(-contains("missing"))
cvp_data = get_ranges(hemo_data = hemo_data,
                      thresholds = cvp_thresholds,
                      hemo_variable = "CVP") %>% 
  select(-contains("missing"))

bivar_bricks = get_ranges_biv(hemo_data, map_thresholds = map_thresholds_bivar,
                              cvp_thresholds = cvp_thresholds)

saveRDS(map_data, here::here("results/data_cabg_manuscript/map_data.rds"))
saveRDS(cvp_data, here::here("results/data_cabg_manuscript/cvp_data.rds"))
saveRDS(bivar_bricks, here::here("results/data_cabg_manuscript/bricks_data.rds"))


set.seed(2123)
# regression for 48h outcome 
# univariate 
map_uni_48h = 
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate",
         outcome = "AKI48")
# adjusted, w cma 
map_adjusted_48h  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted",
         outcome = "AKI48")

# cvp univariate 48h 
cvp_uni_48h = 
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate",
         outcome = "AKI48")

# cvp adjusted 48h with cma 
cvp_adjusted_48h  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, 
                 B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted",
         outcome = "AKI48")



# composite outcome 
map_uni_composite = 
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_akicomposite",
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate",
         outcome = "AKI Composite")

map_adjusted_composite  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_akicomposite",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted",
         outcome = "AKI Composite")

cvp_uni_composite = 
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_akicomposite",
                 adjustment_type = "uni",
                 cma = FALSE, 
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate",
         outcome = "AKI Composite")

cvp_adjusted_composite  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_akicomposite",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted",
         outcome = "AKI Composite")



all_map = bind_rows(map_uni_48h,
                    map_adjusted_48h,
                    map_uni_composite,
                    map_adjusted_composite) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))


all_cvp = bind_rows(cvp_uni_48h,
                    cvp_adjusted_48h,
                    cvp_uni_composite,
                    cvp_adjusted_composite) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

labels = c("Hemodynamic range",
           "Odds ratio",
           "Std. error", 
           "Statistic",
           "P value", 
           "95% CI lower bound",
           "95% CI upper bound",
           "95% bonferroni CI lower bound",
           "95% bonferroni CI upper bound", 
           "FDR adjusted pvalue",
           "Type of model (adjusted or univariate)",
           "Outcome (48h or composite)",
           "95% CMA CI lower bound", 
           "95% CMA CI upper bound"
)

all_map = sjlabelled::set_label(all_map, labels)
all_cvp = sjlabelled::set_label(all_cvp, labels)

saveRDS(all_map, here::here("results/data_cabg_manuscript/map_regressions.rds"))
saveRDS(all_cvp, here::here("results/data_cabg_manuscript/cvp_regressions.rds"))


## bricks 
bricks_uni_48h =
  run_regression_biv(hemo_preds = bivar_bricks,
                     covariate_df = covar_data,
                     outcome = "bin_aki48h",
                     adjustment_type = "uni",
                     cma = FALSE,
                     exponentiate = TRUE) %>%
  mutate(type = "univariate",
         outcome = "AKI 48h")

bricks_uni_composite =
  run_regression_biv(hemo_preds = bivar_bricks,
                     covariate_df = covar_data,
                     outcome = "bin_akicomposite",
                     adjustment_type = "uni",
                     cma = FALSE,
                     exponentiate = TRUE) %>%
  mutate(type = "univariate",
         outcome = "AKI Composite")




bivar_bricks_small = bivar_bricks %>% 
  select(-c("MAP[45,55]CVP(12,14]",
            "MAP[45,55]CVP(14,16]",
            "MAP[45,55]CVP(16,18]",
            "MAP[45,55]CVP(18,20]",
            "MAP(95,105]CVP[0,2]",
            "MAP(95,105]CVP(2,4]",
            "MAP(95,105]CVP(16,18]",
            "MAP(95,105]CVP(18,20]",
            starts_with("MAP(105,115")))

bricks_adjusted_48h_cma = 
  run_regression_biv(hemo_preds = bivar_bricks_small,
                     covariate_df = covar_data,
                     covars_control = primary_and_mediators,
                     outcome = "bin_aki48h",
                     adjustment_type = "multi",
                     cma = TRUE,
                     B = 2,
                     exponentiate = TRUE) %>%
  mutate(type = "adjusted",
         outcome = "AKI48")

# w/o cma to get point estimates 
bricks_adjusted_48h = 
  run_regression_biv(hemo_preds = bivar_bricks,
                     covariate_df = covar_data,
                     outcome = "bin_aki48h",
                     covars_control = primary_and_mediators,
                     adjustment_type = "multi",
                     cma = FALSE,
                     exponentiate = TRUE) %>%
  mutate(type = "adjusted",
         outcome = "AKI48")

bricks_adjusted_composite_cma = 
  run_regression_biv(hemo_preds = bivar_bricks_small,
                     covariate_df = covar_data,
                     covars_control = primary_and_mediators,
                     outcome = "bin_akicomposite",
                     adjustment_type = "multi",
                     cma = TRUE,
                     B = 2,
                     exponentiate = TRUE) %>%
  mutate(type = "adjusted",
         outcome = "AKI Composite")

# w/o cma to get point estimates 
bricks_adjusted_composite = 
  run_regression_biv(hemo_preds = bivar_bricks,
                     covariate_df = covar_data,
                     covars_control = primary_and_mediators,
                     outcome = "bin_akicomposite",
                     adjustment_type = "multi",
                     cma = FALSE,
                     exponentiate = TRUE) %>%
  mutate(type = "adjusted",
         outcome = "AKI Composite")


bricks_greyed_48h = bricks_adjusted_48h %>% 
  filter(term %notin% bricks_adjusted_48h_cma$term)
bricks_greyed_composite = bricks_adjusted_composite %>% 
  filter(term %notin% bricks_adjusted_composite_cma$term)

all_bivariate = bind_rows(bricks_uni_48h,
                          bricks_uni_composite,
                          bricks_greyed_48h,
                          bricks_greyed_composite,
                          bricks_adjusted_48h_cma,
                          bricks_adjusted_composite_cma) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

# map_data %>% left_join(covar_data) %>% 
#   group_by(bin_aki48h) %>%  
#   summarize(across(starts_with("MAP"), ~ sum(.x >5))) %>% 
#   pivot_longer(-bin_aki48h) %>% 
#   arrange(value)

# cvp_data %>% left_join(covar_data) %>% 
#   group_by(bin_aki48h) %>% 
#   summarize(across(starts_with("CVP"), ~ sum(.x >5))) %>% 
#   pivot_longer(-bin_aki48h) %>% 
#   arrange(value)


labels = c("Hemodynamic range",
           "Odds ratio",
           "Std. error", 
           "Statistic",
           "P value", 
           "95% CI lower bound",
           "95% CI upper bound",
           "95% bonferroni CI lower bound",
           "95% bonferroni CI upper bound", 
           "FDR adjusted p-value",
           "Type of model (adjusted or univariate)",
           "Outcome (48h or composite)",
           "95% CMA CI lower bound", 
           "95% CMA CI upper bound"
)

all_bivariate = sjlabelled::set_label(all_bivariate, labels)

saveRDS(all_bivariate, here::here("results/data_cabg_manuscript/bricks_regressions.rds"))

# shell regressions 

range_left = function(x, left, right){
  ifelse(x > left & x <= right, TRUE, FALSE)
}
# group 1: 
shell_df = 
  hemo_data %>% 
  mutate(
    shell = case_when(
      range_left(val_MAP, 95, 115) & between(val_CVP, 0, 8) ~ "group_1",
      (range_left(val_MAP, 75, 95) & between(val_CVP, 0, 8)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 8, 10)) ~ "group_2",
      (range_left(val_MAP, 55, 75) & between(val_CVP, 0, 8)) | 
        (range_left(val_MAP, 65, 85) & range_left(val_CVP, 8, 12)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 10, 12)) ~ "group_3",
      (between(val_MAP, 45, 55) & between(val_CVP, 0, 8)) | 
        (range_left(val_MAP, 55, 65) & range_left(val_CVP, 8, 12)) | 
        (range_left(val_MAP, 65, 115) & range_left(val_CVP, 12, 16)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 16, 18)) ~ "group_4",
      (between(val_MAP, 45, 55) & range_left(val_CVP, 8, 20)) | 
        (range_left(val_MAP, 55, 65) & range_left(val_CVP, 12, 20)) | 
        (range_left(val_MAP, 65, 85) & range_left(val_CVP, 16, 20)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 18, 20)) ~ "group_5"
    )
  ) %>% 
  group_by(id) %>% 
  count(shell, .drop = FALSE) %>% 
  ungroup() %>% 
  drop_na() %>% 
  pivot_wider(names_from =shell, values_from = n, id_cols = id) %>% 
  mutate(across(starts_with("group"), ~ ifelse(is.na(.x), 0, .x)))
saveRDS(shell_df, here::here("results/data_cabg_manuscript/shell_data.rds"))

shell_df_3 = 
  hemo_data %>% 
  mutate(
    shell = case_when(
      range_left(val_MAP, 95, 115) & between(val_CVP, 0, 8) ~ "group_1",
      (range_left(val_MAP, 75, 95) & between(val_CVP, 0, 8)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 8, 10)) ~ "group_2",
      (range_left(val_MAP, 55, 75) & between(val_CVP, 0, 8)) | 
        (range_left(val_MAP, 65, 85) & range_left(val_CVP, 8, 12)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 10, 12)) ~ "group_3",
      (between(val_MAP, 45, 55) & between(val_CVP, 0, 8)) | 
        (range_left(val_MAP, 55, 65) & range_left(val_CVP, 8, 12)) | 
        (range_left(val_MAP, 65, 115) & range_left(val_CVP, 12, 16)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 16, 18)) ~ "group_4",
      (between(val_MAP, 45, 55) & range_left(val_CVP, 8, 20)) | 
        (range_left(val_MAP, 55, 65) & range_left(val_CVP, 12, 20)) | 
        (range_left(val_MAP, 65, 85) & range_left(val_CVP, 16, 20)) | 
        (range_left(val_MAP, 85, 115) & range_left(val_CVP, 18, 20)) ~ "group_5"
    )
  ) %>% 
  mutate(shell_new = 
           case_when(shell %in% c("group_1", "group_2") ~ "zone_1",
                     shell == "group_3" ~ "zone_2",
                     shell %in% c("group_4", "group_5") ~ "zone_3",
                     TRUE ~ NA)) %>% 
  select(-shell) %>% 
  group_by(id) %>% 
  count(shell_new, .drop = FALSE) %>% 
  ungroup() %>% 
  drop_na() %>% 
  pivot_wider(names_from =shell_new, values_from = n, id_cols = id) %>% 
  mutate(across(starts_with("zone"), ~ ifelse(is.na(.x), 0, .x)))
saveRDS(shell_df_3, here::here("results/data_cabg_manuscript/shell_data_threezones.rds"))



regression_df = 
  shell_df %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>% 
  mutate(across(starts_with("group"), ~.x / 5))

run_univariate = function(var, df){
  tempdf = df %>% select(all_of(var), bin_aki48h)
  logistic_reg() %>%
    set_engine("glm") %>%
    fit(bin_aki48h ~ .,
        data = tempdf) %>% 
    tidy(exponentiate = TRUE) %>% 
    filter(grepl("group", term)) 
}

uni_results = map(.x = regression_df %>% select(starts_with("group")) %>% colnames,
                  .f = run_univariate,
                  df = regression_df) %>% 
  bind_rows()



run_multi = function(var, df){
  tempdf = df %>% select(id, all_of(var)) %>% 
    left_join(covar_data %>% select(id, bin_aki48h, all_of(primary_and_mediators))) %>% 
    select(-id) %>% 
    mutate(across(bin_aki48h, as.factor))
  logistic_reg() %>%
    set_engine("glm") %>%
    fit(bin_aki48h ~ .,
        data = tempdf) %>% 
    tidy(exponentiate = TRUE) %>% 
    filter(grepl("group", term)) 
  
}
multi_results = map(.x = regression_df %>% select(starts_with("group")) %>% colnames,
                    .f = run_multi,
                    df = regression_df) %>% 
  bind_rows()

shell_df_covars = 
  shell_df %>% 
  left_join(covar_data %>% select(id, bin_aki48h, all_of(primary_and_mediators))) %>% 
  # select(-val_proctime) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>% 
  mutate(across(starts_with("group"), ~.x / 5)) %>% 
  select(-id) 

shell_df_covars3 = 
  shell_df_3 %>% 
  left_join(covar_data %>% select(id, bin_aki48h, all_of(primary_and_mediators))) %>% 
  # select(-val_proctime) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>% 
  mutate(across(starts_with("zone"), ~.x / 5)) %>% 
  select(-id) 

multi_shell = 
  logistic_reg() %>% 
  set_engine("glm") %>% 
  fit(bin_aki48h ~ ., data = shell_df_covars) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) 

saveRDS(multi_shell, here::here("results", "data_cabg_manuscript", "shell_regressions_5.rds")) 


multi_shell_3 = 
  logistic_reg() %>% 
  set_engine("glm") %>% 
  fit(bin_aki48h ~ ., data = shell_df_covars3) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) 

saveRDS(multi_shell_3, here::here("results", "data_cabg_manuscript", "shell_regressions_3.rds")) 



cor(shell_df_covars3$zone_1, shell_df_covars3$zone_2)
cor(shell_df_covars3$zone_2, shell_df_covars3$zone_3)
cor(shell_df_covars3$zone_1, shell_df_covars3$zone_3)


cor(shell_df_covars$group_1, shell_df_covars$group_2)
cor(shell_df_covars$group_2, shell_df_covars$group_3)
cor(shell_df_covars$group_1, shell_df_covars$group_3)
cor(shell_df_covars$group_1, shell_df_covars$group_4)
cor(shell_df_covars$group_1, shell_df_covars$group_5)
cor(shell_df_covars$group_2, shell_df_covars$group_4)
cor(shell_df_covars$group_2, shell_df_covars$group_5)
cor(shell_df_covars$group_3, shell_df_covars$group_4)
cor(shell_df_covars$group_3, shell_df_covars$group_5)


all_shells = bind_rows(uni_results %>% mutate(type = "univariate"), 
                       multi_results %>% mutate(type = "adjusted"))

saveRDS(all_shells, here::here("results/data_cabg_manuscript/shell_regressions.rds"))


# sensitivity analysis with bigger ranges 
map_thresholds_10 = seq(45, 115, 10)
map_thresholds_15 = seq(45, 120, 15)
map_thresholds_20 = seq(45, 125, 20)


cvp_thresholds_4 = seq(0, 20, 4)
cvp_thresholds_6 = seq(0, 20, 6)
cvp_thresholds_8 = seq(0, 20, 8)

# create time in range regression predictors 
map_data_10 = get_ranges(hemo_data = hemo_data,
                         thresholds = map_thresholds_10,
                         hemo_variable = "MAP") %>% 
  select(-contains("missing"))

map_data_15 = get_ranges(hemo_data = hemo_data,
                         thresholds = map_thresholds_15,
                         hemo_variable = "MAP") %>% 
  select(-contains("missing"))

map_data_20 = get_ranges(hemo_data = hemo_data,
                         thresholds = map_thresholds_20,
                         hemo_variable = "MAP") %>% 
  select(-contains("missing"))

cvp_data_4 = get_ranges(hemo_data = hemo_data,
                        thresholds = cvp_thresholds_4,
                        hemo_variable = "CVP") %>% 
  select(-contains("missing"))

cvp_data_6 = get_ranges(hemo_data = hemo_data,
                        thresholds = cvp_thresholds_6,
                        hemo_variable = "CVP") %>% 
  select(-contains("missing"))

cvp_data_8 = get_ranges(hemo_data = hemo_data,
                        thresholds = cvp_thresholds_8,
                        hemo_variable = "CVP") %>% 
  select(-contains("missing"))

set.seed(2123)
controls = primary_and_mediators[-c(20,22, 24)]

# regression for 48h outcome 
# 
map_adjusted_48h_10  =
  run_regression(hemo_preds = map_data_10,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_10",
         outcome = "AKI48")

map_nomed_48h_10  =
  run_regression(hemo_preds = map_data_10,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_10_nomed",
         outcome = "AKI48")


map_uni_48h_10 = 
  run_regression(hemo_preds = map_data_10,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_10",
         outcome = "AKI48")

map_adjusted_48h_15  =
  run_regression(hemo_preds = map_data_15,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_15",
         outcome = "AKI48")
map_uni_48h_15 = 
  run_regression(hemo_preds = map_data_15,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_15",
         outcome = "AKI48")
map_nomed_48h_15  =
  run_regression(hemo_preds = map_data_15,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_15_nomed",
         outcome = "AKI48")
map_adjusted_48h_20  =
  run_regression(hemo_preds = map_data_20,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_20",
         outcome = "AKI48")
map_uni_48h_20 = 
  run_regression(hemo_preds = map_data_20,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_20",
         outcome = "AKI48")

map_nomed_48h_20  =
  run_regression(hemo_preds = map_data_20,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 2,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_20_nomed",
         outcome = "AKI48")

# cvp adjusted 48h with cma 
cvp_adjusted_48h_4  =
  run_regression(hemo_preds = cvp_data_4,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, 
                 B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_4",
         outcome = "AKI48")
cvp_uni_48h_4 = 
  run_regression(hemo_preds = cvp_data_4,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_4",
         outcome = "AKI48")
cvp_nomed_48h_4  =
  run_regression(hemo_preds = cvp_data_4,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_4_nomed",
         outcome = "AKI48")
cvp_adjusted_48h_6  =
  run_regression(hemo_preds = cvp_data_6,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, 
                 B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_6",
         outcome = "AKI48")
cvp_uni_48h_6 = 
  run_regression(hemo_preds = cvp_data_6,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_6",
         outcome = "AKI48")

cvp_nomed_48h_6  =
  run_regression(hemo_preds = cvp_data_6,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_6_nomed",
         outcome = "AKI48")

cvp_adjusted_48h_8  =
  run_regression(hemo_preds = cvp_data_8,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = TRUE, 
                 B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_8",
         outcome = "AKI48")

cvp_uni_48h_8 = 
  run_regression(hemo_preds = cvp_data_8,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_8",
         outcome = "AKI48")

all_map = bind_rows(map_adjusted_48h_10,
                    map_adjusted_48h_15,
                    map_adjusted_48h_20,
                    map_uni_48h_10,
                    map_uni_48h_15,
                    map_uni_48h_20) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

cvp_nomed_48h_8  =
  run_regression(hemo_preds = cvp_data_8,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 B = 2,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_8_nomed",
         outcome = "AKI48")

all_cvp = bind_rows(cvp_adjusted_48h_4,
                    cvp_adjusted_48h_6,
                    cvp_adjusted_48h_8,
                    cvp_uni_48h_4,
                    cvp_uni_48h_6,
                    cvp_uni_48h_8,
                    cvp_nomed_48h_4,
                    cvp_nomed_48h_6,
                    cvp_nomed_48h_8) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

labels = c("Hemodynamic range",
           "Odds ratio",
           "Std. error", 
           "Statistic",
           "P value", 
           "95% CI lower bound",
           "95% CI upper bound",
           "95% bonferroni CI lower bound",
           "95% bonferroni CI upper bound", 
           "FDR adjusted pvalue",
           "Type of model (adjusted or univariate)",
           "Outcome (48h or composite)",
           "95% CMA CI lower bound", 
           "95% CMA CI upper bound"
)

all_map = sjlabelled::set_label(all_map, labels)
all_cvp = sjlabelled::set_label(all_cvp, labels)

saveRDS(all_map, here::here("results/data_cabg_manuscript/map_regressions_telescoping.rds"))
saveRDS(all_cvp, here::here("results/data_cabg_manuscript/cvp_regressions_telescoping.rds"))



