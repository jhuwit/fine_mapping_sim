# post-hoc sensivity analyses 
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
# bivar_bricks = readRDS(here::here("results/data_cabg_manuscript/bricks_data.rds"))


set.seed(2123)
# analysis 1: by phase of surgery 

map_thresholds = seq(45, 115, 5)
map_thresholds_bivar = seq(45, 115, 10)
cvp_thresholds = seq(0, 20, 2)

# create time in range regression predictors 
map_data_pre = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb == "pre"),
                          thresholds = map_thresholds,
                          hemo_variable = "MAP") %>% 
  select(-contains("missing"))
map_data_intra = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb == "intra"),
                            thresholds = map_thresholds,
                            hemo_variable = "MAP") %>% 
  select(-contains("missing"))
map_data_post = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb == "post"),
                           thresholds = map_thresholds,
                           hemo_variable = "MAP") %>% 
  select(-contains("missing"))

cvp_data_pre = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb == "pre"),
                          thresholds = cvp_thresholds,
                          hemo_variable = "CVP") %>% 
  select(-contains("missing"))
cvp_data_intra = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb == "intra"),
                            thresholds = cvp_thresholds,
                            hemo_variable = "CVP") %>% 
  select(-contains("missing"))
cvp_data_post = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb == "post"),
                           thresholds = cvp_thresholds,
                           hemo_variable = "CVP") %>% 
  select(-contains("missing"))


shell5_data = readRDS(here::here("results/data_cabg_manuscript/shell_data.rds"))
shell3_data = readRDS(here::here("results/data_cabg_manuscript/shell_data_threezones.rds"))

range_left = function(x, left, right){
  ifelse(x > left & x <= right, TRUE, FALSE)
}

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
  group_by(id, cat_cpb) %>% 
  count(shell, .drop = FALSE) %>% 
  ungroup() %>% 
  drop_na() %>% 
  pivot_wider(names_from =shell, values_from = n, id_cols = c(id, cat_cpb)) %>% 
  mutate(across(starts_with("group"), ~ ifelse(is.na(.x), 0, .x))) %>% 
  group_split(cat_cpb)
shell5_pre = shell_df[[3]]  %>% dplyr::select(-cat_cpb)
shell5_intra = shell_df[[1]]  %>% dplyr::select(-cat_cpb)
shell5_post = shell_df[[2]] %>% dplyr::select(-cat_cpb)

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
  group_by(id, cat_cpb) %>% 
  count(shell_new, .drop = FALSE) %>% 
  ungroup() %>% 
  drop_na() %>% 
  pivot_wider(names_from =shell_new, values_from = n, id_cols = c(cat_cpb, id)) %>% 
  mutate(across(starts_with("zone"), ~ ifelse(is.na(.x), 0, .x))) %>% 
  group_split(cat_cpb) 
shell3_pre = shell_df_3[[3]] %>% dplyr::select(-cat_cpb)
shell3_intra = shell_df_3[[1]] %>% dplyr::select(-cat_cpb)
shell3_post = shell_df_3[[2]] %>% dplyr::select(-cat_cpb)



################ map - pre post intra 
map_pre_uni = 
  run_regression(hemo_preds = map_data_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_precpb")
map_pre_adj  =
  run_regression(hemo_preds = map_data_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_precpb")


map_intra_uni = 
  run_regression(hemo_preds = map_data_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_intracpb")

map_intra_adj  =
  run_regression(hemo_preds = map_data_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_intracpb")

map_post_uni = 
  run_regression(hemo_preds = map_data_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_postcpb")
map_post_adj  =
  run_regression(hemo_preds = map_data_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_postcpb")
################ cvp - pre post intra 

cvp_pre_uni = 
  run_regression(hemo_preds = cvp_data_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_precpb")
cvp_pre_adj  =
  run_regression(hemo_preds = cvp_data_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_precpb")
cvp_intra_uni = 
  run_regression(hemo_preds = cvp_data_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_intracpb")
cvp_intra_adj  =
  run_regression(hemo_preds = cvp_data_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_intracpb")

cvp_post_uni = 
  run_regression(hemo_preds = cvp_data_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_postcpb")
cvp_post_adj  =
  run_regression(hemo_preds = cvp_data_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_postcpb")

#### shells 3 
shell3_pre_uni = 
  run_regression(hemo_preds = shell3_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_zone3_precpb")
shell3_pre_adj  =
  run_regression(hemo_preds = shell3_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_zone3_precpb")

shell3_intra_uni = 
  run_regression(hemo_preds = shell3_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_zone3_intracpb")

shell3_intra_adj  =
  run_regression(hemo_preds = shell3_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_zone3_intracpb")

shell3_post_uni = 
  run_regression(hemo_preds = shell3_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_zone3_postcpb")
shell3_post_adj  =
  run_regression(hemo_preds = shell3_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_zone3_postcpb")


shell3_uni = 
  run_regression(hemo_preds = shell3_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_zone3")

shell3_adj = 
  run_regression(hemo_preds = shell3_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_zone3")

# shells 5

shell5_pre_uni = 
  run_regression(hemo_preds = shell5_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_zone5_precpb")
shell5_pre_adj  =
  run_regression(hemo_preds = shell5_pre,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_zone5_precpb")

shell5_intra_uni = 
  run_regression(hemo_preds = shell5_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_zone5_intracpb")

shell5_intra_adj  =
  run_regression(hemo_preds = shell5_intra,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_zone5_intracpb")

shell5_post_uni = 
  run_regression(hemo_preds = shell5_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_zone5_postcpb")
shell5_post_adj  =
  run_regression(hemo_preds = shell5_post,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_zone5_postcpb")


shell5_uni = 
  run_regression(hemo_preds = shell5_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_zone5")

shell5_adj = 
  run_regression(hemo_preds = shell5_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_zone5")

all_map = 
  bind_rows(map_pre_uni,
            map_pre_adj,
            map_intra_uni, 
            map_intra_adj,
            map_post_uni,
            map_post_adj) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term),
         analysis = "MAP")
unique(all_map$type)

all_cvp = 
  bind_rows(cvp_pre_uni,
            cvp_pre_adj,
            cvp_intra_uni, 
            cvp_intra_adj,
            cvp_post_uni,
            cvp_post_adj) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term)) %>% 
  mutate(analysis = "CVP")
unique(all_cvp$type)

all_shells3 = 
  bind_rows(shell3_pre_uni,
            shell3_pre_adj,
            shell3_intra_uni,
            shell3_intra_adj,
            shell3_post_uni,
            shell3_post_adj,
            shell3_adj,
            shell3_uni) %>% 
  mutate(analysis = "shells3")
unique(all_shells3$type)

all_shells5 = 
  bind_rows(shell5_pre_uni,
            shell5_pre_adj,
            shell5_intra_uni,
            shell5_intra_adj,
            shell5_post_uni,
            shell5_post_adj,
            shell5_uni,
            shell5_adj) %>% 
  mutate(analysis = "shells5")
unique(all_shells5$type)

all_res = 
  bind_rows(all_map,
            all_cvp,
            all_shells3,
            all_shells5)
saveRDS(all_res, here::here("results/data_cabg_manuscript/post_hoc_cpb_regressions.rds"))


################# pts with lvef <= 40% 
ids = covar_data %>% 
  filter(val_hdef <= 40) %>% 
  pull(id)

# 209 pts 

map_uni_48h1 = 
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_lvefleq40")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_lvefleq40")

map_uni_48h2 = 
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_lvefg40")

map_adjusted_48h2  =
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_lvefg40")


cvp_uni_48h1 = 
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_lvefleq40")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_lvefleq40")

cvp_uni_48h2 = 
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_lvefg40")

cvp_adjusted_48h2  =
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>%  filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_lvefg40")


shell3_uni_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_lvefleq40")

shell3_adjusted_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_lvefleq40")

shell3_uni_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_lvefg40")

shell3_adjusted_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_lvefg40")


shell5_uni_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_lvefleq40")

shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_lvefleq40")

shell5_uni_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_lvefg40")

shell5_adjusted_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_lvefg40")


res = 
  bind_rows(map_uni_48h1,
            map_adjusted_48h1,
            map_uni_48h2,
            map_adjusted_48h2,
            cvp_uni_48h1,
            cvp_adjusted_48h1,
            cvp_uni_48h2,
            cvp_adjusted_48h2,
            shell3_uni_48h1,
            shell3_adjusted_48h1,
            shell3_uni_48h2,
            shell3_adjusted_48h2,
            shell5_uni_48h1,
            shell5_adjusted_48h1,
            shell5_uni_48h2,
            shell5_adjusted_48h2) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_lvef40_regressions.rds"))

############# pts with egfr < 60
ids = covar_data %>% 
  filter(val_egfr < 60) %>% 
  pull(id)

# 372 pts ; 387 if instead do creatinine >= 1.2
map_uni_48h1 = 
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_egfrgeq60")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_egfrgeq60")

map_uni_48h2 = 
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_egfrl60")

map_adjusted_48h2  =
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_egfrl60")


cvp_uni_48h1 = 
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_egfrgeq60")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_egfrgeq60")

cvp_uni_48h2 = 
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_egfrl60")

cvp_adjusted_48h2  =
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>%  filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_egfrl60")


shell3_uni_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_egfrgeq60")

shell3_adjusted_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_egfrgeq60")

shell3_uni_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_egfrl60")

shell3_adjusted_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_egfrl60")


shell5_uni_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_egfrgeq60")

shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_egfrgeq60")

shell5_uni_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_egfrl60")

shell5_adjusted_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_egfrl60")


res = 
  bind_rows(map_uni_48h1,
            map_adjusted_48h1,
            map_uni_48h2,
            map_adjusted_48h2,
            cvp_uni_48h1,
            cvp_adjusted_48h1,
            cvp_uni_48h2,
            cvp_adjusted_48h2,
            shell3_uni_48h1,
            shell3_adjusted_48h1,
            shell3_uni_48h2,
            shell3_adjusted_48h2,
            shell5_uni_48h1,
            shell5_adjusted_48h1,
            shell5_uni_48h2,
            shell5_adjusted_48h2) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_egfr60_regressions.rds"))

### shock 

ids = covar_data %>% 
  filter(bin_iabp == 1 | bin_emergent == 1) %>% 
  pull(id)

# 1093 pts 
map_uni_48h1 = 
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_shockyes")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_shockyes")

map_uni_48h2 = 
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_shockno")

map_adjusted_48h2  =
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_shockno")


cvp_uni_48h1 = 
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_shockyes")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_shockyes")

cvp_uni_48h2 = 
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_shockno")

cvp_adjusted_48h2  =
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>%  filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_shockno")


shell3_uni_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_shockyes")

shell3_adjusted_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_shockyes")

shell3_uni_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_shockno")

shell3_adjusted_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_shockno")


shell5_uni_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_shockyes")

shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_shockyes")

shell5_uni_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_shockno")

shell5_adjusted_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_shockno")


res = 
  bind_rows(map_uni_48h1,
            map_adjusted_48h1,
            map_uni_48h2,
            map_adjusted_48h2,
            cvp_uni_48h1,
            cvp_adjusted_48h1,
            cvp_uni_48h2,
            cvp_adjusted_48h2,
            shell3_uni_48h1,
            shell3_adjusted_48h1,
            shell3_uni_48h2,
            shell3_adjusted_48h2,
            shell5_uni_48h1,
            shell5_adjusted_48h1,
            shell5_uni_48h2,
            shell5_adjusted_48h2) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_shock_regressions.rds"))


#### use predrenf instead of predmort 
controls = c(primary_and_mediators[-23], "val_predrenf")


map_adjusted_48h1  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_predrenf")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_predrenf")


shell3_adjusted_48h1 = 
  run_regression(hemo_preds = shell3_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_predrenf")


shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_predrenf")


res = 
  bind_rows(map_adjusted_48h1,
            cvp_adjusted_48h1,
            shell3_adjusted_48h1,
            shell5_adjusted_48h1) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_predrenf_regressions.rds"))


# no mediators  
controls = primary_and_mediators[-c(20,22, 24)]


map_adjusted_48h1  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_nomed")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_nomed")


shell3_adjusted_48h1 = 
  run_regression(hemo_preds = shell3_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_nomed")


shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_nomed")


res = 
  bind_rows(map_adjusted_48h1,
            cvp_adjusted_48h1,
            shell3_adjusted_48h1,
            shell5_adjusted_48h1) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_nomediator_regressions.rds"))

### age 
############# pts with age < 65
ids = covar_data %>% 
  filter(val_age < 65) %>% 
  pull(id)

tertiles = quantile(covar_data$val_age, c(1/3, 2/3)); tertiles

t1 = covar_data %>% 
  filter(val_age < 61) %>% 
  pull(id)

t2 = covar_data %>% 
  filter(val_age >= 61 & val_age < 69) %>% 
  pull(id)


t3 = covar_data %>% 
  filter(val_age >= 69) %>% 
  pull(id)

map_uni_48h1 = 
  run_regression(hemo_preds = map_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_aget1")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_aget1")

map_uni_48h2 = 
  run_regression(hemo_preds = map_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_aget2")

map_adjusted_48h2  =
  run_regression(hemo_preds = map_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_aget2")

map_uni_48h3 = 
  run_regression(hemo_preds = map_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_aget3")

map_adjusted_48h3  =
  run_regression(hemo_preds = map_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_aget3")
cvp_uni_48h1 = 
  run_regression(hemo_preds = cvp_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_aget1")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_aget1")

cvp_uni_48h2 = 
  run_regression(hemo_preds = cvp_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_aget2")

cvp_adjusted_48h2  =
  run_regression(hemo_preds = cvp_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_aget2")

cvp_uni_48h3 = 
  run_regression(hemo_preds = cvp_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_aget3")

cvp_adjusted_48h3  =
  run_regression(hemo_preds = cvp_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_aget3")

shell3_uni_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_aget1")

shell3_adjusted_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_aget1")

shell3_uni_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_aget2")

shell3_adjusted_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_aget2")
shell3_uni_48h3 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_aget3")

shell3_adjusted_48h3 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_aget3")

shell5_uni_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_aget1")

shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% t1),
                 covariate_df = covar_data %>% filter(id %in% t1),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_aget1")

shell5_uni_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_aget2")

shell5_adjusted_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% t2),
                 covariate_df = covar_data %>% filter(id %in% t2),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_aget2")

shell5_uni_48h3 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_aget3")

shell5_adjusted_48h3 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% t3),
                 covariate_df = covar_data %>% filter(id %in% t3),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_aget3")

res = 
  bind_rows(map_uni_48h1,
            map_adjusted_48h1,
            map_uni_48h2,
            map_adjusted_48h2,
            map_uni_48h3,
            map_adjusted_48h3,
            cvp_uni_48h1,
            cvp_adjusted_48h1,
            cvp_uni_48h2,
            cvp_adjusted_48h2,
            cvp_uni_48h3,
            cvp_adjusted_48h3,
            shell3_uni_48h1,
            shell3_adjusted_48h1,
            shell3_uni_48h2,
            shell3_adjusted_48h2,
            shell3_uni_48h3,
            shell3_adjusted_48h3,
            shell5_uni_48h1,
            shell5_adjusted_48h1,
            shell5_uni_48h2,
            shell5_adjusted_48h2,
            shell5_uni_48h3,
            shell5_adjusted_48h3) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_agetertiles_regressions.rds"))

### end 
ids = covar_data %>% 
  filter(val_age < 65) %>% 
  pull(id)

map_uni_48h1 = 
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_agel65")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_agel65")

map_uni_48h2 = 
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "univariate_MAP_agegeq65")

map_adjusted_48h2  =
  run_regression(hemo_preds = map_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_agegeq65")


cvp_uni_48h1 = 
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_agel65")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_agel65")

cvp_uni_48h2 = 
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "univariate_CVP_agegeq65")

cvp_adjusted_48h2  =
  run_regression(hemo_preds = cvp_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>%  filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_agegeq65")


shell3_uni_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_agel65")

shell3_adjusted_48h1 = 
  run_regression(hemo_preds = shell3_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_agel65")

shell3_uni_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "univariate_shell3_agegeq65")

shell3_adjusted_48h2 = 
  run_regression(hemo_preds = shell3_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_shell3_agegeq65")


shell5_uni_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_agel65")

shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data %>% filter(id %in% ids),
                 covariate_df = covar_data %>% filter(id %in% ids),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_agel65")

shell5_uni_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "uni",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "univariate_shell5_agegeq65")

shell5_adjusted_48h2 = 
  run_regression(hemo_preds = shell5_data %>% filter(!(id %in% ids)),
                 covariate_df = covar_data %>% filter(!(id %in% ids)),
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_agegeq65")


res = 
  bind_rows(map_uni_48h1,
            map_adjusted_48h1,
            map_uni_48h2,
            map_adjusted_48h2,
            cvp_uni_48h1,
            cvp_adjusted_48h1,
            cvp_uni_48h2,
            cvp_adjusted_48h2,
            shell3_uni_48h1,
            shell3_adjusted_48h1,
            shell3_uni_48h2,
            shell3_adjusted_48h2,
            shell5_uni_48h1,
            shell5_adjusted_48h1,
            shell5_uni_48h2,
            shell5_adjusted_48h2) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_age65_regressions.rds"))

### categorize egfr and lvef here 

covar_data_cat = 
  covar_data %>% 
  mutate(cat_egfr = val_egfr < 60,
         cat_lvef = val_hdef <= 40)

controls = c(primary_and_mediators[-which(primary_and_mediators =="val_creatlst")], "cat_egfr")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_categfr")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_categfr")



shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_categfr")


res = 
  bind_rows(map_adjusted_48h1,
            cvp_adjusted_48h1,
            shell5_adjusted_48h1) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_categfr_regressions.rds"))

controls = c(primary_and_mediators[-which(primary_and_mediators =="val_creatlst")], "val_egfr")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_contegfr")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_contegfr")



shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_contegfr")


res = 
  bind_rows(map_adjusted_48h1,
            cvp_adjusted_48h1,
            shell5_adjusted_48h1) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_contegfr_regressions.rds"))


controls = c(primary_and_mediators[-which(primary_and_mediators =="val_hdef")], "cat_lvef")

map_adjusted_48h1  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_MAP_catlvef")

cvp_adjusted_48h1  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_CVP_catlvef")



shell5_adjusted_48h1 = 
  run_regression(hemo_preds = shell5_data,
                 covariate_df = covar_data_cat,
                 outcome = "bin_aki48h",
                 covars_control = controls,
                 adjustment_type = "multi",
                 cma = FALSE,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_shell5_catlvef")


res = 
  bind_rows(map_adjusted_48h1,
            cvp_adjusted_48h1,
            shell5_adjusted_48h1) %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))
saveRDS(res, here::here("results/data_cabg_manuscript/post_hoc_catlvef_regressions.rds"))