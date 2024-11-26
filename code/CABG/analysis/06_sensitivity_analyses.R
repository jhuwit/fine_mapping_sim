library(here)
library(tidyverse)
library(gt)
options(dplyr.summarise.inform=F)
library(paletteer)
library(ggrepel)
# source(here::here("code/analysis/utilities.R"))
scales::show_col(hue_pal()(3))

col1 = "#7873C0FF"; col2 = "#21B087FF"; col3 = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"

# col1 = "#F8766D"; col2 = "#00BA38"; col3 = "#619CFF"
col_vector = c(col1, col2, col3, col6, col5, col4)
# col_vector = c(col1, col2, col3)

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


map_data = readRDS(here::here("results/data_cabg_manuscript/map_data.rds"))
cvp_data = readRDS(here::here("results/data_cabg_manuscript/cvp_data.rds"))
shell_df = readRDS(here::here("results/data_cabg_manuscript/shell_data.rds"))
shell_df_3 = readRDS(here::here("results/data_cabg_manuscript/shell_data_threezones.rds"))


covar_data = 
  covar_data %>% 
  mutate(bin_ef_leq40 = val_hdef <= 40,
         bin_egfr_l60 = val_egfr < 60)

primary_and_mediators1 = colnames(covar_data)[c(6:18, 20:27, 32, 33, 49, 50)]

primary_and_mediators2 = colnames(covar_data)[c(6:18, 20:26, 32, 33, 49, 50, 51)]

set.seed(2123)
# regression for 48h outcome 
# univariate 

# adjusted, w cma 
map_adjusted_48h  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators1,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_catlvef",
         outcome = "AKI48") %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

map_adjusted_48h2  =
  run_regression(hemo_preds = map_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators2,
                 adjustment_type = "multi",
                 cma = FALSE, B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "MAP") %>%
  mutate(type = "adjusted_catlvefegfr",
         outcome = "AKI48") %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))


# cvp adjusted 48h with cma 
cvp_adjusted_48h  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators1,
                 adjustment_type = "multi",
                 cma = FALSE, 
                 B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_catlvef",
         outcome = "AKI48") %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

cvp_adjusted_48h2  =
  run_regression(hemo_preds = cvp_data,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators2,
                 adjustment_type = "multi",
                 cma = FALSE, 
                 B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "CVP") %>%
  mutate(type = "adjusted_catlvefegfr",
         outcome = "AKI48") %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))


shell5_adjusted_48h  =
  run_regression(hemo_preds = shell_df,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators1,
                 adjustment_type = "multi",
                 cma = FALSE, 
                 B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_catlvef",
         outcome = "AKI48") 
shell5_adjusted_48h2  =
  run_regression(hemo_preds = shell_df,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators2,
                 adjustment_type = "multi",
                 cma = FALSE, 
                 B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "group") %>%
  mutate(type = "adjusted_catlvefegfr",
         outcome = "AKI48") %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))


shell3_adjusted_48h  =
  run_regression(hemo_preds = shell_df_3,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators1,
                 adjustment_type = "multi",
                 cma = FALSE, 
                 B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_catlvef",
         outcome = "AKI48") 
shell3_adjusted_48h2  =
  run_regression(hemo_preds = shell_df_3,
                 covariate_df = covar_data,
                 outcome = "bin_aki48h",
                 covars_control = primary_and_mediators2,
                 adjustment_type = "multi",
                 cma = FALSE, 
                 B = 1000,
                 exponentiate = TRUE,
                 hemo_type = "zone") %>%
  mutate(type = "adjusted_catlvefegfr",
         outcome = "AKI48") %>% 
  mutate(term = sub(".*\\`(.+)\\`.*", "\\1", term))

all = 
  bind_rows(map_adjusted_48h, 
            map_adjusted_48h2, 
            cvp_adjusted_48h,
            cvp_adjusted_48h2,
            shell5_adjusted_48h,
            shell5_adjusted_48h2,
            shell3_adjusted_48h,
            shell3_adjusted_48h2)

plot_j = function(df, new_levels, n_comparisons, xlabs,
                  xname, yname, title, lims){
  colors = col_vector[1:n_comparisons]
  col_vec = shQuote(colors, type = "cmd")
  df %>% 
    mutate(term = factor(term, levels = new_levels)) %>% 
    ggplot(aes(x = term, y = estimate))+
    geom_point(aes(x = term, y = estimate, col = type), position = position_dodge(0.5)) +
    geom_errorbar(width = .35, aes(
      x = term,
      y = estimate,
      ymin = lb,
      ymax = ub,
      col = type
    ),
    position = position_dodge(0.5)) +
    theme_classic() +
    scale_x_discrete(labels = xlabs)+
    scale_y_continuous(limits = lims)+ 
    labs(x = paste(xname), y = paste(yname),
         title = paste(title)) +
    geom_hline(aes(yintercept = 1), col = "black") +
    scale_color_manual(
      values = col_vector[1:n_comparisons],
      name = "95% CI"
    ) +
    # theme(legend.position = pos)
    theme(legend.position = "bottom")
}


temp = all %>% filter(grepl("MAP", term) & !grepl("CVP", term)) 
map_levels = unique(temp$term)
map_labels = sub(".*\\_", "", map_levels)

temp = all %>% filter(grepl("CVP", term) & !grepl("MAP", term)) 
cvp_levels = unique(temp$term)
cvp_labels = sub(".*\\_", "", cvp_levels)


overall_map = readRDS(here::here("results/data_cabg_manuscript/map_regressions.rds"))
overall_cvp = readRDS(here::here("results/data_cabg_manuscript/cvp_regressions.rds"))

pdf(here::here("results/figures_cabg_manuscript/sensitivity_analyses_categorical_covars.pdf"))


temp_df = 
  all %>% 
  filter(grepl("MAP", term)) %>% 
  mutate(type = if_else(type == "adjusted_catlvef", "Cat. LVEF", "Cat. LVEF & Cat EGFR")) %>% 
  bind_rows(overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "Original"))
plot_j(df = temp_df,
       new_levels = map_levels,
       n_comparisons = 3,
       xlabs = map_labels,
       xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
       title = "MAP - Various Adjustments",
       lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  all %>% 
  filter(grepl("CVP", term)) %>% 
  mutate(type = if_else(type == "adjusted_catlvef", "Cat. LVEF", "Cat. LVEF & Cat EGFR")) %>% 
  bind_rows(overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "adjusted_original"))

plot_j(df = temp_df,
       new_levels = cvp_levels,
       n_comparisons = 3,
       xlabs = cvp_labels,
       xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
       title = "CVP - Various Adjustments",
       lims = c(min(temp_df$lb), max(temp_df$ub)))

range_left = function(x, left, right){
  ifelse(x > left & x <= right, TRUE, FALSE)
}
hemo_data = readRDS(here::here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))

plotdf_5shell = hemo_data %>% 
  mutate(
    group = case_when(
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
    ),
    cut_MAP = cut(
      val_MAP,
      breaks = seq(45, 115, 10),
      include.lowest = T
    ),
    cut_CVP = cut(
      val_CVP,
      breaks = seq(0, 20, 2),
      include.lowest = T
    )
  ) %>%
  dplyr::count(cut_MAP, cut_CVP, group, .drop = F) %>% 
  drop_na() %>% group_by(cut_MAP, cut_CVP, group) %>%
  summarize(mean = mean(n)) 

plotdf_3shell = hemo_data %>% 
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
    ),
    shell_new = 
      case_when(shell %in% c("group_1", "group_2") ~ "zone_1",
                shell == "group_3" ~ "zone_2",
                shell %in% c("group_4", "group_5") ~ "zone_3",
                TRUE ~ NA),
    cut_MAP = cut(
      val_MAP,
      breaks = seq(45, 115, 10),
      include.lowest = T
    ),
    cut_CVP = cut(
      val_CVP,
      breaks = seq(0, 20, 2),
      include.lowest = T
    )
  ) %>%
  dplyr::count(cut_MAP, cut_CVP, shell_new, .drop = F) %>% 
  drop_na() %>% group_by(cut_MAP, cut_CVP, shell_new) %>%
  summarize(mean = mean(n)) 

res = readRDS(here::here("results/data_cabg_manuscript/post_hoc_cpb_regressions.rds"))

temp_df = 
  all %>% 
  filter(grepl("group", term)) %>% 
  mutate(type = if_else(type == "adjusted_catlvef", "Cat. LVEF", "Cat. LVEF & Cat EGFR")) %>% 
  bind_rows(res %>% filter(type == "adjusted_zone5") %>% 
              mutate(type = "Original"))

plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = estimate
  )) + facet_wrap(.~type, nrow = 2)+
  geom_tile(col = "black") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio",
    title = "5 Zones, Various Adjustments"
  ) +
  scale_fill_paletteer_c("ggthemes::Red-Blue Diverging", direction = -1)+
  theme(panel.grid = element_blank())


temp_df = 
  all %>% 
  filter(grepl("zone", term)) %>% 
  mutate(type = if_else(type == "adjusted_catlvef", "Cat. LVEF", "Cat. LVEF & Cat EGFR")) %>% 
  bind_rows(res %>% filter(type == "adjusted_zone3") %>% 
              mutate(type = "Original"))

plotdf_3shell %>%
  left_join(temp_df, by = c("shell_new" = "term")) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = estimate
  )) + facet_wrap(.~type, nrow = 2)+
  geom_tile(col = "black") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio",
    title = "3 Zones, Various Adjustments"
  ) +
  scale_fill_paletteer_c("ggthemes::Red-Blue Diverging", direction = -1)+
  theme(panel.grid = element_blank())


dev.off()