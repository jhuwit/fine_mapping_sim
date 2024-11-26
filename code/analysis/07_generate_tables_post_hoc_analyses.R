# process results 
library(tidyverse)
library(gt)
rm(list = ls())
options(dplyr.summarise.inform=F)
library(paletteer)
library(ggrepel)
# source(here::here("code/analysis/utilities.R"))
# scales::show_col(hue_pal()(3))

col1 = "#7873C0FF"; col2 = "#21B087FF"; col3 = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"

# col1 = "#F8766D"; col2 = "#00BA38"; col3 = "#619CFF"
col_vector = c(col1, col2, col3, col6, col5, col4)
# col_vector = c(col1, col2, col3)

shell_res = readRDS(here::here("results/data_cabg_manuscript/shell_regressions.rds"))

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


# step one: j plots separately for pre, post, intra bypass (with overall) 
res = readRDS(here::here("results/data_cabg_manuscript/post_hoc_cpb_regressions.rds"))

temp = res %>% filter(grepl("MAP", term) & !grepl("CVP", term)) 
map_levels = unique(temp$term)
map_labels = sub(".*\\_", "", map_levels)

temp = res %>% filter(grepl("CVP", term) & !grepl("MAP", term)) 
cvp_levels = unique(temp$term)
cvp_labels = sub(".*\\_", "", cvp_levels)


overall_map = readRDS(here::here("results/data_cabg_manuscript/map_regressions.rds"))
overall_cvp = readRDS(here::here("results/data_cabg_manuscript/cvp_regressions.rds"))



### now do analysis for LVEF 

lvef = readRDS(here::here("results/data_cabg_manuscript/post_hoc_lvef40_regressions.rds")) 
egfr = readRDS(here::here("results/data_cabg_manuscript/post_hoc_egfr60_regressions.rds")) 
shock = readRDS(here::here("results/data_cabg_manuscript/post_hoc_shock_regressions.rds")) 
age = readRDS(here::here("results/data_cabg_manuscript/post_hoc_age65_regressions.rds"))
agetert = readRDS(here::here("results/data_cabg_manuscript/post_hoc_agetertiles_regressions.rds"))
cpb = readRDS(here::here("results/data_cabg_manuscript/post_hoc_cpb_regressions.rds"))
mediators =readRDS(here::here("results/data_cabg_manuscript/post_hoc_nomediator_regressions.rds"))


temp_df = 
  mediators %>% 
  filter(grepl("adjusted_MAP_nomed", type)) %>% 
  mutate(type = "Adjusted, no intra-op") %>% 
  bind_rows(overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "Adjusted, intra-op")) %>% 
  bind_rows(overall_map %>% filter(type == "univariate" & outcome == "AKI48") %>% 
              mutate(type = "Univariate")) %>% 
  mutate(type = factor(type, levels = c("Univariate", "Adjusted, no intra-op", "Adjusted, intra-op")),
         term = sub(".*MAP\\_", "", term))

temp_df %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  mutate(across(c(estimate, lb, ub, lb_cma, ub_cma), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"),
         ci_cma = if_else(is.na(lb_cma), "---", paste0("(", lb_cma, ",", ub_cma, ")"))) %>% 
  select("MAP Range" = term, "Odds Ratio" = estimate, "p-value" = p.value, "95 % CI" = ci, "95% CMA CI" = ci_cma, type) %>% 
  group_by(type) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_map_mediators.docx"))
  
temp_df = 
  mediators %>% 
  filter(grepl("adjusted_CVP_nomed", type)) %>% 
  mutate(type = "Adjusted, no intra-op") %>% 
  bind_rows(overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "Adjusted, intra-op")) %>% 
  bind_rows(overall_cvp %>% filter(type == "univariate" & outcome == "AKI48") %>% 
              mutate(type = "Univariate")) %>% 
  mutate(type = factor(type, levels = c("Univariate", "Adjusted, no intra-op", "Adjusted, intra-op")),
         term = sub(".*CVP\\_", "", term))

temp_df %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  mutate(across(c(estimate, lb, ub, lb_cma, ub_cma), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"),
         ci_cma = if_else(is.na(lb_cma), "---", paste0("(", lb_cma, ",", ub_cma, ")"))) %>% 
  select("CVP Range" = term, "Odds Ratio" = estimate, "p-value" = p.value, "95 % CI" = ci, "95% CMA CI" = ci_cma, type) %>% 
  group_by(type) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_cvp_mediators.docx"))



lvef %>% 
  mutate(term = case_when(
    term == "group_1" ~ "zone_1",
    term == "group_2" ~ "zone_2",
    term == "group_3" ~ "zone_3",
    term == "group_4" ~ "zone_4",
    term == "group_5" ~ "zone_5",
    .default = term
  )) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  mutate(type = factor(type,
         labels = c("CVP, Adjusted, LVEF > 40",
                    "CVP, Adjusted, LVEF <= 40",
                    "MAP, Adjusted, LVEF > 40",
                    "MAP, Adjusted, LVEF <= 40",
                    "3 Shells, Adjusted, LVEF > 40",
                    "3 Shells, Adjusted, LVEF <= 40",
                    "5 Shells, Adjusted, LVEF > 40",
                    "5 Shells, Adjusted, LVEF <= 40",
                  "CVP, Univariate, LVEF > 40",
                    "CVP, Univariate, LVEF <= 40",
                    "MAP, Univariate, LVEF > 40",
                    "MAP, Univariate, LVEF <= 40",
                    "3 Shells, Univariate, LVEF > 40",
                    "3 Shells, Univariate, LVEF <= 40",
                    "5 Shells, Univariate, LVEF > 40",
                    "5 Shells, Univariate, LVEF <= 40"
                    ))) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_lvef40.docx"))


egfr %>% 
  mutate(term = case_when(
    term == "group_1" ~ "zone_1",
    term == "group_2" ~ "zone_2",
    term == "group_3" ~ "zone_3",
    term == "group_4" ~ "zone_4",
    term == "group_5" ~ "zone_5",
    .default = term
  )) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3)),) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  mutate(type = factor(type,
                       labels = c("CVP, Adjusted, EGFR >= 60",
                                  "CVP, Adjusted, EGFR < 60",
                                  "MAP, Adjusted, EGFR >= 60",
                                  "MAP, Adjusted, EGFR < 60",
                                  "3 Shells, Adjusted, EGFR >= 60",
                                  "3 Shells, Adjusted, EGFR < 60",
                                  "5 Shells, Adjusted, EGFR >= 60",
                                  "5 Shells, Adjusted, EGFR < 60",
                                  "CVP, Univariate, EGFR >= 60",
                                  "CVP, Univariate, EGFR < 60",
                                  "MAP, Univariate, EGFR >= 60",
                                  "MAP, Univariate, EGFR < 60",
                                  "3 Shells, Univariate, EGFR >= 60",
                                  "3 Shells, Univariate, EGFR < 60",
                                  "5 Shells, Univariate, EGFR >= 60",
                                  "5 Shells, Univariate, EGFR < 60"
                       ))) %>% 
  dplyr::select(-lb, -ub) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_egfr60.docx"))


shock %>% 
  mutate(term = case_when(
    term == "group_1" ~ "zone_1",
    term == "group_2" ~ "zone_2",
    term == "group_3" ~ "zone_3",
    term == "group_4" ~ "zone_4",
    term == "group_5" ~ "zone_5",
    .default = term
  )) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  mutate(type = factor(type,
                       labels = c("CVP, Adjusted, No shock",
                                  "CVP, Adjusted, Shock",
                                  "MAP, Adjusted, No shock",
                                  "MAP, Adjusted, Shock",
                                  "3 Shells, Adjusted, No shock",
                                  "3 Shells, Adjusted, Shock",
                                  "5 Shells, Adjusted, No shock",
                                  "5 Shells, Adjusted, Shock",
                                  "CVP, Univariate, No shock",
                                  "CVP, Univariate, Shock",
                                  "MAP, Univariate, No shock",
                                  "MAP, Univariate, Shock",
                                  "3 Shells, Univariate, No shock",
                                  "3 Shells, Univariate, Shock",
                                  "5 Shells, Univariate, No shock",
                                  "5 Shells, Univariate, Shock"
                       ))) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_shock.docx"))


t1 = cpb %>% 
  filter(grepl("MAP", type)) %>% 
  mutate(phase = sub(".*MAP\\_", "", type),
         adjustment = sub("\\_MAP.*", "", type))  %>% 
  # )) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value, adjustment, type, phase) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  pivot_wider(names_from = phase, values_from = c(estimate, p.value, ci),
              id_cols = c(term, adjustment))


t2 = cpb %>% 
  filter(grepl("CVP", type)) %>% 
  mutate(phase = sub(".*CVP\\_", "", type),
         adjustment = sub("\\_CVP.*", "", type))  %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value, adjustment, type, phase) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  pivot_wider(names_from = phase, values_from = c(estimate, p.value, ci),
              id_cols = c(term, adjustment))


t3 = cpb %>% 
  filter(grepl("zone3", type) & type != "adjusted_zone3") %>% 
  mutate(phase = sub(".*zone3\\_", "", type),
         adjustment = sub("\\_zone3.*", "", type))  %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value, adjustment, type, phase) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  pivot_wider(names_from = phase, values_from = c(estimate, p.value, ci),
              id_cols = c(term, adjustment))


t4 = cpb %>% 
  filter(grepl("zone5", type) & type != "adjusted_zone5") %>% 
  mutate(phase = sub(".*zone5\\_", "", type),
         adjustment = sub("\\_zone5.*", "", type))  %>%
  # mutate(term = case_when(
  #   term == "group_1" ~ "zone_1",
  #   term == "group_2" ~ "zone_2",
  #   term == "group_3" ~ "zone_3",
  #   term == "group_4" ~ "zone_4",
  #   term == "group_5" ~ "zone_5",
  #   .default = term
  # )) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value, adjustment, type, phase) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  pivot_wider(names_from = phase, values_from = c(estimate, p.value, ci),
              id_cols = c(term, adjustment))

t1 %>% 
  bind_rows(t2) %>% 
  bind_rows(t3) %>% 
  bind_rows(t4) %>% 
  group_by(adjustment) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_bycpb.docx"))


age %>% 
  mutate(term = case_when(
    term == "group_1" ~ "zone_1",
    term == "group_2" ~ "zone_2",
    term == "group_3" ~ "zone_3",
    term == "group_4" ~ "zone_4",
    term == "group_5" ~ "zone_5",
    .default = term
  )) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  mutate(type = factor(type,
                       labels = c("CVP, Adjusted, Age >= 65",
                                  "CVP, Adjusted, Age < 65",
                                  "MAP, Adjusted, Age >= 65",
                                  "MAP, Adjusted, Age < 65",
                                  "3 Shells, Adjusted, Age >= 65",
                                  "3 Shells, Adjusted, Age < 65",
                                  "5 Shells, Adjusted, Age >= 65",
                                  "5 Shells, Adjusted, Age < 65",
                                  "CVP, Univariate, Age >= 65",
                                  "CVP, Univariate, Age < 65",
                                  "MAP, Univariate, Age >= 65",
                                  "MAP, Univariate, Age < 65",
                                  "3 Shells, Univariate, Age >= 65",
                                  "3 Shells, Univariate, Age < 65",
                                  "5 Shells, Univariate, Age >= 65",
                                  "5 Shells, Univariate, Age < 65"
                       ))) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_age.docx"))


agetert %>% 
  mutate(term = case_when(
    term == "group_1" ~ "zone_1",
    term == "group_2" ~ "zone_2",
    term == "group_3" ~ "zone_3",
    term == "group_4" ~ "zone_4",
    term == "group_5" ~ "zone_5",
    .default = term
  )) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 3))) %>% 
  dplyr::select(term, estimate, type, lb, ub, p.value) %>% 
  mutate(across(c(estimate, lb, ub), ~round(.x, 3)),
         ci = paste0("(", lb, ",", ub, ")"))  %>% 
  dplyr::select(-lb, -ub) %>% 
  mutate(type = factor(type,
                       labels = c("CVP, Adjusted, Age Tert 1",
                                  "CVP, Adjusted, Age Tert 2",
                                  "CVP, Adjusted, Age Tert 3",
                                  "MAP, Adjusted, Age Tert 1",
                                  "MAP, Adjusted, Age Tert 2",
                                  "MAP, Adjusted, Age Tert 3",
                                  "3 Shells, Adjusted, Age Tert 1",
                                  "3 Shells, Adjusted, Age Tert 2",
                                  "3 Shells, Adjusted, Age Tert 3",
                                  "5 Shells, Adjusted, Age Tert 1",
                                  "5 Shells, Adjusted, Age Tert 2",
                                  "5 Shells, Adjusted, Age Tert 3",
                                  "CVP, Univariate, Age Tert 1",
                                  "CVP, Univariate, Age Tert 2",
                                  "CVP, Univariate, Age Tert 3",
                                  "MAP, Univariate, Age Tert 1",
                                  "MAP, Univariate, Age Tert 2",
                                  "MAP, Univariate, Age Tert 3",
                                  "3 Shells, Univariate, Age Tert 1",
                                  "3 Shells, Univariate, Age Tert 2",
                                  "3 Shells, Univariate, Age Tert 3",
                                  "5 Shells, Univariate, Age Tert 1",
                                  "5 Shells, Univariate, Age Tert 2",
                                  "5 Shells, Univariate, Age Tert 3"
                                  
                       ))) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/regressions_agetert.docx"))



