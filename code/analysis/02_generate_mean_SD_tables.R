library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
set.seed(123)
options(dplyr.summarise.inform=F)

if(!dir.exists(here::here("results/tables_cabg_manuscript"))){
  dir.create(here::here("results/tables_cabg_manuscript"), recursive = TRUE)
}
# cohort = "full" 
cohort = "cabg"
# covariate data 
covar_data = read_csv(here::here("data/analytic/covariates_post_exclusion.csv"),
                      col_types = cols(id = col_character()))
# hemo data 
hemo_data = readRDS(here::here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))

# make categorical variables factors for logistic reg 
covar_data = covar_data %>% 
  mutate(across(starts_with("cat"), as.factor))

# filter to just CABG ids 
# if(cohort == "cabg"){
#   covar_data = covar_data %>% filter(cat_surgcat4 == "1_CABG")
#   hemo_data = hemo_data %>% filter(id %in% covar_data$id)
# }

primary = colnames(covar_data)[c(6:24)]
primary_and_mediators = colnames(covar_data)[c(6:27)]

all_covars <- colnames(covar_data)[c(6:27, 32, 33, 2, 3, 45, 34, 47)]


source(here::here("code/analysis/utilities.R"))

map_thresholds = seq(45, 115, 5)
map_thresholds_bivar = seq(45, 115, 10)
cvp_thresholds = seq(0, 20, 2)

map_data = get_ranges(hemo_data = hemo_data,
                      thresholds = map_thresholds,
                      hemo_variable = "MAP")
cvp_data = get_ranges(hemo_data = hemo_data,
                      thresholds = cvp_thresholds,
                      hemo_variable = "CVP")

bivar_bricks = get_ranges_biv(hemo_data, map_thresholds = map_thresholds_bivar,
                              cvp_thresholds = cvp_thresholds)

table = 
  map_data %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(starts_with("MAP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% pivot_longer(cols = -bin_aki48h) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = c(measure, bin_aki48h), values_from = value) %>% 
  mutate(aki0_meansd = paste0(round(mean_0, 3), " (", round(sd_0, 3), ")"),
         aki1_meansd = paste0(round(mean_1, 3), " (", round(sd_1, 3), ")"),
         aki0_npct = paste0(round(num5_0, 3), " (", round(pct5_0, 3), ")"),
         aki1_npct = paste0(round(num5_1, 3), " (", round(pct5_1, 3), ")")) %>% 
  select(metric, aki0_meansd, aki1_meansd, aki0_npct, aki1_npct) %>% 
  magrittr::set_colnames(c("MAP Range", "No AKI at 48h mean (SD) ", "AKI at 48h mean (SD)",
                           "No AKI at 48h n (%)", "AKI at 48 n (%)")) 

# table %>% 
#   kableExtra::kbl(caption = "Mean and SD of Minutes in Each Range, MAP",
#                   format = "html",
#                   align = "l") %>% 
#   kableExtra::kable_classic(full_width = TRUE)

table %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_mean_sd_48h.docx"))

table2 = 
  cvp_data %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(starts_with("CVP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% pivot_longer(cols = -bin_aki48h) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = c(measure, bin_aki48h), values_from = value) %>% 
  mutate(aki0_meansd = paste0(round(mean_0, 3), " (", round(sd_0, 3), ")"),
         aki1_meansd = paste0(round(mean_1, 3), " (", round(sd_1, 3), ")"),
         aki0_npct = paste0(round(num5_0, 3), " (", round(pct5_0, 3), ")"),
         aki1_npct = paste0(round(num5_1, 3), " (", round(pct5_1, 3), ")")) %>% 
  select(metric, aki0_meansd, aki1_meansd, aki0_npct, aki1_npct) %>% 
  magrittr::set_colnames(c("CVP Range", "No AKI at 48h mean (SD) ", "AKI at 48h mean (SD)",
                           "No AKI at 48h n (%)", "AKI at 48 n (%)")) 


# table2 %>% 
#   kableExtra::kbl(caption = "Mean and SD of Minutes in Each Range, CVP",
#                   format = "html",
#                   align = "l") %>% 
#   kableExtra::kable_classic(full_width = TRUE)

table2 %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_mean_sd_48h.docx"))

table3 = 
  bivar_bricks %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(starts_with("MAP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% pivot_longer(cols = -bin_aki48h) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = sub("\\_.*", "", name)) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = c(measure, bin_aki48h), values_from = value) %>% 
  mutate(aki0_meansd = paste0(round(mean_0, 3), " (", round(sd_0, 3), ")"),
         aki1_meansd = paste0(round(mean_1, 3), " (", round(sd_1, 3), ")"),
         aki0_npct = paste0(round(num5_0, 3), " (", round(pct5_0, 3), ")"),
         aki1_npct = paste0(round(num5_1, 3), " (", round(pct5_1, 3), ")")) %>% 
  select(metric, aki0_meansd, aki1_meansd, aki0_npct, aki1_npct) %>% 
  magrittr::set_colnames(c("Range", "No AKI at 48h mean (SD) ", "AKI at 48h mean (SD)",
                           "No AKI at 48h n (%)", "AKI at 48 n (%)")) 


# table3 %>% 
#   kableExtra::kbl(caption = "Mean and SD of Minutes in Each Range",
#                   format = "html",
#                   align = "l") %>% 
#   kableExtra::kable_classic(full_width = TRUE)

table3 %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_CVP_mean_sd_48h.docx"))

table = 
  map_data %>% 
  left_join(covar_data %>% select(id, bin_akicomposite)) %>% 
  group_by(bin_akicomposite) %>% 
  summarize(across(starts_with("MAP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n())))  %>% pivot_longer(cols = -bin_akicomposite) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = c(measure, bin_akicomposite), values_from = value) %>% 
  mutate(aki0_meansd = paste0(round(mean_0, 3), " (", round(sd_0, 3), ")"),
         aki1_meansd = paste0(round(mean_1, 3), " (", round(sd_1, 3), ")"),
         aki0_npct = paste0(round(num5_0, 3), " (", round(pct5_0, 3), ")"),
         aki1_npct = paste0(round(num5_1, 3), " (", round(pct5_1, 3), ")")) %>% 
  select(metric, aki0_meansd, aki1_meansd, aki0_npct, aki1_npct) %>% 
  magrittr::set_colnames(c("Range", "No AKI mean (SD) ", "AKI mean (SD)",
                           "No AKI n (%)", "AKI at 48 n (%)")) 


table %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_mean_sd_composite.docx"))

table2 = 
  cvp_data %>% 
  left_join(covar_data %>% select(id, bin_akicomposite)) %>% 
  group_by(bin_akicomposite) %>% 
  summarize(across(starts_with("CVP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% pivot_longer(cols = -bin_akicomposite) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = c(measure, bin_akicomposite), values_from = value) %>% 
  mutate(aki0_meansd = paste0(round(mean_0, 3), " (", round(sd_0, 3), ")"),
         aki1_meansd = paste0(round(mean_1, 3), " (", round(sd_1, 3), ")"),
         aki0_npct = paste0(round(num5_0, 3), " (", round(pct5_0, 3), ")"),
         aki1_npct = paste0(round(num5_1, 3), " (", round(pct5_1, 3), ")")) %>% 
  select(metric, aki0_meansd, aki1_meansd, aki0_npct, aki1_npct) %>% 
  magrittr::set_colnames(c("Range", "No AKI mean (SD) ", "AKI mean (SD)",
                           "No AKI n (%)", "AKI at 48 n (%)")) 

cvp_data %>% 
  summarize(across(starts_with("CVP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% 
  mutate(id = "id") %>% 
  pivot_longer(cols = -id) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name, -id) %>% 
  mutate(metric = factor(metric, levels = c("[0,2]", 
                                            "(2,4]",
                                            "(4,6]",
                                            "(6,8]",
                                            "(8,10]",
                                            "(10,12]", 
                                            "(12,14]", "(14,16]",
                                            "(16,18]", "(18,20]"))) %>% 
  filter(measure %in% c("pct5", "mean")) %>% 
  ggplot(aes(x = metric, y = value, fill = measure))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_classic()+
  scale_fill_discrete(name = "", labels = c("Mean minutes", "Proportion"))+
  labs(x = "CVP", y = "") 




table2 %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_mean_sd_composite.docx"))

table3 = 
  bivar_bricks %>% 
  left_join(covar_data %>% select(id, bin_akicomposite)) %>% 
  group_by(bin_akicomposite) %>% 
  summarize(across(starts_with("MAP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% pivot_longer(cols = -bin_akicomposite) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = sub("\\_.*", "", name)) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = c(measure, bin_akicomposite), values_from = value) %>% 
  mutate(aki0_meansd = paste0(round(mean_0, 3), " (", round(sd_0, 3), ")"),
         aki1_meansd = paste0(round(mean_1, 3), " (", round(sd_1, 3), ")"),
         aki0_npct = paste0(round(num5_0, 3), " (", round(pct5_0, 3), ")"),
         aki1_npct = paste0(round(num5_1, 3), " (", round(pct5_1, 3), ")")) %>% 
  select(metric, aki0_meansd, aki1_meansd, aki0_npct, aki1_npct) %>% 
  magrittr::set_colnames(c("Range", "No AKI mean (SD) ", "AKI mean (SD)",
                           "No AKI n (%)", "AKI at 48 n (%)")) 


table3 %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_CVP_mean_sd_composite.docx"))


# quick figure generation 
temp = cvp_data %>% 
  summarize(across(starts_with("CVP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% 
  mutate(id = "id") %>% 
  pivot_longer(cols = -id) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name, -id) %>% 
  mutate(metric = factor(metric, levels = c("[0,2]", 
                                            "(2,4]",
                                            "(4,6]",
                                            "(6,8]",
                                            "(8,10]",
                                            "(10,12]", 
                                            "(12,14]", "(14,16]",
                                            "(16,18]", "(18,20]"))) %>% 
 pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(lb = mean - sd,
         ub = mean + sd)
  
cvp_data %>% 
  summarize(across(starts_with("CVP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% 
  mutate(id = "id") %>% 
  pivot_longer(cols = -id) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name, -id) %>% 
  mutate(metric = factor(metric, levels = c("[0,2]", 
                                            "(2,4]",
                                            "(4,6]",
                                            "(6,8]",
                                            "(8,10]",
                                            "(10,12]", 
                                            "(12,14]", "(14,16]",
                                            "(16,18]", "(18,20]"))) %>% 
  filter(measure %in% c("pct5", "mean")) %>% 
  mutate(value = ifelse(measure == "pct5", value * 100, value)) %>% 
  ggplot(aes(x = metric, y = value, fill = measure))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  scale_fill_brewer(palette = "Dark2", direction = -1, name = "", labels = c("Mean minutes in range", "Percent individuals with at least 5 minutes in range"))+
  labs(x = "CVP", y = "")+
  theme( legend.position = "bottom")


map_data %>% 
  summarize(across(starts_with("MAP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% 
  mutate(id = "id") %>% 
  pivot_longer(cols = -id) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name, -id) %>% 
  mutate(metric = factor(metric, levels = c("[45,50]", 
                                            "(50,55]",
                                            "(55,60]",
                                            "(60,65]",
                                            "(65,70]",
                                            "(70,75]",
                                            "(75,80]",
                                            "(80,85]",
                                            "(85,90]",
                                            "(90,95]",
                                            "(95,100]",
                                            "(100,105]",
                                            "(105,110]",
                                            "(110,115]"))) %>% 
  filter(measure %in% c("pct5", "mean")) %>% 
  mutate(value = ifelse(measure == "pct5", value * 100, value)) %>% 
  ggplot(aes(x = metric, y = value, fill = measure))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  scale_fill_brewer(palette = "Dark2", direction = -1, name = "", labels = c("Mean minutes in range", "Percent individuals with at least 5 minutes in range"))+
  labs(x = "MAP", y = "")+
  theme(legend.position = "bottom")


## tertile tables 
map_65_cvp_12 = 
  hemo_data %>% 
  group_by(id) %>% 
  summarize(
    min_65 = sum(val_MAP < 65, na.rm = TRUE),
    min_12 = sum(val_CVP > 12, na.rm = TRUE)
  )
  
map_tertiles = quantile(map_65_cvp_12$min_65, c(1/3, 2/3)); map_tertiles
cvp_tertiles = quantile(map_65_cvp_12$min_12, c(1/3, 2/3)); cvp_tertiles

map_65_cvp_12 = 
  map_65_cvp_12 %>% 
  mutate(
    map_tertile = case_when(
      min_65 <= 64 ~ "t1",
      min_65 <= 105 ~ "t2", 
      TRUE ~ "t3"
    ),
    cvp_tertile = case_when(
      min_12 <= 21 ~ "t1",
      min_12 <= 64 ~ "t2",
      TRUE ~ "t3"
    )
  )

table_dat = 
  map_65_cvp_12 %>% 
  select(id, map_tertile) %>% 
  left_join(covar_data %>% select(id, all_of(primary_and_mediators), val_proctime)) %>% 
  select(-id)

myvars = colnames(table_dat %>% select(-map_tertile))
factorvars = table_dat %>% select(starts_with("cat"), starts_with("bin")) %>% colnames()
tab = 
  tableone::CreateTableOne(myvars,
                         strata = "map_tertile",
                         factorVars = factorvars, 
                         test = TRUE,
                         data = table_dat)
tabprint = print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(tabprint, here::here("results/tables_cabg_manuscript/map_tertiles_table.csv"))

table_dat = 
  map_65_cvp_12 %>% 
  select(id, cvp_tertile) %>% 
  left_join(covar_data %>% select(id, all_of(primary_and_mediators),val_proctime)) %>% 
  select(-id)

myvars = colnames(table_dat %>% select(-cvp_tertile))
factorvars = table_dat %>% select(starts_with("cat"), starts_with("bin")) %>% colnames()
tab = 
  tableone::CreateTableOne(myvars,
                           strata = "cvp_tertile",
                           factorVars = factorvars, 
                           test = TRUE,
                           data = table_dat)
tabprint = print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(tabprint, here::here("results/tables_cabg_manuscript/cvp_tertiles_table.csv"))

# table one 
table_dat = 
  covar_data %>% 
  select(id, all_of(primary_and_mediators),val_proctime, bin_aki48h) %>% 
  select(-id)

myvars = colnames(table_dat %>% select(-bin_aki48h))
factorvars = table_dat %>% select(starts_with("cat"), starts_with("bin")) %>% colnames()
tab = 
  tableone::CreateTableOne(myvars,
                           strata = "bin_aki48h",
                           factorVars = factorvars, 
                           test = TRUE,
                           data = table_dat)
tabprint = print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(tabprint, here::here("results/tables_cabg_manuscript/table_one_w_pvalues.csv"))


# mean sd tables, not stratified by AKI 

table = 
  map_data %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  summarize(across(starts_with("MAP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% pivot_longer(cols = starts_with("MAP")) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(aki_meansd = paste0(round(mean, 3), " (", round(sd, 3), ")"),
         aki_npct = paste0(round(num5, 3), " (", round(pct5, 3), ")")) %>% 
  select(metric, aki_meansd, aki_npct) %>% 
  magrittr::set_colnames(c("MAP Range", "Mean (SD)", "n (%)")) 

table %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_mean_sd_overall.docx"))

table2 = 
  cvp_data %>% 
  summarize(across(starts_with("CVP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% 
  pivot_longer(cols = starts_with("CVP")) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = strsplit(name, "_")[[1]][2]) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(aki_meansd = paste0(round(mean, 3), " (", round(sd, 3), ")"),
         aki_npct = paste0(round(num5, 3), " (", round(pct5, 3), ")")) %>% 
  select(metric, aki_meansd, aki_npct) %>% 
  magrittr::set_colnames(c("CVP Range", "Mean (SD)", "n (%)")) 

table2 %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_mean_sd_overall.docx"))

table3 = 
  bivar_bricks %>% 
  summarize(across(starts_with("MAP") & !contains("missing"),
                   list(mean = mean,
                        sd = sd,
                        num5 = ~ sum(.x>= 5),
                        pct5 = ~sum(.x >= 5)/n()))) %>% 
  pivot_longer(cols = starts_with("MAP")) %>% 
  rowwise() %>% 
  mutate(measure = sub(".*\\_", "", name),
         metric = sub("\\_.*", "", name)) %>% 
  ungroup() %>% 
  select(-name) %>% 
  pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(aki_meansd = paste0(round(mean, 3), " (", round(sd, 3), ")"),
         aki_npct = paste0(round(num5, 3), " (", round(pct5, 3), ")")) %>% 
  select(metric, aki_meansd, aki_npct) %>% 
  magrittr::set_colnames(c("Range", "Mean (SD)", "n (%)")) 


table3 %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_CVP_mean_sd_overall.docx"))

