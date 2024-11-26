library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
set.seed(123)
options(dplyr.summarise.inform=F)


# cohort = "full" 
cohort = "cabg"
# covariate data 
covar_data = read_csv(here::here("data/analytic/covariates_post_exclusion.csv"))
# hemo data 
hemo_data = readRDS(here::here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))

# make categorical variables factors for logistic reg 
covar_data = covar_data %>% 
  mutate(across(starts_with("cat"), as.factor))

# filter to just CABG ids 
if(cohort == "cabg"){
  covar_data = covar_data %>% filter(cat_surgcat4 == "1_CABG")
  hemo_data = hemo_data %>% filter(id %in% covar_data$id)
}


source(here::here("code/CABG/analysis/utilities.R"))

map_thresholds = seq(45, 115, 5)
cvp_thresholds = seq(0, 20, 2)

map_data = get_ranges(hemo_data = hemo_data,
                      thresholds = map_thresholds,
                      hemo_variable = "MAP")
cvp_data = get_ranges(hemo_data = hemo_data,
                      thresholds = cvp_thresholds,
                      hemo_variable = "CVP")

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
  pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(lb = mean - sd, ub = mean + sd) %>% 
  pivot_longer(cols = c(num5, pct5, mean)) %>% 
  filter(name %in% c("pct5", "mean")) %>% 
  mutate(value = ifelse(name == "pct5", value * 100, value),
         lb = ifelse(name == "mean", lb, NA),
         ub = ifelse(name == "mean", ub, NA)) %>% 
  ggplot(aes(x = metric, y = value, fill = name))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin = lb, ymax = ub), position = position_dodge(),
                size = .5)+
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
  pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(lb = mean - sd, ub = mean + sd) %>% 
  pivot_longer(cols = c(num5, pct5, mean)) %>% 
  filter(name %in% c("pct5", "mean")) %>% 
  mutate(value = ifelse(name == "pct5", value * 100, value),
         lb = ifelse(name == "mean", lb, NA),
         ub = ifelse(name == "mean", ub, NA)) %>% 
  ggplot(aes(x = metric, y = value, fill = name))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin = lb, ymax = ub), position = position_dodge(),
                size = .5)+
  theme_bw()+
  scale_fill_brewer(palette = "Dark2", direction = -1, name = "", labels = c("Mean minutes in range", "Percent individuals with at least 5 minutes in range"))+
  labs(x = "MAP (mmHg)", y = "")+
  theme( legend.position = "bottom")
