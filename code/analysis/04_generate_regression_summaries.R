library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
set.seed(123)
options(dplyr.summarise.inform=F)
covar_data = read_csv(here::here("data/analytic/covariates_post_exclusion.csv"),
                      col_types = cols(id = col_character()))


# read in map results 

map_summaries = readRDS(here::here("results/data_cabg_manuscript/map_regressions.rds"))
cvp_summaries = readRDS(here::here("results/data_cabg_manuscript/cvp_regressions.rds"))
bv_summaries = readRDS(here::here("results/data_cabg_manuscript/bricks_regressions.rds"))

# add in p values 
# se = (temp$ub_cma - temp$lb_cma) / (2*1.96)
# z = temp$estimate / se
# pnorm(abs(z))
# 2*pnorm(-abs(z))
map_summaries = 
  map_summaries %>% 
  rowwise() %>% 
  mutate(cma_p = 2* pnorm(-abs(log(estimate) / 
                                 ((log(ub_cma) - log(lb_cma)) / (2 * 1.96))))) %>% 
  ungroup()

cvp_summaries = 
  cvp_summaries %>% 
  rowwise() %>% 
  mutate(cma_p = 2* pnorm(-abs(log(estimate) / 
                                 ((log(ub_cma) - log(lb_cma)) / (2 * 1.96))))) %>% 
  ungroup()

bv_summaries = 
  bv_summaries %>% 
  rowwise() %>% 
  mutate(cma_p = 2* pnorm(-abs(log(estimate) / 
                                 ((log(ub_cma) - log(lb_cma)) / (2 * 1.96))))) %>% 
  ungroup()

# univariate, primary analysis 
map_summaries %>% 
  filter(type == "univariate" & outcome == "AKI48") %>% 
  select(term, estimate, lb, ub, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("MAP Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_univariate_48h.docx"))

cvp_summaries %>% 
  filter(type == "univariate" & outcome == "AKI48") %>% 
  select(term, estimate, lb, ub, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("CVP Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_univariate_48h.docx"))

bv_summaries %>% 
  filter(type == "univariate" & outcome == "AKI48") %>% 
  select(term, estimate, lb, ub, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("Bivariate Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/bricks_univariate_48h.docx"))

# adjusted 

map_summaries %>% 
  filter(type == "adjusted" & outcome == "AKI48") %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("MAP Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_adjusted_48h.docx"))

cvp_summaries %>% 
  filter(type == "adjusted" & outcome == "AKI48") %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("CVP Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_adjusted_48h.docx"))

bv_summaries %>% 
  filter(type == "adjusted" & outcome == "AKI48") %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/bricks_adjusted_48h.docx"))

# repeat for composite outcome # univariate, primary analysis 
map_summaries %>% 
  filter(type == "univariate" & outcome == "AKI Composite") %>% 
  select(term, estimate, lb, ub, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("MAP Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_univariate_composite.docx"))

cvp_summaries %>% 
  filter(type == "univariate" & outcome == "AKI Composite") %>% 
  select(term, estimate, lb, ub, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("CVP Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_univariate_composite.docx"))

bv_summaries %>% 
  filter(type == "univariate" & outcome == "AKI Composite") %>% 
  select(term, estimate, lb, ub, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("Bivariate Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/bricks_univariate_composite.docx"))

# adjusted 

map_summaries %>% 
  filter(type == "adjusted" & outcome == "AKI Composite") %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("MAP Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_adjusted_composite.docx"))

cvp_summaries %>% 
  filter(type == "adjusted" & outcome == "AKI Composite") %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("CVP Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_adjusted_composite.docx"))

bv_summaries %>% 
  filter(type == "adjusted" & outcome == "AKI Composite") %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/bricks_adjusted_composite.docx"))



