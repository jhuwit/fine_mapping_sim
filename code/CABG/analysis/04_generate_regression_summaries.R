library(here)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
set.seed(123)
options(dplyr.summarise.inform=F)


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



# plot - old stuff  

plot_df <-
  hemo %>%
  group_by(id) %>%
  mutate(
    cut_MAP = cut(
      val_MAP,
      breaks = map_thresholds_small,
      include.lowest = T
    ),
    cut_CVP = cut(
      val_CVP,
      breaks = cvp_thresholds_small,
      include.lowest = T
    )
  ) %>%
  count(cut_MAP, cut_CVP, .drop = F) %>% drop_na() %>% group_by(cut_MAP, cut_CVP) %>%
  summarize(mean = mean(n)) %>% mutate(MAP_CVP = paste0("`MAP", cut_MAP, "CVP", cut_CVP, "`")) %>%
  left_join(bivar_mediators_small, by = c("MAP_CVP" = "term")) %>%
  ungroup() %>% mutate(sig = ifelse(p.value < 0.05, "Yes", "No"))

plot = 
  plot_df %>%
  mutate(label = ifelse(p.value < 0.05, paste0(round(estimate, 2), "*"), paste0(round(estimate, 2)))) %>%
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = estimate,
    label = label
  )) + geom_tile() +
  geom_text() +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio",
    title = "Adjusted"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(),
    legend.text = element_text(
      angle = 45,
      size = 8,
      hjust = 1
    )
  ) +  scale_fill_gradient2(low ="blue",
                            high = "red",
                            mid = "white",
                            midpoint = 1)+ 
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    title = element_text(size = 15),
    legend.text=element_text(size = 11))

plot_df <-
  hemo %>%
  group_by(id) %>%
  mutate(
    cut_MAP = cut(
      val_MAP,
      breaks = map_thresholds_small,
      include.lowest = T
    ),
    cut_CVP = cut(
      val_CVP,
      breaks = cvp_thresholds_small,
      include.lowest = T
    )
  ) %>%
  count(cut_MAP, cut_CVP, .drop = F) %>% drop_na() %>% group_by(cut_MAP, cut_CVP) %>%
  summarize(mean = mean(n)) %>% mutate(MAP_CVP = paste0("`MAP", cut_MAP, "CVP", cut_CVP, "`")) %>%
  left_join(bivar_uni_small, by = c("MAP_CVP" = "term")) %>%
  ungroup() %>% mutate(sig = ifelse(p.value < 0.05, "Yes", "No"))

uni_plot = 
  plot_df %>%
  mutate(label = ifelse(p.value < 0.05, paste0(round(estimate, 2), "*"), paste0(round(estimate, 2)))) %>%
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = estimate,
    label = label
  )) + geom_tile() +
  geom_text() +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio",
    title = "Univariate"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(),
    legend.text = element_text(
      angle = 45,
      size = 8,
      hjust = 1
    )
  ) + 
  scale_fill_gradient2(low ="blue",
                       high = "red",
                       mid = "white",
                       midpoint = 1)+ 
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    title = element_text(size = 15),
    legend.text=element_text(size = 11))

cowplot::plot_grid(uni_plot, plot)



# shell regressions 
shell_res = readRDS(here::here("results", "data_cabg_manuscript", "shell_regressions_5.rds"))


shell_res %>% 
  mutate(ub = estimate - std.error) %>% 
  select(term, estimate, lb = conf.low, ub = conf.high, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("MAP Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/shells_5.docx"))

shell_res = readRDS(here::here("results", "data_cabg_manuscript", "shell_regressions_3.rds"))


shell_res %>% 
  mutate(ub = estimate - std.error) %>% 
  select(term, estimate, lb = conf.low, ub = conf.high, p.value) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("MAP Range", "Odds Ratio", "P Value", "95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/shells_3.docx"))

# next up is coarser shells map_thresholds_10 = seq(45, 115, 10)

map_summaries = readRDS(here::here("results/data_cabg_manuscript/map_regressions_telescoping.rds"))
cvp_summaries = readRDS(here::here("results/data_cabg_manuscript/cvp_regressions_telescoping.rds"))


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


# adjusted 

map_summaries %>% 
  filter(type == "adjusted") %>% 
  # arrange(desc(term)) %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("MAP Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/MAP_adjusted_48h_telescoping.docx"))

cvp_summaries %>% 
  filter(type == "adjusted") %>% 
  select(term, estimate, lb, ub, lb_cma, ub_cma, p.value, cma_p) %>% 
  mutate(across(-term, ~round(.x, 3))) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", p.value)) %>% 
  mutate(cma_p = ifelse(cma_p < 0.001, "<0.001", cma_p)) %>% 
  mutate(CI_reg = paste0("(", lb, ", ", ub, ")"),
         CI_cma =  paste0("(", lb_cma, ", ", ub_cma, ")")) %>% 
  select(-c(starts_with("lb"), starts_with("ub"))) %>% 
  magrittr::set_colnames(c("CVP Range", "Odds Ratio", "P Value", "CMA P Value", "95% CI", "CMA 95% CI")) %>% 
  gt::gt() %>% 
  gt::gtsave(filename = here::here("results/tables_cabg_manuscript/CVP_adjusted_48h_telescoping.docx"))
