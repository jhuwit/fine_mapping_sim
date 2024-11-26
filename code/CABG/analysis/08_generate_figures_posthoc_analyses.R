########## NEW 
# process results 
library(tidyverse)
library(gt)
rm(list = ls())
options(dplyr.summarise.inform=F)
library(paletteer)
library(ggrepel)
# source(here::here("code/analysis/utilities.R"))
# scales::show_col(hue_pal()(3))

col1a = "#7873C0FF"; col2a = "#21B087FF"; col3a = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"
col1 = "#E69F00FF"; col2 = "#56B4E9FF"; col3 = "#009E73FF"; col4 = "#CC79A7FF"

col_vector = c(col1, col2, col3, col6, col5, col4)

shell_res = readRDS(here::here("results/data_cabg_manuscript/shell_regressions.rds"))

range_left = function(x, left, right){
  ifelse(x > left & x <= right, TRUE, FALSE)
}
hemo_data = readRDS(here::here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))

xs = c(0.5, 1.5, 1.5, 2.5, 2.5, 4.5, 4.5, 1.5, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5, 4.5,
       5.5, 5.5)
xends = c(1.5, 1.5, 2.5, 2.5, 4.5, 7.5, 4.5, 1.5, 2.5, 2.5, 7.5, 3.5, 4.5, 4.5, 7.5, 4.5,
          5.5, 7.5)
ys = c(4.5, 4.5, 6.5, 6.5, 8.5, 9.5, 8.5, 0.5, 4.5, 4.5, 6.5, 0.5, 4.5, 5.5, 5.5, 4.5, 
       0.5, 4.5)
yends = c(4.5, 6.5, 6.5, 8.5, 8.5, 9.5, 9.5, 4.5, 4.5, 6.5, 6.5, 4.5, 4.5, 5.5, 5.5, 5.5,
          4.5, 4.5)

dat = tibble(xs, xends, ys, yends)

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

lvef = readRDS(here::here("results/data_cabg_manuscript/post_hoc_lvef40_regressions.rds")) 
egfr = readRDS(here::here("results/data_cabg_manuscript/post_hoc_egfr60_regressions.rds")) 
shock = readRDS(here::here("results/data_cabg_manuscript/post_hoc_shock_regressions.rds")) 
age = readRDS(here::here("results/data_cabg_manuscript/post_hoc_age65_regressions.rds"))
agetert = readRDS(here::here("results/data_cabg_manuscript/post_hoc_agetertiles_regressions.rds"))
mediators =readRDS(here::here("results/data_cabg_manuscript/post_hoc_nomediator_regressions.rds"))
predrenf = readRDS(here::here("results/data_cabg_manuscript/post_hoc_predrenf_regressions.rds"))
catlvef = readRDS(here::here("results/data_cabg_manuscript/post_hoc_catlvef_regressions.rds"))
categfr = readRDS(here::here("results/data_cabg_manuscript/post_hoc_categfr_regressions.rds"))
contegfr = readRDS(here::here("results/data_cabg_manuscript/post_hoc_contegfr_regressions.rds"))

telescoping_map = readRDS(here::here("results/data_cabg_manuscript/map_regressions_telescoping.rds"))
telescoping_cvp = readRDS(here::here("results/data_cabg_manuscript/cvp_regressions_telescoping.rds"))

plot_j_fancy = function(df, new_levels, n_comparisons, xlabs,
                        xname, yname, title, lims, legend_pos = "bottom", 
                        col_vector, dodge = 0.6){
  # paletteer_d("ggthemes::colorblind")
  # col1 = "#7873C0FF"; col2 = "#21B087FF"; col3 = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"
  # col1 = "#E69F00FF"; col2 = "#56B4E9FF"; col3 = "#009E73FF"; col4 = "#CC79A7FF"
  # 
  # col_vector = c(col1, col2, col3, col4)
  df %>% 
    mutate(term = factor(term, levels = new_levels)) %>% 
    ggplot(aes(x = term, y = estimate))+
    geom_point(aes(x = term, y = estimate, col = type, shape = type), position = position_dodge(dodge),
               size = 2) +
    geom_errorbar(width = .5, linewidth = .5, aes(
      x = term,
      y = estimate,
      ymin = lb,
      ymax = ub,
      col = type
    ),
    position = position_dodge(dodge)) +
    scale_shape_manual(values = c(1, 16, 4, 6), name = "95% CI")+
    theme_bw() +
    scale_x_discrete(labels = xlabs)+
    scale_y_continuous(limits = lims)+ 
    labs(x = paste(xname), y = paste(yname),
         title = paste(title)) +
    geom_hline(aes(yintercept = 1), col = "black")+
    scale_color_manual(
      values = col_vector[1:n_comparisons],
      name = "95% CI"
    ) +
    theme(legend.position = legend_pos,
          axis.text.x = element_text(size = 12),
          legend.title = element_blank(),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          # legend.title = element_text(size = 16),
          axis.title.y = element_text(margin = margin(r = 3), size = 14),
          axis.title.x = element_text(margin = margin(t = 10), size = 14))+
    guides(color = guide_legend(nrow = 4))
}

temp_df = 
  res %>% 
  filter(type %in% c("adjusted_zone5", "univariate_zone5")) %>% 
  mutate(type = factor(type, labels = c("adj", "uni")))
temp = temp_df %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(85,95]", 2),
                        rep("(55,65]", 2)),
            cut_CVP = c(rep("[0,2]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(16,18]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p.val"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group))) %>% 
  # mutate(ci = paste0("(", signif(lb, 3), ",", signif(ub, 3), ")")) %>% 
  mutate(ci = paste0("(", round(lb, 2), ",", round(ub, 2), ")")) %>% 
  rename(est = estimate, p = p.value) %>%
  select(zone, est, p, cut_MAP, cut_CVP, type, ci) %>%
  pivot_wider(names_from = type, values_from = c(est, p, ci))

p2 = plotdf_5shell %>%
  left_join(temp_df %>% filter(type == "adj") , by = c("group" = "term")) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, type) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est,
  )) +
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        # strip.text = element_text(size = 14),
        # axis.title = element_text(size = 14),
        legend.position = "none",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1,
        axis.title.x = element_text(margin=margin(t=7,r=0,b=0,l=0)),
        axis.title.y = element_text(margin=margin(t=0,r=7,b=0,l=0)),
        text=element_text(size=16))+
  scale_x_discrete(labels =c("[45-55]","(55-65]","(65-75]","(75-85]","(85-95]","(95-105]","(105-115]"))+
  scale_y_discrete(labels = c("[0-2]","(2,4]","(4-6]","(6-8]", "(8-10]","(10-12]","(12-14]","(14-16]","(16-18]","(18-20]"))+
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =1, linetype ="dashed")  +
  geom_label(data = temp %>% mutate(est = est_adj), aes(x = cut_MAP, y = cut_CVP,
                                                        label = paste0(zone, 
                                                                       "\n Unadj OR (95%CI): ", format(round(est_uni, 2),nsmall=2), " ", ci_uni,
                                                                       "\n Adj OR (95%CI): ", format(round(est_adj, 2),nsmall=2), " ", ci_adj)),
             col = "black", hjust = 0.25,size=4)+scale_color_discrete(name="")

p2
# export 1200 x 800


p2 = plotdf_5shell %>%
  left_join(temp_df %>% filter(type == "adj") , by = c("group" = "term")) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, type) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est,
  )) +
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        # strip.text = element_text(size = 14),
        # axis.title = element_text(size = 14),
        legend.position = "none",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1,
        axis.title.x = element_text(margin=margin(t=7,r=0,b=0,l=0)),
        axis.title.y = element_text(margin=margin(t=0,r=7,b=0,l=0)),
        text=element_text(size=16))+
  scale_x_discrete(labels =c("[45-55]","(55-65]","(65-75]","(75-85]","(85-95]","(95-105]","(105-115]"))+
  scale_y_discrete(labels = c("[0-2]","(2,4]","(4-6]","(6-8]", "(8-10]","(10-12]","(12-14]","(14-16]","(16-18]","(18-20]"))+
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =1, linetype ="dashed")  +
  geom_label(data = temp %>% mutate(est = est_adj), aes(x = cut_MAP, y = cut_CVP,
                                                        label = paste0(zone, 
                                                                       "\n Unadj OR=", format(round(est_uni, 2),nsmall=2), "; ", p_uni,
                                                                       "\n Adj OR=", format(round(est_adj, 2),nsmall=2), "; ", p_adj)),
             col = "black", hjust = 0.25,size=4)+scale_color_discrete(name="")

p2

library(patchwork)
(plotunadjustedci | p2) + plot_annotation(tag_levels = 'a')


# pdf(here::here("results", "figures_cabg_manuscript", "revision_figs_9_6_24.pdf"))
# MAP figure by CPB 

# plots where egfr, lvef, predrenf vs. predmort are categorized vs continuous 
covar_data = readr::read_csv(here::here("data/analytic/covariates_post_exclusion.csv"))

# sample sizes 

covar_data %>% count(val_egfr < 60)
# 621/(1406 + 621)

covar_data %>% count(val_hdef <= 40)
# 281/(281 + 1746)

covar_data %>% count(bin_iabp == 1 | bin_emergent == 1) 
# 139/(139+1888)

# predrenf vs. predmort 
temp_df = 
  overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
  mutate(type = "adjusted_MAP") %>% 
  bind_rows(predrenf %>% filter(type == "adjusted_MAP_predrenf")) %>% 
  mutate(type = factor(type, levels = c("adjusted_MAP", "adjusted_MAP_predrenf"),
                       labels = c("Predmort", "Predrenf")))


plot_j_fancy(df = temp_df,
             new_levels = map_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.2),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
  mutate(type = "adjusted_CVP") %>% 
  bind_rows(predrenf %>% filter(type == "adjusted_CVP_predrenf")) %>% 
  mutate(type = factor(type, levels = c("adjusted_CVP", "adjusted_CVP_predrenf"),
                       labels = c("Predmort", "Predrenf")))


plot_j_fancy(df = temp_df,
             new_levels = cvp_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.8),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))


temp_df = 
  res %>% 
  filter(type == "adjusted_zone5") %>% 
  bind_rows(predrenf %>% filter(type == "adjusted_shell5_predrenf")) %>% 
  mutate(type = factor(type, levels = c("adjusted_zone5", "adjusted_shell5_predrenf"),
                       labels = c("Predmort", "Predrenf")))
temp = temp_df %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p.val"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group))) %>% 
  rename(est = estimate)

plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, type) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est,
  )) + facet_wrap(.~type, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =1, linetype ="dashed")  +
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\n Adj OR=", format(round(est, 2),nsmall=2), "; ", p.value)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")



# j plot 

## cat lvef vs continuous 

temp_df = 
  overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
  mutate(type = "adjusted_MAP") %>% 
  bind_rows(catlvef %>% filter(type == "adjusted_MAP_catlvef")) %>% 
  mutate(type = factor(type, levels = c("adjusted_MAP", "adjusted_MAP_catlvef"),
                       labels = c("Continuous LVEF", "Categorical LVEF (<=40)")))


plot_j_fancy(df = temp_df,
             new_levels = map_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.2),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
  mutate(type = "adjusted_CVP") %>% 
  bind_rows(catlvef %>% filter(type == "adjusted_CVP_catlvef")) %>% 
  mutate(type = factor(type, levels = c("adjusted_CVP", "adjusted_CVP_catlvef"),
                       labels = c("Continuous LVEF", "Categorical LVEF (<=40)")))


plot_j_fancy(df = temp_df,
             new_levels = cvp_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.8),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  res %>% 
  filter(type == "adjusted_zone5") %>% 
  bind_rows(catlvef %>% filter(type == "adjusted_shell5_catlvef")) %>% 
  mutate(type = factor(type, levels = c("adjusted_zone5", "adjusted_shell5_catlvef"),
                       labels = c("Continuous LVEF", "Categorical LVEF (<=40)")))
temp = temp_df %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p.val"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group))) %>% 
  rename(est = estimate)

plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, type) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est,
  )) + facet_wrap(.~type, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =1, linetype ="dashed")  +
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\n Adj OR=", format(round(est, 2),nsmall=2), "; ", p.value)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")

# cat egfr vs. creatlst

temp_df = 
  overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
  mutate(type = "adjusted_MAP") %>% 
  bind_rows(categfr %>% filter(type == "adjusted_MAP_categfr")) %>% 
  mutate(type = factor(type, levels = c("adjusted_MAP", "adjusted_MAP_categfr"),
                       labels = c("Continuous Creatinine", "Categorical eGFR (<60)")))


plot_j_fancy(df = temp_df,
             new_levels = map_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.2),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
  mutate(type = "adjusted_CVP") %>% 
  bind_rows(categfr %>% filter(type == "adjusted_CVP_categfr")) %>% 
  mutate(type = factor(type, levels = c("adjusted_CVP", "adjusted_CVP_categfr"),
                       labels = c("Continuous Creatinine", "Categorical eGFR (<60)")))

plot_j_fancy(df = temp_df,
             new_levels = cvp_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.8),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  res %>% 
  filter(type == "adjusted_zone5") %>% 
  bind_rows(categfr %>% filter(type == "adjusted_shell5_categfr")) %>% 
  mutate(type = factor(type, levels = c("adjusted_zone5", "adjusted_shell5_categfr"),
                       labels = c("Continuous Creatinine", "Categorical eGFR (<60)")))
temp = temp_df %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p.val"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group))) %>% 
  rename(est = estimate)

plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, type) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est,
  )) + facet_wrap(.~type, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =1, linetype ="dashed")  +
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\n Adj OR=", format(round(est, 2),nsmall=2), "; ", p.value)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")

## cont egfr vs categorized egfr

temp_df = 
  contegfr %>% filter(type == "adjusted_MAP_contegfr") %>% 
  bind_rows(categfr %>% filter(type == "adjusted_MAP_categfr")) %>% 
  mutate(type = factor(type, levels = c("adjusted_MAP_contegfr", "adjusted_MAP_categfr"),
                       labels = c("Continuous eGFR", "Categorical eGFR (<60)")))


plot_j_fancy(df = temp_df,
             new_levels = map_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.2),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  contegfr %>% filter(type == "adjusted_CVP_contegfr") %>% 
  mutate(type = "adjusted_CVP") %>% 
  bind_rows(categfr %>% filter(type == "adjusted_CVP_categfr")) %>% 
  mutate(type = factor(type, levels = c("adjusted_CVP", "adjusted_CVP_categfr"),
                       labels = c("Continuous eGFR", "Categorical eGFR (<60)")))

plot_j_fancy(df = temp_df,
             new_levels = cvp_levels,
             n_comparisons = 2,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.8),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1a, col2a),
             lims = c(min(temp_df$lb), max(temp_df$ub)))
temp_df = 
  contegfr %>% 
  filter(type == "adjusted_shell5_contegfr") %>% 
  bind_rows(categfr %>% filter(type == "adjusted_shell5_categfr")) %>% 
  mutate(type = factor(type, levels = c("adjusted_shell5_contegfr", "adjusted_shell5_categfr"),
                       labels = c("Continuous eGFR", "Categorical eGFR (<60)")))
temp = temp_df %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p.val"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group))) %>% 
  rename(est = estimate)

plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, type) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est,
  )) + facet_wrap(.~type, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =1, linetype ="dashed")  +
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\n Adj OR=", format(round(est, 2),nsmall=2), "; ", p.value)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")


temp_df = 
  res %>% 
  filter(grepl("adjusted_MAP", type)) %>% 
  bind_rows(overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "adjusted_MAP")) %>% 
  mutate(type = factor(type, levels = c("adjusted_MAP",
                                        "adjusted_MAP_precpb",
                                        "adjusted_MAP_intracpb",
                                        "adjusted_MAP_postcpb"),
                       labels = c("Entire Surgery", "Pre CPB", "Intra CPB", "Post CPB")))
temp_df = 
  temp_df %>% 
  ungroup() %>%
  mutate(estimate = case_when(
    (term %in% c("MAP_(100,105]", "MAP_(105,110]", "MAP_(110,115]") & 
       type == "Intra CPB") ~ NA_real_,
    (term %in% c("MAP_(110,115]") & type == "Post CPB") ~ NA_real_,
    TRUE ~ estimate),
    ub = case_when(
      (term %in% c("MAP_(100,105]", "MAP_(105,110]", "MAP_(110,115]") & 
         type == "Intra CPB") ~ NA_real_,
      (term %in% c("MAP_(110,115]") & type == "Post CPB") ~ NA_real_,
      TRUE ~ ub),
    lb = case_when(
      (term %in% c("MAP_(100,105]", "MAP_(105,110]", "MAP_(110,115]") & 
         type == "Intra CPB") ~ NA_real_,
      (term %in% c("MAP_(110,115]") & type == "Post CPB") ~ NA_real_,
      TRUE ~ lb))


plot_j_fancy(df = temp_df,
             new_levels = map_levels,
             n_comparisons = 4,
             xlabs = map_labels,
             legend_pos =  c(0.15,0.2),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1, col2, col3, col4),
             lims = c(0.3, 1.7))
# lims = c(min(temp_df$lb), max(temp_df$ub)))

# CVP by CPB 
temp_df = 
  res %>% 
  filter(grepl("adjusted_CVP", type)) %>% 
  bind_rows(overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "adjusted_CVP")) %>% 
  mutate(type = factor(type, levels = c("adjusted_CVP",
                                        "adjusted_CVP_precpb",
                                        "adjusted_CVP_intracpb",
                                        "adjusted_CVP_postcpb"),
                       labels = c("Entire Surgery", "Pre CPB", "Intra CPB", "Post CPB")))


plot_j_fancy(df = temp_df,
             new_levels = cvp_levels,
             n_comparisons = 4,
             xlabs = cvp_labels,
             dodge = 0.5,
             legend_pos =  c(0.22,0.81),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             col_vector = c(col1, col2, col3, col4),
             # lims = c(0, 3.489))
             lims = c(min(temp_df$lb), max(temp_df$ub)))


temp_df = 
  res %>% 
  filter(grepl("zone5", type))

temp = temp_df %>% 
  mutate(type = case_when(type == "adjusted_zone5" ~ "adjusted_zone5_entire",
                          type == "univariate_zone5" ~ "univariate_zone5_entire",
                          TRUE ~ type),
         phase = sub(".*\\_zone5\\_", "", type),
         model = sub("\\_zone5\\_.*", "", type)) %>% 
  select(term, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 4), rep("(75,85]", 4), rep("(55,65]", 4), rep("(95,105]", 4),
                        rep("[45,55]", 4)),
            cut_CVP = c(rep("(2,4]",4), rep("(2,4]",4), rep("(4,6]",4), rep("(14,16]",4),rep("(14,16]", 4))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p_"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(phase = factor(phase, levels = c("entire", "precpb", "intracpb", "postcpb"),
                        labels = c("Entire Surgery", "Pre-CPB", "Intra-CPB", "Post-CPB"))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group)))

# pdf(here::here("results", "figures_cabg_manuscript", "shell_plot_opts.pdf"))


plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  mutate(type = case_when(type == "adjusted_zone5" ~ "adjusted_zone5_entire",
                          type == "univariate_zone5" ~ "univariate_zone5_entire",
                          TRUE ~ type),
         phase = sub(".*\\_zone5\\_", "", type),
         model = sub("\\_zone5\\_.*", "", type)) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  mutate(phase = factor(phase, levels = c("entire", "precpb", "intracpb", "postcpb"),
                        labels = c("Entire Surgery", "Pre-CPB", "Intra-CPB", "Post-CPB"))) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est_adjusted
  )) + facet_wrap(.~phase, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =1, linetype ="dashed") + 
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\nUnadj OR=", format(round(est_univariate, 2),nsmall=2), "; ", p_univariate,
                                             "\n Adj OR=", format(round(est_adjusted, 2),nsmall=2), "; ", p_adjusted)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")




temp_df = 
  mediators %>% 
  filter(grepl("adjusted_MAP_nomed", type)) %>% 
  bind_rows(overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "adjusted_MAP_med")) %>% 
  bind_rows(overall_map %>% filter(type == "univariate" & outcome == "AKI48") %>% 
              mutate(type = "unadjusted")) %>% 
  bind_rows(overall_map %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "cma_MAP_med",
                     ub = ub_cma, lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("unadjusted", 
                                  "adjusted_MAP_nomed",
                                  "adjusted_MAP_med",
                                  "cma_MAP_med"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))
# col1 = "#F8766D"; col2 = "#00BA38"; col3 = "#619CFF"; col4 = "purple4"

plot_j_fancy(df = temp_df,
             new_levels = map_levels,
             n_comparisons = 4,
             legend_pos =  c(0.42,0.81),
             xlabs = map_labels,
             col_vector = c(col1, col2, col3, col4),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(temp_df$lb), max(temp_df$ub)))

#0.7282376 1.3616545
temp_df = 
  mediators %>% 
  filter(grepl("adjusted_CVP_nomed", type)) %>% 
  bind_rows(overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "adjusted_CVP_med")) %>% 
  bind_rows(overall_cvp %>% filter(type == "univariate" & outcome == "AKI48") %>% 
              mutate(type = "unadjusted")) %>% 
  bind_rows(overall_cvp %>% filter(type == "adjusted" & outcome == "AKI48") %>% 
              mutate(type = "cma_CVP_med",
                     ub = ub_cma, lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("unadjusted", 
                                  "adjusted_CVP_nomed",
                                  "adjusted_CVP_med",
                                  "cma_CVP_med"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))

plot_j_fancy(df = temp_df,
             new_levels = cvp_levels,
             n_comparisons = 4,
             legend_pos =  c(0.32,0.81),
             xlabs = cvp_labels,
             dodge = 0.5,
             col_vector = c(col1, col2, col3, col4),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(0.7282376, 1.3616545)) # to match MAP 

# telescoping regressions 
temp_df = 
  telescoping_map %>% 
  filter(grepl("10", type)) %>% 
  bind_rows(telescoping_map %>% 
              filter(type == "adjusted_10") %>% 
              mutate(type = "cma_10",
                     ub = ub_cma,
                     lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("univariate_10", 
                                  "adjusted_10_nomed",
                                  "adjusted_10",
                                  "cma_10"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))
# col1 = "#F8766D"; col2 = "#00BA38"; col3 = "#619CFF"; col4 = "purple4"
map_levels10 = unique(temp_df$term)
map_labels10 = sub(".*\\_", "", map_levels10)

plot_j_fancy(df = temp_df,
             new_levels = map_levels10,
             n_comparisons = 4,
             legend_pos =  c(0.42,0.81),
             xlabs = map_labels10,
             col_vector = c(col1, col2, col3, col4),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  telescoping_map %>% 
  filter(grepl("15", type)) %>% 
  bind_rows(telescoping_map %>% 
              filter(type == "adjusted_15") %>% 
              mutate(type = "cma_15",
                     ub = ub_cma,
                     lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("univariate_15", 
                                  "adjusted_15_nomed",
                                  "adjusted_15",
                                  "cma_15"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))
map_levels15 = unique(temp_df$term)
map_labels15 = sub(".*\\_", "", map_levels15)

plot_j_fancy(df = temp_df,
             new_levels = map_levels15,
             n_comparisons = 4,
             legend_pos =  c(0.32,0.21),
             xlabs = map_labels15,
             col_vector = c(col1, col2, col3, col4),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(temp_df$lb), max(temp_df$ub)))
temp_df = 
  telescoping_map %>% 
  filter(grepl("20", type)) %>% 
  bind_rows(telescoping_map %>% 
              filter(type == "adjusted_20") %>% 
              mutate(type = "cma_20",
                     ub = ub_cma,
                     lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("univariate_20", 
                                  "adjusted_20_nomed",
                                  "adjusted_20",
                                  "cma_20"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))
map_levels20 = unique(temp_df$term)
map_labels20 = sub(".*\\_", "", map_levels20)

plot_j_fancy(df = temp_df,
             new_levels = map_levels20,
             n_comparisons = 4,
             legend_pos =  c(0.32,0.21),
             xlabs = map_labels20,
             col_vector = c(col1, col2, col3, col4),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  telescoping_cvp %>% 
  filter(grepl("4", type)) %>% 
  bind_rows(telescoping_cvp %>% 
              filter(type == "adjusted_4") %>% 
              mutate(type = "cma_4",
                     ub = ub_cma,
                     lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("univariate_4", 
                                  "adjusted_4_nomed",
                                  "adjusted_4",
                                  "cma_4"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))
cvp_levels4 = unique(temp_df$term)
cvp_labels4 = sub(".*\\_", "", cvp_levels4)

plot_j_fancy(df = temp_df,
             new_levels = cvp_levels4,
             n_comparisons = 4,
             legend_pos =  c(0.32,0.8),
             xlabs = cvp_labels4,
             col_vector = c(col1, col2, col3, col4),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(.94, 1.08))

temp_df = 
  telescoping_cvp %>% 
  filter(grepl("8", type)) %>% 
  bind_rows(telescoping_cvp %>% 
              filter(type == "adjusted_8") %>% 
              mutate(type = "cma_8",
                     ub = ub_cma,
                     lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("univariate_8", 
                                  "adjusted_8_nomed",
                                  "adjusted_8",
                                  "cma_8"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))
cvp_levels8 = unique(temp_df$term)
cvp_labels8 = sub(".*\\_", "", cvp_levels8)

plot_j_fancy(df = temp_df, 
             new_levels = cvp_levels8,
             n_comparisons = 8,
             legend_pos =  c(0.32,0.8),
             xlabs = cvp_labels8,
             col_vector = c(col1, col2, col3, col4),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(temp_df$lb), max(temp_df$ub)))

temp_df = 
  telescoping_cvp %>% 
  filter(grepl("6", type)) %>% 
  bind_rows(telescoping_cvp %>% 
              filter(type == "adjusted_6") %>% 
              mutate(type = "cma_6",
                     ub = ub_cma,
                     lb = lb_cma)) %>% 
  mutate(type = factor(type,
                       levels = c("univariate_6", 
                                  "adjusted_6_nomed",
                                  "adjusted_6",
                                  "cma_6"),
                       labels = c("Univariate",
                                  "Adjusted, no intra-op variables",
                                  "Adjusted, full model",
                                  "Adjusted, full model, with CMA")))
cvp_levels6 = unique(temp_df$term)
cvp_labels6 = sub(".*\\_", "", cvp_levels6)

plot_j_fancy(df = temp_df,
             new_levels = cvp_levels6,
             n_comparisons = 6,
             legend_pos =  c(0.32,0.8),
             xlabs = cvp_labels6,
             col_vector = c(col1, col2, col3, col4),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(temp_df$lb), max(temp_df$ub)))

## HERE 

temp_df = res %>% filter(type == "adjusted_zone5") %>% 
  mutate(type = "adjusted_shell5_med") %>% 
  bind_rows(res %>% filter(type == "univariate_zone5") %>% 
              mutate(type = "unadjusted_shell5_med"))  %>% 
  mutate(type = factor(type, levels = c("unadjusted_shell5_med",  
                                        "adjusted_shell5_med"),
                       labels = c("Univariate", "Adjusted, no mediators",
                                  "Adjusted, full model")))

temp = temp_df %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(term, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p_"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(phase = factor(phase, levels = c("lvefleq40", "lvefg40"),
                        labels = c("LVEF <= 40", "LVEF > 40"))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group)))
plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  mutate(phase = factor(phase, levels = c("lvefleq40", "lvefg40"),
                        labels = c("LVEF <= 40", "LVEF > 40"))) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est_adjusted
  )) + facet_wrap(.~phase, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    # palette = function(x) c("#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF", "#FFE1E1"), 
    # breaks=c(0.6,0.7,0.8,0.9,0.97,1.02,1/0.9),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =.9, linetype = "dashed") + 
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\nUnadj OR=", format(round(est_univariate, 2),nsmall=2), "; ", p_univariate,
                                             "\n Adj OR=", format(round(est_adjusted, 2),nsmall=2), "; ", p_adjusted)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")

temp_df = 
  mediators %>% 
  filter(type == "adjusted_shell5_nomed")  %>% 
  bind_rows(res %>% filter(type == "adjusted_zone5") %>% 
              mutate(type = "adjusted_shell5_med")) %>% 
  bind_rows(res %>% filter(type == "univariate_zone5") %>% 
              mutate(type = "unadjusted_shell5_med"))  %>% mutate(type = factor(type, levels = c("unadjusted_shell5_med",
                                                                                                 "adjusted_shell5_nomed",
                                                                                                 "adjusted_shell5_med"),
                                                                                labels = c("Univariate", "Adjusted, no mediators",
                                                                                           "Adjusted, full model")))


# temp = temp_df %>% 
#   mutate(phase = sub(".*\\_shell5\\_", "", type),
#          model = sub("\\_shell5\\_.*", "", type)) %>% 
#   select(term, est = estimate, p = p.value, phase, model) %>% 
#   pivot_wider(names_from = model, values_from = c(p, est)) %>% 
#   arrange(term) %>% 
#   bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
#                         rep("[45,55]", 2)),
#             cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
#   rename(group = term) %>% 
#   mutate(across(contains("p_"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
#   mutate(phase = factor(phase, levels = c("egfrl60", "egfrgeq60"),
#                         labels = c("eGFR < 60", "eGFR >= 60"))) %>% 
#   mutate(zone = paste0("Zone ", sub(".*group\\_", "", group)))
p = plotdf_5shell %>%
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
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF", "#FFE1E1"), 
    breaks=c(0.7,0.8,0.9,0.97,1.02,1/0.9),
    guide = "colorsteps",
    limits = c(0.7, 1.1),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        legend.position =c(.8, .15)) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =.9)



# lvef 

plot_j_fancy(df = lvef %>% filter(grepl("MAP", term) & !grepl("CVP", term) & grepl("adjusted", type)) %>% 
               mutate(type = factor(type, levels = c("adjusted_MAP_lvefleq40",
                                                     "adjusted_MAP_lvefg40"), 
                                    labels = c("LVEF <= 40", "LVEF > 40"))),
             new_levels = map_levels,
             legend_pos =  c(0.32,0.81),
             n_comparisons = 2,
             xlabs = map_labels,
             dodge = 0.5,
             col_vector = c(col1a, col2a),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(lvef$lb), max(lvef$ub)))
plot_j_fancy(df = lvef %>% filter(grepl("CVP", term) & grepl("adjusted", type)) %>% 
               mutate(type = factor(type, levels = c("adjusted_CVP_lvefleq40",
                                                     "adjusted_CVP_lvefg40"), 
                                    labels = c("LVEF <= 40", "LVEF > 40"))),
             new_levels = cvp_levels,
             legend_pos =  c(0.32,0.81),
             n_comparisons = 2,
             xlabs = cvp_labels,
             dodge = 0.5,
             col_vector = c(col1a, col2a),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(0.8, 1.4))


temp_df = 
  lvef %>% 
  filter(grepl("shell5", type))

temp = temp_df %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(term, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p_"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(phase = factor(phase, levels = c("lvefleq40", "lvefg40"),
                        labels = c("LVEF <= 40", "LVEF > 40"))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group)))
plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  mutate(phase = factor(phase, levels = c("lvefleq40", "lvefg40"),
                        labels = c("LVEF <= 40", "LVEF > 40"))) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est_adjusted
  )) + facet_wrap(.~phase, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    # palette = function(x) c("#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF", "#FFE1E1"), 
    # breaks=c(0.6,0.7,0.8,0.9,0.97,1.02,1/0.9),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =.9, linetype = "dashed") + 
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\nUnadj OR=", format(round(est_univariate, 2),nsmall=2), "; ", p_univariate,
                                             "\n Adj OR=", format(round(est_adjusted, 2),nsmall=2), "; ", p_adjusted)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")




plot_j_fancy(df = egfr %>% filter(grepl("MAP", term) & !grepl("CVP", term) & grepl("adjusted", type)) %>% 
               mutate(type = factor(type, levels = c("adjusted_MAP_egfrl60",
                                                     "adjusted_MAP_egfrgeq60"), 
                                    labels = c("eGFR < 60", "eGFR >= 60"))),
             new_levels = map_levels,
             legend_pos =  c(0.32,0.81),
             n_comparisons = 2,
             xlabs = map_labels,
             dodge = 0.5,
             col_vector = c(col1a, col2a),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(min(egfr$lb), max(egfr$ub)))
plot_j_fancy(df = egfr %>% filter(grepl("CVP", term) & grepl("adjusted", type)) %>% 
               mutate(type = factor(type, levels = c("adjusted_CVP_egfrl60",
                                                     "adjusted_CVP_egfrgeq60"), 
                                    labels = c("eGFR < 60", "eGFR >= 60"))),
             new_levels = cvp_levels,
             legend_pos =  c(0.32,0.81),
             n_comparisons = 2,
             xlabs = cvp_labels,
             dodge = 0.5,
             col_vector = c(col1a, col2a),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(0.8, 1.2))

temp_df = 
  egfr %>% 
  filter(grepl("shell5", type))

temp = temp_df %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(term, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p_"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(phase = factor(phase, levels = c("egfrl60", "egfrgeq60"),
                        labels = c("eGFR < 60", "eGFR >= 60"))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group)))
plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  mutate(phase = factor(phase, levels = c("egfrl60", "egfrgeq60"),
                        labels = c("eGFR < 60", "eGFR >= 60"))) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est_adjusted
  )) + facet_wrap(.~phase, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    # palette = function(x) c("#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF", "#FFE1E1"), 
    # breaks=c(0.6,0.7,0.8,0.9,0.97,1.02,1/0.9),
    guide = "colorsteps",
    limits = c(0.6, 1.1),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        # panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =.9, linetype = "dashed") + 
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\nUnadj OR=", format(round(est_univariate, 2),nsmall=2), "; ", p_univariate,
                                             "\n Adj OR=", format(round(est_adjusted, 2),nsmall=2), "; ", p_adjusted)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")




plot_j_fancy(df = shock %>% filter(grepl("MAP", term) & !grepl("CVP", term) &
                                     type %in% c("univariate_MAP_shockyes",
                                                 "univariate_MAP_shockno")) %>% 
               mutate(type = factor(type, levels = c("univariate_MAP_shockyes",
                                                     "univariate_MAP_shockno"), 
                                    labels = c("Shock", "No shock"))),
             new_levels = map_levels,
             legend_pos =  c(0.32,0.81),
             n_comparisons = 2,
             xlabs = map_labels,
             dodge = 0.5,
             col_vector = c(col1a, col2a),
             xname = "MAP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(0.5, 2.2))
plot_j_fancy(df =shock %>% filter(!grepl("MAP", term) & grepl("CVP", term) &
                                    type %in% c("univariate_CVP_shockyes",
                                                "univariate_CVP_shockno")) %>% 
               mutate(type = factor(type, levels = c("univariate_CVP_shockyes",
                                                     "univariate_CVP_shockno"),
                                    labels = c("Shock", "No shock"))),
             new_levels = cvp_levels,
             legend_pos =  c(0.32,0.81),
             n_comparisons = 2,
             xlabs = cvp_labels,
             dodge = 0.5,
             col_vector = c(col1a, col2a),
             xname = "CVP (mmHg) Range", yname = "Odds Ratio (per 5 minutes in range)",
             title = "",
             lims = c(0.8, 1.5))

temp_df = 
  shock %>% 
  filter(grepl("shell5", type) & type != "adjusted_shell5_shockyes")

temp = temp_df %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(term, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p_"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(phase = factor(phase, levels = c("shockyes", "shockno"),
                        labels = c("Shock", "No Shock"))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group)),
         est = if_else(is.na(est_adjusted), est_univariate, est_adjusted)) %>% 
  mutate(lab = if_else(is.na(p_adjusted),
                       paste0(zone, "\nUnadj OR=", format(round(est_univariate, 2),nsmall=2), "; ", p_univariate),
                       paste0(zone, "\nUnadj OR=", format(round(est_univariate, 2),nsmall=2), "; ", p_univariate,
                              "\n Adj OR=", format(round(est_adjusted, 2),nsmall=2), "; ", p_adjusted)))
plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  mutate(phase = factor(phase, levels = c("shockyes", "shockno"),
                        labels = c("Shock", "No Shock")),
         est = if_else(is.na(est_adjusted), est_univariate, est_adjusted)) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est
  )) + facet_wrap(.~phase, nrow = 2)+
  geom_tile(col = "darkgrey") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    # palette = function(x) c("#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF", "#FFE1E1"), 
    # breaks=c(0.6,0.7,0.8,0.9,0.97,1.02,1/0.9),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =.9, linetype = "dashed") + 
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = lab),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")

temp_df = 
  shock %>% 
  filter(grepl("shell5", type) & !grepl("adjusted", type))

temp = temp_df %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(term, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  arrange(term) %>% 
  bind_cols(cut_MAP = c(rep("(95,105]", 2), rep("(75,85]", 2), rep("(55,65]", 2), rep("(95,105]", 2),
                        rep("[45,55]", 2)),
            cut_CVP = c(rep("(2,4]",2), rep("(2,4]",2), rep("(4,6]",2), rep("(14,16]",2),rep("(14,16]", 2))) %>% 
  rename(group = term) %>% 
  mutate(across(contains("p_"), ~ifelse(.x < 0.001, "p<0.001", paste0("p=", round(.x, 3))))) %>% 
  mutate(phase = factor(phase, levels = c("shockyes", "shockno"),
                        labels = c("Shock", "No Shock"))) %>% 
  mutate(zone = paste0("Zone ", sub(".*group\\_", "", group)))

plotdf_5shell %>%
  left_join(temp_df, by = c("group" = "term")) %>% 
  mutate(phase = sub(".*\\_shell5\\_", "", type),
         model = sub("\\_shell5\\_.*", "", type)) %>% 
  select(cut_MAP, cut_CVP, group, mean, est = estimate, p = p.value, phase, model) %>% 
  pivot_wider(names_from = model, values_from = c(p, est)) %>% 
  mutate(phase = factor(phase, levels = c("shockyes", "shockno"),
                        labels = c("Shock", "No Shock"))) %>% 
  ggplot(aes(
    x = cut_MAP,
    y = cut_CVP,
    fill = est_univariate
  )) + facet_wrap(.~phase, nrow = 2)+
  geom_tile(col = "black") +
  labs(
    x = "MAP (mmHg)",
    y = "CVP (mmHg)",
    fill = "Odds Ratio (per 5 minutes in range)",
    title = ""
  ) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) c("#5D82F0", "#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF",
                            "#FFC3C3", "#FF8888" , "#FF6B6B", "#FF4D4D"), 
    breaks=c(0.6, 0.67,0.75, 0.82, 0.89, 0.97, 1.02, 1.04, 1.06, 1.08, 1.12),
    # palette = function(x) c("#93ACF5" ,"#AEC0F7", "#C9D5FA", "#E4EAFC", "#FFFFFF", "#FFE1E1"), 
    # breaks=c(0.6,0.7,0.8,0.9,0.97,1.02,1/0.9),
    guide = "colorsteps",
    limits = c(0.6, 1.12),
    show.limits = TRUE)+
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'),
        legend.title.align = 1) + 
  geom_segment(data = dat, aes(x = xs, xend = xends, y = ys, yend= yends),
               col = "black", inherit.aes = F, linewidth =.9) + 
  geom_label(data = temp, aes(x = cut_MAP, y = cut_CVP,
                              label = paste0(zone, "\nUnadj OR=", format(round(est_univariate, 2),nsmall=2), "; ", p_univariate)),
             col = "black", hjust = 0.25,size=3)+scale_color_discrete(name="")




dev.off()

## supplement bar charts
map_data = readRDS(here::here("results", "data_cabg_manuscript", "map_data.rds"))

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
  theme_bw()+
  labs(x = "MAP (mmHg) Range", y = "")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values=c(col1, col2),labels= c("Mean minutes in range", "Percent individuals with at least 5 minutes in range"),name="")+
  theme(text=element_text(size=17),
        # axis.title.x = element_text(margin=margin(t=8.5,r=0,b=0,l=0)),
        legend.text = element_text(size=18),
        axis.text.y =element_text(size=18),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size=18))+
  geom_bar(position = position_dodge(), stat = "identity",colour="black")+
  geom_errorbar(aes(ymin = value, ymax = ub), position = position_dodge(),
                size = .5)

cvp_data = readRDS(here::here("results", "data_cabg_manuscript", "cvp_data.rds"))


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
  geom_errorbar(aes(ymin = value, ymax = ub), position = position_dodge(width=0.9),
                size = .1,width=0.5)+
  theme_bw()+
  labs(x = "CVP (mmHg) Range", y = "")+
  scale_fill_manual(values=c(col1, col2),labels= c("Mean minutes in range", "Percent individuals with at least 5 minutes in range"),name="")+
  theme(text=element_text(size=17),
        legend.position = "bottom",
        # axis.title.x = element_text(margin=margin(t=8.5,r=0,b=0,l=0)),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size = 12),
        axis.text.y =element_text(size=18),
        axis.title.x = element_text(size=18))+
  geom_bar(position = position_dodge(), stat = "identity",colour="black")+
  geom_errorbar(aes(ymin = value, ymax = ub), position = position_dodge(),
                size = .5)
### telescoping models 

map_telescoping = rmap_regressions_telescoping = 
  readRDS(here::here("results/data_cabg_manuscript/map_regressions_telescoping.rds"))


hemo_data = readRDS(here::here("data/analytic/hemo_timeseries_interp_post_exclusion.rds"))


# source(here::here("code/CABG/analysis/utilities.R"))
get_ranges = function(hemo_data, thresholds, hemo_variable) {
  column_name = paste0("val_", hemo_variable) # get column of hemodynamic vars from hemo dataset
  
  hemo_data %>%
    filter(cat_anes == "intra") %>%
    group_by(id) %>%
    select(id, time, all_of(column_name)) %>%
    mutate(across(
      all_of(column_name),
      ~ cut(.x,
            breaks = thresholds,
            include.lowest = TRUE)
    )) %>%
    count(across(all_of(column_name)), .drop = FALSE) %>%
    pivot_wider(
      values_from = n,
      id_cols = id,
      names_from = all_of(column_name)
    ) %>%
    ungroup() %>%
    rename(missing = 'NA') %>%
    rename_with(~ str_c(hemo_variable, "_", .),-id)
}
map_thresholds = seq(45, 115, 5)

cvp_thresholds = seq(0,20,2)
covar_data = readr::read_csv(here::here("data/analytic/covariates_post_exclusion.csv"))

# create time in range regression predictors 
map_data_intra = get_ranges(hemo_data = hemo_data %>% 
                              filter(cat_cpb == "intra"),
                            thresholds = map_thresholds,
                            hemo_variable = "MAP") %>% 
  select(-contains("missing"))

map_data_pre = get_ranges(hemo_data = hemo_data %>% 
                            filter(cat_cpb == "pre"),
                          thresholds = map_thresholds,
                          hemo_variable = "MAP") %>% 
  select(-contains("missing"))

map_data_post = get_ranges(hemo_data = hemo_data %>% 
                             filter(cat_cpb == "post"),
                           thresholds = map_thresholds,
                           hemo_variable = "MAP") %>% 
  select(-contains("missing"))

shock_id = covar_data %>% 
  filter(bin_iabp == 1 | bin_emergent == 1) %>% 
  pull(id)
map_data_shock = get_ranges(hemo_data = hemo_data  %>% filter(id %in% shock_id),
                            thresholds =  map_thresholds,
                            hemo_variable = "MAP")
cvp_data_shock = get_ranges(hemo_data = hemo_data  %>% filter(id %in% shock_id),
                            thresholds =  cvp_thresholds,
                            hemo_variable = "CVP")

cvp_data_intra = get_ranges(hemo_data = hemo_data %>% 
                              filter(cat_cpb == "intra"),
                            thresholds = cvp_thresholds,
                            hemo_variable = "CVP") %>% 
  select(-contains("missing"))

cvp_data_pre = get_ranges(hemo_data = hemo_data %>% 
                            filter(cat_cpb == "pre"),
                          thresholds = cvp_thresholds,
                          hemo_variable = "CVP") %>% 
  select(-contains("missing"))

cvp_data_post = get_ranges(hemo_data = hemo_data %>% 
                             filter(cat_cpb == "post"),
                           thresholds = cvp_thresholds,
                           hemo_variable = "CVP") %>% 
  select(-contains("missing"))

map_data_intra %>% 
  rename(val_115 = contains("115")) %>% 
  select(val_115) %>% 
  summarize(mean = mean(val_115),
            unique = length(unique(val_115)))
unique(map_data_intra$`MAP_(110,115]`)
# library(tidyverse)

props = 
  map_data_intra %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("MAP"),
                   ~sum(.x > 5)))
# map 100-105, 105-110, 110-115 

props = 
  map_data_pre %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("MAP"),
                   ~sum(.x > 5)))

props = 
  map_data_post %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("MAP"),
                   ~sum(.x > 5)))
# 110-115 post 

## all fine 
props = 
  cvp_data_intra %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("CVP"),
                   ~sum(.x > 5)))


props = 
  cvp_data_pre %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("CVP"),
                   ~sum(.x > 5)))

props = 
  cvp_data_post %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("CVP"),
                   ~sum(.x > 5)))

res %>% 
  filter(grepl("adjusted_MAP_intracpb", type))

props = 
  cvp_data_shock %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("CVP"),
                   ~sum(.x > 5)))
props %>% t()

props = 
  map_data_shock %>% 
  left_join(covar_data %>% select(id, bin_aki48h)) %>% 
  group_by(bin_aki48h) %>% 
  summarize(across(contains("MAP"),
                   ~sum(.x > 5)))
props %>% t()
shock %>% filter(!grepl("adjusted", type)) %>% print(n=Inf)
