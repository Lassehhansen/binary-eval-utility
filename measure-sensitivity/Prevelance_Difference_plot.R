library(tidyverse)

results_simulated_data_extended <- read_csv("simulation_data/results_simulated_data_v6.csv")



calculate_overall_prevalence <- function(prevalence_a, prevalence_b, attribute_ratio) {
  R <- attribute_ratio
  overall_prevalence <- (prevalence_a * R / (R + 1)) + (prevalence_b * 1 / (R + 1))
  return(overall_prevalence)
}


filtered_results_vis <- results_simulated_data_extended %>%
  mutate(
    overall_prevalence = calculate_overall_prevalence(p_pos_space_a, p_pos_space_b, attribute_ratio)
  ) %>% 
  filter(auroc_total < 0.95, auprc_total < 0.95) %>% 
  filter(auroc_total > 0.5, auprc_total > overall_prevalence)

breaks_prevelance <- quantile(filtered_results_vis$overall_prevalence, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Create the age groups
filtered_results_vis$prevelance_bins <- cut(filtered_results_vis$overall_prevalence,
                                       breaks = breaks_prevelance,
                                       labels = c("5-20% Prevelance", "20-35% Prevelance", "35-50% Prevelance"),
                                       include.lowest = TRUE)

filtered_results_vis$prev_diff = filtered_results_vis$p_pos_space_a - filtered_results_vis$p_pos_space_b
filtered_results_vis$prev_diff_ratio = filtered_results_vis$p_pos_space_a/filtered_results_vis$p_pos_space_b
filtered_results_vis$prev_diff_ratio = round(filtered_results_vis$prev_diff_ratio, digits = 2)


filtered_results_vis$prev_diff_factor = as.factor(filtered_results_vis$prev_diff)
filtered_results_vis$prev_diff_ratio_factor = as.factor(filtered_results_vis$prev_diff_ratio)

p_auroc_1 = ggplot(filtered_results_vis, aes(x = auprc_total, y = auroc_total)) +
  # Line for auroc_a
  geom_point() +
  geom_smooth() +
  #geom_violin(fill = "#669bbc") +
  # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  facet_wrap(. ~ attribute_ratio, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute)) +
  labs(title = "AUROC Total Over Different Prevelance Differences",
       subtitle = "Prevelance Diff: Prevelance B - Prevelance A",
       x = "AUPRC",
       y = "AUROC") 

facet_labels_a_attribute <- c(`0.05` = "A - 5% of Population", `0.15` = "A - 15% of Population", `0.25` = "A - 25% of Population")

best_auroc_models_vis <- filtered_results_vis %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auroc_total)) %>%
  slice_head(n = 1000) %>%
  mutate(model_choice = "AUROC",
         model_id = 1:1000)

auroc_summary  = lm(auroc_total ~ ev_pos_a + ev_pos_b + ev_neg_a + ev_neg_b + attribute_ratio +  p_pos_space_a + p_pos_space_b, best_auroc_models_a)
summary(auroc_summary)                                                                                                                                                                             

best_auroc_models_vis$overall_prevalence

p_auroc_1 = ggplot(best_auroc_models_vis, aes(x = prev_diff_factor, y = auroc_total)) +
  # Line for auroc_a
  #geom_point() +
  #geom_smooth(color = "#EE7674") +
  geom_violin(fill = "#669bbc") +
  # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  facet_wrap(. ~ attribute_ratio, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute)) +
  labs(title = "AUROC Total Over Different Prevelance Differences",
       subtitle = "Prevelance Diff: Prevelance B - Prevelance A",
       x = "Prevelance Difference",
       y = "AUROC") 

p_auroc_2 = ggplot(best_auroc_models_vis, aes(x = prev_diff_factor, y = auroc_total)) +
  # Line for auroc_a
  #geom_point() +
  #geom_smooth(color = "#EE7674") +
  geom_violin(fill = "#669bbc") +
  # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  facet_wrap(prevelance_bins ~ attribute_ratio, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute)) +
  labs(title = "AUROC Total Over Different Prevelance Differences",
       subtitle = "Prevelance Diff: Prevelance B - Prevelance A",
       x = "Prevelance Difference",
       y = "AUROC") 
#+ 
  #scale_x_discrete(labels = c("" = "Minority (A)", "b" = "Majority (B)")) #+ 

p_auroc_3 = ggplot(best_auroc_models_vis, aes(x = prev_diff_factor, y = auroc_total)) +
  # Line for auroc_a
  #geom_point() +
  #geom_smooth(color = "#EE7674") +
  geom_violin(fill = "#669bbc") +
  # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  facet_wrap(prevelance_bins ~ ., 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute),
             scales = "free_x") +
  labs(title = "AUROC Total Over Different Prevelance Differences",
       subtitle = "Prevelance Diff: Prevelance B - Prevelance A",
       x = "Prevelance Difference",
       y = "AUROC") +
  scale_y_continuous(limits = c(0.60, 0.95), expand = c(0,0))
#+ 

best_auprc_models_vis <- filtered_results_vis %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auprc_total)) %>%
  slice_head(n = 1000) %>%
  mutate(model_choice = "AUPRC",
         subgroup = "a",
         model_id = 1:1000)

auprc_summary  = lm(auroc_total ~ ev_pos_a + ev_pos_b + ev_neg_a + ev_neg_b + attribute_ratio +  p_pos_space_a + p_pos_space_b, best_auprc_models_a)
summary(auprc_summary)     

p_auprc_1 = ggplot(best_auprc_models_vis, aes(x = prev_diff_factor, y = auprc_total)) +
  # Line for auroc_a
  geom_violin(fill = "#9b2226") +
  geom_smooth() +
    # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  facet_wrap(. ~ attribute_ratio, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute)) +
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
        ) +
  labs(title = "AUPRC Total Over Different Prevelance Differences",
       subtitle = "Prevelance Diff: Prevelance B - Prevelance A",
       x = "Prevelance Difference",
       y = "AUPRC")

p_auprc_2 = ggplot(best_auprc_models_vis, aes(x = prev_diff_factor, y = auprc_total)) +
  # Line for auroc_a
  geom_violin(fill = "#9b2226") +
  geom_smooth() +
  # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  facet_wrap(prevelance_bins ~ attribute_ratio, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute)) +
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  labs(title = "AUPRC Total Over Different Prevelance Differences",
       subtitle = "Prevelance Diff: Prevelance B - Prevelance A",
       x = "Prevelance Difference",
       y = "AUPRC")

p_auprc_3 = ggplot(best_auprc_models_vis, aes(x = prev_diff_factor, y = auprc_total)) +
  geom_violin(fill = "#9b2226") +
  theme_minimal()  + 
  facet_wrap(. ~ prevelance_bins, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute),
             scales = "free_x") +
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  labs(title = "AUPRC Total Over Different Prevelance Differences",
       subtitle = "Prevelance Diff: Prevelance B - Prevelance A",
       x = "Prevelance Difference",
       y = "AUPRC") +
  scale_y_continuous(limits = c(0.60, 0.95), expand = c(0,0))



p_tot1 = ggpubr::ggarrange(p_auroc_1, p_auprc_1,
                  nrow = 2,
                  labels = c("A", "B"))

png(filename = "figures_simulation/disttribution_prev_diff.png",  width = 2200, height = 1700, units = "px", res = 300)

p_tot1

dev.off()

p_tot1

ggsave("figures_simulation/disttribution_prev_diff.eps",   width = 2200, height = 1700, units = "px", dpi = 300)


p_tot2 = ggpubr::ggarrange(p_auroc_2, p_auprc_2,
                  nrow = 2,
                  labels = c("A", "B"))

p_auprc_2

p_tot2

p_tot3 = ggpubr::ggarrange(p_auprc_3, p_auroc_3,
                           nrow = 2,
                           labels = c("A", "B"))
p_tot3

png(filename = "figures_simulation/disttribution_prev_bins.png", width = 2200, height = 1700, units = "px", res = 300)

p_tot3

dev.off()

p_tot3

ggsave("figures_simulation/disttribution_prev_bins.eps",   width = 2200, height = 1700, units = "px", dpi = 300)


auroc_summary <- best_auroc_models_vis %>%
  group_by(prevelance_bins, prev_diff_factor) %>%
  summarise(
    median_auroc = median(auroc_total),
    mean_auroc = mean(auroc_total),
    sd_auroc = sd(auroc_total),
    iqr_auroc = IQR(auroc_total),
    .groups = 'drop'
  )

auprc_summary <- best_auprc_models_vis %>%
  group_by(prevelance_bins, prev_diff_factor) %>%
  summarise(
    median_auroc = median(auprc_total),
    mean_auroc = mean(auprc_total),
    sd_auroc = sd(auprc_total),
    iqr_auroc = IQR(auprc_total),
    .groups = 'drop'
  )
### Models

auprc_model  = lm(auprc_total ~  prev_diff_factor * attribute_ratio, data = best_auprc_models_vis)
auroc_model  = lm(auroc_total ~  prev_diff_factor * attribute_ratio, data = best_auroc_models_vis)
  
auprc_int <- ggeffects::ggpredict(auprc_model, terms = c("prev_diff_factor", "attribute_ratio"))
auroc_int <- ggeffects::ggpredict(auroc_model, terms = c("prev_diff_factor", "attribute_ratio"))

auprc_int_p = ggplot(auprc_int, aes(x, predicted, colour = group)) + 
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 1) +
  labs(x = "",
       y = "Pred. AUPRC",
       fill = "",
       color = "") +
  scale_color_manual(values = c("#669bbc", "#9b2226", "#bc6c25"),
                     labels = c("A - 5% of pop.",
                                "A - 15% of pop.",
                                "A - 25% of pop.")) +
  theme_light() +
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top",
        panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5),
        legend.key.width = unit(1, "cm")) +
  scale_x_discrete(labels = c("a" = "Minority (A)", "b" = "Majority (B)"))

auroc_int_p = ggplot(auroc_int, aes(x, predicted, colour = group)) + 
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 1) +
  labs(x = "",
       y = "Pred. auroc",
       fill = "",
       color = "") +
  scale_color_manual(values = c("#669bbc", "#9b2226", "#bc6c25"),
                     labels = c("A - 5% of pop.",
                                "A - 15% of pop.",
                                "A - 25% of pop.")) +
  theme_light() +
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top",
        panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5),
        legend.key.width = unit(1, "cm")) +
  scale_x_discrete(labels = c("a" = "Minority (A)", "b" = "Majority (B)"))
