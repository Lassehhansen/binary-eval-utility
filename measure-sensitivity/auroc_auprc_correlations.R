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

facet_labels_a_attribute <- c(`0.05` = "A - 5% of Population", `0.15` = "A - 15% of Population", `0.25` = "A - 25% of Population")

best_auroc_models_vis <- filtered_results_vis %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auroc_total)) %>%
  slice_head(n = 1000) %>%
  mutate(model_choice = "AUROC",
         model_id = 1:1000)

auroc_summary  = glm(auroc_total ~ auprc_total*prevelance_bins + ev_pos_b + ev_neg_a + ev_neg_b +  p_pos_space_a + p_pos_space_b, 
                     data = best_auroc_models_vis,
                     quasibinomial(link = "logit"))

auroc_int <- ggeffects::ggpredict(auroc_summary, terms = c("auprc_total", "prevelance_bins"))
plot_auroc_int = plot(auroc_int)

plot_auroc_int = ggplot(auroc_int, aes(x, predicted, color = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey", alpha = 0.3) +
  labs(x = "AUROC - Minority (A)",
       y = "Pred. AUROC Total",
       fill = "",
       color = "AUROC 
Majority (B):",
       title = "") +
  scale_color_manual(values = c("#9b2226", "#bc6c25", "#669bbc")) +
  facet_grid(. ~ facet) +
  theme_light() +
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right",
        panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5),
        legend.key.width = unit(1, "cm")) 






breaks_prevelance <- quantile(best_auroc_models_vis$overall_prevalence, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Create the age groups
best_auroc_models_vis$prevelance_bins <- cut(best_auroc_models_vis$overall_prevalence,
                                             breaks = breaks_prevelance,
                                             labels = c("5-20%", "20-35%", "35-50%"),
                                             include.lowest = TRUE)

best_auroc_models_vis_long = best_auroc_models_vis %>% pivot_longer(cols = c(auroc_a, auroc_b))


p_auroc_subgroup = ggplot(best_auroc_models_vis_long, aes(x = value, y = auroc_total, color = name)) +
  # Line for auroc_a
  #geom_point() +
  geom_smooth() +
  # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  scale_color_manual(values = c("#9b2226", "#bc6c25")) + 
  facet_wrap(attribute_ratio ~ p_pos_space_a, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute)) +
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  labs(title = "AUROC Total ~ AUROC Subgroup",
       y = "AUROC Total",
       x = "AUROC Subgroup")


best_auprc_models_vis <- filtered_results_vis %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auprc_total)) %>%
  slice_head(n = 1000) %>%
  mutate(model_choice = "AUPRC",
         subgroup = "a",
         model_id = 1:1000)

### Models

auprc_model  = lm(auprc_total ~  auprc_a * auprc_b * prevelance_bins, data = best_auprc_models_vis)

auprc_model_glm <- glm(
  auprc_total ~ auprc_a * auprc_b * prevelance_bins,
  quasibinomial(link = "logit"),
  data = best_auprc_models_vis
)

auprc_model_glm <- glm(
  auprc_total ~ auprc_a * auprc_b,
  quasibinomial(link = "logit"),
  data = best_auprc_models_vis %>% filter(overall_prevalence < 0.20,
  )
)

test = best_auprc_models_vis %>% filter(overall_prevalence < 0.20,
                                        auprc_b > 0.85
)

ggplot(test, aes(x = auprc_a, y = auprc_total)) + geom_smooth()

auprc_int <- ggeffects::ggpredict(auprc_model_glm, terms = c("auprc_a[all]", "auprc_b"))
plot_auprc_int = plot(auprc_int)

plotplot_auprc_int = ggplot(auprc_int, aes(x, predicted, color = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey", alpha = 0.3) +
  labs(x = "AUPRC - Minority (A)",
       y = "Pred. AUPRC Total",
       fill = "",
       color = "AUPRC 
Majority (B):",
       title = "") +
  scale_color_manual(values = c("#9b2226", "#bc6c25", "#669bbc")) +
  facet_grid(. ~ facet) +
  #, labeller = labeller(facet = facet_labels_subgroup)) + 
  theme_light() +
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right",
        panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5),
        legend.key.width = unit(1, "cm")) 

breaks_prevelance <- quantile(best_auprc_models_vis$overall_prevalence, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Create the age groups
best_auprc_models_vis$prevelance_bins <- cut(best_auprc_models_vis$overall_prevalence,
                                             breaks = breaks_prevelance,
                                             labels = c("5-20%", "20-35%", "35-50%"),
                                             include.lowest = TRUE)


best_auprc_models_vis_long = best_auprc_models_vis %>% pivot_longer(cols = c(auprc_a, auprc_b))

p_auprc_subgroup = ggplot(best_auprc_models_vis_long, aes(x = auprc_total, y = value, color = name)) +
  # Line for auroc_a
  #geom_point() +
  geom_smooth() +
  # Line for auroc_b
  #geom_smooth(aes(x = auroc_total, y = auroc_b), color = "#bc6c25") +
  # Optional: To improve readability
  theme_minimal()  + 
  scale_color_manual(values = c("#9b2226", "#bc6c25")) + 
  facet_wrap(attribute_ratio ~ p_pos_space_a, 
             labeller = labeller(
               attribute_ratio = facet_labels_a_attribute)) +
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black")
  ) +
  labs(title = "AUPRC Total ~ AUPRC Subgroup",
       y = "AUPRC Total",
       x = "AUPRC Subgroup")

png(filename = "figures_simulation/subgroup_cor.png", width = 3000, height = 2000, units = "px", res = 300)

ggpubr::ggarrange(plot_auroc_int, plot_auprc_int, nrow = 2
)

dev.off()

ggpubr::ggarrange(plot_auroc_int, plot_auprc_int, nrow = 2
)

ggsave("figures_simulation/subgroup_cor.eps",  width = 2500, height = 2000, units = "px", dpi = 300)




