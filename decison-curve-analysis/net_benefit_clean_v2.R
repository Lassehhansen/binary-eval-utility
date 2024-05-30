set.seed(1997)

### Filtering AUROC and AUPRC values with relevance
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pracma) 

results_df = read_csv("simulation_data/results_simulated_data_v6.csv")
# Filter out models that don't meet the criteria

calculate_total_prevalence <- function(prevalence_a, prevalence_b, attribute_ratio) {
  # Calculate the population proportions for subgroup A and B
  population_a = attribute_ratio
  population_b = 1 - attribute_ratio
  
  # Calculate the overall disease prevalence
  total_prevalence = (population_a * prevalence_a) + (population_b * prevalence_b)
  
  return(round(total_prevalence, digits = 2))
}

filtered_results_df <- results_df %>%
  mutate(
    total_prevelance = calculate_total_prevalence(p_pos_space_a, p_pos_space_b, attribute_ratio)
  ) %>% 
  filter(auroc_total < 0.95, auprc_total < 0.95) %>% 
  filter(auroc_total > 0.5, auprc_total > total_prevelance)

net_benefit_data <- function(alpha_pos, beta_pos, alpha_neg, beta_neg, prevalence, total_population) {
  if (alpha_pos <= 0 || beta_pos <= 0 || alpha_neg <= 0 || beta_neg <= 0) {
    stop("Alpha and beta parameters must be positive")
  }
  
  # Create a sequence of threshold values from 0 to tau_max
  tau <- seq(0.0001, 0.9999, length.out = 100)
  
  # Calculate TPR (Sensitivity) and FPR using the CDF (pbeta)
  tpr_values <- sapply(tau, function(t) 1 - pbeta(t, alpha_pos, beta_pos))
  fpr_values <- sapply(tau, function(t) 1 - pbeta(t, alpha_neg, beta_neg))
  
  # Calculate expected counts for TP, FP, TN, FN
  TP_counts <- tpr_values * prevalence * total_population
  FP_counts <- fpr_values * (1 - prevalence) * total_population
  TN_counts <- (1 - fpr_values) * (1 - prevalence) * total_population
  FN_counts <- (1 - tpr_values) * prevalence * total_population
  
  # Calculate Net Benefit, Relative Harm, and Relative Benefit for each threshold
  net_benefits <- sapply(seq_along(tau), function(i) {
    pt <- tau[i]
    w <- pt / (1 - pt)  # Weight calculation
    # Net Benefit calculation
    TP_counts[i] / total_population - FP_counts[i] / total_population * w
  })
  
  relative_harm <- sapply(seq_along(tau), function(i) {
    pt <- tau[i]
    w <- pt / (1 - pt)  # Weight calculation
    # Relative Harm calculation
    w * FP_counts[i] / total_population + FN_counts[i] / total_population
  })
  
  relative_benefit <- sapply(seq_along(tau), function(i) {
    # Relative Benefit calculation
    TP_counts[i] / total_population + TN_counts[i] / total_population
  })
  
  # Create a data frame for plotting or further analysis
  net_benefit_plot_data <- data.frame(
    threshold = tau, 
    net_benefit = net_benefits,
    relative_harm = relative_harm,
    relative_benefit = relative_benefit
  )
  
  return(net_benefit_plot_data)
}

### Best models in each a/b prevelance and attribute_ratio

best_auroc_models_a <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auroc_total)) %>%
  slice_head(n = 20) %>%
  mutate(model_choice = "AUROC",
         subgroup = "a",
         model_id = 1:20)

best_auroc_models_b <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auroc_total)) %>%
  slice_head(n = 20) %>%
  mutate(model_choice = "AUROC",
         subgroup = "b",
         model_id = 1:20)

best_auprc_models_a <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auprc_total)) %>%
  slice_head(n = 20) %>%
  mutate(model_choice = "AUPRC",
         subgroup = "a",
         model_id = 1:20)


best_auprc_models_b <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auprc_total)) %>%
  slice_head(n = 20) %>%
  mutate(model_choice = "AUPRC",
         subgroup = "b",
         model_id = 1:20)

best_auroc_auprc_models = rbind(best_auroc_models_a, best_auroc_models_b, best_auprc_models_a, best_auprc_models_b)

best_test = rbind(best_auroc_models_a, best_auprc_models_a)

best_test$prev_diff = best_test$p_pos_space_a - best_test$p_pos_space_b


best_test %>% group_by(model_choice) %>% summarize(ev_pos_a = mean(ev_pos_a),
                                                   ev_pos_b = mean(ev_pos_b),
                                                   ev_neg_a = mean(ev_neg_a),
                                                   ev_neg_b = mean(ev_neg_b),
                                                   confidence_pos_a = mean(confidence_pos_a),
                                                   confidence_pos_b = mean(confidence_pos_b),
                                                   confidence_neg_a = mean(confidence_neg_a),
                                                   confidence_neg_b = mean(confidence_neg_b),
                                                   p_pos_space_a = mean(p_pos_space_a),
                                                   p_pos_space_b = mean(p_pos_space_b),
                                                   attribute_ratio = mean(attribute_ratio)) %>% view()


# Function to apply cost_threshold_data and return a data frame with additional model information

get_net_benefit <- function(row) {
  # Select the right subgroup based on the 'Subgroup' column
  subgroup_suffix <- ifelse(row$subgroup == "a", "_a", "_b")
  alpha_pos <- row[[paste0("alpha_pos", subgroup_suffix)]]
  beta_pos <- row[[paste0("beta_pos", subgroup_suffix)]]
  alpha_neg <- row[[paste0("alpha_neg", subgroup_suffix)]]
  beta_neg <- row[[paste0("beta_neg", subgroup_suffix)]]
  p_pos_space <- row[[paste0("p_pos_space", subgroup_suffix)]]
  
  # Calculate the expected costs
  expected_costs_df <- net_benefit_data(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, total_population = 100000)
  
  # Add additional columns
  expected_costs_df$model_choice <- row$model_choice
  expected_costs_df$model_id <- row$model_id
  
  expected_costs_df$subgroup <- row$subgroup
  expected_costs_df$p_pos_space_a <- row$p_pos_space_a
  expected_costs_df$p_pos_space_b <- row$p_pos_space_b
  
  expected_costs_df$attribute_ratio <- row$attribute_ratio
  expected_costs_df$alpha_pos_a <- row$alpha_pos_a
  expected_costs_df$alpha_pos_b <- row$alpha_pos_b
  expected_costs_df$alpha_neg_a <- row$alpha_neg_a
  expected_costs_df$alpha_neg_b <- row$alpha_neg_b
  expected_costs_df$beta_pos_a <- row$beta_pos_a
  expected_costs_df$beta_pos_b <- row$beta_pos_b
  expected_costs_df$beta_neg_a <- row$beta_neg_a
  expected_costs_df$beta_neg_b <- row$beta_neg_b
  expected_costs_df$ev_pos_a <-  row$ev_pos_a
  expected_costs_df$ev_pos_b <-  row$ev_pos_b
  expected_costs_df$ev_neg_a <-  row$ev_neg_a
  expected_costs_df$ev_neg_b <-  row$ev_neg_b
  
  expected_costs_df$auroc_total <-  row$auroc_total
  expected_costs_df$auprc_total <-  row$auprc_total
  expected_costs_df$auroc_a <-  row$auroc_a
  expected_costs_df$auroc_b <-  row$auroc_b
  expected_costs_df$auprc_a <-  row$auprc_a
  expected_costs_df$auprc_b <-  row$auprc_b
  
  return(expected_costs_df)
}

# Initialize an empty data frame to store all the results
all_net_benefit <- data.frame()

# Loop over each row of the cost grid and calculate expected costs
for (i in 1:nrow(best_auroc_auprc_models)) {
  row_costs <- get_net_benefit(best_auroc_auprc_models[i, ])
  all_net_benefit <- rbind(all_net_benefit, row_costs)
}

all_net_benefit <- all_net_benefit %>%
  mutate(
    total_prevelance = calculate_total_prevalence(p_pos_space_a, p_pos_space_b, attribute_ratio)
  ) 

breaks_prevelance <- quantile(all_net_benefit$total_prevelance, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Create the age groups
all_net_benefit$prevelance_bins <- cut(all_net_benefit$total_prevelance,
                                       breaks = breaks_prevelance,
                                       labels = c("5-20%", "20-35%", "35-50%"),
                                       include.lowest = TRUE)


all_net_benefit$prev_diff = all_net_benefit$p_pos_space_a - all_net_benefit$p_pos_space_b

all_net_benefit %>% select(p_pos_space_a, p_pos_space_b, prev_diff) %>% view()

all_net_benefit$prev_diff_factor = as.factor(all_net_benefit$prev_diff)

# Create the age groups

all_net_benefit <- all_net_benefit %>% 
  mutate(
    prev_diff_bins =
      ifelse(prev_diff == 0, "Equal Prevelance",
             ifelse(prev_diff < 0, "Prev B > Prev A",
                    ifelse(prev_diff > 0, "Prev A > Prev B", NA))))

all_net_benefit$prev_diff_bins = factor(all_net_benefit$prev_diff_bins ,
                                        levels = c("Prev B > Prev A",
                                                   "Equal Prevelance",
                                                   "Prev A > Prev B"))

test = all_net_benefit %>% 
  filter(threshold < 0.25) %>% 
  group_by(p_pos_space_a, p_pos_space_b, 
           attribute_ratio, model_choice, 
           subgroup) %>% 
  summarize(
    m_net_benefit = mean(net_benefit),
    sd_net_benefit = sd(net_benefit),
    m_relative_harm = mean(relative_harm),
    sd_relative_harm = sd(relative_harm)
  )

test2 = test %>% pivot_wider(., id_cols = c("p_pos_space_a", "p_pos_space_b", 
                                            "attribute_ratio"),
                             names_from = c("model_choice", "subgroup"),
                             values_from = c("m_net_benefit", "sd_net_benefit",
                                             "m_relative_harm", "sd_relative_harm"))

test2$auc_diff_a_nb = test2$m_net_benefit_AUROC_a - m_net_benefit_AUPRC_a
test2$auc_diff_b_nb = test2$m_net_benefit_AUROC_b - m_net_benefit_AUPRC_b
test2$auc_diff_a_rh = test2$m_relative_harm_AUROC_a - m_relative_harm_AUPRC_a
test2$auc_diff_b_rh = test2$m_relative_harm_AUROC_b - m_relative_harm_AUPRC_b

test3 = all_net_benefit %>% 
  filter(threshold < 0.25) %>% 
  group_by(prevelance_bins,
           model_choice, subgroup) %>% 
  summarize(
    m_net_benefit = mean(net_benefit),
    sd_net_benefit = sd(net_benefit),
    m_relative_harm = mean(relative_harm),
    sd_relative_harm = sd(relative_harm)
  )

test4 = test3 %>% pivot_wider(., id_cols = c("prevelance_bins"),
                              names_from = c("model_choice", "subgroup"),
                              values_from = c("m_net_benefit", "sd_net_benefit",
                                              "m_relative_harm", "sd_relative_harm"))

test4$auc_diff_a_nb = test4$m_net_benefit_AUROC_a - test4$m_net_benefit_AUPRC_a
test4$auc_diff_b_nb = test4$m_net_benefit_AUROC_b - test4$m_net_benefit_AUPRC_b
test4$auc_diff_a_rh = test4$m_relative_harm_AUROC_a - test4$m_relative_harm_AUPRC_a
test4$auc_diff_b_rh = test4$m_relative_harm_AUROC_b - test4$m_relative_harm_AUPRC_b


### adding prev diff

all_net_benefit$prev_diff_ratio = all_net_benefit$p_pos_space_a/all_net_benefit$p_pos_space_b
all_net_benefit$prev_diff_ratio = round(all_net_benefit$prev_diff_ratio, digits = 2)
all_net_benefit$prev_diff_ratio_factor = as.factor(all_net_benefit$prev_diff_ratio)

### For which threshold does AUPRC and AUROC have higher net benefits?

interact_nb_threshold  = lm(net_benefit ~  p_pos_space_a + p_pos_space_b +  model_choice * threshold  * subgroup, data = all_net_benefit %>% filter(threshold < 0.9))
interact_harm_threshold  = lm(relative_harm ~  p_pos_space_a + p_pos_space_b + model_choice * threshold  * subgroup, data = all_net_benefit %>% filter(threshold < 0.9))
interact_ben_threshold  = lm(relative_benefit ~  p_pos_space_a + p_pos_space_b + model_choice * threshold  * subgroup, data = all_net_benefit %>% filter(threshold < 0.9))

interact_nb_threshold  = lm(net_benefit ~ model_choice * threshold  * subgroup*prev_diff_factor, data = all_net_benefit %>% filter(threshold < 0.9))


glm(auroc_total ~ auprc_total*prevelance_bins + ev_pos_b + ev_neg_a + ev_neg_b +  p_pos_space_a + p_pos_space_b, 
    data = best_auroc_models_vis,
    quasibinomial(link = "logit"))

interact_nb_threshold  = lm(net_benefit ~  p_pos_space_a + p_pos_space_b +  model_choice * threshold  * subgroup, 
                             data = all_net_benefit %>% filter(threshold < 0.9))

nb_threshold <- ggeffects::ggpredict(interact_nb_threshold, terms = c("threshold[all]", "model_choice", "subgroup", "prev_diff_factor"))
harm_threshold <- ggeffects::ggpredict(interact_harm_threshold, terms = c("threshold", "model_choice", "subgroup"))
ben_threshold <- ggeffects::ggpredict(interact_ben_threshold, terms = c("threshold", "model_choice", "subgroup"))

### Plot net benefit threshold

facet_labels_subgroup <- c(`a` = "Minority (A)", `b` = "Majority (B)")

p_nb_threshold = ggplot(nb_threshold, aes(x, predicted, color = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey", alpha = 0.3) +
  labs(x = "",
       y = "Pred. Net Benefit",
       fill = "",
       color = "",
       title = "") +
  scale_color_manual(values = c(AUPRC = "#9b2226", 
                                AUROC = "#669bbc")) +
  facet_grid(panel ~ facet, labeller = labeller(facet = facet_labels_subgroup)) + 
  theme_light() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        plot.margin = margin(t = 5, l = 5, r = 5, b = -10, unit = "pt"),
        legend.position = "right",
        #panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)) +
  scale_x_continuous(limits = c(0, 0.9),
                     expand = c(0,0),
                     breaks = c(0.04761905, 0.1666667, 0.33, 0.50, 0.65, 0.80),
                     labels = c("5% 
(1:20)","17%, 
(1:5)", "", "50% 
(1:1)", "", "80% 
(1:0.25)")) + scale_y_continuous(labels = scales::label_number(accuracy = 0.01))  # format to 2 decimal places


p_harm_threshold = ggplot(harm_threshold, aes(x, predicted, color = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey", alpha = 0.3) +
  labs(x = "Threshold",
       y = "Pred. Relative Harm",
       fill = "",
       color = "",
       title = "") +
  scale_color_manual(values = c(AUPRC = "#9b2226", 
                                AUROC = "#669bbc")) +
  facet_grid(. ~ facet, labeller = labeller(facet = facet_labels_subgroup)) + 
  theme_light() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "right",
        plot.margin = margin(t = -10, l = 5, r = 5, b = 2, unit = "pt"),
        #panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5)
        #legend.key.width = unit(1, "cm")
  ) +
  scale_x_continuous(limits = c(0, 0.9),
                     expand = c(0,0),
                     breaks = c(0.04761905, 0.1666667, 0.33, 0.50, 0.65, 0.80),
                     labels = c("5% 
(1:20)","17%, 
(1:5)", "", "50% 
(1:1)", "", "80% 
(1:0.25)")) + scale_y_continuous(labels = scales::label_number(accuracy = 0.01))  # format to 2 decimal places


int_p_threhold_tot = ggpubr::ggarrange(p_nb_threshold, p_harm_threshold, labels = c("A", "B"),
                                       common.legend = T, nrow = 2,
                                       align = "hv")

png(filename = "figures_simulation/utility/threshold_utility_nb.png", width = 2300, height = 1600, units = "px", res = 300)

int_p_threhold_tot

dev.off()

int_p_threhold_tot

ggsave("figures_simulation/utility/threshold_utility_nb.eps",  width = 2300, height = 1600, units = "px", dpi = 300)


### Model choice and subgroup prediction when thresh < 0.20


nb_model_int  = lm(net_benefit ~  p_pos_space_a + p_pos_space_b + model_choice  * subgroup, data = all_net_benefit %>% filter(threshold < 0.20))
harm_model_int  = lm(relative_harm ~  p_pos_space_a + p_pos_space_b + model_choice  * subgroup, data = all_net_benefit %>% filter(threshold < 0.20))
ben_model_int  = lm(relative_benefit ~  p_pos_space_a + p_pos_space_b + model_choice  * subgroup, data = all_net_benefit %>% filter(threshold < 0.20))

nb_int <- ggeffects::ggpredict(nb_model_int, terms = c("subgroup", "model_choice"))
harm_int <- ggeffects::ggpredict(harm_model_int, terms = c("subgroup", "model_choice"))
ben_int <- ggeffects::ggpredict(ben_model_int, terms = c("subgroup", "model_choice"))


# Convert point plots to bar plots for predicted net benefit (nb_int)
nb_int_p <- ggplot(nb_int, aes(x = x, y = predicted, fill = group)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0.8), 
                width = 0.2, 
                size = 0.5) +
  labs(x = "",
       y = "Pred. Net Benefit",
       fill = "",
       color = "") +
  scale_fill_manual(values = c(AUPRC = "#9b2226", 
                               AUROC = "#669bbc")) +
  facet_grid(. ~ x, labeller = labeller(x = facet_labels_subgroup)) + 
  theme_light() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "right",
        panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5),
        legend.key.width = unit(1, "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))  # format to 2 decimal places

# Convert point plots to bar plots for predicted relative harm (harm_int)
harm_int_p <- ggplot(harm_int, aes(x = x, y = predicted, fill = group)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0.8), 
                width = 0.2, 
                size = 0.5) +
  labs(x = "",
       y = "Pred. Relative Harm",
       fill = "",
       color = "") +
  scale_fill_manual(values = c(AUPRC = "#9b2226", 
                               AUROC = "#669bbc")) +
  facet_grid(. ~ x, labeller = labeller(x = facet_labels_subgroup)) + 
  theme_light() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "right",
        panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5),
        legend.key.width = unit(1, "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))  # format to 2 decimal places

int_p = ggpubr::ggarrange(nb_int_p, harm_int_p, labels = c("A", "B"),
                          common.legend = T, nrow = 2)

png(filename = "figures_simulation/utility/mod_sub_utility_nb.png", width = 2100, height = 1400, units = "px", res = 300)

int_p

dev.off()

int_p

ggsave("figures_simulation/utility/mod_sub_utility_nb.eps", width = 2100, height = 1400, units = "px", dpi = 300)




### Prevelance bins

int_nb_prev_bins  = lm(net_benefit ~  model_choice * prevelance_bins  * subgroup, data = all_net_benefit %>% filter(threshold < 0.20))
int_harm_prev_bins  = lm(relative_harm ~  model_choice * prevelance_bins  * subgroup, data = all_net_benefit %>% filter(threshold < 0.20))

nb_prev_bins <- ggeffects::ggpredict(int_nb_prev_bins, terms = c("prevelance_bins", "model_choice", "subgroup"))
harm_prev_bins <- ggeffects::ggpredict(int_harm_prev_bins, terms = c("prevelance_bins", "model_choice", "subgroup"))


nb_prev_bins_p <- ggplot(nb_prev_bins, aes(x = x, y = predicted, fill = group)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0.8), 
                width = 0.2, 
                size = 1) +
  labs(x = "",
       y = "Pred. Net Benefit",
       fill = "",
       color = "") +
  scale_fill_manual(values = c(AUPRC = "#9b2226", 
                               AUROC = "#669bbc")) +
  facet_grid(. ~ facet, labeller = labeller(facet = facet_labels_subgroup)) + 
  theme_light() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "right",
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))  # format to 2 decimal places


# Convert point plots to bar plots for predicted relative harm
harm_prev_bins_p <- ggplot(harm_prev_bins, aes(x = x, y = predicted, fill = group)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0.8), 
                width = 0.2, 
                size = 1) +
  labs(x = "",
       y = "Pred. Relative Harm",
       fill = "",
       color = "") +
  scale_fill_manual(values = c(AUPRC = "#9b2226", 
                               AUROC = "#669bbc")) +
  facet_grid(. ~ facet, labeller = labeller(facet = facet_labels_subgroup)) + 
  theme_light() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "right",
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))  # format to 2 decimal places

prev_bins_tot_p = ggpubr::ggarrange(nb_prev_bins_p, harm_prev_bins_p, labels = c("A", "B"),
                                    common.legend = T, nrow = 2)

alltot = ggpubr::ggarrange(p_nb_threshold, p_harm_threshold, nb_prev_bins_p, harm_prev_bins_p, 
                           common.legend = T, 
                           legend = "bottom",
                           ncol = 2, nrow = 2,
                           heights = c(1.1, 1),
                           labels = c("A", "B",
                                      "C", "D"))
png(filename = "figures_simulation/utility/prev_bins_utility_nb.png", width = 2100, height = 1500, units = "px", res = 300)

prev_bins_tot_p

dev.off()

prev_bins_tot_p

ggsave("figures_simulation/utility/prev_bins_utility_nb.eps", width = 2100, height = 1500, units = "px", dpi = 300)


