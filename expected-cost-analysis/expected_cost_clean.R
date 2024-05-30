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

### Defining different types of expected cost functions:

cost_threshold_data <- function(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, cost_fp, cost_fn, n) {
  tau <- round(seq(0, 1, length.out = 100), digits = 2)
  
  # Ensure all parameters are positive
  if (alpha_pos <= 0 || beta_pos <= 0 || alpha_neg <= 0 || beta_neg <= 0) {
    stop("Alpha and beta parameters must be positive")
  }
  
  tpr_values <- sapply(tau, function(t) 1 - pbeta(t, alpha_pos, beta_pos))
  fpr_values <- sapply(tau, function(t) 1 - pbeta(t, alpha_neg, beta_neg))
  
  expected_fps <- n * (1 - p_pos_space) * fpr_values
  expected_fns <- n * p_pos_space * (1 - tpr_values)
  
  expected_costs <- (cost_fp * expected_fps + cost_fn * expected_fns) / n
  
  # Find the minimum cost threshold
  min_cost_index <- which.min(expected_costs)
  min_cost_threshold <- tau[min_cost_index]
  min_cost_value <- expected_costs[min_cost_index]
  
  # Repeat the minimum cost threshold and value for each tau to match the length
  min_cost_thresholds <- rep(min_cost_threshold, 100)
  min_cost_values <- rep(min_cost_value, 100)
  
  # Create a data frame for plotting or further analysis
  cost_dat <- data.frame(
    threshold = tau, 
    expected_cost = expected_costs, 
    min_cost_threshold = min_cost_thresholds, 
    min_cost_value = min_cost_values
  )
  
  return(cost_dat)
}


### Best models in each a/b prevelance and attribute_ratio

best_auroc_models_a <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auroc_total)) %>%
  slice_head(n = 10) %>%
  mutate(model_choice = "AUROC",
         subgroup = "a",
         model_id = 1:10)

best_auroc_models_b <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auroc_total)) %>%
  slice_head(n = 10) %>%
  mutate(model_choice = "AUROC",
         subgroup = "b",
         model_id = 1:10)

best_auprc_models_a <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auprc_total)) %>%
  slice_head(n = 10) %>%
  mutate(model_choice = "AUPRC",
         subgroup = "a",
         model_id = 1:10)


best_auprc_models_b <- filtered_results_df %>%
  group_by(p_pos_space_a, p_pos_space_b, attribute_ratio) %>%
  arrange(desc(auprc_total)) %>%
  slice_head(n = 10) %>%
  mutate(model_choice = "AUPRC",
         subgroup = "b",
         model_id = 1:10)

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

### mergin with relative cost ratio

# Prepare a grid for cost_fn and cost_fp combinations
mesh_density <- 8
cost_fn_space <- 1

#tau <- seq(0.1, 0.9, by = 0.1)  # Thresholds

# Calculate the weights
#cost_fp_space <- tau / (1 - tau)
cost_fp_space <- c(0.04, 0.05, 0.10, 0.20, 1.00, 4.00)
#cost_fp_space <- round(pracma::logspace(-1, 0.5, mesh_density), digits = 2)

# Create an empty data frame to store the cost grid
cost_grid <- expand.grid(cost_fn = cost_fn_space, cost_fp = cost_fp_space)

# Left join the best model data onto the cost grid
cost_grid_with_models <- cost_grid %>%
  cross_join(best_auroc_auprc_models)


get_cost_thresh <- function(row, n) {
  # Ensure you are extracting and passing the right suffix based variables
  subgroup_suffix <- ifelse(row$subgroup == "a", "_a", "_b")
  alpha_pos <- row[[paste0("alpha_pos", subgroup_suffix)]]
  beta_pos <- row[[paste0("beta_pos", subgroup_suffix)]]
  alpha_neg <- row[[paste0("alpha_neg", subgroup_suffix)]]
  beta_neg <- row[[paste0("beta_neg", subgroup_suffix)]]
  p_pos_space <- row[[paste0("p_pos_space", subgroup_suffix)]]
  
  # Calculate the expected costs
  expected_costs_df <- cost_threshold_data(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, row$cost_fp, row$cost_fn, n)
  # Add additional model-specific columns
  expected_costs_df$model_choice <- row$model_choice
  expected_costs_df$subgroup <- row$subgroup
  expected_costs_df$cost_fp <- round(row$cost_fp, digits = 2)
  expected_costs_df$cost_fn <- round(row$cost_fn, digits = 2)
  
  # Add all relevant parameters to the dataframe
  for (col in c("p_pos_space_a", "p_pos_space_b", "attribute_ratio", "alpha_pos_a", "alpha_pos_b", "alpha_neg_a", "alpha_neg_b", "beta_pos_a", "beta_pos_b", "beta_neg_a", "beta_neg_b", "ev_pos_a", "ev_pos_b", "ev_neg_a", "ev_neg_b")) {
    expected_costs_df[[col]] <- row[[col]]
  }
  
  return(expected_costs_df)
}

# Correcting the Loop Call
all_expected_costs__new <- data.frame()

for (i in 1:nrow(cost_grid_with_models)) {
  row_costs <- get_cost_thresh(cost_grid_with_models[i, ], n = 100000)
  all_expected_costs__new <- rbind(all_expected_costs__new, row_costs)
}

all_expected_costs__new <- all_expected_costs__new %>%
  mutate(
    total_prevelance = calculate_total_prevalence(p_pos_space_a, p_pos_space_b, attribute_ratio)
  ) 

write_csv(all_expected_costs__new, "expected_costs_v6_new.csv")

all_expected_costs__new = read_csv("expected_costs_v6_new.csv")

breaks_prevelance <- quantile(all_expected_costs__new$total_prevelance, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Create the age groups
all_expected_costs__new$prevelance_bins <- cut(all_expected_costs__new$total_prevelance,
                                       breaks = breaks_prevelance,
                                       labels = c("5-20%", "20-35%", "35-50%"),
                                       include.lowest = TRUE)


all_expected_costs__new$prev_diff = all_expected_costs__new$p_pos_space_a - all_expected_costs__new$p_pos_space_b

all_expected_costs__new$prev_diff_factor = as.factor(all_expected_costs__new$prev_diff)



### testing difference

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
test2$auc_diff_a_nb = test2$m_net_benefit_AUROC_a - test2$m_net_benefit_AUPRC_a
test2$auc_diff_b_nb = test2$m_net_benefit_AUROC_b - test2$m_net_benefit_AUPRC_b
test2$auc_diff_a_rh = test2$m_relative_harm_AUROC_a - test2$m_relative_harm_AUPRC_a
test2$auc_diff_b_rh = test2$m_relative_harm_AUROC_b - test2$m_relative_harm_AUPRC_b
### Cost weight relationship

all_expected_costs__new$weight_relationship = 1/all_expected_costs__new$cost_fp
all_expected_costs__new$weight_relationship_round = round(all_expected_costs__new$weight_relationship, digits = 2)

facet_labels_cost <- c(`0.04` = "1:25", `0.05` = "1:20", `0.1` = "1:10",
                       `0.2` = "1:5", `1` = "1:1", `4` = "1:0.25"
)

cost_fp_space <- c(0.04, 0.05, 0.10, 0.20, 1.00, 4.00)


### Plot net benefit threshold

ep_thresh_mod  = glm(expected_cost ~ threshold * model_choice  * subgroup*cost_fp, data = all_expected_costs__new %>% 
                       filter(weight_relationship_round %in% c(20, 5, 1, 0.25)),
                      family=Gamma(link = "inverse"))


ep_thresh_int <- ggeffects::ggpredict(ep_thresh_mod, terms = c("threshold[all]", "model_choice", "subgroup", "cost_fp"))

facet_labels_subgroup <- c(`a` = "Minority (A)", `b` = "Majority (B)")
#facet_labels_cost <- c(`0.25` = "Cost FP: 0.25", `1` = "Cost FP: 1.00", `4` = "Cost FP: 4.00")


p_ep_fp_cost = ggplot(ep_thresh_int, aes(x, predicted, color = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey", alpha = 0.3) +
  labs(x = "Threshold",
       y = "Pred. Expected Cost",
       fill = "",
       color = "",
       title = "") +
  # scale_color_manual(values = c(AUPRC = "#669bbc", 
  #                               AUROC = "#9b2226")) +
  scale_color_manual(values = c(AUPRC = "#669bbc", 
                                AUROC = "#9b2226")) +
  facet_grid(panel ~ facet, labeller = labeller(facet = facet_labels_subgroup,
                                                panel = facet_labels_cost
                                                ),
             scales="free_y") + 
  theme_light() +
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        
        
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        
        strip.text.y  = element_text(angle = 0, size = 14, color = "black"),
        strip.text.x = element_text(angle = 0, size = 14, color = "black"),
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        panel.border = element_rect(fill = NA, color = "black", size = 1),  # Make sure the size is consistent
        legend.position = "none",
        
        #strip.text.x = element_text(size = 14, color = "black"),
        #strip.text.y = element_text(size = 14, color = "black"),
        
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.spacing.x = unit(1, "lines")
        #legend.key.width = unit(1, "cm")
        ) +
  scale_y_continuous(n.breaks = 4)

p_ep_fp_cost

png(filename = "figures_simulation/cost/cost_fp_utility_ec.png", width = 2300, height = 1600, units = "px", res = 300)

p_ep_fp_cost

dev.off()

p_ep_fp_cost

ggsave("figures_simulation/cost/cost_fp_utility_ec.eps",  width = 2300, height = 1600, units = "px", dpi = 300)


ep_prev_bin_mod  = glm(expected_cost ~ prevelance_bins * model_choice* subgroup, data = all_expected_costs__new %>% filter(cost_fp < 1,
                                                                                                                   threshold < 0.20),
                     family=Gamma(link = "inverse"))

ep_prev_bin_mod  = glm(expected_cost ~ prevelance_bins * model_choice* subgroup*cost_fp, data = all_expected_costs__new %>% filter(weight_relationship_round %in% c(20, 5, 1, 0.25),
                                                                                                                           threshold < 0.20),
                       family=Gamma(link = "inverse"))

ep_prev_bin_p_int <- ggeffects::ggpredict(ep_prev_bin_mod, terms = c("prevelance_bins", "model_choice", "subgroup", "cost_fp"))


facet_labels_subgroup <- c(`a` = "Minority (A)", `b` = "Majority (B)")


p_ep_prev_bin = ggplot(ep_prev_bin_p_int, aes(x, predicted, color = group)) + 
  geom_point() +
  #geom_bar(stat = "identity", color = "black") +
  #geom_bar() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 1) +
  labs(x = "Prevelance Bin",
       y = "Pred. Expected Cost",
       fill = "",
       color = "",
       title = "") +
  # scale_color_manual(values = c(AUPRC = "#669bbc", 
  #                               AUROC = "#9b2226")) +
  scale_color_manual(values = c(AUPRC = "#669bbc", 
                                AUROC = "#9b2226")) +
  #facet_grid(. ~ facet, labeller = labeller(facet = facet_labels_subgroup)) + 
  facet_grid(panel ~ facet, labeller = labeller(facet = facet_labels_subgroup,
                                                panel = facet_labels_cost),
             scales = "free") + 
  theme_minimal() +
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        strip.text.y  = element_text(angle = 0, size = 14, color = "black"),
        #strip.text = element_text(angle = 0, size = 11, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
        
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right",
        panel.spacing.x = unit(1, "lines"),
        #strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        plot.title = element_text(size = 16, color = "black", hjust = 0.5),
        legend.key.width = unit(1, "cm")) +
  scale_y_continuous(n.breaks = 3)


png(filename = "figures_simulation/cost/prev_bin_utility_ec.png", width = 1900, height = 1200, units = "px", res = 300)

p_ep_prev_bin

dev.off()

p_ep_prev_bin

ggsave("figures_simulation/cost/prev_bin_utility_ec.eps",  width = 1900, height = 1200, units = "px", dpi = 300)


### THIS IS THE ONE I AM USING ATM

ggpubr::ggarrange(p_ep_fp_cost, p_ep_prev_bin, common.legend = T,
                  legend = "bottom", labels = c("A", "B"))

ep_prev_diff_mod  = glm(expected_cost ~ prev_diff * model_choice* subgroup, data = all_expected_costs__new %>% filter(cost_fp < 1,
                                                                                                                           threshold < 0.20),
                       family=Gamma(link = "inverse"))

ep_prev_diff_p_int <- ggeffects::ggpredict(ep_prev_diff_mod, terms = c("prev_diff", "model_choice", "subgroup"))


facet_labels_subgroup <- c(`a` = "Minority (A)", `b` = "Majority (B)")


p_ep_prev_diff = ggplot(ep_prev_diff_p_int, aes(x, predicted, color = group)) + 
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 1) +
  labs(x = "Prevelance Bin",
       y = "Pred. Expected Cost",
       fill = "",
       color = "",
       title = "") +
  # scale_color_manual(values = c(AUPRC = "#669bbc", 
  #                               AUROC = "#9b2226")) +
  scale_color_manual(values = c(AUPRC = "#669bbc", 
                                AUROC = "#EE7674")) +
  facet_grid(. ~ facet, labeller = labeller(facet = facet_labels_subgroup)) + 
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


png(filename = "figures_simulation/cost/prev_diff_utility_ec.png", width = 1900, height = 1200, units = "px", res = 300)

p_ep_prev_diff

dev.off()

p_ep_prev_diff

ggsave("figures_simulation/cost/prev_diff_utility_ec.eps",  width = 1900, height = 1200, units = "px", dpi = 300)





