# Plots and Tables

setwd('C:/Users/timo/Documents/Master/UIBK/Master_Thesis/R/data')

library(tidyverse)
library(rfPermute)
library(randomForest)
library(ggplot2)
library(patchwork)
library(sf)
library(dplyr)

# Principle Component Analysis ####

# Read all PCA files
EP_PCA <-readRDS('rf_maps/PCA_EP.rds')
IF_PCA <- readRDS('rf_maps/PCA_IF.rds')
LM100_PCA <- readRDS('rf_maps/PCA_LM100.rds')
LM250_PCA <- readRDS('rf_maps/PCA_LM250.rds')
LM500_PCA <- readRDS('rf_maps/PCA_LM500.rds')#

## Tables  ####

### Summary (eigenvalues + variance) ####

pca_summary_EP <- data.frame(
  PC = paste0("PC", 1:length(EP_PCA$model$sdev)),
  Eigenvalue = round(EP_PCA$model$sdev^2, 3),
  Variance = round((EP_PCA$model$sdev^2/sum(EP_PCA$model$sdev^2))*100, 1),
  Cumulative = round(cumsum(EP_PCA$model$sdev^2/sum(EP_PCA$model$sdev^2))*100, 1)
)

write.csv(pca_summary_EP, "plots/PCA_summary_EP.csv", row.names = FALSE)

pca_summary_IF <- data.frame(
  PC = paste0("PC", 1:length(IF_PCA$model$sdev)),
  Eigenvalue = round(IF_PCA$model$sdev^2, 3),
  Variance = round((IF_PCA$model$sdev^2/sum(IF_PCA$model$sdev^2))*100, 1),
  Cumulative = round(cumsum(IF_PCA$model$sdev^2/sum(IF_PCA$model$sdev^2))*100, 1)
)

write.csv(pca_summary_IF, "plots/PCA_summary_IF.csv", row.names = FALSE)

pca_summary_LM100 <- data.frame(
  PC = paste0("PC", 1:length(LM100_PCA$model$sdev)),
  Eigenvalue = round(LM100_PCA$model$sdev^2, 3),
  Variance = round((LM100_PCA$model$sdev^2/sum(LM100_PCA$model$sdev^2))*100, 1),
  Cumulative = round(cumsum(LM100_PCA$model$sdev^2/sum(LM100_PCA$model$sdev^2))*100, 1)
)

write.csv(pca_summary_LM100, "plots/PCA_summary_LM100.csv", row.names = FALSE)

pca_summary_LM250 <- data.frame(
  PC = paste0("PC", 1:length(LM250_PCA$model$sdev)),
  Eigenvalue = round(LM250_PCA$model$sdev^2, 3),
  Variance = round((LM250_PCA$model$sdev^2/sum(LM250_PCA$model$sdev^2))*100, 1),
  Cumulative = round(cumsum(LM250_PCA$model$sdev^2/sum(LM250_PCA$model$sdev^2))*100, 1)
)

write.csv(pca_summary_LM250, "plots/PCA_summary_LM250.csv", row.names = FALSE)

pca_summary_LM500 <- data.frame(
  PC = paste0("PC", 1:length(LM500_PCA$model$sdev)),
  Eigenvalue = round(LM500_PCA$model$sdev^2, 3),
  Variance = round((LM500_PCA$model$sdev^2/sum(LM500_PCA$model$sdev^2))*100, 1),
  Cumulative = round(cumsum(LM500_PCA$model$sdev^2/sum(LM500_PCA$model$sdev^2))*100, 1)
)

write.csv(pca_summary_LM500, "plots/PCA_summary_LM500.csv", row.names = FALSE)

### Loadings ####

loadings_EP <- EP_PCA$model$loadings[,1:5, drop=FALSE]
varimax_loadings_EP <- varimax(loadings_EP)$loadings[1:nrow(loadings_EP), 1:5]

loadings_EP_df <- data.frame(
  Variable = rownames(varimax_loadings_EP),
  PC1 = round(varimax_loadings_EP[,1], 3),
  PC2 = round(varimax_loadings_EP[,2], 3), 
  PC3 = round(varimax_loadings_EP[,3], 3),
  PC4 = round(varimax_loadings_EP[,4], 3),
  PC5 = round(varimax_loadings_EP[,5], 3)
)

write.csv(loadings_EP_df, "plots/PCA_loadings_EP.csv", row.names = FALSE)

loadings_IF <- IF_PCA$model$loadings[,1:5, drop=FALSE]
varimax_loadings_IF <- varimax(loadings_IF)$loadings[1:nrow(loadings_IF), 1:5]

loadings_IF_df <- data.frame(
  Variable = rownames(varimax_loadings_IF),
  PC1 = round(varimax_loadings_IF[,1], 3),
  PC2 = round(varimax_loadings_IF[,2], 3), 
  PC3 = round(varimax_loadings_IF[,3], 3),
  PC4 = round(varimax_loadings_IF[,4], 3),
  PC5 = round(varimax_loadings_IF[,5], 3)
)

write.csv(loadings_IF_df, "plots/PCA_loadings_IF.csv", row.names = FALSE)

loadings_LM100 <- LM100_PCA$model$loadings[,1:3, drop=FALSE]
varimax_loadings_LM100 <- varimax(loadings_LM100)$loadings[1:nrow(loadings_LM100), 1:3]

loadings_LM100_df <- data.frame(
  Variable = rownames(varimax_loadings_LM100),
  PC1 = round(varimax_loadings_LM100[,1], 3),
  PC2 = round(varimax_loadings_LM100[,2], 3), 
  PC3 = round(varimax_loadings_LM100[,3], 3)
)

write.csv(loadings_LM100_df, "plots/PCA_loadings_LM100.csv", row.names = FALSE)

loadings_LM250 <- LM250_PCA$model$loadings[,1:3, drop=FALSE]
varimax_loadings_LM250 <- varimax(loadings_LM250)$loadings[1:nrow(loadings_LM250), 1:3]

loadings_LM250_df <- data.frame(
  Variable = rownames(varimax_loadings_LM250),
  PC1 = round(varimax_loadings_LM250[,1], 3),
  PC2 = round(varimax_loadings_LM250[,2], 3), 
  PC3 = round(varimax_loadings_LM250[,3], 3)
)

write.csv(loadings_LM250_df, "plots/PCA_loadings_LM250.csv", row.names = FALSE)

loadings_LM500 <- LM500_PCA$model$loadings[,1:3, drop=FALSE]
varimax_loadings_LM500 <- varimax(loadings_LM500)$loadings[1:nrow(loadings_LM500), 1:3]

loadings_LM500_df <- data.frame(
  Variable = rownames(varimax_loadings_LM500),
  PC1 = round(varimax_loadings_LM500[,1], 3),
  PC2 = round(varimax_loadings_LM500[,2], 3), 
  PC3 = round(varimax_loadings_LM500[,3], 3)
)

write.csv(loadings_LM500_df, "plots/PCA_loadings_LM500.csv", row.names = FALSE)

## Plots ####

# Environmental predictors (climate & topography)
loadings_matrix_EP <- EP_PCA$model$loadings[, 1:5]
varimax_results_EP <- varimax(loadings_matrix_EP)
rotated_loadings_EP <- varimax_results_EP$loadings

rotated_long_EP <- rotated_loadings_EP[1:nrow(rotated_loadings_EP), 1:5] %>%
  as.data.frame(row.names = rownames(rotated_loadings_EP)) %>%
  rownames_to_column("Variable") %>%
  rename(`EP_PC1` = Comp.1, `EP_PC2` = Comp.2, 
         `EP_PC3` = Comp.3, `EP_PC4` = Comp.4, `EP_PC5` = Comp.5) %>%
  pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
  mutate(PC = factor(PC, levels = c("EP_PC1", "EP_PC2", "EP_PC3", "EP_PC4", "EP_PC5"))) %>% 
  mutate(
    Label = ifelse(abs(Loading) >= 0.25, 
                   round(Loading, 2), NA_character_)
  )


Varimax_EP <- ggplot(rotated_long_EP, aes(PC, reorder(Variable, desc(Loading)), fill = Loading)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Loading") +
  labs(x = 'Principal Components',
       y = 'Climate and topography predictors') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Distance-based predictors

loadings_matrix_IF <- IF_PCA$model$loadings[, 1:5]
varimax_results_IF <- varimax(loadings_matrix_IF)
rotated_loadings_IF <- varimax_results_IF$loadings

rotated_long_IF <- rotated_loadings_IF[1:nrow(rotated_loadings_IF), 1:5] %>%
  as.data.frame(row.names = rownames(rotated_loadings_IF)) %>%
  rownames_to_column("Variable") %>%
  rename(`IF_PC1` = Comp.1, `IF_PC2` = Comp.2, 
         `IF_PC3` = Comp.3, `IF_PC4` = Comp.4, `IF_PC5` = Comp.5) %>%
  pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
  mutate(PC = factor(PC, levels = c("IF_PC1", "IF_PC2", "IF_PC3", "IF_PC4", "IF_PC5"))) 

rotated_long_IF <- rotated_long_IF %>%
  mutate(
    Label = ifelse(abs(Loading) >= 0.25, 
                   round(Loading, 2), NA_character_)
  )

Varimax_IF <- ggplot(rotated_long_IF, aes(PC, reorder(Variable, desc(Loading)), fill = Loading)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Loading") +
  labs(x = 'Principal Components',
       y = 'Distance-based predictors') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Landscape Metrics 100m

loadings_matrix_LM100 <- LM100_PCA$model$loadings[, 1:3]
varimax_results_LM100 <- varimax(loadings_matrix_LM100)
rotated_loadings_LM100 <- varimax_results_LM100$loadings

rotated_long_LM100 <- rotated_loadings_LM100[1:nrow(rotated_loadings_LM100), 1:3] %>%
  as.data.frame(row.names = rownames(rotated_loadings_LM100)) %>%
  rownames_to_column("Variable") %>%
  rename(`LM100_PC1` = Comp.1, `LM100_PC2` = Comp.2, 
         `LM100_PC3` = Comp.3) %>%
  pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
  mutate(PC = factor(PC, levels = c('LM100_PC1', 'LM100_PC2', 'LM100_PC3'))) 

rotated_long_LM100 <- rotated_long_LM100 %>%
  mutate(
    Label = ifelse(abs(Loading) >= 0.25, 
                   round(Loading, 2), NA_character_)
  )

Varimax_LM100 <- ggplot(rotated_long_LM100, aes(PC, reorder(Variable, desc(Loading)), fill = Loading)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Loading") +
  labs(x = 'Principal Components',
       y = 'Landscape metrics') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Landscape Metrics 250m

loadings_matrix_LM250 <- LM250_PCA$model$loadings[, 1:3]
varimax_results_LM250 <- varimax(loadings_matrix_LM250)
rotated_loadings_LM250 <- varimax_results_LM250$loadings

rotated_long_LM250 <- rotated_loadings_LM250[1:nrow(rotated_loadings_LM250), 1:3] %>%
  as.data.frame(row.names = rownames(rotated_loadings_LM250)) %>%
  rownames_to_column("Variable") %>%
  rename(`LM250_PC1` = Comp.1, `LM250_PC2` = Comp.2, 
         `LM250_PC3` = Comp.3) %>%
  pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
  mutate(PC = factor(PC, levels = c('LM250_PC1', 'LM250_PC2', 'LM250_PC3'))) 

rotated_long_LM250 <- rotated_long_LM250 %>%
  mutate(
    Label = ifelse(abs(Loading) >= 0.25, 
                   round(Loading, 2), NA_character_)
  )

Varimax_LM250 <- ggplot(rotated_long_LM250, aes(PC, reorder(Variable, desc(Loading)), fill = Loading)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Loading") +
  labs(x = 'Principal Components',
       y = 'Landscape metrics') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Landscape Metrics 500m

loadings_matrix_LM500 <- LM500_PCA$model$loadings[, 1:3]
varimax_results_LM500 <- varimax(loadings_matrix_LM500)
rotated_loadings_LM500 <- varimax_results_LM500$loadings

rotated_long_LM500 <- rotated_loadings_LM500[1:nrow(rotated_loadings_LM500), 1:3] %>%
  as.data.frame(row.names = rownames(rotated_loadings_LM500)) %>%
  rownames_to_column("Variable") %>%
  rename(`LM500_PC1` = Comp.1, `LM500_PC2` = Comp.2, 
         `LM500_PC3` = Comp.3) %>%
  pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
  mutate(PC = factor(PC, levels = c('LM500_PC1', 'LM500_PC2', 'LM500_PC3'))) 

rotated_long_LM500 <- rotated_long_LM500 %>%
  mutate(
    Label = ifelse(abs(Loading) >= 0.25, 
                   round(Loading, 2), NA_character_)
  )

Varimax_LM500 <- ggplot(rotated_long_LM500, aes(PC, reorder(Variable, desc(Loading)), fill = Loading)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Loading") +
  labs(x = 'Principal Components',
       y = 'Landscape metrics') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Variable Importance and Partial Dependence ####

# Read model data

occurrence_data_extract <- read.table('predictor_variables/predictor_variables_resampled/occurence_data_extract.txt')

training_data <- occurrence_data_extract[occurrence_data_extract$folds < 4,] # folds 1-3 (~ 60%)
training_data$occ <- as.factor(training_data$occ)
training_data$LULC <- as.factor(training_data$LULC)


# Baseline model: Environmental predictors, LULC at point-scale

baseline_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC) %>% 
  na.omit()

# Landscape metrics model: Standard model + Landscape Metrics 250m radius

LM250_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM250'))) %>% 
  na.omit()

# HSI-integrated model: Standard model + Landscape Metrics & HSI in 250m radius

LM_HSI250_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM250')), HSI_250) %>% 
  na.omit()

# Fine-scale model 100m: Standard model + Landscape Metrics & HSI in 100m radius

LM_HSI100_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM100')), HSI_100) %>% 
  na.omit()

# Larger-scale model 500m: Standard model + Landscape Metrics & HSI in 500m radius

LM_HSI500_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM500')), HSI_500) %>% 
  na.omit()


## Baseline model ####

rp_Baseline <- readRDS('rf_maps_new/rp_Baseline.rds')

Baseline_Importance <- importance(rp_Baseline)

# Main importance scores
main_imp_baseline <- Baseline_Importance %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  select(Variable, `0`, `1`, MeanDecreaseAccuracy) %>%
  pivot_longer(-Variable, names_to = "Class", values_to = "Importance") %>%
  mutate(Class = case_when(
    Class == "0" ~ "Absence",
    Class == "1" ~ "Presence", 
    Class == "MeanDecreaseAccuracy" ~ "Overall"
  ))

# P-values
pvals_baseline <- Baseline_Importance %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  select(Variable, `0.pval`, `1.pval`, `MeanDecreaseAccuracy.pval`) %>%
  pivot_longer(-Variable, names_to = "Class_raw", values_to = "pval") %>%
  mutate(Class = case_when(
    Class_raw == "0.pval" ~ "Absence",
    Class_raw == "1.pval" ~ "Presence",
    Class_raw == "MeanDecreaseAccuracy.pval" ~ "Overall"
  ),
  Sig = case_when(
    pval < 0.01 ~ "**", 
    pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Global order by Accuracy
accuracy_order_baseline <- Baseline_Importance %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  pull(Variable)

# Join & plot
imp_plot_baseline <- main_imp_baseline %>%
  left_join(pvals_baseline %>% select(Variable, Class, Sig), by = c("Variable", "Class")) %>% 
  mutate(Variable = factor(Variable, levels = rev(accuracy_order_baseline)),
         Class = factor(Class, levels = c("Overall", "Absence", "Presence")))

ggplot(imp_plot_baseline, aes(Variable, Importance, fill = Class)) +
  geom_col(color = "black", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = Sig), hjust = -0.5, fontface = "bold", size = 3) +
  coord_flip() +
  facet_wrap(~ Class, scales = "free_y", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Predictors", y = "Mean Decrease Accuracy") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size = 13))


pdp_EP_1_baseline <- partialPlot(rp_Baseline$rf, baseline_data, EP_PC1, which.class = 1)
pdp_EP_2_baseline <- partialPlot(rp_Baseline$rf, baseline_data, EP_PC2, which.class = 1)
pdp_EP_3_baseline <- partialPlot(rp_Baseline$rf, baseline_data, EP_PC3, which.class = 1)
pdp_LULC_baseline <- partialPlot(rp_Baseline$rf, baseline_data, LULC, which.class = 1)

PD_EP1_baseline <- ggplot(data.frame(x = pdp_EP_1_baseline$x, y = pdp_EP_1_baseline$y), aes(x = x, y = y)) +
  geom_line(size = 1.5, color = "darkred") +
  labs(x = "EP_PC1", y = "Partial Dependence (log-odds)") +
  theme_bw(base_size = 12)

PD_EP2_baseline <- ggplot(data.frame(x = pdp_EP_2_baseline$x, y = pdp_EP_2_baseline$y), aes(x = x, y = y)) +
  geom_line(size = 1.5, color = "darkred") +
  labs(x = "EP_PC2", y = "Partial Dependence (log-odds)") +
  theme_bw(base_size = 12)

PD_EP3_baseline <- ggplot(data.frame(x = pdp_EP_3_baseline$x, y = pdp_EP_3_baseline$y), aes(x = x, y = y)) +
  geom_line(size = 1.5, color = "darkred") +
  labs(x = "EP_PC3", y = "Partial Dependence (log-odds)") +
  theme_bw(base_size = 12)

pdp_LULC_df <- data.frame(x = as.factor(pdp_LULC_baseline$x), y = pdp_LULC_baseline$y)
pdp_LULC_df$x_label <- factor(pdp_LULC_df$x, levels = 1:14, 
                              labels = c('Artificial surfaces', 'Sparse settlement areas', 'Arable land', 'Permanent_crops',
                                         'Intensively used grasslands', 'Extensively used grasslands',
                                         'Natural grasslands', 'Coniferous forest', 'Deciduous forest',
                                         'Mixed forest', 'Alpine open land', 'Water bodies',
                                         'Wetland', 'Krummholz'))

PD_LULC_baseline <- ggplot(pdp_LULC_df, aes(x = x_label, y = y)) +
  geom_col(fill = "steelblue", width = 0.7) +
  labs(x = "LULC", y = "Partial Dependence (log-odds)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

(PD_EP1_baseline | PD_EP2_baseline | Varimax_EP) /
  (PD_EP3_baseline | PD_LULC_baseline | plot_spacer())

## Landscape metrics model ####

rp_LM250 <- readRDS('rf_maps_new/rp_LM_250m.rds')

LM250_Importance <- importance(rp_LM250)

# Main importance scores
main_imp_LM <- LM250_Importance %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  select(Variable, `0`, `1`, MeanDecreaseAccuracy) %>%
  pivot_longer(-Variable, names_to = "Class", values_to = "Importance") %>%
  mutate(Class = case_when(
    Class == "0" ~ "Absence",
    Class == "1" ~ "Presence", 
    Class == "MeanDecreaseAccuracy" ~ "Overall"
  ))

# P-values
pvals_LM <- LM250_Importance %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  select(Variable, `0.pval`, `1.pval`, `MeanDecreaseAccuracy.pval`) %>%
  pivot_longer(-Variable, names_to = "Class_raw", values_to = "pval") %>%
  mutate(Class = case_when(
    Class_raw == "0.pval" ~ "Absence",
    Class_raw == "1.pval" ~ "Presence",
    Class_raw == "MeanDecreaseAccuracy.pval" ~ "Overall"
  ),
  Sig = case_when(
    pval < 0.01 ~ "**", 
    pval < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Global order by Accuracy
accuracy_order_LM <- LM250_Importance %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  pull(Variable)

# Join & plot
imp_plot_LM <- main_imp_LM %>%
  left_join(pvals_LM %>% select(Variable, Class, Sig), by = c("Variable", "Class")) %>% 
  mutate(Variable = factor(Variable, levels = rev(accuracy_order_LM)),
         Class = factor(Class, levels = c("Overall", "Absence", "Presence")))

ggplot(imp_plot_LM, aes(Variable, Importance, fill = Class)) +
  geom_col(color = "black", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = Sig), hjust = -0.5, fontface = "bold", size = 2.5) +
  coord_flip() +
  facet_wrap(~ Class, scales = "free_y", ncol = 3) +
  labs(x = "Predictors", y = "Mean Decrease Accuracy") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

pdp_LM250_1 <- partialPlot(rp_LM250$rf, LM250_data, LM250_PC1, which.class = 1)
pdp_LM250_2 <- partialPlot(rp_LM250$rf, LM250_data, LM250_PC2, which.class = 1)
pdp_LM250_3 <- partialPlot(rp_LM250$rf, LM250_data, LM250_PC3, which.class = 1)

PD_LM250_1 <- ggplot(data.frame(x = pdp_LM250_1$x, y = pdp_LM250_1$y), aes(x = x, y = y)) +
  geom_line(size = 1.5, color = "darkred") +
  labs(x = "LM250_PC1", y = "Partial Dependence (log-odds)") +
  theme_bw(base_size = 12)

PD_LM250_2 <- ggplot(data.frame(x = pdp_LM250_2$x, y = pdp_LM250_2$y), aes(x = x, y = y)) +
  geom_line(size = 1.5, color = "darkred") +
  labs(x = "LM250_PC2", y = "Partial Dependence (log-odds)") +
  theme_bw(base_size = 12)

PD_LM250_3 <- ggplot(data.frame(x = pdp_LM250_3$x, y = pdp_LM250_3$y), aes(x = x, y = y)) +
  geom_line(size = 1.5, color = "darkred") +
  labs(x = "LM250_PC3", y = "Partial Dependence (log-odds)") +
  theme_bw(base_size = 12)

(PD_LM250_1 | PD_LM250_2)/
  (PD_LM250_3 | Varimax_LM250)

# Out-of-bag accuracy ####

acc_df <- read.csv('plots/oob_results.csv', sep = ',')

ggplot(acc_df, aes(x = Model, y = Accuracy)) +
  geom_col(fill = "grey70", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = LCI_95, ymax = UCI_95),
                width = 0.2, linewidth = 0.8) +
  labs(y = "OOB Accuracy (%)", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(82, 85))

# Quantile Regression ####

results_table <- read.csv('plots/quantile_regression_results.csv', sep = ',')

# Prepare slopes long format

slopes_H1_2 <- results_table %>% 
  select(tau, Baseline_slope, LM_250m_slope, LM_HSI_250m_slope) %>% 
  pivot_longer(cols = -tau, names_to = "model", values_to = "slope") %>% 
  mutate(model = sub("_slope$", "", model))

slopes_H3 <- results_table %>% 
  select(tau, LM_HSI_100m_slope, LM_HSI_250m_slope, LM_HSI_500m_slope) %>% 
  pivot_longer(cols = -tau, names_to = "model", values_to = "slope") %>% 
  mutate(model = sub("_slope$", "", model))

# Prepare R1 long format

R1_H1_2 <- results_table %>% 
  select(tau, Baseline_R1, LM_250m_R1, LM_HSI_250m_R1) %>% 
  pivot_longer(cols = -tau, names_to = "model", values_to = "R1") %>% 
  mutate(model = sub("_R1$", "", model))

R1_H3 <- results_table %>% 
  select(tau, LM_HSI_100m_R1, LM_HSI_250m_R1, LM_HSI_500m_R1) %>% 
  pivot_longer(cols = -tau, names_to = "model", values_to = "R1") %>% 
  mutate(model = sub("_R1$", "", model))

# Tau vs Slope H1/2
H1_2_slope <- ggplot(slopes_H1_2, aes(x = tau, y = slope, color = model, shape = model)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(x = expression(tau ~"(Quantile)"), y = "Slope",
       color = "Model", shape = "Model") +
  scale_color_manual(values = c("Baseline" = "#2E86AB", 
                                "LM_250m" = "#A23B72", 
                                "LM_HSI_250m" = "#F18F01"),
                     labels = c("Baseline" = "Baseline model", 
                                "LM_250m" = "Landscape metrics model", 
                                "LM_HSI_250m" = "HSI-integrated model")) +
  scale_shape_manual(values = c("Baseline" = 16, "LM_250m" = 17, "LM_HSI_250m" = 15),
                     labels = c("Baseline" = "Baseline model", 
                                "LM_250m" = "Landscape metrics model", 
                                "LM_HSI_250m" = "HSI-integrated model")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))


# Tau vs Slope H3
H3_slope <- ggplot(slopes_H3, aes(x = tau, y = slope, color = model, shape = model)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(x = expression(tau ~"(Quantile)"), y = "Slope",
       color = "Model", shape = "Model") +
  scale_color_manual(values = c("LM_HSI_100m" = "#006BA4", 
                                "LM_HSI_250m" = "#FF8E00", 
                                "LM_HSI_500m" = "#C52E60"),
                     labels = c("LM_HSI_100m" = "Fine-scale model (100m)", 
                                "LM_HSI_250m" = "Reference model (250m)", 
                                "LM_HSI_500m" = "Larger-scale model (500m)")) +
  scale_shape_manual(values = c("LM_HSI_100m" = 16, "LM_HSI_250m" = 17, "LM_HSI_500m" = 15),
                     labels = c("LM_HSI_100m" = "Fine-scale model (100m)", 
                                "LM_HSI_250m" = "Reference model (250m)", 
                                "LM_HSI_500m" = "Larger-scale model (500m)")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))


# Tau vs R1 H1/2
H1_2_R1 <- ggplot(R1_H1_2, aes(x = tau, y = R1, color = model, shape = model)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(x = expression(tau ~"(Quantile)"), y = "Goodness-of-fit R1",
       color = "Model", shape = "Model") +
  scale_color_manual(values = c("Baseline" = "#2E86AB", 
                                "LM_250m" = "#A23B72", 
                                "LM_HSI_250m" = "#F18F01"),
                     labels = c("Baseline" = "Baseline model", 
                                "LM_250m" = "Landscape metrics model", 
                                "LM_HSI_250m" = "HSI-integrated model")) +
  scale_shape_manual(values = c("Baseline" = 16, "LM_250m" = 17, "LM_HSI_250m" = 15),
                     labels = c("Baseline" = "Baseline model", 
                                "LM_250m" = "Landscape metrics model", 
                                "LM_HSI_250m" = "HSI-integrated model")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))


# Tau vs R1 H3
H3_R1 <- ggplot(R1_H3, aes(x = tau, y = R1, color = model, shape = model)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(x = expression(tau ~"(Quantile)"), y = "Goodness-of-fit R1",
       color = "Model", shape = "Model") +
  scale_color_manual(values = c("LM_HSI_100m" = "#006BA4", 
                                "LM_HSI_250m" = "#FF8E00", 
                                "LM_HSI_500m" = "#C52E60"),
                     labels = c("LM_HSI_100m" = "Fine-scale model (100m)", 
                                "LM_HSI_250m" = "Reference model (250m)", 
                                "LM_HSI_500m" = "Larger-scale model (500m)")) +
  scale_shape_manual(values = c("LM_HSI_100m" = 16, "LM_HSI_250m" = 17, "LM_HSI_500m" = 15),
                     labels = c("LM_HSI_100m" = "Fine-scale model (100m)", 
                                "LM_HSI_250m" = "Reference model (250m)", 
                                "LM_HSI_500m" = "Larger-scale model (500m)")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))


(H1_2_slope | H1_2_R1) +
  plot_layout(nrow = 2) & 
  theme(legend.position = "bottom")

(H3_slope | H3_R1) +
  plot_layout(nrow = 2) & 
  theme(legend.position = "bottom")

# Overview Map ####
districts <- st_read('Jagdreviere_nach_1983/Jagdreviere_nach_1983.shp')
hunting_points <- st_read('Red_Deer_Winter.gpkg')

districts <- st_transform(districts, st_crs(hunting_points))

districts$n_records <- lengths(st_intersects(districts, hunting_points))

# Density (records per km²)
districts$area_km2 <- units::drop_units(st_area(districts) / 1e6)
districts$density <- districts$n_records / districts$area_km2

plot(st_geometry(districts), col = hcl.colors(length(unique(districts$density)), "viridis"))

st_write(districts, "hunting_districts_density.shp")
