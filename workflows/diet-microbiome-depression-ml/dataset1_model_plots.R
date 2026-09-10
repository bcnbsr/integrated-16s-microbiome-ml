#!/usr/bin/env Rscript

# ML reporting plots for the selected ASV or OTU branch.
# Linear models use signed coefficients; tree models use feature importance.

required_packages <- c("dplyr", "ggplot2", "pROC")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

library(dplyr)
library(ggplot2)
library(pROC)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "dataset1")
plot_dir <- if (length(args) >= 2) args[[2]] else file.path(results_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

classification_file <- file.path(results_dir, "classification_results_all_splits.csv")
prediction_file <- file.path(results_dir, "model_predictions_all_splits.csv")
attribution_file <- file.path(results_dir, "model_feature_attributions_all_splits.csv")
required_files <- c(classification_file, prediction_file, attribution_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Required ML output file(s) not found. Rerun the ML-analysis cell before plotting:\n",
    paste(missing_files, collapse = "\n")
  )
}

require_columns <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(data_name, " is missing required column(s): ", paste(missing, collapse = ", "))
  }
}

theme_thesis_clean <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_plot_pair <- function(plot_object, file_stem, width, height) {
  ggsave(file.path(plot_dir, paste0(file_stem, ".pdf")), plot_object, width = width, height = height, units = "in")
  ggsave(file.path(plot_dir, paste0(file_stem, ".png")), plot_object, width = width, height = height, units = "in", dpi = 600)
}

short_taxon_label <- function(feature_name, feature_group) {
  if (identical(feature_group, "Diet")) return(feature_name)
  parts <- trimws(strsplit(as.character(feature_name), ";", fixed = TRUE)[[1]])
  genus <- sub("^g__", "", parts[grepl("^g__", parts)])
  genus <- genus[nzchar(genus)]
  if (length(genus) > 0) return(tail(genus, 1))
  named_parts <- sub("^[a-z]__", "", parts)
  named_parts <- named_parts[nzchar(named_parts)]
  if (length(named_parts) > 0) return(tail(named_parts, 1))
  "Unclassified taxon"
}

classification_results <- read.csv(classification_file, check.names = FALSE, stringsAsFactors = FALSE)
model_predictions <- read.csv(prediction_file, check.names = FALSE, stringsAsFactors = FALSE)
model_attributions <- read.csv(attribution_file, check.names = FALSE, stringsAsFactors = FALSE)

require_columns(classification_results, c("Split", "Input features", "Model", "ROC-AUC"), "Classification-results file")
require_columns(model_predictions, c("Split", "Model", "y_true", "y_pred", "score"), "Model-predictions file")
require_columns(model_attributions, c("Split", "Input features", "Model", "Attribution type", "Feature ID", "Feature name", "Feature group", "Attribution", "Rank by absolute attribution"), "Model-attribution file")

model_order <- c("LR", "SVM", "RF", "BRF")
model_colours <- c("LR" = "#4C78A8", "SVM" = "#F58518", "RF" = "#54A24B", "BRF" = "#B279A2")
split_colours <- c("1" = "#4C78A8", "2" = "#F58518", "3" = "#54A24B", "4" = "#E45756", "5" = "#B279A2")

# Plot 1: all classifiers' ROC-AUC distributions.
auc_data <- classification_results %>%
  filter(`Input features` == "Diet + taxa") %>%
  mutate(Model = factor(Model, levels = model_order))
if (nrow(auc_data) == 0) stop("No combined diet-taxa classification results were found.")
auc_plot <- ggplot(auc_data, aes(x = Model, y = `ROC-AUC`, fill = Model)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.6, alpha = 0.75) +
  geom_jitter(aes(colour = Model), width = 0.08, size = 2.5, alpha = 0.9) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", colour = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.5, colour = "grey40") +
  scale_fill_manual(values = model_colours) + scale_colour_manual(values = model_colours) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Model ROC-AUC Across Five Repeated Stratified Splits", x = "Classifier", y = "ROC-AUC") +
  theme_thesis_clean() + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")
save_plot_pair(auc_plot, "01_Model_AUC_Across_Five_Repeated_Splits", 7, 5)

# Plot 2: one mean ROC panel for each classifier.
fpr_grid <- seq(0, 1, length.out = 501)
roc_by_split <- model_predictions %>%
  mutate(Split = as.integer(Split)) %>%
  group_by(Model, Split) %>%
  group_modify(~ {
    roc_object <- pROC::roc(.x$y_true, .x$score, levels = c(0, 1), direction = "<", quiet = TRUE)
    points <- data.frame(FPR = 1 - roc_object$specificities, TPR = roc_object$sensitivities) %>%
      group_by(FPR) %>% summarise(TPR = max(TPR), .groups = "drop") %>% arrange(FPR)
    data.frame(FPR = fpr_grid, TPR = approx(points$FPR, points$TPR, xout = fpr_grid, rule = 2)$y)
  }) %>% ungroup()
mean_roc <- roc_by_split %>% group_by(Model, FPR) %>% summarise(TPR = mean(TPR), .groups = "drop")
roc_plot <- ggplot(mean_roc, aes(x = FPR, y = TPR, colour = Model)) +
  geom_step(linewidth = 1.05) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = model_colours, breaks = model_order) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(~ Model, ncol = 2) +
  labs(title = "Mean ROC Curves Across Five Repeated Splits", x = "False positive rate", y = "True positive rate") +
  theme_thesis_clean() + theme(legend.position = "none")
save_plot_pair(roc_plot, "02_Mean_ROC_Curves_All_Models", 8.5, 7)

# Plot 3: pooled confusion matrix for each classifier.
confusion_rows <- lapply(model_order, function(method) {
  observed <- model_predictions %>%
    filter(Model == method) %>%
    mutate(True_label = ifelse(y_true == 1, "Depressive", "Healthy"), Predicted_label = ifelse(y_pred == 1, "Depressive", "Healthy")) %>%
    count(True_label, Predicted_label, name = "n")
  grid <- expand.grid(True_label = c("Healthy", "Depressive"), Predicted_label = c("Healthy", "Depressive"), stringsAsFactors = FALSE)
  merged <- merge(grid, observed, by = c("True_label", "Predicted_label"), all.x = TRUE)
  merged$n[is.na(merged$n)] <- 0
  merged %>% group_by(True_label) %>% mutate(row_percent = 100 * n / sum(n), label = paste0(n, "\n", sprintf("%.1f", row_percent), "%"), Model = method) %>% ungroup()
}) %>% bind_rows() %>% mutate(Model = factor(Model, levels = model_order), True_label = factor(True_label, levels = c("Healthy", "Depressive")), Predicted_label = factor(Predicted_label, levels = c("Healthy", "Depressive")))
confusion_plot <- ggplot(confusion_rows, aes(x = Predicted_label, y = True_label, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.9) + geom_text(aes(label = label), size = 4, fontface = "bold") +
  scale_fill_gradient(low = "#F7FBFF", high = "#4C78A8") + facet_wrap(~ Model, ncol = 2) +
  labs(title = "Pooled Confusion Matrices Across Five Repeated Splits", x = "Predicted label", y = "True label", fill = "Test instances") +
  theme_thesis_clean() + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))
save_plot_pair(confusion_plot, "03_Pooled_Confusion_Matrices_All_Models", 8.5, 7)

make_attribution_plot <- function(method, attribution_type, file_stem) {
  data <- model_attributions %>%
    filter(Model == method, `Input features` == "Diet + taxa", `Attribution type` == attribution_type) %>%
    mutate(Split = as.integer(Split), In_top20_this_split = `Rank by absolute attribution` <= 20)
  frequency <- data %>% group_by(`Feature ID`, `Feature name`, `Feature group`) %>%
    summarise(Splits_in_top20 = sum(In_top20_this_split), .groups = "drop") %>% filter(Splits_in_top20 >= 2)
  if (nrow(frequency) == 0) {
    warning("No ", method, " features were in the top 20 for at least two splits; no attribution plot was created.")
    return(invisible(NULL))
  }
  plot_data <- data %>% inner_join(frequency, by = c("Feature ID", "Feature name", "Feature group")) %>%
    filter(In_top20_this_split) %>%
    rowwise() %>% mutate(Short_feature_name = short_taxon_label(`Feature name`, `Feature group`)) %>% ungroup()
  label_table <- plot_data %>% distinct(`Feature ID`, Short_feature_name, Splits_in_top20, `Feature group`) %>%
    arrange(desc(Splits_in_top20), `Feature group`, Short_feature_name) %>%
    mutate(Feature_label = make.unique(paste0(Short_feature_name, "\n(", Splits_in_top20, "/5)")), Feature_y = rev(seq_len(n())))
  plot_data <- plot_data %>% inner_join(label_table %>% select(`Feature ID`, Feature_label, Feature_y), by = "Feature ID") %>%
    mutate(Split_offset = (Split - 3) * 0.13, Plot_y = Feature_y + Split_offset)
  ranges <- plot_data %>% group_by(Feature_label, Feature_y) %>% summarise(xmin = min(Attribution), xmax = max(Attribution), .groups = "drop")
  display_labels <- label_table %>% arrange(Feature_y)
  quantity_label <- if (attribution_type == "coefficient") paste(method, "coefficient") else paste(method, "feature importance")
  plot_title <- if (attribution_type == "coefficient") paste(method, "Coefficients Across Five Repeated Splits") else paste(method, "Feature Importances Across Five Repeated Splits")
  attribution_plot <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey45") +
    geom_segment(data = ranges, aes(x = xmin, xend = xmax, y = Feature_y, yend = Feature_y), linewidth = 0.8, colour = "grey55") +
    geom_point(data = plot_data, aes(x = Attribution, y = Plot_y, fill = factor(Split)), shape = 21, size = 5.8, colour = "black", stroke = 0.45) +
    geom_text(data = plot_data, aes(x = Attribution, y = Plot_y, label = Split), size = 2.8, fontface = "bold") +
    scale_fill_manual(values = split_colours, breaks = as.character(1:5)) +
    scale_y_continuous(breaks = display_labels$Feature_y, labels = display_labels$Feature_label) +
    labs(title = plot_title, x = quantity_label, y = "Feature", fill = "Split") +
    theme_thesis_clean() + theme(axis.text.y = element_text(size = 10), legend.position = "bottom")
  save_plot_pair(attribution_plot, file_stem, 10, max(6.5, 0.42 * nrow(label_table) + 2.5))
}

make_attribution_plot("SVM", "coefficient", "04_SVM_Coefficients_Across_Five_Repeated_Splits")
make_attribution_plot("LR", "coefficient", "05_LR_Coefficients_Across_Five_Repeated_Splits")
make_attribution_plot("RF", "feature_importance", "06_RF_Feature_Importances_Across_Five_Repeated_Splits")
make_attribution_plot("BRF", "feature_importance", "07_BRF_Feature_Importances_Across_Five_Repeated_Splits")

cat("ML plots saved in:\n", plot_dir, "\n")
