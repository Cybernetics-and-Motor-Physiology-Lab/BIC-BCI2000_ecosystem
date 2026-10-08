# Long-term stimulation spike-count analysis
# Negative-binomial models with longitudinal adjustment and recording-day clustering
# Primary model: negative-binomial mixed model (glmmTMB)
#
# Required CSV columns:
#   condition, spikeCount, logExposure, spikeRate, daysSinceImplant, dayID
#
# The model uses raw spike counts as the response and logExposure as an offset,
# where exposure = recording duration (min) x number of valid channels.

# -----------------------------
# 1. Packages
# -----------------------------
required_packages <- c(
  "MASS", "glmmTMB", "ggplot2", "dplyr", "readr", "stringr", "splines"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(MASS)
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(splines)

# -----------------------------
# 2. User settings
# -----------------------------
input_file <- "sub-c04-longtermStim.csv"
output_prefix <- "sub-c04-longtermStim-NBGLMM"

# Baseline is the reference category. Remove unused levels automatically.
condition_levels <- c(
  "Baseline", "Stim50Hz", "Stim100Hz", "Stim125Hz"
) #, "Stim133Hz", "Stim200Hz"

# Multiple-comparison correction for condition-vs-baseline coefficients.
p_adjust_method <- "holm"

# MATLAB default colors for the first four plotted conditions.
condition_colors <- c(
  "Baseline"  = rgb(0.0000, 0.4470, 0.7410),
  "Stim50Hz"  = rgb(0.8500, 0.3250, 0.0980),
  "Stim100Hz" = rgb(0.9290, 0.6940, 0.1250),
  "Stim125Hz" = rgb(0.4940, 0.1840, 0.5560)
)
# "Stim133Hz" = "#77AC30",
# "Stim200Hz" = "#4DBEEE"

# -----------------------------
# 3. Load and prepare data
# -----------------------------
df <- read_csv(input_file, show_col_types = FALSE)

required_columns <- c(
  "condition", "spikeCount", "logExposure", "spikeRate",
  "daysSinceImplant", "dayID"
)

missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

df <- df %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    spikeCount = as.integer(spikeCount),
    logExposure = as.numeric(logExposure),
    spikeRate = as.numeric(spikeRate),
    daysSinceImplant = as.numeric(daysSinceImplant),
    dayID = factor(dayID)
  ) %>%
  filter(
    !is.na(condition),
    !is.na(spikeCount),
    !is.na(logExposure),
    !is.na(spikeRate),
    !is.na(daysSinceImplant),
    !is.na(dayID)
  )

df$condition <- droplevels(df$condition)
df$condition <- relevel(df$condition, ref = "Baseline")

# -----------------------------
# 4. Define longitudinal time automatically
# -----------------------------
# The earliest recording in this long-term stimulation dataset is treated as
# the start of the stimulation analysis period (t = 0).
stimStartDay <- min(df$daysSinceImplant, na.rm = TRUE)

df <- df %>%
  mutate(
    daysSinceStimStart = daysSinceImplant - stimStartDay,
    weeksSinceStimStart = daysSinceStimStart / 7
  )

cat("\nLong-term stimulation analysis start:\n")
cat("  Day since implant: ", stimStartDay, "\n", sep = "")
cat("  Observations: ", nrow(df), "\n", sep = "")
cat(
  "  Follow-up: ",
  round(max(df$weeksSinceStimStart, na.rm = TRUE), 1),
  " weeks\n\n",
  sep = ""
)

# -----------------------------
# 5. Fit candidate models
# -----------------------------
# M0: stimulation condition only
m0 <- glm.nb(
  spikeCount ~ condition + offset(logExposure),
  data = df
)

# M1: linear longitudinal adjustment
m1 <- glm.nb(
  spikeCount ~ condition + weeksSinceStimStart + offset(logExposure),
  data = df
)

# M2: nonlinear longitudinal adjustment (sensitivity analysis)
m2 <- glm.nb(
  spikeCount ~ condition + ns(weeksSinceStimStart, df = 3) + offset(logExposure),
  data = df
)

# M3: PRIMARY MODEL
# Linear longitudinal trend + random intercept for recording day.
m3 <- glmmTMB(
  spikeCount ~
    condition +
    weeksSinceStimStart +
    offset(logExposure) +
    (1 | dayID),
  family = nbinom2,
  data = df
)

# M4: nonlinear-time mixed model (sensitivity analysis)
m4 <- glmmTMB(
  spikeCount ~
    condition +
    ns(weeksSinceStimStart, df = 3) +
    offset(logExposure) +
    (1 | dayID),
  family = nbinom2,
  data = df
)

# -----------------------------
# 6. Compare model specifications
# -----------------------------
model_comparison <- AIC(m0, m1, m2, m3, m4)
model_comparison$deltaAIC <- model_comparison$AIC - min(model_comparison$AIC)

cat("\nAIC model comparison:\n")
print(model_comparison)

cat("\nLikelihood-ratio test: M0 vs M1 (linear time adjustment)\n")
print(anova(m0, m1, test = "LRT"))

cat("\nMixed-model comparison: M3 linear time vs M4 nonlinear time\n")
print(anova(m3, m4))

write.csv(
  model_comparison,
  paste0(output_prefix, "_model_comparison_AIC.csv"),
  row.names = TRUE
)

# -----------------------------
# 7. Primary model results
# -----------------------------
nb_model <- m3
model_summary <- summary(nb_model)
print(model_summary)

coef_table <- as.data.frame(model_summary$coefficients$cond)
coef_table$Variable <- rownames(coef_table)
rownames(coef_table) <- NULL

coef_table <- coef_table %>%
  mutate(
    IRR = exp(Estimate),
    Lower95 = exp(Estimate - 1.96 * `Std. Error`),
    Upper95 = exp(Estimate + 1.96 * `Std. Error`),
    pValue = `Pr(>|z|)`
  )

# Condition effects relative to Baseline.
condition_results <- coef_table %>%
  filter(str_detect(Variable, "^condition")) %>%
  mutate(
    condition = str_replace(Variable, "^condition", ""),
    p_adj = p.adjust(pValue, method = p_adjust_method),
    sigLabel = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

cat("\nPrimary-model condition effects relative to Baseline:\n")
print(condition_results)

cat("\nLongitudinal effect:\n")
print(coef_table %>% filter(Variable == "weeksSinceStimStart"))

write.csv(
  coef_table,
  paste0(output_prefix, "_all_coefficients.csv"),
  row.names = FALSE
)

write.csv(
  condition_results,
  paste0(output_prefix, "_condition_effects_vs_baseline.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Prepare annotations for plotting
# -----------------------------
plot_results <- condition_results %>%
  mutate(
    condition = factor(condition, levels = levels(df$condition)),
    label_full = sprintf(
      "IRR = %.2f\n95%% CI %.2f-%.2f\n%s",
      IRR,
      Lower95,
      Upper95,
      sigLabel
    )
  )

# -----------------------------
# 9. Boxplot + individual observations
# -----------------------------
# Plotting-only normalization:
# spikeCount / log(number of valid channels x recording duration)
df$normalizedSpikeCount <- df$spikeCount / df$logExposure

y_max <- max(df$normalizedSpikeCount, na.rm = TRUE)

# Keep annotations inside the requested 0-to-maximum data range.
plot_results$y <- y_max * 0.96

# Keep only colors corresponding to conditions actually present.
plot_colors <- condition_colors[
  names(condition_colors) %in% levels(df$condition)
]

p <- ggplot(
  df,
  aes(
    x = condition,
    y = normalizedSpikeCount,
    fill = condition
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    alpha = 0.65,
    color = "black"
  ) +
  geom_jitter(
    aes(color = condition),
    width = 0.20,
    height = 0,
    size = 2,
    alpha = 0.30
  ) +
  geom_text(
    data = plot_results,
    aes(
      x = condition,
      y = y,
      label = label_full
    ),
    inherit.aes = FALSE,
    size = 3.8,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = plot_colors) +
  scale_color_manual(values = plot_colors) +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = pretty(c(0, y_max), n = 10),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "",
    y = "Normalized spike count"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    
    panel.grid.major.y = element_line(
      colour = "grey80",
      linewidth = 0.4
    ),
    
    panel.grid.minor.y = element_line(
      colour = "grey92",
      linewidth = 0.25
    ),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    plot.margin = margin(10, 20, 10, 10)
  )

print(p)

ggsave(
  paste0(output_prefix, "_boxplot_scatter_IRR.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  paste0(output_prefix, "_boxplot_scatter_IRR.svg"),
  p,
  width = 8,
  height = 6
)


# -----------------------------
# 10. Longitudinal diagnostic plot
# -----------------------------

p_time <- ggplot(
  df,
  aes(
    x = weeksSinceStimStart,
    y = normalizedSpikeCount,
    color = condition
  )
) +
  geom_point(
    alpha = 0.65,
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  scale_color_manual(values = plot_colors) +
  scale_y_continuous(
    limits = c(0, max(df$normalizedSpikeCount, na.rm = TRUE)),
    breaks = pretty(
      c(0, max(df$normalizedSpikeCount, na.rm = TRUE)),
      n = 10
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Weeks since long-term stimulation start",
    y = "Normalized spike count",
    color = "Condition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.35
    ),
    
    panel.grid.minor.y = element_line(
      colour = "grey93",
      linewidth = 0.2
    )
  )

print(p_time)

ggsave(
  paste0(output_prefix, "_longitudinal_diagnostic.png"),
  p_time,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  paste0(output_prefix, "_longitudinal_diagnostic.svg"),
  p_time,
  width = 9,
  height = 6
)

cat("\nAnalysis complete. Primary model: M3\n")
cat("Primary formula:\n")
print(formula(m3))
