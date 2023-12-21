
library(tidyverse)
library(gt)
library(panelView)
library(fixest)
library(data.table)
library(purrr)
library(broom)
library(modelsummary)
library(purrr)
library(gtsummary)
library(interflex)
library(binsreg)


vcf <- read.csv("data/vcf_data_complete.csv")

gfc <- read.csv("data/gfc_dta.csv")

#summary table
vcf_final <- vcf |>
  select(forest_index, green_index, built_index, sch, cover_1990)
colnames(vcf_final) <- c(
  "Forest cover index (0-100)",
  "Non-forest green index (0-100)",
  "Non-green index (0-100)",
  "Scheduled Status",
  "Forest Cover in 1990 (Ex-Ante)"
)
vcf_final |>
  summarise_each(funs(mean = mean,
                      sd = sd,
                      n = n()))  |>
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_([^_]*)"
  ) |>
  gt() |>
  fmt_number(columns = c(mean, sd), decimals = 2) 




gfc_final <- gfc |>
  select(def_ha, sch, pref_mean)
colnames(gfc_final) <- c(
  "Deforested Area (Hectares)",
  "Scheduled Status",
  "Ex-ante forest cover in 2000 (ex-ante)"
)
gfc_final  |>
  summarise_each(funs(mean = mean,
                      sd = sd,
                      n = n())) |>
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_([^_]*)"
  ) |>
  gt() |>
  fmt_number(columns = c(mean, sd), decimals = 2) 


#panel view figure
state_status <- vcf %>%
  filter(year >= 1990) %>%
  group_by(state, year) %>%
  summarize(out = 1, treat = max(D))


panel_figure <- panelView::panelview(
  out ~ treat,
  data = as.data.frame(state_status),
  index = c("state", "year"),
  axis.lab.gap = c(1, 0),
  by.timing = TRUE,
  legendOff = TRUE,
  background = "white"
)
panel_figure



#main table
# %% GFC regressions
m00 = feols(def_ha ~ D |
              village + year, cluster = "block", gfc |> filter(pref == 1))
# 2wFE + state year FEs
m01 = feols(def_ha ~ D | village + styear,
            cluster = "block",
            gfc |>
              filter(gfc3 == TRUE & pref == 1))
m003 = feols(def_ha ~ D |
               village + village[t] + styear,
             cluster = "block",
             gfc |>
               filter(gfc3 == TRUE & pref == 1))

#vcf
vcf0 <- vcf %>%
  filter(year >= 1990)

vcf1 <- vcf0 %>%
  mutate(t = year - 1995,
         ex_ante_med = quantile(cover_1990, 0.5))

# Filter data and add never_treated column in one step
vcf2 <- vcf1 %>%
  filter(cover_1990 > ex_ante_med) %>%
  group_by(cellid) %>%
  mutate(never_treated = max(D) == 0)
controls_pre1 <- vcf2 |>
  filter(never_treated == 1 &
           year < first_pesa_exposure) %>%
  summarize(mean(forest_index))

treat_pre1 <- vcf2 |>
  filter(never_treated == 0 &
           year < first_pesa_exposure) |>
  summarise(mean(forest_index))
m0 = feols(forest_index ~ D |
             cellid + year,
           data = vcf2,
           cluster = "blk")
m1 = feols(forest_index ~ D |
             cellid + styear,
           data = vcf2,
           cluster = "blk")
m3 = feols(forest_index ~ D |
             cellid[t] + styear,
           data = vcf2,
           cluster = "blk")

col1 <-
  c("", "", "", "2001-2017", "GFC", "", "8.809", "12.30")
col2 <-
  c("", "", "", "2001-2017", "GFC", "", "8.809", "12.30")
col3 <-
  c("", "X", "", "2001-2017", "GFC", "", "8.809", "12.30")
col4 <-
  c("", "", "", "1995-2017", "VCF", "", "0.0800", "0.1300")
col5 <-
  c("", "", "", "1995-2017", "VCF", "", "0.0800", "0.1300")
col6 <-
  c("", "", "X", "1995-2017", "VCF", "", "0.0800", "0.1300")
col0 <-
  c(
    "time trends",
    "t (Pixel)",
    "t (Village)",
    "timespan",
    "Dataset",
    "summary statistics",
    "Mean Y (Non-Sch)",
    "Mean Y (Sch)"
  )
row <-
  data.frame(col0, col1, col2, col3, col4, col5, col6)
gm <-
  c("nobs",
    "R2",
    "RMSE",
    "Std.Errors",
    "FE: styear",
    "FE: cellid",
    "FE: village",
    "FE: year")
mod <- modelsummary(
  list(m00, m01, m003, m0, m1, m3),
  stars = TRUE,
  output = "gt",
  gof_map = c(
    "r.squared",
    "nobs",
    "rmse",
    "Std.Errors",
    "FE: styear",
    "FE: styear",
    "FE: village",
    "FE: year"
  ),
  add_rows = row,
  coef_rename = "PESA"
) |>
  tab_footnote(footnote = "Notes: Standard errors are clustered at the block level and reported in parentheses.")

#main diagram
h <- c(1:10)
fitter = function(cutoff) {
  dat <- gfc |>
    filter(pref_bin >= cutoff)
  m = feols(def_ha ~ D |
              village[t] + styear, cluster = ~ block, dat)
  tidy(m)
}

fitter2 = function(cutoff) {
  dat <- gfc |>
    filter(pref_bin >= cutoff)
  m = feols(def_ha ~ D |
              village[t] + styear, cluster = ~ block, dat)
  return(m)
}

modelsummary(map(h, fitter2))
#GFC plot

cutoff_res = map(1:10, fitter) |>
  list_rbind(names_to = "h")

rob_fit_gfc <-
  ggplot(cutoff_res, aes(h, estimate)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymax = estimate + 1.96 * std.error,
      ymin = estimate - 1.96 * std.error
    ),
    alpha = 1,
    width = 0
  ) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = 'GFC',
    y = "Effect (Deforested Area)",
    x = "Inclusion Threshold Decile of Forest Cover in 2000 ",
    caption = ''
  ) +
  theme(
    axis.text = element_text(size = 22),
    text = element_text(size = 23),
    plot.title = element_text(size = 26),
    legend.position = 'none',
    panel.grid = element_line(color = "gray")
  ) +
  theme_classic()
rob_fit_gfc
ggsave("gfc.png", height = 10, width = 20)
#vcf plot

vcf4 <-
  vcf |>  mutate(pref_bin = ntile(cover_1990, 10),
                 t = vcf$year - 1995)

fitter3 = function(cut) {
  m = feols(
    forest_index ~ D | cellid[t] + styear,
    data = vcf4 |> filter(pref_bin >= cut),
    cluster = ~ blk
  )
  tidy(m)
}
fitter4 = function(cut) {
  m = feols(
    forest_index ~ D | cellid[t] + styear,
    data = vcf4 |> filter(cover_1990 >= cut),
    cluster = ~ blk
  )
  return(m)
}
cutmods <- map(h, fitter3) |>
  list_rbind(names_to = "h")

rob_fit_vcf = ggplot(cutmods, aes(h, estimate)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymax = estimate + 1.96 * std.error,
      ymin = estimate - 1.96 * std.error
    ),
    alpha = 1,
    width = 0
  ) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = 1:10) +
  scale_colour_brewer(palette = 'Set1') +
  labs(
    title = 'VCF',
    y = "Effect (Forest Index)",
    x = "Inclusion Threshold Decile of Forest Cover in 1990",
    caption = ''
  ) +
  theme(
    axis.text = element_text(size = 22),
    text = element_text(size = 23),
    plot.title = element_text(size = 26),
    legend.position = 'none'
  ) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_classic()
rob_fit_vcf





#appendix
#FRA
h <- c(1:10)
vcff <- vcf |>
  filter(year >= 1990) |>
  mutate(
    pref_bin = ntile(cover_1990, 10),
    D_fra = sch * (year >= 2008),
    time_fra = year - 2008
  )
fitter6 <- function(cut) {
  m = feols(
    forest_index ~ D_fra | cellid + styear,
    data = vcff |> filter(pref_bin >= cut),
    cluster = ~ blk
  )
  tidy(m)
}

cutmods1 <-  map(h, fitter6) |>
  list_rbind(names_to = "h")

# %%
FRA_fit_vcf <-
  ggplot(cutmods1, aes(h, estimate)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymax = estimate + 1.96 * std.error,
      ymin = estimate - 1.96 * std.error
    ),
    alpha = 1,
    width = 0
  ) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = 1:10) +
  scale_colour_brewer(palette = 'Set1') +
  labs(
    title = 'VCF',
    y = "Effect (Forest Index)",
    x = "Inclusion Threshold Decile of Forest Cover in 1990",
    caption = ''
  ) +
  theme(
    axis.text = element_text(size = 22),
    text = element_text(size = 23),
    plot.title = element_text(size = 26),
    legend.position = 'none'
  ) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_classic()


FRA_fit_vcf



#figure 8 mines

vcf4 <- vcf |>
  filter(year >= 1990) |>
  mutate(
    pref_bin = ntile(cover_1990, 10),
    time_fra = year - 2008,
    D_fra = sch * (year >= 2008)
  )


xcreg <- vcf4 |>
  filter(year == first_pesa_exposure - 1 &
           min_dist_to_mine <= 1) |>
  mutate(del_forest = cover_1990 - forest_index,
         dist2 = min_dist_to_mine * 100)



gen <- binsreg(
  y = xcreg$del_forest,
  x = xcreg$dist2,
  nbins = 30,
  polyreg = 2
)
gen$bins_plot +
  labs(y = 'Decrease in Forest Index \nRelative to 1990 Baseline',
       x = 'Distance (km)') +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 5)
  )


vcf5 <- vcf %>%
  # Convert specified columns to factors
  mutate(
    id_f = factor(cellid),
    sty_f = factor(styear),
    year_f = factor(year),
    D_f = factor(D),
    distance = min_dist_to_mine * 110
  )
intsamp <- vcf5 |>
  filter(distance <= 100)
mining_all = interflex(
  Y = "forest_index",
  D = "D_f",
  X = "distance",
  FE = c("id_f", "sty_f"),
  cl = "id_f",
  data = intsamp,
  theme.bw = T,
  Xlabel = 'Distance',
  Dlabel = "PESA",
  Ylabel = '\nForest Index (VCF)\n',
  cutoffs = seq(0, 100, 10),
  estimator = 'binning',
  CI = TRUE,
  cex.lab = 1.2,
  cex.axis = 1.2,
)


mining_all
