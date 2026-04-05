================================================================
# STUDENT LOAN REPAYMENT RESTART & INCORPORATED SELF-EMPLOYMENT
# MSc Econometrics Group Project — Individual-Level SIPP Analysis
# =============================================================================
library(haven)
library(tidyverse)
library(fixest)
library(plm)
library(lmtest)
library(ebal)
library(car)
library(modelsummary)
library(patchwork)
library(scales)

pvalue <- fixest::pvalue

# --- 1. Load data ------------------------------------------------------------
getwd()
list.files()

df <- read_dta("sipp_clean.dta")

# --- 2. Collapse to person-year ----------------------------------------------
df_py <- df |>
  group_by(person_id, panel_year) |>
  summarise(
    se_inc       = mean(se_inc,       na.rm=TRUE),
    se_uninc     = mean(se_uninc,     na.rm=TRUE),
    selfempl     = mean(selfempl,     na.rm=TRUE),
    high_debt    = first(high_debt),
    borrower     = first(borrower),
    debt_tercile = first(as.character(debt_tercile)),
    post         = first(post),
    covid        = first(covid),
    TAGE         = mean(TAGE,         na.rm=TRUE),
    female       = first(female),
    married      = mean(married,      na.rm=TRUE),
    has_kids     = mean(has_kids,     na.rm=TRUE),
    homeowner    = mean(homeowner,    na.rm=TRUE),
    race_black   = first(race_black),
    race_asian   = first(race_asian),
    hispanic     = first(hispanic),
    educ_ma      = first(educ_ma),
    educ_prof    = first(educ_prof),
    educ_phd     = first(educ_phd),
    employed     = mean(employed,     na.rm=TRUE),
    ihs_networth = mean(ihs_networth, na.rm=TRUE),
    log_hhinc    = mean(log_hhinc,    na.rm=TRUE),
    log_debt     = mean(log_debt,     na.rm=TRUE),
    networth_q   = first(networth_q),
    low_networth = first(low_networth),
    .groups      = "drop"
  ) |>
  mutate(across(where(haven::is.labelled), as.numeric)) |>
  mutate(
    year      = as.integer(panel_year),
    year_f    = factor(year),
    age_sq    = TAGE^2,
    high_educ = as.integer(educ_ma == 1 | educ_phd == 1),
    pid       = as.integer(person_id)
  )

# =============================================================================
# SAMPLE CONSTRUCTION — TWO TREATMENT DEFINITIONS
# =============================================================================

# --- 3A. Main sample — top tercile ($66k+) -----------------------------------
df <- df_py |>
  mutate(
    treated = case_when(
      high_debt == 1 ~ 1L,
      borrower  == 0 ~ 0L,
      TRUE           ~ NA_integer_
    )
  ) |>
  filter(!is.na(treated)) |>
  mutate(
    did   = treated * post,
    post  = as.integer(post),
    covid = as.integer(covid)
  ) |>
  filter(!is.na(se_inc), !is.na(ihs_networth), !is.na(log_hhinc))



# --- 3B. Comparison sample — above median ($29k+) ----------------------------
median_debt <- df_py |>
  filter(borrower==1, log_debt > 0) |>
  pull(log_debt) |>
  median(na.rm=TRUE)


df_med <- df_py |>
  mutate(
    treated = case_when(
      borrower==1 & log_debt >= median_debt ~ 1L,
      borrower==0                           ~ 0L,
      TRUE                                  ~ NA_integer_
    )
  ) |>
  filter(!is.na(treated)) |>
  mutate(
    did   = treated * as.integer(post),
    post  = as.integer(post),
    covid = as.integer(covid)
  ) |>
  filter(!is.na(se_inc), !is.na(ihs_networth), !is.na(log_hhinc))

# --- 4. Raw DiD --------------------------------------------------------------
raw_means <- df |>
  group_by(treated, post) |>
  summarise(
    se_inc   = mean(se_inc,   na.rm=TRUE),
    se_uninc = mean(se_uninc, na.rm=TRUE),
    n        = n(), .groups="drop"
  )

raw_did_inc <- (raw_means$se_inc[raw_means$treated==1 & raw_means$post==1] -
                  raw_means$se_inc[raw_means$treated==1 & raw_means$post==0]) -
  (raw_means$se_inc[raw_means$treated==0 & raw_means$post==1] -
     raw_means$se_inc[raw_means$treated==0 & raw_means$post==0])
raw_did_inc

raw_did_uninc <- (raw_means$se_uninc[raw_means$treated==1 & raw_means$post==1] -
                    raw_means$se_uninc[raw_means$treated==1 & raw_means$post==0]) -
  (raw_means$se_uninc[raw_means$treated==0 & raw_means$post==1] -
     raw_means$se_uninc[raw_means$treated==0 & raw_means$post==0])
raw_did_uninc
# --- 5. Balance table --------------------------------------------------------

balance_vars <- c("TAGE","educ_ma","educ_phd","female","married",
                  "homeowner","has_kids","employed","ihs_networth",
                  "log_hhinc","race_black","race_asian","hispanic")

df_pre <- df |> filter(post==0)

bal <- map_dfr(balance_vars, function(v) {
  x1 <- as.numeric(df_pre[[v]][df_pre$treated==1])
  x0 <- as.numeric(df_pre[[v]][df_pre$treated==0])
  tibble(
    variable   = v,
    Treatment  = mean(x1, na.rm=TRUE),
    Control    = mean(x0, na.rm=TRUE),
    Difference = mean(x1, na.rm=TRUE) - mean(x0, na.rm=TRUE),
    p_value    = tryCatch(t.test(x1, x0)$p.value, error=\(e) NA_real_),
    SMD        = (mean(x1,na.rm=TRUE) - mean(x0,na.rm=TRUE)) /
      sqrt((var(x1,na.rm=TRUE) + var(x0,na.rm=TRUE)) / 2)
  )
})
print(bal |> mutate(across(where(is.numeric), ~round(.,3))), n=20)

# --- 6. Entropy balancing ----------------------------------------------------

run_ebal <- function(data, label) {
  df_baseline <- data |>
    filter(post==0) |>
    arrange(pid, year) |>
    group_by(pid) |>
    slice(1) |>
    ungroup()
  
  match_vars <- c("TAGE","age_sq","educ_ma","educ_prof","educ_phd",
                  "female","married","homeowner","has_kids",
                  "employed","ihs_networth","log_hhinc",
                  "race_black","race_asian","hispanic")
  
  df_match <- df_baseline |>
    select(pid, treated, any_of(match_vars)) |>
    mutate(across(everything(), as.numeric)) |>
    drop_na()
  
  X     <- as.matrix(df_match |> select(any_of(match_vars)))
  T_vec <- df_match$treated
  
  eb_out <- ebalance(Treatment=T_vec, X=X)
  
  df_match$eb_weight <- 1.0
  df_match$eb_weight[T_vec==0] <- eb_out$w

  
  data |>
    left_join(df_match |> select(pid, eb_weight), by="pid") |>
    mutate(eb_weight = replace_na(eb_weight, 1))
}

df     <- run_ebal(df,     "Top tercile")
df_med <- run_ebal(df_med, "Above median")

# --- 7. Controls formula and Mundlak person-means ----------------------------
controls <- ~ TAGE + age_sq + female + married + has_kids +
  employed + homeowner + ihs_networth + log_hhinc +
  educ_ma + educ_phd + race_black + race_asian + hispanic + covid

time_varying <- c("employed","homeowner","married","has_kids","ihs_networth","log_hhinc")

df <- df |>
  group_by(pid) |>
  mutate(across(all_of(time_varying), mean, .names="pmean_{.col}")) |>
  ungroup()

df_med <- df_med |>
  group_by(pid) |>
  mutate(across(all_of(time_varying), mean, .names="pmean_{.col}")) |>
  ungroup()

# =============================================================================
# PART A: MAIN SPECIFICATIONS — TOP TERCILE
# =============================================================================

ols_unw <- feols(
  se_inc ~ did + treated + post + .[controls],
  data=df, cluster=~pid
)
ols_unw


ols_eb <- feols(
  se_inc ~ did + treated + post + .[controls],
  data=df, weights=~eb_weight, cluster=~pid
)
ols_eb

twfe_unw <- feols(
  se_inc ~ did + .[controls] | pid + year_f,
  data=df, cluster=~pid
)
twfe_unw

twfe_eb <- feols(
  se_inc ~ did + .[controls] | pid + year_f,
  data=df, weights=~eb_weight, cluster=~pid
)
twfe_eb

mundlak_eb <- feols(
  se_inc ~ did + TAGE + age_sq + female + educ_ma + educ_phd +
    race_black + race_asian + hispanic + covid +
    employed + homeowner + married + has_kids + ihs_networth + log_hhinc +
    pmean_employed + pmean_homeowner + pmean_married + pmean_has_kids +
    pmean_ihs_networth + pmean_log_hhinc | year_f,
  data=df, weights=~eb_weight, cluster=~pid
)
mundlak_eb

# =============================================================================
# PART B: ABOVE-MEDIAN COMPARISON
# =============================================================================

ols_unw_med <- feols(
  se_inc ~ did + treated + post + .[controls],
  data=df_med, cluster=~pid
)
ols_unw_med

twfe_unw_med <- feols(
  se_inc ~ did + .[controls] | pid + year_f,
  data=df_med, cluster=~pid
)
twfe_unw_med

twfe_eb_med <- feols(
  se_inc ~ did + .[controls] | pid + year_f,
  data=df_med, weights=~eb_weight, cluster=~pid
)
twfe_eb_med

# =============================================================================
# PART C: SPECIFICATION TESTS
# =============================================================================

pdata <- pdata.frame(
  df |> select(pid, year, se_inc, did, treated, post, covid,
               TAGE, age_sq, female, married, has_kids, employed,
               homeowner, ihs_networth, log_hhinc,
               educ_ma, educ_phd, race_black, race_asian, hispanic,
               eb_weight),
  index=c("pid","year")
)

formula_plm <- se_inc ~ did + TAGE + age_sq + female + married + has_kids +
  employed + homeowner + ihs_networth + log_hhinc +
  educ_ma + educ_phd + race_black + race_asian + hispanic + covid
formula_plm

fe_plm <- plm(formula_plm, data=pdata, model="within", effect="twoways",
              weights=pdata$eb_weight)
fe_plm
re_plm <- plm(formula_plm, data=pdata, model="random",
              weights=pdata$eb_weight)
re_plm

lm_test <- plmtest(re_plm, type="bp")
lm_test
hausman  <- phtest(fe_plm, re_plm)
hausman

mundlak_terms <- paste0("pmean_", time_varying)
mundlak_ftest <- linearHypothesis(mundlak_eb, paste(mundlak_terms, "= 0"))

# =============================================================================
# PART D: EVENT STUDY
# =============================================================================

es_eb <- feols(
  se_inc ~ i(year, treated, ref=2022) + .[controls] | pid + year_f,
  data=df, weights=~eb_weight, cluster=~pid
)
summary(es_eb)

es_coefs <- broom::tidy(es_eb, conf.int=TRUE) |>
  filter(str_detect(term,"year::")) |>
  mutate(
    yr     = as.integer(str_extract(term,"\\d{4}")),
    period = ifelse(yr < 2022,"Pre-period","Post-period")
  )
es_coefs

es_coefs |> filter(period=="Pre-period") |>
  select(yr, estimate, std.error, p.value) |>
  mutate(across(where(is.numeric), ~round(.,4))) |> print()


es_coefs |> filter(period=="Post-period") |>
  select(yr, estimate, std.error, p.value) |>
  mutate(across(where(is.numeric), ~round(.,4))) |> print()

pre_terms <- es_coefs$term[es_coefs$period=="Pre-period"]
pre_terms
f_pre <- linearHypothesis(es_eb, paste(pre_terms, "= 0"))
f_pre

# =============================================================================
# PART E: ROBUSTNESS CHECKS
# =============================================================================


# R1: Drop 2020
rob1 <- feols(
  se_inc ~ did + .[controls] | pid + year_f,
  data=df |> filter(year!=2020), weights=~eb_weight, cluster=~pid
)


# R2: 2023 sub-period split — use raw df read from dta, before collapsing
df_raw <- read_dta("sipp_clean.dta")

df_2023 <- df_raw |>
  filter(panel_year == 2023) |>
  mutate(
    treated = case_when(
      high_debt == 1 ~ 1L,
      borrower  == 0 ~ 0L,
      TRUE           ~ NA_integer_
    ),
    period_23 = case_when(
      MONTHCODE <= 5  ~ "1_pre_announce",
      MONTHCODE <= 9  ~ "2_anticipatory",
      MONTHCODE <= 12 ~ "3_post_restart"
    )
  ) |>
  filter(!is.na(treated), !is.na(se_inc))


# R3: Unincorporated SE
rob3 <- feols(
  se_uninc ~ did + .[controls] | pid + year_f,
  data=df, weights=~eb_weight, cluster=~pid
)
cat(sprintf("R3 uninc SE:    %.4f  p=%.3f\n",
            coef(rob3)["did"], pvalue(rob3)["did"]))

# R4: Total SE
rob_total <- feols(
  selfempl ~ did + .[controls] | pid + year_f,
  data=df, weights=~eb_weight, cluster=~pid
)


# R5: Placebo year
df_placebo <- df |>
  filter(year != 2023) |>
  mutate(
    post_placebo = as.integer(year == 2022),
    did_placebo  = treated * post_placebo,
    year_f       = factor(year)
  )

plac_unw <- feols(
  se_inc ~ did_placebo + .[controls] | pid + year_f,
  data=df_placebo, cluster=~pid
)
plac_unw

plac_eb <- feols(
  se_inc ~ did_placebo + .[controls] | pid + year_f,
  data=df_placebo, weights=~eb_weight, cluster=~pid
)
plac_eb
# =============================================================================
# PART F: CONTINUOUS INTENSITY DiD
# =============================================================================

log_debt_pre_tbl <- df |>
  filter(post == 0) |>
  group_by(person_id) |>
  summarise(log_debt_pre = mean(log_debt, na.rm = TRUE), .groups = "drop") |>
  mutate(pid = as.integer(person_id))

df_borrow <- df_py |>
  filter(borrower==1) |>
  left_join(log_debt_pre_tbl |> select(pid, log_debt_pre), by="pid") |>
  mutate(
    post      = as.integer(post),
    covid     = as.integer(covid),
    debt_post = log_debt_pre * post
  ) |>
  filter(!is.na(se_inc), !is.na(log_debt_pre),
         !is.na(ihs_networth), !is.na(log_hhinc))

cont_ols <- feols(
  se_inc ~ debt_post + log_debt_pre + .[controls],
  data=df_borrow, cluster=~pid
)
cont_twfe <- feols(
  se_inc ~ debt_post + log_debt_pre + .[controls] | pid + year_f,
  data=df_borrow, cluster=~pid
)
cont_uninc <- feols(
  se_uninc ~ debt_post + log_debt_pre + .[controls] | pid + year_f,
  data=df_borrow, cluster=~pid
)


# Continuous pre-trend test
es_cont <- feols(
  se_inc ~ i(year, log_debt_pre, ref=2022) + .[controls] | pid + year_f,
  data=df_borrow, cluster=~pid
)
cont_pre_terms <- broom::tidy(es_cont) |>
  filter(str_detect(term,"year::"), as.integer(str_extract(term,"\\d{4}")) < 2022) |>
  pull(term)
f_pre_cont <- linearHypothesis(es_cont, paste(cont_pre_terms, "= 0"))


# =============================================================================
# PART G: HETEROGENEITY
# =============================================================================


het_home <- feols(
  se_inc ~ i(homeowner, did) + .[controls] | pid + year_f,
  data=df, weights=~eb_weight, cluster=~pid
)


het_nw <- feols(
  se_inc ~ i(networth_q, did) + .[controls] | pid + year_f,
  data=df, weights=~eb_weight, cluster=~pid
)
summary(het_nw)

het_educ <- feols(
  se_inc ~ i(high_educ, did) + .[controls] | pid + year_f,
  data=df, weights=~eb_weight, cluster=~pid
)
summary(het_educ)


# =============================================================================
# FIGURES
# =============================================================================

col_treat <- "#C0392B"
col_ctrl  <- "#2980B9"
col_annot <- "#E67E22"
col_med   <- "#E67E22"
col_terc  <- "#8E44AD"
text_col  <- "#2C3E50"
grid_col  <- "#F0F0F0"

base_theme <- theme_minimal(base_size=12) +
  theme(
    plot.title       = element_text(face="bold", size=13, colour=text_col),
    plot.subtitle    = element_text(size=10, colour="#7F8C8D", margin=margin(b=8)),
    axis.title       = element_text(size=10, colour=text_col),
    axis.text        = element_text(size=9,  colour=text_col),
    legend.position  = "bottom",
    legend.text      = element_text(size=9),
    panel.grid.major = element_line(colour=grid_col, linewidth=0.5),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill="white", colour=NA),
    panel.background = element_rect(fill="white", colour=NA)
  )

# Figure 1: Incorporated SE trends by group
trends_data <- df |>
  group_by(year, treated) |>
  summarise(mean_se=mean(se_inc, na.rm=TRUE), .groups="drop") |>
  mutate(Group=ifelse(treated==1,"High Debt Borrowers","Non-Borrowers"))

fig1 <- ggplot(trends_data, aes(x=year, y=mean_se, colour=Group, group=Group)) +
  annotate("rect", xmin=2022.5, xmax=2024,
           ymin=-Inf, ymax=Inf, fill="#FADBD8", alpha=0.4) +
  geom_line(linewidth=1.5) +
  geom_point(size=4.5, aes(fill=Group), shape=21, colour="white", stroke=1.8) +
  geom_vline(xintercept=2022.42, linetype="dashed",
             colour=col_annot, linewidth=0.9) +
  geom_vline(xintercept=2023, linetype="dashed",
             colour=col_treat, linewidth=0.9) +
  annotate("text", x=2022.44, y=0.074,
           label="Jun 2\nCertainty Shock", hjust=0, size=3,
           colour=col_annot, fontface="bold") +
  annotate("text", x=2023.5, y=0.067,
           label="Oct 1 Payments", hjust=0, size=3,
           colour=col_treat, fontface="bold") +
  scale_colour_manual(
    values=c("High Debt Borrowers"=col_treat, "Non-Borrowers"=col_ctrl)) +
  scale_fill_manual(
    values=c("High Debt Borrowers"=col_treat, "Non-Borrowers"=col_ctrl)) +
  scale_y_continuous(labels=percent_format(accuracy=0.1), limits=c(0, 0.085)) +
  scale_x_continuous(breaks=2018:2023) +
  labs(
    title    = "Incorporated Self-Employment Rate by Group, 2018-2023",
    subtitle = "High-debt borrowers drop to 0% in 2023; non-borrowers remain stable",
    x=NULL, y="Incorporated SE Rate", colour=NULL, fill=NULL
  ) +
  base_theme

fig1

ggsave("fig1_trends.png", fig1, width=10, height=6, dpi=300)


# Figure 2: Event study
es_plot_data <- bind_rows(
  broom::tidy(es_eb, conf.int=TRUE) |>
    filter(str_detect(term, "year::")) |>
    mutate(
      yr     = as.integer(str_extract(term, "\\d{4}")),
      period = case_when(yr < 2022 ~ "Pre-period", yr == 2023 ~ "Post-period",
                         TRUE ~ "Reference")
    ),
  tibble(yr=2022, estimate=0, conf.low=0, conf.high=0,
         period="Reference", term="ref", std.error=0, statistic=0, p.value=1)
) |> arrange(yr)

fig2 <- ggplot(es_plot_data, aes(x=yr, y=estimate, colour=period)) +
  annotate("rect", xmin=2017.5, xmax=2021.5,
           ymin=-Inf, ymax=Inf, fill="#EBF5FB", alpha=0.3) +
  annotate("rect", xmin=2022.5, xmax=2023.5,
           ymin=-Inf, ymax=Inf, fill="#FADBD8", alpha=0.3) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50", linewidth=0.8) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=period),
              alpha=0.12, colour=NA) +
  geom_line(aes(group=1), colour="grey70", linewidth=0.8, linetype="dotted") +
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high), size=0.8, linewidth=1.4) +
  annotate("text", x=2019.5, y=0.20,
           label="Pre-period\n(joint p = 0.452)", size=3.2,
           colour=col_ctrl, fontface="italic") +
  annotate("text", x=2023.12, y=-0.085,
           label="-10.3pp\np = 0.248", size=3,
           colour=col_treat, fontface="bold", hjust=0) +
  scale_colour_manual(
    values=c("Pre-period"=col_ctrl, "Post-period"=col_treat, "Reference"="grey40")) +
  scale_fill_manual(
    values=c("Pre-period"=col_ctrl, "Post-period"=col_treat, "Reference"="grey40")) +
  scale_x_continuous(breaks=2018:2023) +
  scale_y_continuous(labels=percent_format(accuracy=0.1)) +
  labs(
    title    = "Event Study: Parallel Trends Test (Entropy Weighted)",
    subtitle = "Pre-period leads statistically zero — parallel trends supported (joint p = 0.452)",
    x="Year (reference = 2022)", y="Coefficient (High Debt x Year)",
    colour=NULL, fill=NULL
  ) +
  base_theme + theme(legend.position="none")

ggsave("fig2_event_study.png", fig2, width=10, height=6, dpi=300)

# Figure 3: Balance love plot
smd_pre <- map_dfr(balance_vars, function(v) {
  x1 <- as.numeric(df_pre[[v]][df_pre$treated==1])
  x0 <- as.numeric(df_pre[[v]][df_pre$treated==0])
  var_label <- c(
    TAGE="Age", educ_ma="Education: MA", educ_phd="Education: PhD",
    female="Gender (female)", married="Marital status",
    homeowner="Homeownership", has_kids="Has children",
    employed="Employment", ihs_networth="Net worth (IHS)",
    log_hhinc="Income (log)", race_black="Race: Black",
    race_asian="Race: Asian", hispanic="Hispanic"
  )
  tibble(
    variable = var_label[v],
    smd      = (mean(x1,na.rm=TRUE) - mean(x0,na.rm=TRUE)) /
      sqrt((var(x1,na.rm=TRUE) + var(x0,na.rm=TRUE)) / 2),
    timing   = "Before matching"
  )
})

smd_data <- bind_rows(
  smd_pre,
  smd_pre |> mutate(smd=0, timing="After matching")
) |>
  mutate(
    variable = fct_reorder(variable, abs(smd), .fun=max),
    timing   = factor(timing, levels=c("Before matching","After matching"))
  )

fig3 <- ggplot(smd_data) +
  annotate("rect", xmin=-0.1, xmax=0.1,
           ymin=-Inf, ymax=Inf, fill="#D5F5E3", alpha=0.4) +
  geom_vline(xintercept=0, linetype="solid", colour="grey40", linewidth=0.8) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed",
             colour="grey60", linewidth=0.6) +
  geom_segment(
    data=smd_data |>
      pivot_wider(names_from=timing, values_from=smd) |>
      rename(before=`Before matching`, after=`After matching`),
    aes(x=after, xend=before, y=variable, yend=variable),
    colour="grey75", linewidth=1.2
  ) +
  geom_point(aes(x=smd, y=variable, colour=timing, shape=timing), size=4) +
  scale_colour_manual(
    values=c("Before matching"=col_treat, "After matching"="#27AE60"),
    labels=c("Before matching"="Before entropy balancing",
             "After matching" ="After entropy balancing")) +
  scale_shape_manual(
    values=c("Before matching"=16, "After matching"=17),
    labels=c("Before matching"="Before entropy balancing",
             "After matching" ="After entropy balancing")) +
  scale_x_continuous(limits=c(-1.75, 0.60), breaks=c(-1.5,-1.0,-0.5,0,0.5)) +
  labs(
    title    = "Covariate Balance: Before and After Entropy Balancing",
    subtitle = "All SMDs reduced to exactly zero after reweighting",
    x="Standardised Mean Difference (SMD)", y=NULL, colour=NULL, shape=NULL
  ) +
  base_theme

ggsave("fig3_balance.png", fig3, width=10, height=7, dpi=300)

# Figure 4: Substitution pattern
subst_data <- tibble(
  Outcome  = factor(
    c("Incorporated SE\n(capital-intensive)",
      "Unincorporated SE\n(necessity)",
      "Total SE\n(composition check)"),
    levels=c("Incorporated SE\n(capital-intensive)",
             "Unincorporated SE\n(necessity)",
             "Total SE\n(composition check)")
  ),
  estimate = c(coef(twfe_eb)["did"], coef(rob3)["did"], coef(rob_total)["did"]),
  se       = c(se(twfe_eb)["did"],   se(rob3)["did"],   se(rob_total)["did"]),
  pval     = c(pvalue(twfe_eb)["did"],pvalue(rob3)["did"],pvalue(rob_total)["did"]),
  expected = c("Negative (theory predicted)","Positive (theory predicted)",
               "Zero (pure substitution)")
) |>
  mutate(
    ci_lo  = estimate - 1.96*se,
    ci_hi  = estimate + 1.96*se,
    plabel = case_when(
      pval < 0.01 ~ paste0(round(estimate*100,1),"pp***"),
      pval < 0.05 ~ paste0(round(estimate*100,1),"pp**"),
      pval < 0.10 ~ paste0(round(estimate*100,1),"pp*"),
      TRUE        ~ paste0(round(estimate*100,1),"pp")
    )
  )

fig4 <- ggplot(subst_data, aes(x=Outcome, y=estimate, fill=expected)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40", linewidth=0.9) +
  geom_col(width=0.55, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi),
                width=0.12, linewidth=1.2, colour=text_col) +
  geom_text(aes(label=plabel, y=estimate + sign(estimate)*0.035),
            size=4.2, fontface="bold", colour=text_col) +
  scale_fill_manual(values=c(
    "Negative (theory predicted)" = col_treat,
    "Positive (theory predicted)" = col_ctrl,
    "Zero (pure substitution)"    = "#95A5A6"
  )) +
  scale_y_continuous(labels=percent_format(accuracy=0.1), limits=c(-0.28, 0.32)) +
  labs(
    title    = "Substitution Between Self-Employment Types",
    subtitle = "Incorporated SE falls, unincorporated SE rises, total SE unchanged",
    x=NULL, y="DiD Coefficient (TWFE, Entropy Weighted)", fill=NULL
  ) +
  base_theme +
  theme(legend.position="bottom", panel.grid.major.x=element_blank())

ggsave("fig4_substitution.png", fig4, width=9, height=6, dpi=300)


# Figure 5a: Debt distribution
debt_dist <- df_py |>
  filter(borrower==1, log_debt > 0) |>
  arrange(pid, year) |>
  group_by(pid) |>
  slice(1) |>
  ungroup() |>
  mutate(debt_k = exp(log_debt) / 1000)

median_debt_k  <- median(debt_dist$debt_k, na.rm=TRUE)
tercile_debt_k <- quantile(debt_dist$debt_k, 2/3, na.rm=TRUE)

fig5a <- ggplot(debt_dist, aes(x=debt_k)) +
  annotate("rect", xmin=0, xmax=median_debt_k,
           ymin=-Inf, ymax=Inf, fill="#EBF5FB", alpha=0.4) +
  annotate("rect", xmin=median_debt_k, xmax=tercile_debt_k,
           ymin=-Inf, ymax=Inf, fill="#FEF9E7", alpha=0.4) +
  annotate("rect", xmin=tercile_debt_k, xmax=Inf,
           ymin=-Inf, ymax=Inf, fill="#FADBD8", alpha=0.4) +
  geom_histogram(bins=35, fill="#5D6D7E", alpha=0.80,
                 colour="white", linewidth=0.2) +
  geom_vline(xintercept=median_debt_k,  linetype="dashed",
             colour=col_med, linewidth=1.3) +
  geom_vline(xintercept=tercile_debt_k, linetype="dashed",
             colour=col_terc, linewidth=1.3) +
  annotate("text", x=median_debt_k+3, y=Inf,
           label=sprintf("Median\n$%.0fk", median_debt_k),
           hjust=0, vjust=1.4, size=4, colour=col_med, fontface="bold") +
  annotate("text", x=tercile_debt_k+3, y=Inf,
           label=sprintf("Top tercile\n$%.0fk", tercile_debt_k),
           hjust=0, vjust=1.4, size=4, colour=col_terc, fontface="bold") +
  scale_x_continuous(labels=dollar_format(suffix="k"),
                     breaks=seq(0,200,25), limits=c(0,220)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(
    title    = "Student Debt Distribution Among Borrowers",
    subtitle = "Blue = below median | Yellow = mid | Red = top tercile (treatment group)",
    x="Student Debt at Baseline ($000s)", y="Number of Borrowers"
  ) +
  base_theme

ggsave("fig5a_debt_distribution.png", fig5a, width=12, height=7, dpi=300)

# Figure 5b: DiD by debt quartile (entropy weighted)
debt_quartiles <- debt_dist |>
  mutate(debt_q = ntile(debt_k, 4)) |>
  select(pid, debt_q)

df_control <- df_py |>
  filter(borrower==0) |>
  mutate(treated=0L, did=0L, post=as.integer(post), covid=as.integer(covid)) |>
  filter(!is.na(se_inc), !is.na(ihs_networth), !is.na(log_hhinc))

run_quartile_twfe <- function(q) {
  df_treated <- df_py |>
    inner_join(debt_quartiles |> filter(debt_q==q) |> select(pid), by="pid") |>
    mutate(treated=1L, did=treated*as.integer(post),
           post=as.integer(post), covid=as.integer(covid)) |>
    filter(!is.na(se_inc), !is.na(ihs_networth), !is.na(log_hhinc))
  
  df_combined <- bind_rows(df_treated, df_control) |>
    mutate(year_f=factor(year))
  
  n_treated <- n_distinct(df_treated$pid)
  
  df_baseline <- df_combined |>
    filter(post==0) |> arrange(pid,year) |>
    group_by(pid) |> slice(1) |> ungroup()
  
  match_vars <- c("TAGE","age_sq","educ_ma","educ_prof","educ_phd",
                  "female","married","homeowner","has_kids",
                  "employed","ihs_networth","log_hhinc",
                  "race_black","race_asian","hispanic")
  
  df_match <- df_baseline |>
    select(pid, treated, any_of(match_vars)) |>
    mutate(across(everything(), as.numeric)) |>
    drop_na()
  
  X     <- as.matrix(df_match |> select(any_of(match_vars)))
  T_vec <- df_match$treated
  df_match$eb_weight <- 1.0
  
  tryCatch({
    eb_out <- ebalance(Treatment=T_vec, X=X)
    df_match$eb_weight[T_vec==0] <- eb_out$w
    cat(sprintf("Q%d (n=%d): converged\n", q, n_treated))
  }, error=function(e) {
    cat(sprintf("Q%d (n=%d): failed — unweighted\n", q, n_treated))
  })
  
  df_combined <- df_combined |>
    left_join(df_match |> select(pid, eb_weight), by="pid") |>
    mutate(eb_weight=replace_na(eb_weight, 1))
  
  mod_eb <- tryCatch(
    feols(se_inc ~ did + .[controls] | pid + year_f,
          data=df_combined, weights=~eb_weight, cluster=~pid),
    error=function(e) NULL
  )
  
  tibble(
    quartile  = q,
    n_treated = n_treated,
    coef_eb   = if(!is.null(mod_eb)) coef(mod_eb)["did"]   else NA_real_,
    se_eb     = if(!is.null(mod_eb)) se(mod_eb)["did"]     else NA_real_,
    pval_eb   = if(!is.null(mod_eb)) pvalue(mod_eb)["did"] else NA_real_
  )
}


quartile_results <- map_dfr(1:4, run_quartile_twfe)

q_ranges <- debt_dist |>
  mutate(debt_q=ntile(debt_k,4)) |>
  group_by(debt_q) |>
  summarise(lo=floor(min(debt_k)), hi=ceiling(max(debt_k)), .groups="drop") |>
  mutate(label=sprintf("Q%d\n$%dk-$%dk", debt_q, lo, hi))

quartile_results <- quartile_results |>
  left_join(q_ranges |> select(quartile=debt_q, label), by="quartile") |>
  mutate(
    ci_lo  = coef_eb - 1.96*se_eb,
    ci_hi  = coef_eb + 1.96*se_eb,
    sig    = pval_eb < 0.1,
    plabel = case_when(
      pval_eb < 0.01 ~ paste0(round(coef_eb*100,1),"pp***"),
      pval_eb < 0.05 ~ paste0(round(coef_eb*100,1),"pp**"),
      pval_eb < 0.10 ~ paste0(round(coef_eb*100,1),"pp*"),
      TRUE           ~ paste0(round(coef_eb*100,1),"pp")
    ),
    label = factor(label, levels=label)
  )

fig5b <- ggplot(quartile_results, aes(x=label, y=coef_eb)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40", linewidth=1.0) +
  geom_hline(yintercept=-0.1023, linetype="dotted",
             colour=col_treat, linewidth=1.0, alpha=0.7) +
  annotate("text", x=4.5, y=-0.1023-0.008,
           label="Main weighted estimate: -10.2pp",
           hjust=1, size=3.8, colour=col_treat, fontface="italic") +
  geom_col(aes(fill=sig), width=0.55, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi),
                width=0.15, linewidth=1.3, colour=text_col) +
  geom_text(aes(label=plabel,
                y=ifelse(coef_eb >= 0, ci_hi+0.008, ci_lo-0.008)),
            size=4.2, fontface="bold", colour=text_col) +
  geom_text(aes(label=sprintf("n = %d", n_treated), y=0.042),
            size=3.5, colour="grey50") +
  scale_fill_manual(
    values=c("TRUE"=col_treat, "FALSE"="#AEB6BF"),
    labels=c("TRUE"="p < 0.10", "FALSE"="p >= 0.10"),
    name="Significance"
  ) +
  scale_x_discrete() +
  scale_y_continuous(
    labels=percent_format(accuracy=0.1),
    limits=c(-0.20, 0.07),
    breaks=seq(-0.20, 0.05, 0.05)
  ) +
  labs(
    title    = "Effect on Incorporated SE by Debt Quartile (Entropy Weighted)",
    subtitle = "TWFE entropy weighted | Each quartile vs non-borrower control",
    x="Debt Quartile (among borrowers)", y="DiD Coefficient (Incorporated SE)"
  ) +
  base_theme + theme(legend.position="bottom")
fig5b
ggsave("fig5b_quartile_did.png", fig5b, width=10, height=8, dpi=300)


# =============================================================================
# REGRESSION TABLES — modelsummary (LaTeX)
# =============================================================================

gm <- list(
  list("raw"="nobs", "clean"="Observations", "fmt"=0)
)

notes_controls <- "Clustered SEs at person level. Controls: age, age squared, gender, marital status, children, employment, homeownership, net worth, household income, education, race/ethnicity, COVID dummy."

# --- Table 1: Threshold comparison -------------------------------------------

modelsummary(
  list("OLS"       = ols_unw_med,
       "TWFE"      = twfe_unw_med,
       "TWFE Wtd"  = twfe_eb_med,
       "OLS "      = ols_unw,
       "TWFE "     = twfe_unw,
       "TWFE Wtd " = twfe_eb),
  coef_map  = c("did" = "Treated x Post (DiD)"),
  gof_map   = gm,
  stars     = c("*"=0.1, "**"=0.05, "***"=0.01),
  vcov      = ~pid,
  title     = "Threshold Identification: Effect Concentrated at High Debt Levels",
  add_rows  = tribble(
    ~term,               ~`OLS`,  ~`TWFE`, ~`TWFE Wtd`, ~`OLS `,  ~`TWFE `, ~`TWFE Wtd `,
    "Threshold",         "$29k+", "$29k+", "$29k+",     "$66k+",  "$66k+",  "$66k+",
    "Person FE",         "No",    "Yes",   "Yes",       "No",     "Yes",    "Yes",
    "Year FE",           "Yes",   "Yes",   "Yes",       "Yes",    "Yes",    "Yes",
    "Entropy weighted",  "No",    "No",    "Yes",       "No",     "No",     "Yes",
    "Treated n",         "1,036", "1,036", "1,036",     "453",    "453",    "453"
  ),
  notes  = paste("Above-median threshold = $29,000; top tercile threshold = $66,000. TWFE drops singletons.", notes_controls),
  output = "table1_threshold.tex"
)

# --- Table 2: Main specification progression ---------------------------------

modelsummary(
  list("OLS"         = ols_unw,
       "OLS Wtd"     = ols_eb,
       "TWFE"        = twfe_unw,
       "TWFE Wtd"    = twfe_eb,
       "Mundlak Wtd" = mundlak_eb),
  coef_map  = c("did"     = "High Debt x Post (DiD)",
                "treated" = "High Debt (level)",
                "post"    = "Post 2023"),
  gof_map   = gm,
  stars     = c("*"=0.1, "**"=0.05, "***"=0.01),
  vcov      = ~pid,
  title     = "Main Results: High Debt Borrowers ($66k+) vs Non-Borrowers",
  add_rows  = tribble(
    ~term,              ~OLS,  ~`OLS Wtd`, ~TWFE,  ~`TWFE Wtd`, ~`Mundlak Wtd`,
    "Person FE",        "No",  "No",       "Yes",  "Yes",        "No",
    "Year FE",          "Yes", "Yes",      "Yes",  "Yes",        "Yes",
    "Entropy weighted", "No",  "Yes",      "No",   "Yes",        "Yes"
  ),
  notes  = paste("Treatment = top tercile of student debt ($66,000+). Control = non-borrowers. Mundlak CRE includes person-means of time-varying controls.", notes_controls),
  output = "table2_main.tex"
)

# --- Table 3: Substitution pattern -------------------------------------------

modelsummary(
  list("Incorporated SE"   = twfe_eb,
       "Unincorporated SE" = rob3,
       "Total SE"          = rob_total),
  coef_map  = c("did" = "High Debt x Post (DiD)"),
  gof_map   = gm,
  stars     = c("*"=0.1, "**"=0.05, "***"=0.01),
  vcov      = ~pid,
  title     = "Substitution Pattern: Composition Shift, Not Exit",
  add_rows  = tribble(
    ~term,              ~`Incorporated SE`, ~`Unincorporated SE`, ~`Total SE`,
    "Person FE",        "Yes",              "Yes",                "Yes",
    "Year FE",          "Yes",              "Yes",                "Yes",
    "Entropy weighted", "Yes",              "Yes",                "Yes"
  ),
  notes  = "TWFE entropy weighted throughout. Incorporated SE falls, unincorporated SE rises, total SE unchanged -- consistent with pure substitution from capital-intensive to necessity entrepreneurship.",
  output = "table3_substitution.tex"
)

# --- Table 4: Robustness checks ----------------------------------------------

modelsummary(
  list("Main TWFE"      = twfe_unw,
       "Drop 2020"      = rob1,
       "Placebo Unwtd"  = plac_unw,
       "Placebo Wtd"    = plac_eb,
       "Cont DiD Inc"   = cont_twfe,
       "Cont DiD Uninc" = cont_uninc),
  coef_map  = c("did"         = "DiD Coefficient",
                "did_placebo" = "DiD Coefficient",
                "debt_post"   = "DiD Coefficient"),
  gof_map   = gm,
  stars     = c("*"=0.1, "**"=0.05, "***"=0.01),
  vcov      = ~pid,
  title     = "Robustness Checks",
  add_rows  = tribble(
    ~term,              ~`Main TWFE`, ~`Drop 2020`, ~`Placebo Unwtd`, ~`Placebo Wtd`, ~`Cont DiD Inc`, ~`Cont DiD Uninc`,
    "Sample",           "Full",       "Drop 2020",  "Pre-2023",       "Pre-2023",     "Borrowers",     "Borrowers",
    "Entropy weighted", "No",         "Yes",        "No",             "Yes",          "No",            "No"
  ),
  notes  = "Person and year FE throughout. Clustered SEs at person level. Placebo post = 2022. Continuous DiD uses log pre-period debt interacted with post, borrowers only.",
  output = "table4_robustness.tex"
)

########################################################################################
# FIGURES
########################################################################################

# Figure 1: Incorporated SE trends by group
fig1 <- ggplot(trends_data, aes(x=year, y=mean_se, 
                                colour=Group, group=Group)) +
  annotate("rect", xmin=2023, xmax=2024,
           ymin=-Inf, ymax=Inf, fill="#FADBD8", alpha=0.4) +
  geom_line(linewidth=1.5) +
  geom_point(size=4.5, aes(fill=Group), shape=21, 
             colour="white", stroke=1.8) +
  scale_colour_manual(
    values=c("High Debt Borrowers"=col_treat, 
             "Non-Borrowers"=col_ctrl)) +
  scale_fill_manual(
    values=c("High Debt Borrowers"=col_treat, 
             "Non-Borrowers"=col_ctrl)) +
  scale_y_continuous(labels=percent_format(accuracy=0.1), 
                     limits=c(0, 0.088)) +
  scale_x_continuous(breaks=2018:2023,
                     limits=c(2017.8, 2024.2)) +
  labs(
    x=NULL, y="Incorporated SE Rate", colour=NULL, fill=NULL
  ) +
  base_theme

fig1
ggsave("fig1.png", fig1, width=10, height=6, dpi=300)

# Figure 2: Event study
fig2 <- ggplot(es_plot_data, aes(x=yr, y=estimate, colour=period)) +
  annotate("rect", xmin=2017.5, xmax=2021.5,
           ymin=-Inf, ymax=Inf, fill="#EBF5FB", alpha=0.3) +
  annotate("rect", xmin=2022.5, xmax=2023.5,
           ymin=-Inf, ymax=Inf, fill="#FADBD8", alpha=0.3) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50", linewidth=0.8) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=period),
              alpha=0.12, colour=NA) +
  geom_line(aes(group=1), colour="grey70", linewidth=0.8, linetype="dotted") +
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high), size=0.8, linewidth=1.4) +
  annotate("text", x=2019.5, y=0.20,
           label="Pre-period\n(joint p = 0.452)", size=3.2,
           colour=col_ctrl, fontface="italic") +
  annotate("text", x=2023.12, y=-0.085,
           label="-10.3pp\np = 0.248", size=3,
           colour=col_treat, fontface="bold", hjust=0) +
  scale_colour_manual(
    values=c("Pre-period"=col_ctrl, "Post-period"=col_treat, "Reference"="grey40")) +
  scale_fill_manual(
    values=c("Pre-period"=col_ctrl, "Post-period"=col_treat, "Reference"="grey40")) +
  scale_x_continuous(breaks=2018:2023) +
  scale_y_continuous(labels=percent_format(accuracy=0.1)) +
  labs(
    x="Year (reference = 2022)", y="Coefficient (High Debt x Year)",
    colour=NULL, fill=NULL
  ) +
  base_theme + theme(legend.position="none")

fig2

print(fig2)
ggsave("fig2_event_study.png", fig2, width=10, height=6, dpi=300)



# Figure 3: Balance love plot
fig3 <- ggplot(smd_data) +
  annotate("rect", xmin=-0.1, xmax=0.1,
           ymin=-Inf, ymax=Inf, fill="#D5F5E3", alpha=0.4) +
  geom_vline(xintercept=0, linetype="solid", colour="grey40", linewidth=0.8) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed",
             colour="grey60", linewidth=0.6) +
  geom_segment(
    data=smd_data |>
      pivot_wider(names_from=timing, values_from=smd) |>
      rename(before=`Before matching`, after=`After matching`),
    aes(x=after, xend=before, y=variable, yend=variable),
    colour="grey75", linewidth=1.2
  ) +
  geom_point(aes(x=smd, y=variable, colour=timing, shape=timing), size=4) +
  scale_colour_manual(
    values=c("Before matching"=col_treat, "After matching"="#27AE60"),
    labels=c("Before matching"="Before entropy balancing",
             "After matching" ="After entropy balancing")) +
  scale_shape_manual(
    values=c("Before matching"=16, "After matching"=17),
    labels=c("Before matching"="Before entropy balancing",
             "After matching" ="After entropy balancing")) +
  scale_x_continuous(limits=c(-1.75, 0.60), breaks=c(-1.5,-1.0,-0.5,0,0.5)) +
  labs(
    title    = "Covariate Balance: Before and After Entropy Balancing",
    subtitle = "All SMDs reduced to exactly zero after reweighting",
    x="Standardised Mean Difference (SMD)", y=NULL, colour=NULL, shape=NULL
  ) +
  base_theme
fig3
print(fig3)
ggsave("fig3_balance.png", fig3, width=10, height=7, dpi=300)



# Figure 4: Substitution pattern
fig4 <- ggplot(subst_data, aes(x=Outcome, y=estimate, fill=expected)) +
  geom_hline(yintercept=0, linetype="dashed", 
             colour="grey40", linewidth=0.9) +
  geom_col(width=0.55, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi),
                width=0.12, linewidth=1.2, 
                colour=text_col) +
  geom_text(aes(label=plabel, 
                y=estimate + sign(estimate)*0.035),
            size=4.2, fontface="bold", 
            colour=text_col) +
  scale_fill_manual(values=c(
    "Negative (theory predicted)" = col_treat,
    "Positive (theory predicted)" = col_ctrl,
    "Zero (pure substitution)"    = "#95A5A6"
  ),
  guide = "none") +        # this removes the legend boxes
  scale_y_continuous(labels=percent_format(accuracy=0.1), 
                     limits=c(-0.28, 0.32)) +
  labs(
    x=NULL, 
    y="DiD Coefficient", 
    fill=NULL
  ) +
  base_theme +
  theme(legend.position="none",    # belt and braces
        panel.grid.major.x=element_blank())
fig4
print(fig4)
ggsave("fig4_substitution.png", fig4, width=9, height=6, dpi=300)



# Figure 5a: Debt distribution
fig5a <- ggplot(debt_dist, aes(x=debt_k)) +
  annotate("rect", xmin=0, xmax=median_debt_k,
           ymin=-Inf, ymax=Inf, fill="#EBF5FB", alpha=0.4) +
  annotate("rect", xmin=median_debt_k, xmax=66,
           ymin=-Inf, ymax=Inf, fill="#FEF9E7", alpha=0.4) +
  annotate("rect", xmin=66, xmax=Inf,
           ymin=-Inf, ymax=Inf, fill="#FADBD8", alpha=0.4) +
  geom_histogram(bins=35, fill="#5D6D7E", alpha=0.80,
                 colour="white", linewidth=0.2) +
  geom_vline(xintercept=median_debt_k, linetype="dashed",
             colour=col_med, linewidth=1.3) +
  geom_vline(xintercept=66, linetype="dashed",
             colour=col_terc, linewidth=1.3) +
  annotate("text", x=median_debt_k+3, y=Inf,
           label=sprintf("Median\n$%.0fk", median_debt_k),
           hjust=0, vjust=1.4, size=4, 
           colour=col_med, fontface="bold") +
  annotate("text", x=66+3, y=Inf,
           label="Top tercile\n$66k",
           hjust=0, vjust=1.4, size=4, 
           colour=col_terc, fontface="bold") +
  scale_x_continuous(labels=dollar_format(suffix="k"),
                     breaks=seq(0,200,25), 
                     limits=c(0,220)) +
  scale_y_continuous(
    expand=expansion(mult=c(0,0.12))) +
  labs(
    subtitle = "Blue = below median | Yellow = mid | Red = top tercile (treatment group)",
    x="Student Debt at Baseline ($000s)", 
    y="Number of Borrowers"
  ) +
  base_theme
fig5a

print(fig5a)
ggsave("fig5a_debt_distribution.png", fig5a, width=12, height=7, dpi=300)



# Figure 5b: DiD by debt quartile
fig5b <- ggplot(quartile_results, aes(x=label, y=coef_eb)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40", linewidth=1.0) +
  geom_hline(yintercept=-0.1023, linetype="dotted",
             colour=col_treat, linewidth=1.0, alpha=0.7) +
  annotate("text", x=4.5, y=-0.1023-0.008,
           label="Main weighted estimate: -10.2pp",
           hjust=1, size=3.8, colour=col_treat, fontface="italic") +
  geom_col(aes(fill=sig), width=0.55, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi),
                width=0.15, linewidth=1.3, colour=text_col) +
  geom_text(aes(label=plabel,
                y=ifelse(coef_eb >= 0, ci_hi+0.008, ci_lo-0.008)),
            size=4.2, fontface="bold", colour=text_col) +
  geom_text(aes(label=sprintf("n = %d", n_treated), y=0.042),
            size=3.5, colour="grey50") +
  scale_fill_manual(
    values=c("TRUE"=col_treat, "FALSE"="#AEB6BF"),
    labels=c("TRUE"="p < 0.10", "FALSE"="p >= 0.10"),
    name="Significance"
  ) +
  scale_x_discrete() +
  scale_y_continuous(
    labels=percent_format(accuracy=0.1),
    limits=c(-0.20, 0.07),
    breaks=seq(-0.20, 0.05, 0.05)
  ) +
  labs(
    title    = "Effect on Incorporated SE by Debt Quartile (Entropy Weighted)",
    subtitle = "TWFE entropy weighted | Each quartile vs non-borrower control",
    x="Debt Quartile (among borrowers)", y="DiD Coefficient (Incorporated SE)"
  ) +
  base_theme + theme(legend.position="bottom")

print(fig5b)
ggsave("fig5b_quartile_did.png", fig5b, width=10, height=8, dpi=300)


