#Seminar paper – Universität Regensburg
# Author: Aditya Singh Negi
# - Data expected in ./R_Data
# - Figures saved to ../latexGraphics
# - Tables saved to ../latexTables

## Clear Working Space
rm(list = ls())

## Set Working Directory: run this script from its code folder
get_script_dir <- function(){
  if (requireNamespace('rstudioapi', quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error=function(e) '')
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
setwd(get_script_dir())


# 0) Packages
libs <- c("dplyr","tidyr","stringr","lubridate","readr","fixest","ggplot2","broom","readxl","xtable")
to_install <- libs[!sapply(libs, requireNamespace, quietly = TRUE)]
if(length(to_install) > 0) install.packages(to_install)
invisible(lapply(libs, library, character.only = TRUE))

#Template folders
# Input data should be placed in ./R_Data
# Figures saved to ../latexGraphics and tables to ../latexTables
dir.create('R_Data', showWarnings = FALSE, recursive = TRUE)
dir.create('../latexGraphics', showWarnings = FALSE, recursive = TRUE)
dir.create('../latexTables', showWarnings = FALSE, recursive = TRUE)

TABLE_DIR <- '../latexTables'
GRAPHICS_DIR <- '../latexGraphics'


# 1) Paths

DATA_FILE <- file.path('R_Data','panel_emp_monthly_with_wdi_macro.csv')

OUT_DIR <- 'Seminar_output'  # local results folder
# Template output folders are TABLE_DIR and GRAPHICS_DIR
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

# Helpers: write text outputs
write_txt_safe <- function(df, path){
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  message("Saved: ", path)
}

write_capture_txt <- function(obj, path){
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  capture.output(obj, file = path)
  message("Saved: ", path)
}


# 2) Load final merged dataset
#'/Users/adityanegi/Desktop/Genesis Dataset/My Dataset/panel_emp_monthly_with_wdi_macro.csv'
stopifnot(file.exists(DATA_FILE))

df <- read.csv(DATA_FILE, stringsAsFactors = FALSE) %>%
  mutate(
    date = as.Date(date),
    code = as.factor(code),
    year = as.integer(year(date)),
    moy  = lubridate::month(date),
    # global monthly time index (useful for trends)
    t = as.integer((year(date) - min(year(date), na.rm = TRUE)) * 12 + month(date))
  )

# Core event date
war_start <- as.Date("2022-03-01")

df <- df %>%
  mutate(
    war = as.integer(date >= war_start),
    # relative month index for event study
    m_id   = year(date) * 12 + month(date),
    m_war0 = year(war_start) * 12 + month(war_start),
    t_war  = m_id - m_war0
  )

# Ensure core outcomes exist (if missing, create from raw)
if(!"ln_prod" %in% names(df) && "prod_index" %in% names(df)){
  df <- df %>% mutate(ln_prod = ifelse(prod_index > 0, log(prod_index), NA_real_))
}
if(!"ln_emp_persons" %in% names(df) && "emp_persons" %in% names(df)){
  df <- df %>% mutate(ln_emp_persons = ifelse(emp_persons > 0, log(emp_persons), NA_real_))
}
if(!"ln_emp_hours" %in% names(df) && "emp_hours_worked" %in% names(df)){
  df <- df %>% mutate(ln_emp_hours = ifelse(emp_hours_worked > 0, log(emp_hours_worked), NA_real_))
}

# Ensure energy_std exists (if not, create it)
if(!"energy_std" %in% names(df) && "energy_share" %in% names(df)){
  df <- df %>% mutate(energy_std = as.numeric(scale(energy_share)))
}

# Basic sample restriction (manufacturing WZ 10-33) if div2 exists
if("div2" %in% names(df)){
  df <- df %>% filter(between(div2, 10, 33))
}


# 3) Main WAR-only FE models

m_prod_war <- feols(
  ln_prod ~ war:energy_std | code + date,
  data = df,
  cluster = ~code
)
m_emp_war <- feols(
  ln_emp_persons ~ war:energy_std | code + date,
  data = df,
  cluster = ~code
)
m_hours_war <- feols(
  ln_emp_hours ~ war:energy_std | code + date,
  data = df,
  cluster = ~code
)

# Save regression table (text)
tab_main <- etable(m_prod_war, m_emp_war, m_hours_war, tex = FALSE)
write_capture_txt(tab_main, file.path(TABLE_DIR, 'table_main_war_only.txt'))

# Save tidy coef table (text)
tidy_main <- bind_rows(
  broom::tidy(m_prod_war)  %>% mutate(model = "m_prod_war"),
  broom::tidy(m_emp_war)   %>% mutate(model = "m_emp_war"),
  broom::tidy(m_hours_war) %>% mutate(model = "m_hours_war")
)
write_txt_safe(tidy_main, file.path(TABLE_DIR, 'main_models_tidy.txt'))


# 4) Event-study DiD (energy heterogeneity)

WIN_PRE  <- -24
WIN_POST <-  36

df_es <- df %>% filter(t_war >= WIN_PRE, t_war <= WIN_POST)

es_war_prod <- feols(
  ln_prod ~ i(t_war, energy_std, ref = -1) | code + date,
  data = df_es,
  cluster = ~code
)
es_war_emp <- feols(
  ln_emp_persons ~ i(t_war, energy_std, ref = -1) | code + date,
  data = df_es,
  cluster = ~code
)
es_war_hours <- feols(
  ln_emp_hours ~ i(t_war, energy_std, ref = -1) | code + date,
  data = df_es,
  cluster = ~code
)

# Save plots as PDF (instead of PNG)
pdf(file.path(GRAPHICS_DIR, 'eventstudy_war_ln_prod.pdf'), width = 10, height = 7)
iplot(es_war_prod,
      main = "WAR event study: ln(prod) x energy intensity",
      xlab = "Months relative to war start (t=0 at 2022-03)",
      ylab = "Effect per 1 SD energy intensity")
abline(v = 0, lty = 2)
dev.off()

pdf(file.path(GRAPHICS_DIR, 'eventstudy_war_ln_emp.pdf'), width = 10, height = 7)
iplot(es_war_emp,
      main = "WAR event study: ln(employment) x energy intensity",
      xlab = "Months relative to war start (t=0 at 2022-03)",
      ylab = "Effect per 1 SD energy intensity")
abline(v = 0, lty = 2)
dev.off()

pdf(file.path(GRAPHICS_DIR, 'eventstudy_war_ln_hours.pdf'), width = 10, height = 7)
iplot(es_war_hours,
      main = "WAR event study: ln(hours) x energy intensity",
      xlab = "Months relative to war start (t=0 at 2022-03)",
      ylab = "Effect per 1 SD energy intensity")
abline(v = 0, lty = 2)
dev.off()


# 5) Parallel-trends / joint tests (ROBUST) 

parse_wald_out <- function(w){
  if(is.null(w)) return(list(stat=NA_real_, df1=NA_real_, df2=NA_real_, p=NA_real_))
  if(inherits(w, "htest")){
    return(list(
      stat = as.numeric(unname(w$statistic)),
      df1  = NA_real_,
      df2  = NA_real_,
      p    = as.numeric(w$p.value)
    ))
  }
  if(is.list(w)){
    stat <- suppressWarnings(as.numeric(w$stat %||% w$statistic))
    df1  <- suppressWarnings(as.numeric(w$df1))
    df2  <- suppressWarnings(as.numeric(w$df2))
    p    <- suppressWarnings(as.numeric(w$p %||% w$p.value))
    return(list(stat=stat, df1=df1, df2=df2, p=p))
  }
  if(is.atomic(w) && is.numeric(w)){
    if(!is.null(names(w))){
      stat <- suppressWarnings(as.numeric(w["stat"]))
      df1  <- suppressWarnings(as.numeric(w["df1"]))
      df2  <- suppressWarnings(as.numeric(w["df2"]))
      p    <- suppressWarnings(as.numeric(w["p"]))
      if(any(is.na(c(stat,df1,df2,p))) && length(w) >= 4){
        stat <- w[1]; df1 <- w[2]; df2 <- w[3]; p <- w[4]
      }
      return(list(stat=stat, df1=df1, df2=df2, p=p))
    } else {
      if(length(w) >= 4){
        return(list(stat=w[1], df1=w[2], df2=w[3], p=w[4]))
      } else {
        return(list(stat=NA_real_, df1=NA_real_, df2=NA_real_, p=NA_real_))
      }
    }
  }
  list(stat=NA_real_, df1=NA_real_, df2=NA_real_, p=NA_real_)
}

`%||%` <- function(a, b) if(!is.null(a)) a else b

coef_in_window <- function(model, k_min, k_max, x="energy_std", var="t_war"){
  cn <- names(coef(model))
  pattern <- paste0("^", var, "::(-?\\d+):", x, "$")
  keep <- cn[grepl(pattern, cn)]
  k <- suppressWarnings(as.integer(sub(paste0("^", var, "::"), "", sub(paste0(":", x, "$"), "", keep))))
  keep[!is.na(k) & k >= k_min & k <= k_max]
}

joint_test_window_fixest <- function(model, k_min, k_max, x="energy_std", var="t_war"){
  keep <- coef_in_window(model, k_min, k_max, x=x, var=var)
  if(length(keep) == 0){
    return(data.frame(window=paste0(k_min,"..",k_max), stat=NA, df1=NA, df2=NA, p_value=NA, n_coef=0))
  }
  hyps <- paste0(keep, " = 0")
  w <- tryCatch(fixest::wald(model, hyps), error=function(e) NULL)
  ww <- parse_wald_out(w)
  data.frame(
    window = paste0(k_min, "..", k_max),
    stat = ww$stat,
    df1 = ww$df1,
    df2 = ww$df2,
    p_value = ww$p,
    n_coef = length(keep)
  )
}

# Pre and Post tests (baseline ES)
pre_prod  <- joint_test_window_fixest(es_war_prod,  -24, -2)
pre_emp   <- joint_test_window_fixest(es_war_emp,   -24, -2)
pre_hours <- joint_test_window_fixest(es_war_hours, -24, -2)

post_prod  <- joint_test_window_fixest(es_war_prod,  0, 36)
post_emp   <- joint_test_window_fixest(es_war_emp,   0, 36)
post_hours <- joint_test_window_fixest(es_war_hours, 0, 36)

wald_table_baseline <- bind_rows(
  pre_prod   %>% mutate(outcome="ln_prod",        test="pre"),
  pre_emp    %>% mutate(outcome="ln_emp_persons", test="pre"),
  pre_hours  %>% mutate(outcome="ln_emp_hours",   test="pre"),
  post_prod  %>% mutate(outcome="ln_prod",        test="post"),
  post_emp   %>% mutate(outcome="ln_emp_persons", test="post"),
  post_hours %>% mutate(outcome="ln_emp_hours",   test="post")
)

write_txt_safe(wald_table_baseline, file.path(TABLE_DIR, 'wald_joint_tests_baseline.txt'))


# 6) Fix for significant pre-trends: sector-specific linear trends 
#    Correct varying-slopes syntax: code[t]

df_es <- df_es %>%
  mutate(
    code_fe = as.factor(as.character(code)),
    date_fe = as.factor(as.character(date)),
    t = as.numeric(t)   # slope variable MUST be numeric
  )

es_trend_prod <- feols(
  ln_prod ~ i(t_war, energy_std, ref = -1) | code_fe + date_fe + code_fe[t],
  data = df_es,
  cluster = ~code_fe
)

es_trend_emp <- feols(
  ln_emp_persons ~ i(t_war, energy_std, ref = -1) | code_fe + date_fe + code_fe[t],
  data = df_es,
  cluster = ~code_fe
)

es_trend_hours <- feols(
  ln_emp_hours ~ i(t_war, energy_std, ref = -1) | code_fe + date_fe + code_fe[t],
  data = df_es,
  cluster = ~code_fe
)

# Joint tests (trend-adjusted)
pre_prod_T  <- joint_test_window_fixest(es_trend_prod,  -24, -2)
pre_emp_T   <- joint_test_window_fixest(es_trend_emp,   -24, -2)
pre_hours_T <- joint_test_window_fixest(es_trend_hours, -24, -2)

post_prod_T  <- joint_test_window_fixest(es_trend_prod,  0, 36)
post_emp_T   <- joint_test_window_fixest(es_trend_emp,   0, 36)
post_hours_T <- joint_test_window_fixest(es_trend_hours, 0, 36)

wald_table_trends <- bind_rows(
  pre_prod_T   %>% mutate(outcome="ln_prod",        test="pre"),
  pre_emp_T    %>% mutate(outcome="ln_emp_persons", test="pre"),
  pre_hours_T  %>% mutate(outcome="ln_emp_hours",   test="pre"),
  post_prod_T  %>% mutate(outcome="ln_prod",        test="post"),
  post_emp_T   %>% mutate(outcome="ln_emp_persons", test="post"),
  post_hours_T %>% mutate(outcome="ln_emp_hours",   test="post")
)

write_txt_safe(wald_table_trends, file.path(TABLE_DIR, 'wald_joint_tests_with_sector_trends.txt'))

# Trend-adjusted plots as PDF
pdf(file.path(GRAPHICS_DIR, 'eventstudy_war_ln_prod_with_trends.pdf'), width = 10, height = 7)
iplot(es_trend_prod,
      main="WAR event study (with sector trends): ln(prod) x energy intensity",
      xlab="Months relative to war start",
      ylab="Effect per 1 SD energy intensity")
abline(v=0, lty=2)
dev.off()

pdf(file.path(GRAPHICS_DIR, 'eventstudy_war_ln_emp_with_trends.pdf'), width = 10, height = 7)
iplot(es_trend_emp,
      main="WAR event study (with sector trends): ln(employment) x energy intensity",
      xlab="Months relative to war start",
      ylab="Effect per 1 SD energy intensity")
abline(v=0, lty=2)
dev.off()

pdf(file.path(GRAPHICS_DIR, 'eventstudy_war_ln_hours_with_trends.pdf'), width = 10, height = 7)
iplot(es_trend_hours,
      main="WAR event study (with sector trends): ln(hours) x energy intensity",
      xlab="Months relative to war start",
      ylab="Effect per 1 SD energy intensity")
abline(v=0, lty=2)
dev.off()


# 7) Macro war patterns (GDP, investment, capital) from merged WDI columns

macro_vars <- c("gdp_growth_pct","gdp_real_2015usd","gfcf_real_2015usd","gfcf_share_gdp",
                "K_pim","ln_gdp","ln_gfcf","ln_K_pim")
macro_vars <- macro_vars[macro_vars %in% names(df)]

macro_year <- df %>%
  distinct(year, .keep_all = TRUE) %>%
  select(year, all_of(macro_vars)) %>%
  arrange(year) %>%
  mutate(
    war_y = as.integer(year >= 2022),
    covid_y = as.integer(year >= 2020)
  )

# Save macro data as TEXT
write_txt_safe(macro_year, file.path(TABLE_DIR, 'macro_series_by_year.txt'))

# Plot GDP growth and investment share if available (PDF)
if("gdp_growth_pct" %in% names(macro_year)){
  p <- ggplot(macro_year, aes(x=year, y=gdp_growth_pct)) +
    geom_line() +
    geom_vline(xintercept=2022, linetype="dashed") +
    labs(title="Germany GDP growth (WDI), war break at 2022", x="Year", y="GDP growth (%)")
  ggsave(file.path(GRAPHICS_DIR, 'macro_gdp_growth.pdf'), p, width=9, height=5)
}
if("gfcf_share_gdp" %in% names(macro_year)){
  p <- ggplot(macro_year, aes(x=year, y=gfcf_share_gdp)) +
    geom_line() +
    geom_vline(xintercept=2022, linetype="dashed") +
    labs(title="Germany GFCF share of GDP (WDI), war break at 2022", x="Year", y="GFCF share (%)")
  ggsave(file.path(GRAPHICS_DIR, 'macro_gfcf_share.pdf'), p, width=9, height=5)
}
if("ln_K_pim" %in% names(macro_year)){
  p <- ggplot(macro_year, aes(x=year, y=ln_K_pim)) +
    geom_line() +
    geom_vline(xintercept=2022, linetype="dashed") +
    labs(title="Germany capital stock (PIM from WDI GFCF), war break at 2022", x="Year", y="log(K)")
  ggsave(file.path(GRAPHICS_DIR, 'macro_capital_lnK.pdf'), p, width=9, height=5)
}

# Simple (descriptive) break regressions – tiny N, interpret cautiously
macro_reg_results <- list()

if("ln_gdp" %in% names(macro_year)){
  macro_reg_results$gdp <- summary(feols(ln_gdp ~ war_y + covid_y, data=macro_year))
  write_capture_txt(
    macro_reg_results$gdp,
    file.path(TABLE_DIR, 'macro_reg_ln_gdp.txt')
  )
}

if("ln_gfcf" %in% names(macro_year)){
  macro_reg_results$gfcf <- summary(feols(ln_gfcf ~ war_y + covid_y, data=macro_year))
  write_capture_txt(
    macro_reg_results$gfcf,
    file.path(TABLE_DIR, 'macro_reg_ln_gfcf.txt')
  )
}

if("ln_K_pim" %in% names(macro_year)){
  macro_reg_results$K <- summary(feols(ln_K_pim ~ war_y + covid_y, data=macro_year))
  write_capture_txt(
    macro_reg_results$K,
    file.path(TABLE_DIR, 'macro_reg_ln_K.txt')
  )
}


# Pre vs war averages
macro_summary <- macro_year %>%
  mutate(period = ifelse(year <= 2021, "prewar_2015_2021", "war_2022_2025")) %>%
  group_by(period) %>%
  summarise(across(all_of(intersect(c("gdp_growth_pct","ln_gdp","ln_gfcf","gfcf_share_gdp","ln_K_pim"), names(.))),
                   ~mean(.x, na.rm=TRUE)),
            .groups="drop")

write_txt_safe(macro_summary, file.path(TABLE_DIR, 'macro_prewar_vs_war_summary.txt'))


# 8) Save model objects for reproducibility

saveRDS(
  list(
    m_prod_war=m_prod_war, m_emp_war=m_emp_war, m_hours_war=m_hours_war,
    es_war_prod=es_war_prod, es_war_emp=es_war_emp, es_war_hours=es_war_hours,
    es_trend_prod=es_trend_prod, es_trend_emp=es_trend_emp, es_trend_hours=es_trend_hours,
    wald_baseline=wald_table_baseline, wald_trends=wald_table_trends,
    macro_year=macro_year
  ),
  file = file.path(TABLE_DIR, 'Seminar_models_and_outputs.rds')
)

message("\nDONE. All outputs saved in: ", OUT_DIR)


# 9) German-Russia Export shock (Exports = sheet 1, Imports = sheet 2)
TRADE_FILE <- file.path('R_Data','Ger_Export_Russ.xlsx')      
war_start <- as.Date("2022-03-01")
parse_num <- function(x){
  x <- as.character(x)
  x <- na_if(x, ":")
  readr::parse_number(x)
}

# Read exports (Sheet 1)
stopifnot(file.exists(TRADE_FILE))
exp_raw <- readxl::read_excel(TRADE_FILE, sheet = 1)

exports_ru <- exp_raw %>%
  filter(str_detect(as.character(Time), "^\\d{4}-\\d{2}$")) %>%
  transmute(
    date = as.Date(paste0(Time, "-01")),
    exports_ru = parse_num(`Value`)
  )
# Imports 
imp_raw <- readxl::read_excel(TRADE_FILE, sheet = 2)

imports_ru <- imp_raw %>%
  filter(str_detect(as.character(Time), "^\\d{4}-\\d{2}$")) %>%
  transmute(
    date = as.Date(paste0(Time, "-01")),
    imports_ru = parse_num(`Value`)
  )
# Merge exports & imports
trade_df <- exports_ru %>%
  left_join(imports_ru, by = "date") %>%
  arrange(date) %>%
  mutate(net_ru = exports_ru - imports_ru)

# Standardise using pre-war period
trade_pre <- trade_df %>%
  filter(date >= as.Date("2015-01-01"),
         date <= as.Date("2019-12-01"))

mu_ru <- mean(trade_pre$net_ru, na.rm = TRUE)
sd_ru <- sd(trade_pre$net_ru, na.rm = TRUE)

trade_df <- trade_df %>%
  mutate(ru_net_shock = (net_ru - mu_ru) / sd_ru)

#Merge into sector panel 
df <- df %>%
  left_join(trade_df %>% select(date, ru_net_shock), by = "date")

message("Russia trade shock merged. Missing values: ", sum(is.na(df$ru_net_shock)))

# 9B) Trade shock interactions (Table)
#     ru_net_shock × (energy intensity, foreign exposure)

if(!"exp_foreign" %in% names(df)){
  
  # Try to build it from turnover variables if available
  if(all(c("emp_turnover_th_eur","emp_turn_for_th_eur") %in% names(df))){
    
    df <- df %>%
      mutate(share_for = emp_turn_for_th_eur / emp_turnover_th_eur)
    
    exp_map <- df %>%
      filter(date < war_start) %>%
      group_by(code) %>%
      summarise(exp_foreign = mean(share_for, na.rm = TRUE), .groups = "drop")
    
    df <- df %>% left_join(exp_map, by = "code")
    
  } else if("share_for" %in% names(df)){
    
    exp_map <- df %>%
      filter(date < war_start) %>%
      group_by(code) %>%
      summarise(exp_foreign = mean(share_for, na.rm = TRUE), .groups = "drop")
    
    df <- df %>% left_join(exp_map, by = "code")
    
  } else {
    warning("Cannot construct exp_foreign: missing turnover/share_for variables.")
  }
}

#Estimate the trade-shock interaction models

if(all(c("ru_net_shock","energy_std","exp_foreign") %in% names(df))){
  
  df <- df %>% mutate(t_num = as.numeric(t))
  
  m_prod_trade <- feols(
    ln_prod ~ ru_net_shock:energy_std + ru_net_shock:exp_foreign | code + date + code[t_num],
    data = df,
    cluster = ~code
  )
  
  m_emp_trade <- feols(
    ln_emp_persons ~ ru_net_shock:energy_std + ru_net_shock:exp_foreign | code + date + code[t_num],
    data = df,
    cluster = ~code
  )
  
  m_hours_trade <- feols(
    ln_emp_hours ~ ru_net_shock:energy_std + ru_net_shock:exp_foreign | code + date + code[t_num],
    data = df,
    cluster = ~code
  )
  
# 3) Export Table 3
  tab_trade <- etable(
    m_prod_trade, m_emp_trade, m_hours_trade,
    tex = FALSE,
    dict = c(
      "ru_net_shock:energy_std"  = "ru_shock × energy_std",
      "ru_net_shock:exp_foreign" = "ru_shock × exp_foreign"
    )
  )
  
  write_capture_txt(tab_trade, file.path(TABLE_DIR, 'table_trade_shock_interactions.txt'))
  
# Also export a tidy coefficient table (useful for writing θ1, θ2 values)
  tidy_trade <- bind_rows(
    broom::tidy(m_prod_trade)  %>% mutate(model="prod"),
    broom::tidy(m_emp_trade)   %>% mutate(model="emp"),
    broom::tidy(m_hours_trade) %>% mutate(model="hours")
  )
  write_txt_safe(tidy_trade, file.path(TABLE_DIR, 'table_trade_shock_interactions_tidy.txt'))
  
} else {
  warning("Trade interaction models skipped: need ru_net_shock, energy_std, and exp_foreign.")
}



# Plot: aggregate exports to Russia over time (saved to latexGraphics)
if(requireNamespace('readxl', quietly = TRUE) && file.exists(TRADE_FILE)) {
  stopifnot(file.exists(TRADE_FILE))
exp_raw <- readxl::read_excel(TRADE_FILE, sheet = 1)
  exports_ru <- exp_raw %>%
    mutate(Time = as.character(Time)) %>%
    filter(stringr::str_detect(Time, '^\\d{4}-\\d{2}$')) %>%
    transmute(date = as.Date(paste0(Time, '-01')), exports_ru = parse_num(Value)) %>%
    arrange(date)

  p_exp <- ggplot(exports_ru, aes(x = date, y = exports_ru)) +
    geom_line() +
    geom_vline(xintercept = war_start, linetype = 'dashed') +
    labs(title = 'Aggregate German exports to Russia (monthly)',
         subtitle = 'Vertical line: March 2022 (war start)',
         x = 'Date', y = 'Exports (value)')

  ggsave(filename = file.path(GRAPHICS_DIR, 'figure_exports_russia_over_time.pdf'),
         plot = p_exp, width = 10, height = 5)
  message('Saved: ', file.path(GRAPHICS_DIR, 'figure_exports_russia_over_time.pdf'))
}

getwd()
list.files("Seminar_output", recursive = TRUE)
list.files("../latexTables", recursive = TRUE)
list.files("../latexGraphics", recursive = TRUE)
