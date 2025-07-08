



# ---- libs ----


library("data.table")
library("tictoc")
library("dplyr")
library("tidyr")
library("readr")
library("tibble")
library("stringr")
library("lubridate")
library("stringr")
library("ggplot2") 
library("ggthemes") 



# ---- lets_get_it ----



### read in raw data

use_test_dat <- FALSE
(file_to_use <- ifelse(use_test_dat, "ukb_test.tab", "ukb_real.tab"))
# files are not in this repo - securely stored elsewhere
file_to_use <- off_world_file(file_to_use)

datx <- fread(file = file_to_use, header = TRUE)


### note: processing/factor_encodings.R below was run over raw csv data provided
### by biobank it contains the variable level and label conversions via
### factor()/ordered() as well as date conversions via as.Date() etc

# source("processing/factor_encodings.R")


### now we select the columns we need which are named in a numeric biobank way

# speed up variable exploration with variable search and summary print
find_col <- function(search_str, tabulate = TRUE, fixed = TRUE, this_dat = ukbb_raw) {
  dat_cns <- colnames(this_dat)
  found_cols <- dat_cns[grepl(search_str, dat_cns, fixed = fixed)]
  if (!tabulate) {
    cat("\nVariables found: '", paste(found_cols, collapse = "', '"), "'\n", sep = "")
    return(invisible(found_cols))
  }
  if (tabulate & (length(found_cols) > 0)) {
    for (i in 1:length(found_cols)) {
      cat("\n\ntable of variable: '", found_cols[i], "'\n", sep = "")
      print(table(this_dat[[found_cols[i]]], useNA = "ifany"))
    }
  }
}


find_col_grep_datx <- function(search_str, tabulate = TRUE, fixed = TRUE) {
  find_col(search_str = search_str, tabulate = tabulate, fixed = fixed, this_dat = datx) 
}



kp_col <- list()

kp_col[['id']] <- find_col_grep_datx('^f\\.eid')    # id
kp_col[['sex']] <- find_col_grep_datx('^f\\.31\\.')    # Sex
kp_col[['yob']] <- find_col_grep_datx('^f\\.34\\.')    # Year of birth
kp_col[['mob']] <- find_col_grep_datx('^f\\.52\\.')    # Month of birth
kp_col[['dte_assess_ctre']] <- find_col_grep_datx('^f\\.53\\.')    # Date of attending assessment centre
kp_col[['day_walk_10p']] <- find_col_grep_datx('^f\\.864\\.')    # Number of days/week walked 10+ minutes
kp_col[['alc']] <- find_col_grep_datx('^f\\.1558\\.')    # Alcohol use
kp_col[['soc_iso']] <- find_col_grep_datx('^f\\.2020\\.')    # Social isolation (do you feel lonely?)
kp_col[['hearing']] <- find_col_grep_datx('^f\\.2247\\.')    # Hearing problems (do you have difficulty with your hearing?)
kp_col[['edu']] <- find_col_grep_datx('^f\\.6138\\.')    # Highest educational attainment
kp_col[['smoke']] <- find_col_grep_datx('^f\\.20116\\.')    # Smoking status
kp_col[['pairs_match']] <- find_col_grep_datx('^f\\.20132\\.')    # pairs matching: number of incorrect matches in round
kp_col[['dtm_pairs_match']] <- find_col_grep_datx('^f\\.20134\\.')    # When pairs test completed
kp_col[['dtm_fluid_intel']] <- find_col_grep_datx('^f\\.20135\\.')    # When fluid intelligence test completed
kp_col[['dtm_trail_make']] <- find_col_grep_datx('^f\\.20136\\.')    # When trail making test completed
kp_col[['dtm_symbol_digit']] <- find_col_grep_datx('^f\\.20137\\.')    # When symbol digit substitution test completed
kp_col[['dtm_num_mem']] <- find_col_grep_datx('^f\\.20138\\.')    # When numeric memory test completed 
kp_col[['trail_make_a']] <- find_col_grep_datx('^f\\.20156\\.')    # Trail making: duration to complete Trails A (numeric path)
kp_col[['trail_make_b']] <- find_col_grep_datx('^f\\.20157\\.')    # trail making: duration to complete Trails B (alphanumeric path)
kp_col[['symbol_digit_correct']] <- find_col_grep_datx('^f\\.20159\\.')    # Symbol digit substitutuion: number of symbol digit matches made correctly
kp_col[['fluid_intel']] <- find_col_grep_datx('^f\\.20191\\.')    # Fluid intelligence: fluid intelligence score global
kp_col[['symbol_digit_attempt']] <- find_col_grep_datx('^f\\.20195\\.')    # Symbol digit substitution: number of symbol digit matches attempted
kp_col[['dur_symbol_digit']] <- find_col_grep_datx('^f\\.20230\\.')    # Symbol digit substitution: duration to entering symbol choice
kp_col[['num_mem']] <- find_col_grep_datx('^f\\.20240\\.')    # Numeric memory: maximum digits remembered correctly
kp_col[['valid_fluid_intel']] <- find_col_grep_datx('^f\\.20242\\.')    # Fluid intelligence completion status
kp_col[['valid_pairs_match']] <- find_col_grep_datx('^f\\.20244\\.')    # Pairs matching completion status
kp_col[['valid_symbol_digit']] <- find_col_grep_datx('^f\\.20245\\.')    # Symbol digit substition completion status
kp_col[['valid_trail_make']] <- find_col_grep_datx('^f\\.20246\\.')    # Trail making completion status
kp_col[['ethn']] <- find_col_grep_datx('^f\\.21000\\.')    # Ethnicity
kp_col[['bmi']] <- find_col_grep_datx('^f\\.21001\\.')    # BMI
kp_col[['age_assess_ctre']] <- find_col_grep_datx('^f\\.21003\\.')    # Age when attended assessment centre
kp_col[['sl']] <- find_col_grep_datx('^f\\.40046\\.')    # Time in sleep
kp_col[['sb']] <- find_col_grep_datx('^f\\.40047\\.')    # Time in SB
kp_col[['lpa']] <- find_col_grep_datx('^f\\.40048\\.')    # Time in LPA
kp_col[['mvpa']] <- find_col_grep_datx('^f\\.40049\\.')    # Time in MVPA
kp_col[['dte_dem_all_cause']] <- find_col_grep_datx('^f\\.42018\\.')    # Date of all cause dementia report
kp_col[['src_dem_all_cause']] <- find_col_grep_datx('^f\\.42019\\.')    # Source of all cause dementia report
kp_col[['dte_alz']] <- find_col_grep_datx('^f\\.42020\\.')    # Date of alzheimer's disease report
kp_col[['src_alz']] <- find_col_grep_datx('^f\\.42021\\.')    # Source of alzheimer's disease report
kp_col[['dte_dem_vasc']] <- find_col_grep_datx('^f\\.42022\\.')    # date that vascular dementia diagnosis first reported
kp_col[['src_dem_vasc']] <- find_col_grep_datx('^f\\.42023\\.')    # Source of vascular dementia report
kp_col[['accel_start']] <- find_col_grep_datx('^f\\.90010\\.')    # Start time of accelerometer wear
kp_col[['valid_wear']] <- find_col_grep_datx('^f\\.90015\\.')    # Valid wear time (data quality)
kp_col[['data_prob_ind']] <- find_col_grep_datx('^f\\.90002\\.')    # Data problem indicator
kp_col[['depres']] <- find_col_grep_datx('^f\\.4598\\.')    # depression
kp_col[['hypten']] <- find_col_grep_datx('^f\\.6150\\.')    # hypertention
kp_col[['diabet']] <- find_col_grep_datx('^f\\.2443\\.')  # diabetes  
kp_col[['day_walk_10p']] <- find_col_grep_datx('^f\\.864\\.')  # https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=864



### this is subsetting columns to wanted only

kp_cols <- unlist(kp_cols)
names(kp_cols) <- NULL
kp_cols

col_nms <- colnames(datx)
(valid_wear_col <- col_nms[grepl("90015", col_nms)])
(rows_to_peek <- if(use_test_dat) { 1:9 } else { 1:30 })

# peaky blinders: having a look at the data momentarily
datx[rows_to_peek, ..kp_cols]
datx[, .N, by = c(valid_wear_col)]


# subset using data.table syntax (package deserves more recognition)
ukb_raw <- datx[, ..kp_cols]
nrow(ukb_raw)
ncol(ukb_raw)


format(object.size(ukb_raw), "GB")

rm(list = "datx") # remove initial import now we have a subsetted version
gc() # garbage clean-up: frees unused memory



### easier (time-wise) to do column renaming after the data is subsetted

dat_dic <- NULL
for (j in 1:length(kp_cols)) { # j <- 10
  from_ <- kp_cols[[j]]
  this_nc <- length(from_)
  to_ <- rep(names(kp_cols)[j], this_nc)
  if (this_nc > 1) {
    
    postfix_nos1 <- str_extract(from_, "([0-9]+)\\.([0-9]+)$", group = 1)
    max_char_postfix_nos1 <- max(nchar(postfix_nos1))
    postfix_nos1 <- 
      sprintf(
        paste0("%0", max_char_postfix_nos1, ".0f"), 
        as.integer(postfix_nos1)
      )
    
    postfix_nos2 <- str_extract(from_, "([0-9]+)\\.([0-9]+)$", group = 2)
    max_char_postfix_nos2 <- max(nchar(postfix_nos2))
    postfix_nos2 <- 
      sprintf(
        paste0("%0", max_char_postfix_nos2, ".0f"), 
        as.integer(postfix_nos2)
      )
    
    postfix_nos <- paste0("_", postfix_nos1, "_", postfix_nos2)
    to_ <- paste0(to_, postfix_nos)
  }
  
  names(to_) <- from_
  
  dat_dic <- c(dat_dic, to_)
  
}


dat_dic


ukb_raw <- ukb_raw %>% rename(., all_of(dat_dic))



format(object.size(ukb_raw), "GB") # "0.4 Gb"

gc() # garbage clean-up: frees unused memory


### can convert to "tibble" for tidyverse-ness instead of "data.frame" or "data.table"
tic()
ukb_raw <- as_tibble(ukb_raw)
toc() # takes < 1 sec

format(object.size(ukb_raw), "GB") # still "0.4 Gb"


### Remove duplicate rows

ukb_raw %>%
  group_by(id) %>% 
  summarise(n = n()) %>%
  ungroup() %>% 
  dplyr::filter(n > 1)  %>%
  inner_join(
    ukb_raw,
    .,
    "id"
  )

nrow(ukb_raw)
ukb_raw <-
  ukb_raw %>%
  group_by(id) %>% 
  # get the last row ("most recent entry") for each id
  dplyr::filter(row_number() == n()) %>%
  ungroup() 
nrow(ukb_raw)





# ---- remove_optouts ----

to_rm <- read_csv(off_world_file("opted_out_participants.csv"), col_names = FALSE)

(old_n <- nrow(ukb_raw))
ukb_raw <-
  anti_join(
    ukb_raw,
    to_rm,
    c("id" = "X1")
  )
(new_n <- nrow(ukb_raw))
cat("---- ids removed:", old_n - new_n, "\n\n")




# ---- explore ----

ukb_raw_cns <- colnames(ukb_raw)
# (ukb_raw_cns)

### note this is our include variable
# initialised as 0 = included
# values < 0 are to be exclude
ukb_raw[["exclude_ind"]] <- 0
with(ukb_raw, table(exclude_ind, useNA = "ifany"))




which_valid_lvls <- function(x, vlvls) {
  case_when(
    x %in% vlvls ~ TRUE,
    TRUE         ~ FALSE # catch other cases including NAs
  )
}

create_new_valid_var <- function(base_var, vlvls) {
  new_var <- 
    which_valid_lvls(
      ukb_raw[[base_var]], 
      vlvls = vlvls
    )
  cat("summary table:")
  print(table(new_var, useNA = "ifany"))
  return(new_var)
}



# ---- valid_wear ----

# valid wear col and summarise
with(ukb_raw, table(valid_wear, exclude_ind, useNA = "ifany"))

ukb_raw <-
  ukb_raw %>%
  mutate(
    exclude_ind =
      case_when(
        exclude_ind != 0   ~ exclude_ind, # don't change already non-zero values 
        is.na(valid_wear)  ~ 1,           # gotta catch 'em all (NA pokemons)
        valid_wear == "No" ~ 1,
        TRUE               ~ 0            # else keep
      )
  )

with(ukb_raw, table(exclude_ind, useNA = "ifany"))


# ---- data_prob_ind ----

# valid wear col and summarise
with(ukb_raw, table(data_prob_ind, exclude_ind, useNA = "ifany"))

rm_strs <- 
  c(
    "Dataset previously flagged as unreliable, now believed valid",
    "Unreliable due to unexpectedly small size",
    "Unreliable due to unexpectedly large size"
  )

ukb_raw <-
  ukb_raw %>%
  mutate(
    exclude_ind =
      case_when(
        exclude_ind != 0           ~ exclude_ind, # don't change already non-zero values 
        data_prob_ind %in% rm_strs ~  2,
        TRUE                       ~  0            # else keep (these are the NAs)
      )
  )

with(ukb_raw, table(exclude_ind, useNA = "ifany"))


# ---- cog_outcome_var_wrangle ----

## ---- num_mem ----


find_col("dtm_num_mem", tabulate = FALSE)

# see date range
range(ukb_raw[["dtm_num_mem"]], na.rm = TRUE)
# table of valid dates in num_mem
with(ukb_raw, table(!is.na(dtm_num_mem)))


ukb_raw[["valid_num_mem"]] <- !is.na(ukb_raw[["dtm_num_mem"]])
with(ukb_raw, table(valid_num_mem))


## ---- pairs_match ----

find_col("valid_pairs_match", tabulate = TRUE)
# table of pairs_match vs valid dates in num_mem
with(ukb_raw, table(valid_pairs_match, useNA = "ifany"))
with(ukb_raw, table(valid_pairs_match, valid_num_mem, useNA = "ifany"))


ukb_raw[["valid_pairs_match"]] <-
  create_new_valid_var("valid_pairs_match", "Completed") # Completed with pause not included

with(ukb_raw, table(valid_pairs_match, useNA = "ifany"))

colnames(ukb_raw)






## ---- fluid_intel ----

find_col("valid_fluid_intel", tabulate = TRUE)
# table of fluid_intel vs valid dates in num_mem
with(ukb_raw, table(valid_fluid_intel, useNA = "ifany"))
with(ukb_raw, table(valid_fluid_intel, valid_num_mem, useNA = "ifany"))


ukb_raw[["valid_fluid_intel"]] <-
  create_new_valid_var("valid_fluid_intel", "Completed")


with(ukb_raw, table(valid_fluid_intel, useNA = "ifany"))



## ---- trail_make ----

find_col("valid_trail_make", tabulate = TRUE)
# table of trail_make vs valid dates in num_mem
with(ukb_raw, table(valid_trail_make, useNA = "ifany"))
with(ukb_raw, table(valid_trail_make, valid_num_mem, useNA = "ifany"))


ukb_raw[["valid_trail_make"]] <-
  create_new_valid_var("valid_trail_make", "Completed")


with(ukb_raw, table(valid_trail_make, useNA = "ifany"))



## ---- symbol_digit ----

find_col("valid_symbol_digit", tabulate = FALSE)
# table of symbol_digit vs valid dates in num_mem
with(ukb_raw, table(valid_symbol_digit, useNA = "ifany"))
with(ukb_raw, table(valid_symbol_digit, valid_num_mem, useNA = "ifany"))


ukb_raw[["valid_symbol_digit"]] <-
  create_new_valid_var("valid_symbol_digit", "Completed")

with(ukb_raw, table(valid_symbol_digit, useNA = "ifany"))



# score summaries 

# 20240: Maximum digits remembered correctly
# Longest number correctly recalled during the numeric memory test. 
# A value of -1 is recorded if the participant chose to abandon the 
# test before completing the first round. 
find_col("20240")
find_col("20138", tabulate = FALSE)
range(ukb_raw[["dtm_num_mem"]], na.rm = TRUE)
find_col("f.52.") # month of birth
find_col("f.34.") # year of birth

# valid test scores
ukb_raw %>%
  dplyr::filter(valid_num_mem) %>%
  with(., table(f.20240.0.0, useNA = "ifany"))

ukb_raw %>%
  dplyr::filter(valid_num_mem) %>%
  pull(dtm_num_mem) %>%
  range(., na.rm = FALSE)





# ---- intersection_of_cog_test_completion ----


ukb_raw <-
  ukb_raw %>%
  mutate(
    cog_complete_4_of_5 =
      ((valid_num_mem + valid_pairs_match + valid_fluid_intel + 
          valid_trail_make + valid_symbol_digit) >= 4)
  )


with(ukb_raw, table(cog_complete_4_of_5, useNA = "ifany"))



# these are the cog test valid columns
cog_valid_cols <- 
  c(
    "valid_num_mem",
    "valid_pairs_match",
    "valid_fluid_intel",
    "valid_trail_make",
    "valid_symbol_digit"
  )


# peak
ukb_raw %>%
  select(all_of(cog_valid_cols))


# count of records with valid entry (TRUE/FALSE) combinations of each variable
(summ_valid <-
    ukb_raw %>%
    select(all_of(cog_valid_cols)) %>%
    # mutate(across(all_of(names(cog_cols1)), is_valid_dttm)) %>%
    # mutate(across(all_of(names(cog_cols2)), is_valid_fct)) %>%
    mutate(nvalid = rowSums(across(everything()))) %>%
    group_by(across(everything())) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(nvalid), desc(n)) %>%
    print(., n = nrow(.)))


summ_valid %>%
  dplyr::filter(nvalid >= 4) %>%
  pull(n) %>%
  sum(.)


with(ukb_raw, table(cog_complete_4_of_5, useNA = "ifany"))


### NOTE: limiting to 4+ valid cog records

with(ukb_raw, table(cog_complete_4_of_5, exclude_ind, useNA = "ifany"))

ukb_raw <-
  ukb_raw %>%
  mutate(
    exclude_ind =
      case_when(
        exclude_ind != 0           ~ exclude_ind, # don't change already non-zero values 
        cog_complete_4_of_5 ~ 0,
        TRUE                ~ 3            # else exclude (<= 3 cog complete)
      )
  )

with(ukb_raw, table(exclude_ind, useNA = "ifany"))



# ---- cog_test_dates ----



# 7590	20134-0.0	118496	Time	When pairs test completed
# 7591	20135-0.0	123581	Time	When fluid intelligence test completed
# 7592	20136-0.0	120427	Time	When trail making test completed
# 7593	20137-0.0	120143	Time	When symbol digit substitution test completed
# 7594	20138-0.0	111032	Time	When numeric memory test completed
date_cols <- 
  c(
    "dtm_pairs_match",
    "dtm_fluid_intel",
    "dtm_trail_make",
    "dtm_symbol_digit",
    "dtm_num_mem"
  )

dates_dat_tmp <-
  ukb_raw %>%
  select(id, all_of(date_cols)) %>%
  mutate(across(all_of(date_cols), as_date)) %>%
  mutate(
    max_d = pmax(
      dtm_pairs_match, dtm_fluid_intel, dtm_trail_make, 
      dtm_symbol_digit, dtm_num_mem, na.rm = TRUE
    ),
    min_d = pmin(
      dtm_pairs_match, dtm_fluid_intel, dtm_trail_make, 
      dtm_symbol_digit, dtm_num_mem, na.rm = TRUE
    ),
    max_diff = (min_d %--% max_d) / days(1)
  )

dates_dat_tmp
hist(dates_dat_tmp$max_diff)



# ---- assess_centre_dates ----

# find assess centre variables
assess_cols <- find_col("dte_assess_ctre", tabulate = FALSE)
length(assess_cols) # four of them

# look at date ranges
range(ukb_raw[[assess_cols[[1]]]], na.rm = TRUE)
range(ukb_raw[[assess_cols[[2]]]], na.rm = TRUE)
range(ukb_raw[[assess_cols[[3]]]], na.rm = TRUE)
range(ukb_raw[[assess_cols[[4]]]], na.rm = TRUE)

# peak
ukb_raw %>%
  select(all_of(assess_cols))

# count of records with valid entry (TRUE/FALSE) combinations of each variable
ukb_raw %>%
  select(all_of(assess_cols)) %>%
  mutate(across(everything(), ~!is.na(.))) %>%
  mutate(nvalid = rowSums(across(everything()))) %>%
  group_by(across(everything())) %>%
  summarise( n = n(), .groups = "drop") %>%
  arrange(desc(nvalid), desc(n))

# ---- dob_derivation ----


head(ukb_raw[["mob"]])
levels(ukb_raw[["mob"]]) <- substr(levels(ukb_raw[["mob"]]), 1, 3)
head(ukb_raw[["mob"]])
table(ukb_raw[["mob"]], useNA = "ifany")

# approx dob table
ukb_raw %>%
  select(yob, mob) %>%
  with(., table(yob, mob, useNA = "ifany"))

# set approx dob
ukb_raw <-
  ukb_raw %>%
  mutate(
    approx_dob_str = 
      paste0(yob, "-", sprintf("%02.0f", as.integer(mob)), "-15"),
    approx_dob = as_date(approx_dob_str)
  )

min(ukb_raw[["approx_dob_str"]])
min(ukb_raw[["approx_dob"]])
max(ukb_raw[["approx_dob_str"]])
max(ukb_raw[["approx_dob"]])

ukb_raw %>% sample_n(20) %>% select(yob, mob, approx_dob_str, approx_dob)
range(ukb_raw[["approx_dob"]])

age_at_nummem <-
  ukb_raw %>%
  # dtm_num_mem is nummem date
  mutate(age_at_nummem = (approx_dob %--% dtm_num_mem) / years(1)) %>%
  pull(age_at_nummem)

table(is.na(age_at_nummem), useNA = "ifany")

hist(age_at_nummem)

# ---- date_tests ----

# 7590	20134-0.0	118496	Time	When pairs test completed
# 7591	20135-0.0	123581	Time	When fluid intelligence test completed
# 7592	20136-0.0	120427	Time	When trail making test completed
# 7593	20137-0.0	120143	Time	When symbol digit substitution test completed
# 7594	20138-0.0	111032	Time	When numeric memory test completed
date_cols <- 
  c(
    "dte_assess_ctre_0_0",   # assess centre first visit
    "accel_start", # Start time of wear
    "dtm_pairs_match",
    "dtm_fluid_intel",
    "dtm_trail_make",
    "dtm_symbol_digit",
    "dtm_num_mem"
  )

dates_dat <-
  ukb_raw %>%
  select(id, all_of(date_cols)) %>%
  mutate(across(all_of(date_cols), as_date))

# this gets max date of all date columns
### NOTE: if diagnosis date > max(date cols) then meets exclude criteria
dates_dat <-
  dates_dat %>%
  pivot_longer(cols = -id, values_to = "date") %>%
  dplyr::filter(!is.na(date)) %>%
  arrange(id, desc(date)) 

dates_dat %>%
  ggplot(aes(date, fill = name)) +
  geom_histogram(alpha = 0.5, col = NA, bins = 90) +
  # geom_density(alpha = 0.5, col = NA) +
  facet_grid(name ~ ., scales = "free_y") +
  theme_bw() +
  scale_fill_colorblind()

dates_dat_accel <-
  dates_dat %>%
  dplyr::filter(name == "accel_start") %>%
  select(id, date_accel = date)


# want 0, no dups
dates_dat_accel %>%
  group_by(id) %>% 
  summarise(n = n()) %>%
  ungroup() %>% 
  dplyr::filter(n > 1) 


dates_dat_cog_strt <-
  dates_dat %>%
  dplyr::filter(!(name %in% c("dte_assess_ctre_0_0", "accel_start"))) %>%
  arrange(id, date) %>%
  group_by(id) %>% 
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  select(id, first_cog_date = date, first_cog_date_src = name)


dates_dat_max <-
  dates_dat %>%
  # allow all dates
  # dplyr::filter(!(name %in% c("dte_assess_ctre_0_0", "accel_start"))) %>%
  arrange(id, date) %>%
  group_by(id) %>% 
  dplyr::filter(row_number() == n()) %>%
  ungroup() %>%
  select(id, max_date = date, max_date_src = name)

nrow(ukb_raw)
nrow(dates_dat_accel)
nrow(dates_dat_cog_strt)
nrow(dates_dat_max)

dates_dat <-
  left_join(
    ukb_raw %>% select(id),
    dates_dat_accel,
    "id"
  ) %>%
  left_join(
    .,
    dates_dat_cog_strt,
    "id"
  ) %>%
  left_join(
    .,
    dates_dat_max,
    "id"
  ) 
nrow(dates_dat)

dates_dat

# anyone without any dates?
dates_dat %>%
  dplyr::filter(is.na(date_accel), is.na(first_cog_date), is.na(max_date)) 


with(dates_dat, table(is.na(date_accel), is.na(first_cog_date), useNA = "ifany"))

time_diff_dat <-
  dates_dat %>%
  mutate(
    accel_cog_diff = 
      if_else(
        is.na(date_accel) | is.na(first_cog_date),
        as.numeric(NA),
        (date_accel %--% first_cog_date) / days(1)
      )
  )

with(time_diff_dat, table(is.na(accel_cog_diff), useNA = "ifany"))

time_diff_dat <-
  time_diff_dat %>%
  select(-matches("^max")) %>%
  dplyr::filter(!is.na(accel_cog_diff))


time_diff_dat %>%
  mutate(yr_accel = as.character(year(date_accel))) %>%
  ggplot(aes(accel_cog_diff, fill = yr_accel, group = yr_accel)) +
  # geom_histogram(alpha = 0.5, col = NA, bins = 90) +
  geom_density(alpha = 0.5, col = NA) +
  # facet_grid( yr_accel ~ ., scales = "free_y") +
  theme_bw() +
  scale_fill_colorblind()


nrow(time_diff_dat)
nrow(ukb_raw)
ukb_raw <-
  left_join(
    ukb_raw ,
    time_diff_dat,
    "id"
  ) 
nrow(ukb_raw)

with(ukb_raw, table(is.na(accel_cog_diff), useNA = "ifany"))



### CALCULATE age at first cognitive completion date
# as this is the "start date" of the predictor variables in the time course


dates_dat_cog_strt %>%
  dplyr::filter(is.na(first_cog_date))



with(ukb_raw, table(is.na(first_cog_date), useNA = "ifany"))


ukb_raw <-
  ukb_raw %>%
  mutate(
    age_first_cog = 
      if_else(
        is.na(first_cog_date),
        as.numeric(NA),
        (approx_dob %--% first_cog_date) / years(1)
      )
  ) 

ukb_raw %>%
  select(id, approx_dob, first_cog_date, age_first_cog)


with(ukb_raw, table(is.na(age_first_cog), useNA = "ifany"))

with(ukb_raw, hist(age_first_cog))

ukb_raw %>%
  dplyr::filter(!is.na(age_first_cog)) %>%
  ggplot(aes(age_first_cog, fill = sex)) +
  geom_histogram(alpha = 0.5, col = NA, bins = 90) +
  # geom_density(alpha = 0.5, col = NA) +
  facet_grid(sex ~ ., scales = "free_y") +
  theme_bw() +
  scale_fill_colorblind()







# ---- icd10_processing1 ----


# F02 dementia
# F03
# F04
# F05
# G10-G14 atrophies primarily affecting the central nervous system 
# G10-G14
# G10-G14
# G10-G14
# G25 Other extrapyramidal and movement disorders
# G20-G26 Extrapyramidal and movement disorders
# G31 Other degenerative diseases of nervous system, not elsewhere classified
# G32 Other degenerative disorders of nervous system in diseases classified elsewhere
# G30-G32 Other degenerative diseases of the nervous system
# G35-G37 Demyelinating diseases of the central nervous system
# H54 Blindness and low vision


# files are not in this repo - securely stored elsewhere
icd10_f <- off_world_file("raw/icd10_raw_pivoted_longer.rds") 
# one dated icd10 per id per row:
ukb_raw_icd10 <- read_rds(file = icd10_f)



ukb_raw_icd10 <-
  ukb_raw_icd10 %>%
  mutate(icd10_block = str_sub(icd10, 1, 3))

with(ukb_raw_icd10, table(icd10_block, useNA = "ifany")) %>%
  as.data.frame() %>%
  as_tibble()


icd10_exclude <- 
  read_csv(
    "icd10_of_exclusion.csv", 
    col_types = cols(icd10_block = col_character())
  )

icd10_exclude <-
  icd10_exclude %>%
  distinct(icd10_block)

icd10_exclude


icd10_exclude <-
  ukb_raw_icd10 %>%
  dplyr::filter(type != "seco") %>%
  inner_join(
    .,
    icd10_exclude,
    "icd10_block"
  ) %>%
  arrange(id, type, seq, date)

nrow(icd10_exclude)
icd10_exclude

knitr::kable(data.frame(icd10codes = sort(unique(icd10_exclude$icd10_block))))


dates_max_to_join <- 
  dates_dat %>%
  select(id, max_date, max_date_src) %>%
  dplyr::filter(!is.na(max_date)) 
dates_max_to_join

incl_dat <- 
  ukb_raw %>%
  dplyr::filter(exclude_ind == 0) %>%
  select(id) 

nrow(incl_dat)
incl_dat <-
  left_join(
    incl_dat,
    dates_max_to_join,
    "id"
  )
nrow(incl_dat)


nrow(icd10_exclude)
icd10_exclude <-
  icd10_exclude %>%
  left_join(
    .,
    incl_dat,
    "id"
  ) 
nrow(icd10_exclude)

icd10_exclude <-
  icd10_exclude %>%
  dplyr::filter(!is.na(max_date))
nrow(icd10_exclude)

icd10_exclude <-
  icd10_exclude %>%
  mutate(
    after_period =
      if_else(
        is.na(max_date),
        NA,
        max_date < date 
      )
  )

nrow(icd10_exclude)
icd10_exclude
with(icd10_exclude, table(after_period, useNA = "ifany"))

icd10_exclude <-
  icd10_exclude %>%
  dplyr::filter(
    after_period == FALSE
  )

nrow(icd10_exclude)

icd10_exclude <-
  icd10_exclude %>%
  distinct(id) %>%
  mutate(diag_excl = 1)


nrow(icd10_exclude)



nrow(ukb_raw)
ukb_raw <-
  ukb_raw %>%
  left_join(
    .,
    icd10_exclude,
    "id"
  )
nrow(ukb_raw)


ukb_raw <-
  ukb_raw %>%
  mutate(
    exclude_ind =
      case_when(
        exclude_ind != 0 ~ exclude_ind,
        !is.na(diag_excl) ~ 4,
        TRUE ~ exclude_ind
      )
  )




# ---- icd10_processing2 ----

# F33	Depression (recurrent depressive episode)
# I10	History of hypertension (diagnosis of essential hypertension)
# E11	Type 2 Diabetes
# S06	Traumatic brain injury (intracranial injury)

covar_icds <- c("F33", "I10", "E11", "S06")



# F33	Depression (recurrent depressive episode)
# I10	History of hypertension (diagnosis of essential hypertension)
# E11	Type 2 Diabetes
# S06	Traumatic brain injury (intracranial injury)

rename_map <- c(
  "id"     = "id",
  "depres" = "F33", 
  "hypten" = "I10", 
  "t2diab" = "E11", 
  "brninj" = "S06"
)


covars_diag_dat <-
  ukb_raw_icd10 %>%
  dplyr::filter(type != "seco", icd10_block %in% covar_icds)

covars_diag_dat


dates_max_to_join <- 
  dates_dat %>%
  select(id, max_date, max_date_src) %>%
  dplyr::filter(!is.na(max_date)) 

dates_max_to_join


nrow(covars_diag_dat)
covars_diag_dat <-
  left_join(
    covars_diag_dat,
    dates_max_to_join,
    "id"
  )
nrow(covars_diag_dat)

covars_diag_dat

covars_diag_dat <-
  covars_diag_dat %>%
  mutate(
    before_end =
      if_else(
        is.na(max_date),
        NA,
        max_date > date 
      )
  )

with(covars_diag_dat, table(before_end, useNA = "ifany"))


covars_diag_dat <-
  covars_diag_dat %>%
  dplyr::filter(before_end)

nrow(covars_diag_dat)

covars_diag_dat <-
  covars_diag_dat %>%
  arrange(id, icd10_block, date) %>%
  group_by(id, icd10_block) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()

nrow(covars_diag_dat)

covars_diag_dat %>%
  group_by(icd10_block) %>%
  summarise(n = n(), .groups = "drop")


covars_diag_dat <-
  covars_diag_dat %>%
  select(id, icd10b = icd10_block) %>%
  mutate(ind = 1) %>%
  pivot_wider(id_cols = id, names_from = icd10b, values_from = ind, values_fill = 0)

covars_diag_dat



covars_diag_dat <-
  covars_diag_dat %>%
  select(all_of(rename_map))

covars_diag_dat



nrow(ukb_raw)
ukb_raw <-
  ukb_raw %>%
  left_join(
    .,
    covars_diag_dat,
    "id"
  )
nrow(ukb_raw)


ukb_raw %>%
  select(all_of(names(rename_map)))


sum_no_na <- function(x) sum(x, na.rm = TRUE)
(diag_cols <- names(rename_map)[-1])

## covars_diag_dat %>%
##   group_by(icd10_block) %>%
##   summarise(n = n(), .groups = "drop")
# icd10_block     n
# <chr>       <int>
# 1 E11         12239
# 2 F33           656
# 3 I10         53163
# 4 S06           739


ukb_raw %>%
  select(all_of(diag_cols)) %>%
  summarise(across(everything(), sum_no_na))

mk_na_0 <- function(x) { x[is.na(x)] <- 0; return(x) }
mk_na_0(c(1,0,NA, 0)) # test

ukb_raw %>%
  select(all_of(names(rename_map)))

ukb_raw <-
  ukb_raw %>%
  mutate(across(all_of(diag_cols), mk_na_0)) 

ukb_raw %>%
  select(all_of(names(rename_map)))

ukb_raw %>%
  select(all_of(diag_cols)) %>%
  summarise(across(everything(), sum_no_na))

# ---- inc_excl_summ ----

with(ukb_raw, table(exclude_ind, useNA = "ifany")) %>%
  as.data.frame(.) %>%
  knitr::kable(.)



# ---- write_data ----


colnames(ukb_raw)

write_rds(ukb_raw, file = off_word_file("raw/ukb_raw.rds"))

