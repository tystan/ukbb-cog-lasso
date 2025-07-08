




# ---- libs ----

# CoDA packages
library("compositions")
library("robCompositions")
library("zCompositions")
library("codaredistlm")

# tidyverse and adjascent
library("dplyr")
library("tidyr")
library("readr")
library("forcats")
library("purrr")
library("ggplot2")
library("ggthemes")

library("tictoc") # how long between the tic and toc
library("mice")   # a quiet little friend that helps with missingness

library("gtsummary") # summary stats on whole dataset
library("officer")   # summary stats on whole dataset redux




# --- constants ----


# files are not in this repo - securely stored elsewhere
(raw_f <- off_world_file("raw/ukb_raw.rds"))




# ---- ilr_funcs -----

sanitise_ilrs <- function(x) {
  
  if ("rmult" %in% class(x)) {
    class(x) <- NULL # remove "rcomp" class, will either result in numeric vector or matrix
    attr(x, "orig") <- NULL # remove original composition info (issue with indexes)
    if ("numeric" %in% class(x)) { # if vector turn into 1 row matrix
      x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
    }
  }
  return(x)
  
}

poly2 <- function(x, just_names = FALSE) {
  
  # make sure is matrix
  x <- sanitise_ilrs(x)
  
  n <- ncol(x) 
  cnames <- colnames(x)
  
  if (is.null(cnames)) {
    cnames <- paste0("c", 1L:n)
  }
  
  # get all tuples of (j,k) where j <= k
  tups <- subset(expand.grid(j = 1:n, k = 1:n), j <= k)
  tups <- tups[order(tups$j, tups$k), ] # make sure consistent ordering
  j <- tups$j
  k <- tups$k
  
  # drop = FALSE is to make sure 1 row matrices don't become vectors
  sq_out <- x[, j, drop = FALSE] * x[, k, drop = FALSE] 
  colnames(sq_out) <- paste0(cnames[j], ":", cnames[k])
  
  if (just_names) {
    return(colnames(sq_out))
  } else {
    return(sq_out)
  }
  
}

# ---- read ----


tic()
ukb_raw <- read_rds(raw_f)
toc()

format(object.size(ukb_raw), "MB") 

gc() # garbage clean-up: frees unused memory

# strange variable formatting
ukb_raw[["sex"]] <- factor(ukb_raw[["sex"]], ordered = FALSE)
summary(ukb_raw[, "sex"])
class(ukb_raw[["sex"]])


# ---- exclusions ----



# remove participants with exclusion indicator not equal to 0
ukb_raw <- subset(ukb_raw, exclude_ind == 0)

# also remove participants with reported -2 for "Number of days/week 
# walked 10+ minutes" (prior to predictor/outcome collections)
with(ukb_raw, table(day_walk_10p_0_0, useNA = "ifany"))
tibble(dte = ukb_raw[["dte_assess_ctre_0_0"]]) %>% 
  ggplot(aes(x = dte)) + geom_histogram()

with(
  subset(ukb_raw, day_walk_10p_0_0 > -2), 
  table(day_walk_10p_1_0, useNA = "ifany")
)
tibble(dte = ukb_raw[["dte_assess_ctre_1_0"]]) %>% 
  ggplot(aes(x = dte)) + geom_histogram()

with(
  subset(ukb_raw, (day_walk_10p_0_0 > -2) & (day_walk_10p_1_0 > -2)), 
  table(day_walk_10p_2_0, useNA = "ifany")
)
tibble(dte = ukb_raw[["dte_assess_ctre_2_0"]]) %>% 
  ggplot(aes(x = dte)) + geom_histogram()
subset(
  ukb_raw, 
  (day_walk_10p_0_0 > -2) & (day_walk_10p_1_0 > -2) & (day_walk_10p_2_0 == -2)
)

with(
  subset(
    ukb_raw, 
    (day_walk_10p_0_0 > -2) & (day_walk_10p_1_0 > -2) & (day_walk_10p_2_0 > -2)
  ), 
  table(day_walk_10p_3_0, useNA = "ifany")
)
tibble(dte = ukb_raw[["dte_assess_ctre_3_0"]]) %>% 
  ggplot(aes(x = dte)) + geom_histogram()

(n_old <- nrow(ukb_raw))
ukb_raw <- 
  subset(
    ukb_raw, 
    !((day_walk_10p_0_0 %in% -2) | (day_walk_10p_1_0 %in% -2))
  )
nrow(ukb_raw)
n_old - nrow(ukb_raw)

colnames(ukb_raw)

# relevel heart problems
to_del <-
  ukb_raw %>%
  rowwise() %>%
  mutate(
    hypertension =
      if_else(
        any(c_across(starts_with("hypten_0")) %in% "High blood pressure"), 
        "yes", 
        if_else(
          any(c_across(starts_with("hypten_0")) %in% "Prefer not to answer"), 
          "unknown", 
          "no")
      )
  ) %>%
  ungroup()

table(to_del$hypten_0_0)
table(to_del$hypertension)
table(to_del$hypertension, ukb_raw$hypten)
with(to_del, table(hypten_0_0, hypertension, useNA = "ifany"))

rm(to_del); gc()


ukb_raw %>%
  select(starts_with("hypten_0")) %>%
  arrange(across(everything())) %>%
  head() # tail()

with(ukb_raw, table(hypten_0_0, useNA = "ifany"))
with(ukb_raw, table(hypten_0_1, useNA = "ifany"))
with(ukb_raw, table(hypten_0_2, useNA = "ifany"))
with(ukb_raw, table(hypten_0_3, useNA = "ifany"))

with(ukb_raw, levels(hypten_0_1))
with(ukb_raw, levels(hypten_0_2))
with(ukb_raw, levels(hypten_0_3))
with(ukb_raw, levels(hypten_0_3))


# remembering the level ordering is:
lvl_str <- with(ukb_raw, levels(hypten_0_1))
lvl_str[1] <- "Prefer not to answer" # equivalent with NAs
lvl_str[2] <- "None of the above" # equivalent with no

paste(lvl_str, "->", 1:length(lvl_str))

# Takes a vector of any type, assuming its within a given range, and converts
# to a different defined set.

set.seed(1234); (test_vec <- sample(lvl_str, size = 10, replace = TRUE))
# conv_from_to(test_vec, from = lvl_str, to = 1:length(lvl_str))
as.integer(factor(as.character(test_vec), levels = lvl_str))

# test
(test_df <-
    ukb_raw %>%
    select(starts_with("hypten_0"))  %>%
    mutate(
      across(
        starts_with("hypten_0"), 
        ~ factor(as.character(.x), levels = lvl_str)
      )
    )
) # 
# at least one NA in 2nd to 4th col
# %>% dplyr::filter(!(is.na(hypten_0_1) & is.na(hypten_0_2) & is.na(hypten_0_3))))


ifelse(
  is.na(test_df$hypten_0_0[1:10]), 
  NA_integer_, 
  as.integer(as.numeric(test_df$hypten_0_0[1:10]))
)
to_int <- function(x) ifelse(is.na(x), NA_integer_, as.integer(as.numeric(x)))
to_int(test_df$hypten_0_0[1:10])



set.seed(1234)
test_df %>% 
  dplyr::filter(!(is.na(hypten_0_1) & is.na(hypten_0_2) & is.na(hypten_0_3))) %>%
  sample_n(15) %>%
  mutate(across(everything(), to_int))

# fancy boi version
rowwise_max <- 
  function(x) {
    x %>%
      mutate(across(everything(), to_int), tmp = 0) %>%
      rowwise() %>% 
      mutate(max_row_wise_fboi = max(c_across(everything()), na.rm = TRUE)) %>%
      ungroup() %>%
      select(-tmp) %>%
      mutate(
        max_row_wise_fboi = 
          if_else(
            max_row_wise_fboi == 0L, 
            NA_integer_, 
            max_row_wise_fboi
          )
      )
  }

test_df <-
  test_df %>% 
  rowwise_max()

test_df <-
  test_df %>%      
  mutate(
    max_row_wise = 
      pmax(hypten_0_0, hypten_0_1, hypten_0_2, hypten_0_3, na.rm = TRUE)
  )


(tbl_chk <- with(test_df, table(max_row_wise, max_row_wise_fboi, useNA = "ifany")))
# all diag? yup
stopifnot(sum(diag(tbl_chk)) == sum(tbl_chk))



max_rowwise_vec <-
  ukb_raw %>%      
  select(starts_with("hypten_0")) %>%
  mutate(across(everything(), ~ factor(as.character(.x), levels = lvl_str))) %>%
  rowwise_max() %>%
  mutate(
    # NAs ---> 1 which is the unknown category
    max_row_wise_fboi = 
      if_else(is.na(max_row_wise_fboi), 1L, max_row_wise_fboi),
    max_row_wise = 
      factor(max_row_wise_fboi, levels = 1:length(lvl_str), labels = lvl_str)
  )

max_rowwise_vec

max_rowwise_vec <- pull(max_rowwise_vec, max_row_wise)

ukb_raw <-
  ukb_raw %>%
  mutate(
    hypertension_fct = max_rowwise_vec,
    hypertension =
      case_when(
        hypertension_fct == lvl_str[1] ~ "unknown",
        hypertension_fct == lvl_str[2] ~ "no",
        hypertension_fct %in% lvl_str[3:length(lvl_str)] ~ "yes",
        TRUE ~ "!!ERROR!!"
      )
  )

# correct logic results 
with(ukb_raw, table(hypten_0_0, hypertension, useNA = "ifany"))


# rename/relevel depression


ukb_raw <-
  ukb_raw %>%
  mutate(
    depression =
      case_when(
        depres_0_0 %in% c('Prefer not to answer', 'Do not know') ~ "unknown",
        depres_0_0 == "Yes" ~ "yes",
        depres_0_0 == "No" ~ "no", 
        TRUE ~ "unknown"
      )
  )
table(ukb_raw$depression)
table(ukb_raw$depres_0_0)

with(ukb_raw, table(depres_0_0, depression, useNA = "ifany"))


# relevel diabetes


ukb_raw <-
  ukb_raw %>%
  mutate(
    diabetes =
      case_when(
        diabet_0_0 %in% c('Prefer not to answer', 'Do not know') ~ "unknown",
        diabet_0_0 == "Yes" ~ "yes",
        diabet_0_0 == "No" ~ "no", 
        TRUE ~ "unknown"
      )
    )

table(ukb_raw$diabetes)
table(ukb_raw$diabet_0_0)


with(ukb_raw, table(diabet_0_0, diabetes, useNA = "ifany"))


#in the final dataset, include "hypertension", "diabetes", "depression"






# ---- Maddis cog data prep ----


ukb_raw_cns <- colnames(ukb_raw)
head(ukb_raw_cns)

covars <- c(
  "id", "sex", "age_first_cog", 
  "bmi_0_0", 'alc_0_0', 'soc_iso_0_0', 'hearing_0_0', 'smoke_0_0', 
  'day_walk_10p_0_0', 
  'edu_0_0', 'edu_0_1', 'edu_0_2', 'edu_0_3', 'edu_0_4', 'edu_0_5', 
  'depression', 'hypertension', 'diabetes', 'brninj', 
  'ethn_0_0', 'accel_cog_diff'
)
covars <- ukb_raw[covars]

cog <- c(
  'id', 
  'pairs_match_0_0', 'trail_make_a', 'trail_make_b', 'symbol_digit_correct', 
  'fluid_intel', 'num_mem', 'valid_fluid_intel', 
  'valid_pairs_match', 'valid_symbol_digit', 'valid_trail_make'
)
cogdat <- ukb_raw[cog]

#explore cogdat

tbl_summary(covars, statistic= list(all_continuous() ~ "{mean} ({sd})"))

# pairs match number of incorrect ranges from 0-50
# trails A time ranges from 14 to 734 seconds
# trails b time ranges from 21 to 747 seconds
# symbol digit correct ranges from 0 to 103
# fluid intelligence ranges from 0 to 14
# numeric memory ranges from 2 to 11


class(ukb_raw[["sex"]])
class(covars[["sex"]])


# ---- reverse score ----

#first, check the value of pairs match
#record a couple of the values and ID numbers for comparison after reverse scoring
cogdat[1:10, "pairs_match_0_0"] 

#reverse score pairs match
cogdat$pairs_match_0_0 <- cogdat$pairs_match_0_0*(-1)

#sanity check - look at the values for the same ID numbers and make sure the reversal worked
cogdat[1:10, "pairs_match_0_0"]
hist(cogdat$pairs_match_0_0)

#first, check the value of trails A
#record a couple of the values and ID numbers for comparison after reverse scoring
cogdat[1:10, "trail_make_a"]

#reverse score trails A
cogdat$trail_make_a <- cogdat$trail_make_a*(-1)

#sanity check - look at the values for the same ID numbers and make sure the reversal worked
cogdat[1:10, "trail_make_a"]
hist(cogdat$trail_make_a)

#first, check the value of trails B
#record a couple of the values and ID numbers for comparison after reverse scoring
cogdat[1:10, "trail_make_b"] 

#reverse score trails B
cogdat$trail_make_b <- cogdat$trail_make_b*(-1)

#sanity check - look at the values for the same ID numbers and make sure the reversal worked
cogdat[1:10, "trail_make_b"]
hist(cogdat$trail_make_b)


# ---- transform skewed data ----

hist(cogdat$trail_make_a)
hist(cogdat$pairs_match_0_0);
hist(cogdat$trail_make_b)
hist(cogdat$symbol_digit_correct)
hist(cogdat$fluid_intel)
hist(cogdat$num_mem)

# truncate at 300 for trails
maddi_trunc <- function(x, thresh = 300) {
  x[x > thresh] <- thresh
  return(x)
}

#trails b

cogdat %>%
  ggplot(., aes(x = -log(-trail_make_b), fill = valid_trail_make)) +
  geom_histogram() +
  facet_wrap(~valid_trail_make, ncol = 1, scales = "free_y")

hist(cogdat$trail_make_b)
cogdat$trail_make_b <- -log(maddi_trunc(-cogdat$trail_make_b))
hist(cogdat$trail_make_b)

#pairs match

cogdat %>%
  ggplot(., aes(x = -(-pairs_match_0_0 + 1), fill = valid_pairs_match)) +
  geom_histogram() +
  facet_wrap(~valid_pairs_match, ncol = 1, scales = "free_y")

table(cogdat$valid_pairs_match)
quantile(
  cogdat$pairs_match_0_0, 
  c(0, 0.01, 0.025, 0.05, 0.01, 0.25, 0.5, 1), 
  na.rm = TRUE
)
cogdat$pairs_match_0_0[cogdat$pairs_match_0_0 < -7] <- -7
hist(cogdat$pairs_match_0_0)


#trails a
cogdat %>%
  ggplot(., aes(x = -log(-trail_make_a), fill = valid_trail_make)) +
  geom_histogram(bins = 20) +
  facet_wrap(~valid_trail_make, ncol = 1, scales = "free_y")


hist(cogdat$trail_make_a)
cogdat$trail_make_a[cogdat$trail_make_a < -300] <- -300
hist(cogdat$trail_make_a)
cogdat$trail_make_a <- -log(-cogdat$trail_make_a)
hist(cogdat$trail_make_a)


#symbol digit

cogdat %>%
  ggplot(., aes(x = -(-symbol_digit_correct), fill = valid_symbol_digit)) +
  geom_histogram() +
  facet_wrap(~valid_symbol_digit, ncol = 1, scales = "free_y")

table(cogdat$valid_symbol_digit)
quantile(
  cogdat$symbol_digit_correct, 
  1 - c(0, 0.010, 0.025, 0.050, 0.010, 0.250, 0.500, 1), 
  na.rm = TRUE
)
hist(cogdat$symbol_digit_correct)
cogdat$symbol_digit_correct[cogdat$symbol_digit_correct > 31] <- 31
hist(cogdat$symbol_digit_correct)

#need to join subsetted data before z-scoring so we can use age variable

class(covars[["sex"]])

nrow(cogdat)
nrow(covars)
agecog <- left_join(cogdat, covars, by="id")
colnames(agecog)
agecog <- rename(agecog, "age" = "age_first_cog")
hist(agecog$age)

# create age group categorical var

zscore <- function(x) {
  z <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  return(z)
}


#create age categories in agecog dataset
agecog$age_cat <- 
  cut(
    agecog$age, c(35, 65, 110), 
    include.lowest = TRUE, 
    right = FALSE
  )
table(agecog$age_cat)

#create mean/SD cols 
m_sd_agecog <-
  agecog %>%
  select(
    age_cat, sex, 
    pairs_match_0_0, trail_make_a,trail_make_b, 
    symbol_digit_correct, fluid_intel, num_mem
  ) %>%
  pivot_longer(cols = -c(age_cat, sex)) %>%
  dplyr::filter(!is.na(value)) %>%
  group_by(age_cat ,sex   , name) %>%
  summarise(
    m_val = mean(value),
    sd_val = sd(value),
    .groups = "drop"
  ) %>%
  arrange(name, sex, age_cat)

#pivot agecog
agecog_lng <-
  agecog %>%
  select(
    id, age_cat, sex, pairs_match_0_0, 
    trail_make_a,trail_make_b, symbol_digit_correct, 
    fluid_intel, num_mem
  ) %>%
  pivot_longer(cols = -c(id, age_cat, sex)) 

agecog[, "sex"]
summary(agecog_lng[, "sex"])
class(agecog_lng[["sex"]])

#join agecog (long version) and dataset with means and SDs
nrow(agecog_lng)
agecog_lng <-
  left_join(
    agecog_lng,
    m_sd_agecog,
    c("age_cat", "sex", "name")
  )
nrow(agecog_lng)

#create z-score col
agecog_lng <-
  agecog_lng %>%
  mutate( 
    z = 
      if_else(
        is.na(value),
        as.numeric(NA),
        (value - m_val) / sd_val
      )
  )
# dplyr::filter(!is.na(z)) 


# create global composite
agecog_lng_glob <-
  agecog_lng %>%
  dplyr::filter(!is.na(z)) %>%
  group_by(id) %>%
  summarise( 
    globalcog = mean(z),
    n = n()
  )

hist(agecog_lng_glob$globalcog)

#create memory composite

agecog_lng_mem <-
  agecog_lng %>%
  dplyr::filter(!is.na(z), name %in% c("num_mem", "pairs_match_0_0")) %>%
  group_by(id) %>%
  summarise( 
    memory = mean(value),
    n_mem = n()
  )
hist(agecog_lng_mem$memory)

#create processing speed composite

agecog_lng_ps <-
  agecog_lng %>%
  dplyr::filter(!is.na(z), name %in% c("symbol_digit", "trail_make_a")) %>%
  group_by(id) %>%
  summarise( 
    procspeed = mean(value),
    n_mem = n()
  )
hist(agecog_lng_ps$procspeed)

#filter out those with <4 cog tests from global cog data
nrow(agecog_lng_glob)
agecog_lng_glob <-
  agecog_lng_glob %>%
  dplyr::filter(n >= 4)
nrow(agecog_lng_glob)


#join global cog, memory and proc speed 
finaldat1 <- left_join(agecog_lng_glob, agecog_lng_mem, by="id")
head(finaldat1)
finaldat2 <- 
  left_join(
    finaldat1, 
    agecog_lng_ps %>% select(-n_mem), 
    by="id"
  )
head(finaldat2)

#join other cog variables to composites
final_dat <- left_join(cogdat, finaldat2, by="id")
colnames(final_dat)
cogcols <- 
  c('id', 'globalcog', 'memory', 'procspeed', 'fluid_intel', 'trail_make_b')
final_dat <- final_dat[, cogcols]
head(final_dat)

#join covariates to cogdat
final_dat <- left_join(final_dat, covars, by="id")
colnames(final_dat)

#create reasoning variable (renaming fluid intel)
final_dat$reasoning <- final_dat$fluid_intel
hist(final_dat$reasoning)

#create executive function variable (renaming trails b)
final_dat$execfunc <- final_dat$trail_make_b
hist(final_dat$execfunc)

#check covars + comp dat

tbl_summary(final_dat, statistic= list(all_continuous() ~ "{mean} ({sd})"))

# ---- create new categorical variables for covars ----

#check sex
table(final_dat$sex, useNA = "ifany") #all good
class(final_dat$sex)

#check age
hist(final_dat$age_first_cog) #all good
class(final_dat$age_first_cog)
final_dat$age_first_cog <- as.numeric(final_dat$age_first_cog)

final_dat$age_cat <- cut(final_dat$age_first_cog, c(35, 65, 110), include.lowest = TRUE, right = FALSE)
table(final_dat$age_cat, useNA = "ifany")
head(final_dat$age_cat)



#check bmi
hist(final_dat$bmi_0_0)
final_dat$bmi_0_0 <- as.numeric(final_dat$bmi_0_0)

final_dat <- 
  final_dat %>%
  mutate(bmicat = case_when(
    bmi_0_0 <= 25 ~  "not_overweight",
    bmi_0_0 >  25 ~ "overweight",
    TRUE ~ NA
  ))

table(final_dat$bmicat, useNA = "ifany")
final_dat$bmicat <- as.factor(final_dat$bmicat)
head(final_dat$bmicat)

final_dat <- 
  final_dat %>% 
  rename(bmi = bmi_0_0)

#check alcohol
table(final_dat$alc_0_0)

final_dat <- 
  final_dat %>%
  mutate(
    alcohol = 
      case_when(
        alc_0_0 %in% 
          c("One to three times a month", "Special occasions only", "Once or twice a week", "Never") ~ 
          "sometimes/never", 
        alc_0_0 %in% c("Daily or almost daily", "Three to four times a week") ~ 
          "often/very often",
        alc_0_0 == "Prefer not to answer" ~ 
          "unknown",
        TRUE ~ "unknown"
      )
  )

table(final_dat$alcohol, useNA = "ifany")
head(final_dat$alcohol)
final_dat$alcohol <- as.factor(final_dat$alcohol)
final_dat$alcohol <- relevel(final_dat$alcohol, ref = "sometimes/never")
head(final_dat$alcohol)

#social isolation
table(final_dat$soc_iso_0_0)

final_dat <- 
  final_dat %>%
  mutate(
    isolation = 
      case_when(
        soc_iso_0_0 %in% c("Do not know", "Prefer not to answer") ~ "unknown", 
        soc_iso_0_0 == "No" ~ "no",
        soc_iso_0_0 == "Yes" ~ "yes",
        TRUE ~ "unknown"
      )
  )

table(final_dat$isolation, useNA = "ifany")
final_dat$isolation <- factor(final_dat$isolation, levels = c("no", "yes", "unknown"))
head(final_dat$isolation)
table(final_dat$isolation, useNA = "ifany")

#hearing loss 
table(final_dat$hearing_0_0, useNA = "ifany")

final_dat <- 
  final_dat %>%
  mutate(hearing = case_when(
    hearing_0_0 %in% c("Yes", "I am completely deaf") ~ "yes", 
    hearing_0_0 == "No" ~ "no",
    hearing_0_0 %in% c("Prefer not to answer", "Do not know") ~ "unknown",
    TRUE ~ "unknown"
  )
  )

table(final_dat$hearing, useNA = "ifany")
final_dat$hearing <- factor(final_dat$hearing, levels = c("no", "yes", "unknown"))
table(final_dat$hearing, useNA = "ifany")

#smoking
table(final_dat$smoke_0_0, useNA = "ifany")
final_dat <- rename(final_dat, "smoking" = "smoke_0_0")
table(final_dat$smoking, useNA = "ifany")

final_dat <- 
  final_dat %>%
  mutate(smoking = case_when(
    smoking == "Prefer not to answer" ~ "unknown",
    smoking == "Never" ~ "never",
    smoking == "Previous" ~ "previous", 
    smoking == "Current" ~ "current",
    TRUE ~ "unknown"
  )
  )

table(final_dat$smoking, useNA = "ifany")
final_dat$smoking <- factor(final_dat$smoking, levels = c("never", "previous", "current", "unknown"))
table(final_dat$smoking, useNA = "ifany")

#education
table(final_dat$edu_0_0, useNA = "ifany") %>% knitr::kable(.)


with(final_dat, mean(globalcog, na.rm = TRUE))

# Define the ordered levels (ranked education in order)
qual_levels <- c("unknown", "high school", "cert III/diploma", "other professional qual", "college/university")

relevel_edu <- function(dat, col_nm, new_col_nm) {
  
  col_vec <- dat[[col_nm]]
  
  new_col <-
    case_when(
      col_vec %in% c("A levels/AS levels or equivalent", "O levels/GCSEs or equivalent", "CSEs or equivalent") ~ qual_levels[2], 
      col_vec == "NVQ or HND or HNC or equivalent" ~ qual_levels[3],
      col_vec == "Other professional qualifications eg: nursing, teaching" ~ qual_levels[4],
      col_vec == "College or University degree" ~ qual_levels[5],
      col_vec %in% c("None of the above", "Prefer not to answer") ~ qual_levels[1],
      TRUE ~ qual_levels[1]
    )
  
  # Convert education variables to factors with ordered levels
  dat[[new_col_nm]] <- factor(new_col, levels = qual_levels, ordered = TRUE)
  
  return(dat)
  
}

final_dat <- 
  final_dat %>%
  mutate(education1 = case_when(
    edu_0_0 %in% c("A levels/AS levels or equivalent", "O levels/GCSEs or equivalent", "CSEs or equivalent") ~ "high school", 
    edu_0_0 == "Other professional qualifications eg: nursing, teaching" ~ "other professional qual",
    edu_0_0 %in% c("None of the above", "Prefer not to answer") ~ "unknown",
    edu_0_0 == "College or University degree" ~ "college/university",
    edu_0_0 == "NVQ or HND or HNC or equivalent" ~ "cert III/diploma",
    TRUE ~ "unknown"
  )
  )

final_dat$education1 <- factor(final_dat$education1, levels = qual_levels, ordered = TRUE)

table(final_dat$education1, useNA = "ifany") %>% knitr::kable(.)

final_dat <- relevel_edu(final_dat, "edu_0_0", "education1a") 

final_dat %>% 
  select(education1, education1a) %>%
  dplyr::filter(education1 != education1a)

# now column created by function is OK, remove original and rename new
final_dat <-
  final_dat %>% 
  select(-education1) %>%
  rename(education1 = education1a)

table(final_dat$education1, useNA = "ifany")

final_dat <- relevel_edu(final_dat, "edu_0_1", "education2") 
table(final_dat$education2, useNA = "ifany")

final_dat <- relevel_edu(final_dat, "edu_0_2", "education3") 
table(final_dat$education3, useNA = "ifany")

final_dat <- relevel_edu(final_dat, "edu_0_3", "education4") 
table(final_dat$education4, useNA = "ifany")

final_dat <- relevel_edu(final_dat, "edu_0_4", "education5") 
table(final_dat$education5, useNA = "ifany")

final_dat <- relevel_edu(final_dat, "edu_0_5", "education6") 
table(final_dat$education6, useNA = "ifany")



# Now, you can use the pmax function as before to find the highest qualification
final_dat <- 
  final_dat %>%
  mutate(highestqual = pmax(education1, education2, education3, education4, education5, education6, na.rm = TRUE))

final_dat %>%
  select(highestqual, education1, education2, education3, education4, education5, education6)

table(final_dat$highestqual, useNA = "ifany")
class(final_dat$highestqual)
final_dat$highestqual <- factor(final_dat$highestqual, ordered = FALSE) 
final_dat$highestqual <- relevel(final_dat$highestqual, ref = "high school")
table(final_dat$highestqual, useNA = "ifany")
class(final_dat$highestqual)


#TBI history

table(final_dat$brninj, useNA = "ifany")

# ethnicity


# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000
# from website
#
# Category Count
# *	White 14393 
#    -	British 472098
#    -	Irish 13992
#    -	Any other white background 17023
# *	Mixed151 
#    -	White and Black Caribbean 649
#    -	White and Black African 449
#    -	White and Asian 863
#    -	Any other mixed background 1079
# *	Asian or Asian British 228 
#    -	Indian 6108
#    -	Pakistani 1903
#    -	Bangladeshi 240
#    -	Any other Asian background 1853
# *	Black or Black British 164 
#    -	Caribbean 4606
#    -	African 3460
#    -	Any other Black background 126
# *	Chinese 1690
# *	Other ethnic group 4764
# *	Do not know 226
# *	Prefer not to answer 1778


levels(final_dat$ethn_0_0)

# create more coarse categorisations but hoping the potentially problematic
# variable(s) is helpful in highlighting disparity as a proxy/confounder
# with environmental/societal factors
white_cat <- c(
  "White", "British", "Irish", "Any other white background"
)
mixed_cat <- c(
  "Mixed",
  "White and Black Caribbean", "White and Black African", 
  "White and Asian", "Any other mixed background"
)
asian_cat <- c(
  "Asian or Asian British",
  "Indian", "Pakistani", "Bangladeshi", "Any other Asian background"
)
black_cat <- c(
  "Black or Black British",
  "Caribbean", "African", "Any other Black background"
)
chine_cat <- c("Chinese")
# not our language, would prefer not to "other" anyone:
other_cat <- c("Other ethnic group") 
unkow_cat <- c("Prefer not to answer", "Do not know")


final_dat <-
  final_dat %>%
  mutate(
    ethn_char = as.character(ethn_0_0),
    ethn_toplvl =
      case_when(
        is.na(ethn_char)         ~ "Unkown",
        ethn_char %in% white_cat ~ "White",
        ethn_char %in% mixed_cat ~ "Mixed",
        ethn_char %in% asian_cat ~ "Asian", # "or Asian British",
        ethn_char %in% black_cat ~ "Black", # " or Black British",
        ethn_char %in% chine_cat ~ "Chinese",
        ethn_char %in% other_cat ~ "Not specified ethnic group",
        ethn_char %in% unkow_cat ~ "Unkown",
        TRUE ~ "ERROR!!!" # we've missed categories if gets to this
      )
  )


with(
  final_dat,
  table(
    ethn_0_0, 
    ethn_toplvl,
    useNA = "ifany"
  )
)


table(final_dat$ethn_toplvl, useNA = "ifany") # no missing
# class(final_dat$ethn_toplvl)
# final_dat$ethn_toplvl <- as.character(final_dat$ethn_toplvl)
head(final_dat$ethn_toplvl)

final_dat <- 
  final_dat %>%
  mutate(ethnicity = case_when(
    ethn_toplvl == "White" ~ "White",
    ethn_toplvl %in% 
      c("Asian", "Black", "Chinese", "Mixed", "Not specified ethnic group", "Unkown") ~ 
      "Nonwhite", 
    TRUE ~ NA
  ))

table(final_dat$ethnicity, useNA = "ifany")
final_dat$ethnicity <- factor(final_dat$ethnicity, levels = c("White", "Nonwhite"))
table(final_dat$ethnicity, useNA = "ifany")

#hypertension


colnames(final_dat)
#filter dataset to what will be needed for modelling
cols <- c(
  'id', 'globalcog', 'memory', 'procspeed', 'execfunc', 'reasoning',
  'alcohol', 'isolation', 'hearing', 'bmi',
  'highestqual', 'age_first_cog', 'sex', 'smoking',
  'diabetes', 'depression', 'ethnicity', 'brninj', 'hypertension', 'accel_cog_diff'
)

data <- final_dat[, cols]
head(data)

table(apply(data, 1, function(x) sum(is.na(x)) > 0))
summary(data)

md.pattern(data, rotate.names = TRUE)

# remove predictors with NAs
nrow(data)
bmi_nas <- is.na(data$bmi)
data <- data[!bmi_nas, ]
nrow(data)


# now at least one cog outcome var present for all rows
md.pattern(data, rotate.names = TRUE)


# ---- create time-use compositions ----

comps0 <- c('sleep', 'sb', 'lpa', 'mvpa')
comps <- c('id', comps0)

compdat <- ukb_raw[, comps]
nrow(compdat)
nrow(data)
data <- left_join(data, compdat, by = "id")
nrow(data)
colnames(data)
table(is.na(data$sleep))
table(is.na(data$sb))
table(is.na(data$lpa))
table(is.na(data$mvpa))

comp_chk <- 
  data %>% 
  select(all_of(comps0))


hist(1440 * rowSums(comp_chk))

table(1440 * rowSums(comp_chk), useNA = "ifany")


data %>% 
  select(all_of(comps)) %>%
  pivot_longer(cols = -id) %>%
  mutate(name = fct_inorder(name)) %>%
  ggplot(., aes( x = value, fill = name)) +
  geom_histogram(bins = 40) +
  facet_wrap( ~name) + # , scale = "free_y") +
  theme_bw() +
  scale_fill_tableau()

ggsave("fig/raw_comp_dist.png", width = 8, height = 5)


#check for zeroes in composition
min(data$sleep) #no zeroes
min(data$sb) #no zeroes
min(data$lpa) #zeroes present
min(data$mvpa) #zeroes present

#nest compositional variables in "comp"
data$comp <- data %>% dplyr::select(sleep, sb, lpa, mvpa)
data$comp <- as.data.frame(data$comp)
class(data$comp)

#this should give one value which is the smallest non-zero value nested composition
(zero_thesh <- min(data$comp[data$comp > 0], na.rm=TRUE)) 
(detect_lim <- 7 * zero_thesh) # as 7 day average

#replace XXXX here with minimum non-zero value detected previously
data$icomp = zCompositions::lrEM(data$comp, label=0, dl=rep(detect_lim, 4)) 
#in case R doesn’t know it’s a composition – I don’t think lrEM produces compositional variables
data$icomp=acomp(data$icomp)


test_obj <- impRZilr(data$comp, dl=rep(7 * zero_thesh, 4), eps=0.001, method = "lm")

test_obj$x <- acomp(test_obj$x)

test_obj$wind
test_obj$iter
# diffs_ii <- data$id[which(rowSums(test_obj$wind) > 0)]
diffs_ii <- which(rowSums(test_obj$wind) > 0)


ids_w_imputed_dat <- data$id[which(rowSums(test_obj$wind) > 0)]

#pull imputed values to separate columns again
data$sleep2 = data$icomp[,1]
data$sb2 = data$icomp[,2]
data$lpa2 = data$icomp[,3]
data$mvpa2 =  data$icomp[,4]

#check imputation worked
min(data$mvpa2) #all good
min(data$lpa2) #all good

#create new dataframes called cols_pred (predictors), cols_outcome, and cols_covar

#first rename age, bmi, ethnicity
colnames(data)
data <- rename(data, "age" = "age_first_cog")

rename_tu <-
  c(
    "sleep"="sleep2",
    "sb" = "sb2",
    "lpa"="lpa2",
    "mvpa" = "mvpa2"
  )

data <- 
  data %>%
  select(-all_of(names(rename_tu)))

colnames(data) 

#we will rename the columns in finaldat for future simplicity
data <- 
  data %>%
  rename(all_of(rename_tu))

colnames(data) 

cols_tu <- names(rename_tu)

cols_outcome <- 
  c("globalcog", "memory", "procspeed", "execfunc", "reasoning")

pred_vars_cat <- 
  sort(c(
    'sex', 'highestqual', 'hearing', 
    'isolation', 'alcohol', 'smoking', 
    'ethnicity', 'depression', 'diabetes', 
    'hypertension', 'brninj'
  ))

pred_vars_cts <- c('accel_cog_diff', 'age', 'bmi', 'ilr')

cols_pred_all <- c(pred_vars_cat, pred_vars_cts)

cols_covar <- cols_pred_all[!(cols_pred_all %in% "ilr")]


cols_covar[!(cols_covar %in% colnames(data))]
colnames(data)[!(colnames(data) %in% cols_covar)]

cols_want <- c("id", cols_outcome, cols_tu, cols_covar)

data <- data[, cols_want]
colnames(data)

table(abs(rowSums(data[, cols_tu]) - 1) < 1e-12, useNA = "ifany")

#apply closure to 1440
data[, cols_tu] <- clo(data[, cols_tu], total = 1440)

table(abs(rowSums(data[, cols_tu]) - 1440) < 1e-12, useNA = "ifany")

cols_tu #check order

colnames(data)

ilr_nms <- paste0("ilr", 1:3)

#build the sequential binary partition matrix with 4 compositional parts (therefore, 3 ilrs)
sbp4 <- matrix(c(
  +1,  0,  0,
  -1, +1,  0,
  -1, -1, +1,
  -1, -1, -1),
  byrow=TRUE, 
  ncol=3,
  dimnames = list(cols_tu, ilr_nms))

sbp4

psi4 <- compositions::gsi.buildilrBase(sbp4)

ilrs <- compositions::ilr(data[, cols_tu], V = psi4)


head(ilrs) #check that the renaming of column titles worked - looks good!

#join the ilrs to the end of the main dataset
data <- cbind(data[, c("id", cols_outcome, cols_covar, cols_tu)], ilrs)
head(data)



# ---- finalise_dataset ----


#define matrix of predictor variables

prep_dat <- 
  data %>% 
  select(all_of(c("id", cols_outcome, ilr_nms, cols_covar))) %>%
  select(id, everything())

nrow(prep_dat)
rows_to_rm <- na.omit(prep_dat %>% dplyr::select(-all_of(cols_outcome[-1])))
str(rows_to_rm)
rows_to_rm <- unname(attr(rows_to_rm, "na.action"))
prep_dat <- prep_dat[-rows_to_rm, ]
nrow(prep_dat)

prep_dat$ilr <- with(prep_dat, cbind(ilr1, ilr2, ilr3))
head(prep_dat)

prep_dat <- prep_dat %>% dplyr::select(-all_of(ilr_nms))
head(prep_dat)
colnames(prep_dat)


### testing:
# colnames(prep_dat)[ grepl("bmi", colnames(prep_dat))]
# head(poly(prep_dat[, "bmi", drop = TRUE], 2, raw = TRUE))
# head(poly2(prep_dat[, "bmi", drop = FALSE]))
# poly2(as.matrix(prep_dat[, "bmi", drop = FALSE]), just_names = TRUE)
# class(as.matrix(prep_dat[, "bmi", drop = FALSE]))
# class(prep_dat[, "ilr", drop = FALSE])
# class(prep_dat[, "ilr", drop = TRUE])
# poly2(prep_dat[, "ilr", drop = TRUE], just_names = TRUE)
# poly2(prep_dat[, "ilr", drop = TRUE])
# colnames(as.matrix(prep_dat[["ilr"]]))
# colnames(as.matrix(prep_dat[, "bmi", drop = FALSE]))
# poly2(as.matrix(prep_dat[, "ilr", drop = FALSE]), just_names = TRUE)

# this creates the R formula with the logic to treat ilrs, other cts variables 
# and categorical variables differently
n_int <- length(pred_vars_cts)
int_terms <- NULL
for (i in 1:n_int) {
  if (pred_vars_cts[i] == "ilr") {
    prep_dat$ilr_sq <- poly2(prep_dat[, pred_vars_cts[i], drop = TRUE])
    int_terms <- 
      c(
        int_terms,
        "ilr_sq"
      )
    
  } else {
    int_terms <- 
      c(
        int_terms,
        poly2(as.matrix(prep_dat[, pred_vars_cts[i], drop = FALSE]), just_names = TRUE)
      )
  }
}
int_terms
(int_terms <- sort(int_terms))

class(prep_dat[, "ilr"])
class(prep_dat[, "ilr_sq"])

form_complex <-
  paste(
    " ~ ",
    "(", paste(c(pred_vars_cat, pred_vars_cts), collapse = " + "), ")^2",
    "+", paste(int_terms, collapse = " + ")
  )

form_complex <- as.formula(form_complex)


form_complex
tibble(prep_dat)
levels(prep_dat$highestqual)

# ---- plot_outcomes ----


prep_dat %>%
  select(all_of(c("id", cols_outcome))) %>%
  pivot_longer(cols = -id) %>%
  mutate(name = fct_inorder(name)) %>%
  na.omit(.) %>%
  group_by(name) %>%
  summarise(
    wdths = (max(value) - min(value)) / 40,
    .groups = "drop"
  )

bn_wdth_fn <- function(x) {
  ifelse(
    ((max(x) - min(x)) / 40) < 0.4, 
    (max(x) - min(x)) / 40,
    1
  )
}



prep_dat %>%
  select(all_of(c("id", cols_outcome))) %>%
  pivot_longer(cols = -id) %>%
  mutate(name = fct_inorder(name)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = value, fill = name)) +
  geom_histogram(binwidth = bn_wdth_fn) +
  facet_wrap(~ name, ncol = 2, scales = "free_x") +
  theme_bw() +
  scale_fill_tableau()


ggsave("fig/cog_outcome_dist.png", width = 8, height = 5)


# ---- individual_outcome_datasets ----

cols_outcome[1]
(n_outc <- length(cols_outcome))

data_filt <- 
  prep_dat %>% 
  select(-all_of(cols_outcome))

terms(form_complex)
attr(terms(form_complex), "variables")

nrow(data_filt)
data_filt <- na.omit(data_filt)
nrow(data_filt)



y_list_unstd <- vector(mode = "list", length = n_outc)
names(y_list_unstd) <- cols_outcome
for (i in 1:n_outc) {
  y_list_unstd[[i]] <- prep_dat[[cols_outcome[i]]]
}

glimpse(y_list_unstd)


std_vec <- function(x, na_rm = TRUE) {
  out_obj <- list(x_std = NULL, m_s_vals = NULL)
  m_x <- mean(x, na.rm = na_rm)
  s_x <- sd(x, na.rm = na_rm)
  out_obj$x_std <- (x - m_x) / s_x
  out_obj$m_s_vals <- c("m" = m_x, "s" = s_x)
  return(out_obj)
}
# test
v <- 1:10
(v_std <- std_vec(v))
stopifnot(all((v_std$x_std * v_std$s + v_std$m) == v))

std_by_col <- function(lst) {
  n_l <- length(lst)
  l_nms <- names(lst)
  # lst_std <- lst
  
  lst_std_obj <- lapply(lst, std_vec)
  lst_std <- lapply(lst_std_obj, function(x) x$x_std)
  lst_msd <- lapply(lst_std_obj, function(x) x$m_s_vals)
  names(lst_std) <- names(lst_msd) <- l_nms
  
  return(list(dat = lst_std, m_s_vals = lst_msd))
  
}
# test
vu_std_obj <- std_by_col(list(v = 1:10, u = -(10:1)))
vu_std_obj$dat
vu_std_obj$m_s_vals


# actual
col_std_obj <- std_by_col(y_list_unstd)
y_list <- col_std_obj$dat
lapply(y_list, head)
(y_list_m_s <- col_std_obj$m_s_vals)


bind_cols(y_list) %>%
  mutate(id = 1:n()) %>%
  select(all_of(c("id", cols_outcome))) %>%
  pivot_longer(cols = -id) %>%
  mutate(name = fct_inorder(name)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = value, fill = name)) +
  geom_histogram() + # binwidth = bn_wdth_fn) +
  facet_wrap(~ name, ncol = 2, scales = "free_x") +
  theme_bw() +
  scale_fill_tableau()


ggsave("fig/cog_outcome_dist_stdised.png", width = 8, height = 5)





# ---- final clean of data ----

data_filt


get_lvls <- function(x) {
  # x <- pull(x)
  x_lvls <- NULL
  if ("character" %in% class(x)) {
    x_lvls <- levels(factor(x)) # is char so convert to factor
  } else {
    x_lvls <- levels(x)         # is already factor, extract lvls
  }
  return(x_lvls)
}


this_xlev <- lapply(data_filt, get_lvls)
this_xlev


data_filt_fctised <- data_filt
xcnames <- colnames(data_filt_fctised)
for (j in 1:ncol(data_filt_fctised)) {
  if ("character" %in% class(data_filt_fctised[[j]])) {
    cname_j <- xcnames[j]
    data_filt_fctised[[j]] <- 
      factor(
        data_filt_fctised[[j]], 
        levels = this_xlev[[cname_j]]
      )
  }
}

as_tibble(data_filt) # some char columns
as_tibble(data_filt_fctised) # all char columns are now factors

data_filt <- data_filt_fctised

# ---- make design matrix ----

head(x <- model.matrix(form_complex, data = data_filt, xlev = this_xlev))
colnames(x)
str(x)

x_dsgn_meta <-
  tibble(
    var_and_lvl = colnames(x),
    var_src = attr(x, "assign")
  ) %>%
  mutate(col_no = 1:(n() + 0)) %>%
  select(col_no, everything())
x_dsgn_meta %>% print(., n = nrow(.))


x_src_meta <-
  tibble(
    var = attr(terms(form_complex, simplify = TRUE), "term.labels"),
    var_order = attr(terms(form_complex, simplify = TRUE), "order")
  ) %>%
  mutate(var_src = 1:n()) %>%
  select(var_src, everything())
x_src_meta %>% print(., n = nrow(.))

x_dsgn_meta <-
  left_join(
    x_dsgn_meta,
    x_src_meta,
    "var_src"
  )
x_dsgn_meta <-
  x_dsgn_meta %>% 
  mutate(var_order = if_else(is.na(var_order), 0L, var_order))

x_dsgn_meta %>% print(., n = nrow(.))

### testing
# attr(terms(form_complex), "variables")
# str(terms(form_complex))
# attr(terms(form_complex), "intercept")
# str(terms(form_complex, simplify = TRUE))
# attr(terms(form_complex), "term.labels")
# attr(terms(form_complex), "order")
# rowSums(attr(terms(form_complex, simplify = TRUE), "factors"))
# attr(x, "contrasts")$highestqual 
# all.vars(terms(form_complex))
# all.vars((form_complex))
# terms(form_complex, simplify = TRUE)
# attr(terms(form_complex, simplify = TRUE), "term.labels")
# attr(terms(form_complex, simplify = TRUE), "order")


# ---- check_for_non-0_column_counts_inder_thresh ----

# need to have, say, at least XX non-zero categorical interactions otherwise,
# there is literally no information the interaction provides
min_non_zero_entries <- 10


colsms <- unlist(lapply(as.data.frame(x), function(x) sum(abs(x))))
names(colsms)[colsms == 0]

to_rm_cols <- which(colsms == 0)
to_rm_cols

sumabsvals <- apply(x, 2, function(x) sum(abs(x) > 1e-12)) # 0

sumabsvals[sumabsvals < min_non_zero_entries]
length(sumabsvals[sumabsvals < min_non_zero_entries])
length(sumabsvals)
sort(sumabsvals[sumabsvals < min_non_zero_entries])

to_rm_df <-
  tibble(
    col_no = 1:length(sumabsvals),
    var_and_lvl_chk = names(sumabsvals),
    non_0_cnts = sumabsvals,
    col_inc = if_else(sumabsvals < min_non_zero_entries, 0L, 1L)
  )
to_rm_df %>% print(., n = nrow(.))

x_dsgn_meta <-
  left_join(
    x_dsgn_meta,
    to_rm_df,
    "col_no"
  )

if (with(x_dsgn_meta, any(var_and_lvl == var_and_lvl_chk))) {
  # checking col names ok, rm redundant col
  x_dsgn_meta <- x_dsgn_meta %>% select(-var_and_lvl_chk)
} else {
  stop("join has non-conforming columns, need to check logic")
}

x_dsgn_meta %>% print(., n = nrow(.))

# use formula object to identify the intercept col if exists
int_col_no <- attr(terms(form_complex), "intercept")
if (int_col_no > 0) {
  x_dsgn_meta$col_inc[int_col_no] <- 0L
}
x_dsgn_meta
sum(x_dsgn_meta$col_inc == 1L) # " == 1L" redundant but for readability
sum(x_dsgn_meta$col_inc == 0L)

# manually make ilr_sq order == 2
x_dsgn_meta %>% filter(row_number() %in% 20:40) %>% print(., n = nrow(.))
x_dsgn_meta <- 
  x_dsgn_meta %>% 
  mutate(var_order = if_else(var == "ilr_sq", 2L, var_order))
x_dsgn_meta %>% filter(row_number() %in% 20:40) %>% print(., n = nrow(.))


cols_rm_names <- names(sort(sumabsvals[sumabsvals < min_non_zero_entries]))

sum(colnames(x) %in% cols_rm_names)

nrow(x)
ncol(x)
colnames(x)
head(x[, 1:10])


# ---- check_all_outcomes_for_non-0_column_counts ----

cat(
  "\n\n---- checking for each outcome variable whether all",
  "variables still have >=", 
  min_non_zero_entries, 
  "non-0 entries ----\n\n"
)
for (i in 1:n_outc) {
  na_i <- is.na(y_list[[i]])
  cat("\n\n#:::::\n\noutcome", cols_outcome[i], "has", sum(na_i), "NA values.\n\n")
  sumabsvals_i <- apply(x[!na_i, ], 2, function(x) sum(abs(x) > 1e-10))
  col_names_thresh_i <- sort(sumabsvals_i[sumabsvals_i < min_non_zero_entries])
  df_additional_col_i <- 
    tibble(
      cnm = names(col_names_thresh_i), 
      non_0_entries = col_names_thresh_i
    ) %>%
    anti_join(
      .,
      x_dsgn_meta %>% filter(col_inc == 0L),
      c("cnm" = "var_and_lvl")
    )
  
  if (nrow(df_additional_col_i) > 0) {
    cat("These are the additional columns below the threshold to consider:\n")
    print(df_additional_col_i)
  } else {
    cat("No additional columns to consider. Great!")
  }
  cat("\n\n#:::::\n\n")
}




### testing

# lapply(as.data.frame(x), class)[grepl("sex", colnames(x))]
colnames(x)[ grepl("sex", colnames(x))]
summary(x[, "sexMale"])


ncol(x)

# ncol(x)

sum(colnames(x) %in% cols_rm_names)
all(cols_rm_names %in% x_dsgn_meta$var_and_lvl[x_dsgn_meta$col_inc == 0L])
(x_dsgn_meta$var_and_lvl[x_dsgn_meta$col_inc == 0L] %in% cols_rm_names)

x_dsgn <- x # redundant: [,  !(colnames(x) %in% cols_rm_names), drop = FALSE]



colnames(x_dsgn)
ncol(x_dsgn) #375 (+ intercept + 19 to rm) possible coefficients
nrow(x_dsgn) #~53k


x_df <- data_filt

# find purely predictor columns without redundancy
# ironically, this is close to redundant itself
strt_col <- which(colnames(x_df) == "id") + 1
colnames(x_df)[strt_col - 1] # id col
x_df <- x_df[, strt_col:ncol(x_df)]




head(model.frame(form_complex, data = x_df))
head(x_df)

all(colnames(x_df) %in% colnames(model.frame(form_complex, data = x_df)))
all(colnames(model.frame(form_complex, data = x_df)) %in% colnames(x_df))


colnames(x_df)
ncol(x_df) 
nrow(x_df) #~53k
head(x_df)
sapply(x_df, function(x) class(x)[1])
sapply(x_df %>% select(starts_with("ilr")), function(x) class(x[1, 1]))
(this_xlev <- sapply(x_df, get_lvls))



# ---- save_data ----

### NB:
# files are not in this repo - securely stored elsewhere


# note: 
# only need to save 
#   | x_dsgn, x_df, form_complex, cols_rm_names and xlvl_lst ONCE.
# y are dependent on cog outcome and 
#   | are stored the outcomes in a list() of length 5 for each


# outcome variables:
(y_rds <- off_world_file("dat/y_lst.rds"))
write_rds(y_list, file = y_rds) 

(y_m_s_rds <- off_world_file("dat/y_orig_m_s.rds"))
write_rds(y_list_m_s, file = y_m_s_rds) 

# unstandardised version
(y_list_unstd_rds <- off_world_file("dat/y_list_unstd.rds"))
write_rds(y_list_unstd, file = y_list_unstd_rds) 

# design matrix (`model.matrix`/numeric matrix with dummy vars encoding factors)
# with all possible interactions (less removed columns)
write_rds(x_dsgn, file = off_world_file("dat/x_dsgn.rds")) 

# dataset with all predictor variables + squared ilrs (as a data.frame with factors)
write_rds(x_df, file = off_world_file("dat/x_df.rds")) 

# model formula
write_rds(form_complex, file = off_world_file("dat/mod_form_cplx.rds")) 

# info to use in modifying the design matrix when required, including
# columns to remove if using the mod_form_cplx model formula and x_df 
# in model call
write_rds(x_dsgn_meta, file = off_world_file("dat/x_dsgn_meta.rds"))  
# i.e., 
# (previous approximation of saved vector of names)
cols_to_rm <- x_dsgn_meta$var_and_lvl[x_dsgn_meta$col_inc == 0L][-1] # remove "(Int)"
as.formula(
  paste(
    "~ ",
    as.character(form_complex)[2],
    paste0("- `", cols_to_rm, "`",  collapse = " ")
  )
)


# levels of categorical variables
write_rds(this_xlev, file = off_world_file("dat/xlvl_lst.rds")) 






