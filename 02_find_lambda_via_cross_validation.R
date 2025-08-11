


# ---- libs ----

# tidyverse and adjacent
library("dplyr")
library("tidyr")
library("readr")
library("forcats")
library("purrr")
library("furrr") # the "future" version of purrr (https://www.futureverse.org/)
library("stringr")
library("ggplot2")
library("ggthemes")
library("scales")
library("tidyheatmaps") # thanks Maddi!

library("knitr")
library("tictoc")
library("foreach") # alternative loopin', and allows parallel computation too
library("arrow")   # save data light and fast

# Lasso reg
# library("glmnet") # need grp lasso which isn't available here

# group lasso
# install.packages("grpreg")
library("grpreg")

# make sure the wanted versions are the global env go-to
select <- dplyr::select # thanks Dot!
filter <- dplyr::filter # thanks Dot!



# ---- import -----


y_lst           <- read_rds(file = off_world_file("dat/y_lst.rds"))
x_df            <- read_rds(file = off_world_file("dat/x_df.rds"))
x_dsgn          <- read_rds(file = off_world_file("dat/x_dsgn.rds"))
x_dsgn_meta     <- read_rds(file = off_world_file("dat/x_dsgn_meta.rds"))
mod_form_cplx   <- read_rds(file = off_world_file("dat/mod_form_cplx.rds"))
xlvl_lst        <- read_rds(file = off_world_file("dat/xlvl_lst.rds"))

# having a peak
ncol(x_dsgn)
nrow(x_dsgn)
colnames(x_dsgn)
head(x_dsgn[, 1:10])
mod_form_cplx

# names of vars to include in the analysis
x_dsgn_meta$var_and_lvl[x_dsgn_meta$col_inc == 0L]
stopifnot(nrow(x_dsgn_meta) == ncol(x_dsgn))

# final design matrix
X <- x_dsgn[, x_dsgn_meta$col_inc == 1L]
nrow(X)
ncol(X)
head(X)


table(is.na(X), useNA = "ifany")
lapply(y_lst, function(x) table(is.na(x), useNA = "ifany"))

# create groups for group lasso (the groups are levels within variables)
v_grps <- pull(x_dsgn_meta[x_dsgn_meta$col_inc == 1L, "var"])
v_grps


# testing
set.seed(1234)
(n_x <- nrow(X))
(n_k <- ceiling(n_x / 10))
n_x - 9 * n_k
kfold_ind <- c(rep(1:9, each = n_k), rep(10, n_x - 9 * n_k))
kfold_ind <- sample(kfold_ind)
table(kfold_ind, useNA = "ifany")
# stop if: x and kfold_ind don't have the same length
stopifnot(nrow(X) == length(kfold_ind))
head(kfold_ind)



# ---- gaussian_outcomes ----


# wrapper function to fit 10-fold cross-validation of grp lassos
# has some checks for complete records and let user know what is being deleted
fit_cv_grplasso <- function(y, X_ = X, grps_ = v_grps, folds_ = kfold_ind, eps_ = 1e-4) {
  
  # stop if: x and y don't have the same number of obs
  stopifnot(length(y) == nrow(X_))
  rm_idx <- is.na(y) | (rowSums(is.na(X_)) > 0)

  if (sum(rm_idx) > 0) {
    message(
      "Removing the following indexes from 'y', 'X_' and 'folds_':\n", 
      paste(which(rm_idx), collapse = ", ")
    )
  }
  
  return(cv.grpreg(
    # keep matrix structure no matter what, 
    # even if there's a fire (I'm not calling him dad)
    X_[!rm_idx, , drop = FALSE],
    y[!rm_idx], 
    group = grps_, 
    fold = folds_[!rm_idx], 
    trace = TRUE,
    seed = 9876,
    penalty = "grLasso",
    eps = eps_
    # returnX = TRUE
  ) )
}



# create generic lists to hold models then insert
(mod_nms <- paste0("cv_grplasso_", names(y_lst)))
# make a tibble that holds models in elements (and other data types comming up)
(mod_df <- tibble(mod = mod_nms, y = y_lst))

 
# perform k-fold cross-validation to find optimal lambda value


### example usage
# try global cog first
set.seed(4566)
y_i <- y_lst[[1]]
sub_i <- sort(sample(length(y_i), 1000))
table(is.na(y_i[sub_i]))
tic()
example_mod_1 <- fit_cv_grplasso(y_i[sub_i], X_ = X[sub_i,], folds_ = kfold_ind[sub_i])
toc()  

# now an outcome with NAs to comlete the testinf
set.seed(4566)
y_i <- y_lst[[3]]
sub_i <- sort(sample(length(y_i), 1000))
table(is.na(y_i[sub_i]))
tic()
example_mod_3 <- fit_cv_grplasso(y_i[sub_i], X_ = X[sub_i,], folds_ = kfold_ind[sub_i])
toc()  





# ---- start_future_life ----

### the following is only useful if a speed up is required via parallel computation
# (helpful for us as multiple outcomes, multiple datasets (folds), multiple lambdas)


# set available memory for subprocesses
options(future.globals.maxSize = 1e3 * 1024 ^ 2) # = 1 GB


# furrr parallel workers/cores setup
# "multicore" might be better on non-Windows machines
(multisesh_threads <- min(parallelly::availableCores() - 1, 5))
plan("multisession", workers = multisesh_threads)  

### test parallel works
# test code from https://furrr.futureverse.org/
# sequential
tic()
dev_null <- map(c(2, 2, 2, 2, 2), ~Sys.sleep(.x))
toc() # ~ 10 + (some overhead )sec
# parallel: should be a fifth of the time of sequential
tic()
dev_null <- future_map(c(2, 2, 2, 2, 2), ~Sys.sleep(.x))
toc() # ~ 2 + (more overhead) sec



# ---- live_that_future_life ----

# takes ~ 11 min on R9 5900X machine 
tic()
mod_df <-
  mod_df %>%
  mutate(
    grp_lasso_obj =
      future_map(
        .x = y,
        .f = fit_cv_grplasso # \(x) fit_cv_grplasso(y = x)
      )
  )
toc()

mod_df


# ---- end_future_life ----

plan("sequential") # back to normal programming. double meaning, happy with that


# ---- mod_diagnostics_extract ----




# test extraction of elements
str(mod_df$grp_lasso_obj[[1]])
mod_df$grp_lasso_obj[[3]]$lambda.min



# find optimal lambdaz valuez that minimizes test MSE
# R^2 values over lambdaz valuez too
pdf("fig/grplasso_cv_error_5_outc.pdf", width = 10, height = 12)
par(mfrow = c(5, 2))
for (i in seq_along(mod_df$grp_lasso_obj)) {
  plot(mod_df$grp_lasso_obj[[i]], sub = mod_df$mod[i])
  plot(mod_df$grp_lasso_obj[[i]], sub = mod_df$mod[i], type = "rsq")
}
par(mfrow = c(1, 1))
dev.off()




get_grp_lasso_cv_results <- function(mod_obj) {
  tibble(
    lambda_argmin_mse = mod_obj$lambda.min,
    ln_lambda = log(mod_obj$lambda),   
    lambda = mod_obj$lambda,   
    cvmse =  mod_obj$cve,      
    cvmse_se = mod_obj$cvse  
  )
}


mod_df <-
  mod_df %>%
  mutate(
    cv_df = 
      map(
        .x = grp_lasso_obj,
        .f = get_grp_lasso_cv_results
      )
  )

mod_df

mod_df <-
  mod_df %>%
  mutate(
    mod = 
      case_when(
        mod == "cv_grplasso_globalcog" ~ "Global cognition",
        mod == "cv_grplasso_memory"    ~ "Memory",
        mod == "cv_grplasso_procspeed" ~ "Processing speed",
        mod == "cv_grplasso_execfunc"  ~ "Executive function",
        mod == "cv_grplasso_reasoning" ~ "Reasoning",
        TRUE               ~ "ERROR!!"
      ),
    mod = fct_inorder(mod)
  ) 
  
plt_mod_df <-
  mod_df %>%
  select(mod, cv_df)  %>%
  unnest(cols = cv_df)


plt_mod_df2 <-
  mod_df %>%
  select(mod, grp_lasso_obj  )%>%
  mutate(summ_stats_obj = map(.x = grp_lasso_obj, .f = summary))

glimpse(plt_mod_df2$summ_stats_obj[[1]])

get_summ_stats_tib <- function(grp_lasso_summary_obj) {
  min_ii <- grp_lasso_summary_obj$min
  tibble(
    min_lambda = grp_lasso_summary_obj$lambda[min_ii],
    min_ln_lambda = log(grp_lasso_summary_obj$lambda[min_ii]),
    cve = grp_lasso_summary_obj$cve[min_ii],
    r2 = grp_lasso_summary_obj$r.squared[min_ii],
    sigma = grp_lasso_summary_obj$sigma[min_ii],
    snr = grp_lasso_summary_obj$snr[min_ii],
    ngroups = grp_lasso_summary_obj$ngroups[min_ii],
    nvars = grp_lasso_summary_obj$nvars[min_ii],
    max_cve = max(grp_lasso_summary_obj$cve),
    first_ln_lam = log(min(grp_lasso_summary_obj$lambda))
  )
}

plt_mod_df2 <-
  plt_mod_df2 %>% 
  mutate(summ_stats = map(.x = summ_stats_obj, .f = get_summ_stats_tib)) %>% 
  select(mod, summ_stats  )%>%
  unnest(cols = summ_stats)

plt_mod_df2 %>% 
  select(-max_cve, -first_ln_lam) %>% 
  knitr::kable(., digits = c(0, 4, 1, 3, 3, 3, 3, 0, 0))

# |mod                | min_lambda| min_ln_lambda|   cve|    r2| sigma|   snr| ngroups| nvars|
# |:------------------|----------:|-------------:|-----:|-----:|-----:|-----:|-------:|-----:|
# |Global cognition   |     0.0022|          -6.1| 0.882| 0.118| 0.939| 0.134|      76|   279|
# |Memory             |     0.0054|          -5.2| 0.952| 0.048| 0.976| 0.050|      41|   128|
# |Processing speed   |     0.0030|          -5.8| 0.871| 0.129| 0.933| 0.148|      67|   240|
# |Executive function |     0.0036|          -5.6| 0.798| 0.202| 0.893| 0.254|      56|   193|
# |Reasoning          |     0.0025|          -6.0| 0.851| 0.149| 0.922| 0.175|      63|   239|




plt_mod_df  %>%
  # pivot_longer(
  #   cols = c(cvmse, cvmse_se)
  # ) %>%
  ggplot(aes(x = ln_lambda , y = cvmse, group = mod, col = mod)) +
  geom_vline(
    data =
      plt_mod_df %>% 
      distinct(mod, lambda_argmin_mse) %>% 
      mutate(ln_lambda_argmin_mse = log(lambda_argmin_mse), name = "cvmse"), 
    aes(xintercept = ln_lambda_argmin_mse, col = mod),
    linewidth = 1/2,
    alpha = 1
  ) +
  geom_line(alpha = 4/5) +
  geom_point(alpha = 2/5) +
  geom_ribbon( 
    aes(fill = mod, ymin = cvmse - 2 * cvmse_se, ymax = cvmse + 2 * cvmse_se), 
    alpha = 1/4,
    col = NA
  ) +
  # scale_y_log10() + 
  geom_label(
    data = plt_mod_df2,
    aes(
      x = min_ln_lambda, 
      y = max_cve, 
      col = mod, 
      label = sprintf("R^2  ==  %01.3f", r2)
    ),
    parse = TRUE,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_colour_tableau()  +
  scale_fill_tableau() +
  # scale_x_continuous(breaks = c(-9, -6, -3), labels = exp(c(-9, -6, -3))) +
  facet_wrap(~ mod, ncol = 1, scales = "free_y") +
  labs(
    col = "Outcome\nvariable",
    fill = "Outcome\nvariable",
    y = "Ten-fold cross-validation\nmean squared error with 95% confidence intervals",
    x = expression(Logged~value~of~regularisation~penalty~term:~ln(lambda))
  )


ggsave(filename = "fig/grplasso_cv_error_common_axes.pdf", width = 10, height = 8)




# ---- mod_coefs_extract ----


get_grp_lasso_coefs <- function(mod_obj) {
  # extract coefficients at a single value of lambda
  coefs_mat <- as.matrix(coef(mod_obj, s = mod_obj$lambda.min))
  tibble(
    coef_nm = row.names(coefs_mat), 
    coef_val = coefs_mat[, 1]
  )

}

tic()
mod_df <-
  mod_df %>%
  mutate(
    coef_df = 
      map(
        .x = grp_lasso_obj,
        .f = get_grp_lasso_coefs
      )
  )
toc()

mod_df


# save models and associated data for further diagnostics etc
write_rds(file = off_world_file("res/grplasso_cv_purrr_tbl.rds"))


heat_mod_df <-
  mod_df %>%
  select(mod, coef_df) %>%
  unnest(cols = coef_df)

heat_mod_df


nrow(heat_mod_df)
heat_mod_df <-
  heat_mod_df %>%
  inner_join(
    .,
    x_dsgn_meta %>%
      mutate(
        var = if_else(is.na(var), var_and_lvl, var),
        var_order = if_else(is.na(var_order), 0L, var_order) ,
      ) %>% 
      rename(coef_var_grp = var),
    c("coef_nm" = "var_and_lvl")
  )
nrow(heat_mod_df)


# save model coefficients for future use  
# (as it is a simple matrix we can use arrow::write_parquet)
write_parquet(heat_mod_df, sink = off_world_file("res/grplasso_cv_coef.parquet"))



