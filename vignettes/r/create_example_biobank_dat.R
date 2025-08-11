

# ---- libs ----

library("compositions")
library("simplexity")
library("mvtnorm")
library("tidyheatmaps") 

library("dplyr")
library("tidyr")
library("readr")
library("forcats")
library("readr")
library("purrr")
library("stringr")
library("ggplot2")
library("ggthemes")
library("knitr")




# ---- consts ----

(is_interactv <- interactive())
select <- dplyr::select 
filter <- dplyr::filter


# ---- helper_functions ----

source("vignettes/r/aux_functions.R")


# ---- real_data ----


### NOTE:
# there is none here however, the file 
# /vignettes/fig/cor_heatmap_real.pdf
# contains a figure containing the correlation values found in the real data

# Also, the files:
# /vignettes/dat/mu_bb.rds, and
# /vignettes/dat/sigma_bb.rds
#  contain the mean vector and var-covar matrix from the biobank data




# ---- gen_fake ----

# how many observations do we want? 100k sounds good
n_gen <- 1e+5 

# the mean-variance of the real data (and assumed multivariate Gaussian)
(mu_bb <- read_rds("vignettes/dat/mu_bb.rds"))
(sigma_bb <- read_rds("vignettes/dat/sigma_bb.rds"))

# These are cumulative proportion of records in each category
# (for converting continuous values to discrete while keep cor structure)
sex_quants <- c(
  "Female" = 0.5671448, 
  "Male" = 1.0000000
)
smk_quants <- c(
  "never" = 0.5754566, 
  "previous" = 0.9328081, 
  "current" = 0.9977760, 
  "unknown" = 1.0000000
)
    

# randomly generate
set.seed(12349876)
fake_bb <- rmvnorm(n = n_gen, mean = mu_bb, sigma = sigma_bb)


colnames(fake_bb) <- names(mu_bb)
fake_bb <-
  fake_bb %>%
  as_tibble()

fake_bb


### discretise factors

(sex_cts_qs <- quantile(fake_bb$sex, sex_quants[-length(sex_quants)]))
(smk_cts_qs <- quantile(fake_bb$smoking, smk_quants[-length(smk_quants)]))


quantile(fake_bb$sex, seq(0, 1, 0.25))
quantile(fake_bb$smoking, seq(0, 1, 0.25))


fake_bb$sex <- 
  cut(
    fake_bb$sex, 
    breaks = c(-Inf, sex_cts_qs, +Inf), 
    labels = names(sex_quants)
  )
fake_bb$smoking <- 
  cut(
    fake_bb$smoking, 
    breaks = c(-Inf, smk_cts_qs, +Inf), 
    labels = names(smk_quants)
  )

levels(fake_bb$sex)
levels(fake_bb$smoking)
table(fake_bb$sex, fake_bb$smoking, useNA = "ifany")
table(fake_bb$sex, useNA = "ifany")
table(fake_bb$smoking, useNA = "ifany")



# ---- plot_cors ----

plt_fake_df <-
  fake_bb %>%  
  mutate(
    sex = rank(sex),
    smoking = rank(smoking),
  ) %>%
  cor(.) %>%
  as.data.frame(.) %>%
  mutate(
    rwnm = rownames(.),
  ) %>%
  pivot_longer(
    cols = -c(rwnm),
    names_to = "clnm",
    values_to = "cor"
  )

plt_fake_df <-
  plt_fake_df %>%
  mutate(
    cor = if_else(rwnm == clnm, as.numeric(NA), cor) ,
    is_outc_rw = if_else(rwnm == "glb_cog", "Yas", "Nah"),
    is_outc_cl = if_else(clnm == "glb_cog", "Yas", "Nah")
  )

(max_abs <- max(abs(plt_fake_df$cor), na.rm = TRUE))

cor_plot <-
  tidyheatmap(
    df = plt_fake_df,
    rows = rwnm,
    columns = clnm,
    values = cor,
    colors = scl_col, 
    scale = "none",
    color_na = "white",
    color_legend_min = -max_abs,
    color_legend_max = max_abs,
    show_rownames = TRUE,
    gaps_row = is_outc_rw,
    gaps_col = is_outc_cl,
    annotation_row = c(is_outc_rw),
    annotation_col = c(is_outc_cl),
    annotation_colors = ann_colors,
    cellwidth = 50,
    cellheight = 50,
    color_legend_n = 11,
    fontsize_row = 8,
    fontsize_number = 8,
    display_numbers = TRUE,
    number_format = "%0.2f",
    number_color = "white",
    filename = "vignettes/fig/cor_heatmap_fake.pdf"
  )

cor_plot


# ---- finalise_fake ----


ilr_df <-
  fake_bb %>%
  select(starts_with("ilr"))

colnames(ilr_df) <- 
  str_replace(colnames(ilr_df), "^ilr\\.", "")


tu_labs <- c("sleep", "sb", "lpa", "mvpa")
ilr_labs <- colnames(ilr_df)

# this was the original ilr basis used, so will use to back transform
(sbp4_b0 <-
  matrix(
    c(
      +1,  0,  0,
      -1, +1,  0,
      -1, -1, +1,
      -1, -1, -1
    ),
    byrow = TRUE,
    ncol = 3,
    dimnames = list(tu_labs, ilr_labs)
  ))
(psi4_b0 <- compositions::gsi.buildilrBase(sbp4_b0))



tu_df <-
  ilr_df %>%
  as.data.frame(.) %>%
  ilrInv(., V = psi4_b0) %>%
  as.data.frame(.)

tu_df <- 1440 * tu_df
head(tu_df)

fake_bb_timeuse <-
  fake_bb %>%
  select(-starts_with("ilr")) %>%
  bind_cols(., tu_df)


fake_bb_timeuse


# ---- save_fake_data ----

write_rds(fake_bb_timeuse, file = "vignettes/dat/fake_bb.rds")





