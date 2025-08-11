

# ---- consts ----


scl_col <- c("#C34B6B","lightyellow2", "#5071DC") #, "dodgerblue")
scl_col_alt <- c("#C34B6B","grey20", "#5071DC") 
is_outc_col = c("Yas" = "turquoise", "Nah" = "orange")

ann_colors <- list(
  is_outc_rw = is_outc_col,
  is_outc_cl = is_outc_col
)

# ---- funcs ----

get_rotate_mat <- function(from_v, to_v) {
  return(t(to_v) %*% from_v)
}

euclid_norm <- function(x) {
  sqrt(sum(x^2))
}

rotate_ilrs <- function(ilrs_from, R) {
  # keep as matrix data
  as.matrix(ilrs_from) %*% t(R)
  # as.data.frame(as.matrix(ilrs_from) %*% t(R))
}



orthog_scale_grp <- function(dat, grp_cns = NULL) {
  
  n_c <- ncol(dat)
  cns <- colnames(dat)
  if (!is.null(grp_cns)) {
    cns <- grp_cns
  } else if (is.null(cns)) {
    cns <- paste0("v", 1:n_c)
  }
  
  dat <- as.matrix(dat)
  (m_grp <- colMeans(dat))
  (v_grp <- var(dat))
  cis <- covar_inverse_sqrt(v_grp)
  
  osg_dat <- t(t(dat) - m_grp) %*% cis
  colnames(osg_dat) <- cns
  
  return(list(
    osg_dat = osg_dat,
    m_grp = m_grp,
    v_grp = v_grp
    # cis = cis
  ))
  
}


orthog_scale_grp_inv <- function(dat_sc, m_grp, v_grp, grp_cns = NULL) {
  
  n_c <- ncol(dat_sc)
  
  cns <- colnames(dat_sc)
  if (!is.null(grp_cns)) {
    cns <- grp_cns
  } else if (is.null(cns)) {
    cns <- paste0("v", 1:n_c)
  }
  
  dat_sc <- as.matrix(dat_sc)
  # as.vector(matrix(0, nrow = 3), mode = "double")
  # as.vector(matrix(0, nrow = 3), mode = "numeric")
  stopifnot(is.matrix(v_grp))
  m_grp <- as.vector(m_grp, mode = "double")
  
  # v_grp_sqrt <- expm::sqrtm(v_grp)
  cis <- covar_inverse_sqrt(v_grp)
  
  dat <- t(t(dat_sc %*% cis %*% v_grp) + m_grp)
  colnames(dat) <- cns
  
  return(dat)
  
}


# ---- ilr_funcs -----

sanitise_ilrs <- function(x) {
  
  if ("rmult" %in% class(x)) {
    class(x) <- NULL # remove "rcomp" class, will either result in numeric vector or matrix
    attr(x, "orig") <- NULL # remove original composition info (issue with indexes)
    attr(x, "V") <- NULL # clear orthonormal basis attr
    if ("numeric" %in% class(x)) { # if vector turn into 1 row matrix
      x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
    }
  }
  stopifnot(is.matrix(x))
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
