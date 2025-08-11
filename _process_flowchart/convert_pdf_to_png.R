
# install.packages("pdftools")
library("pdftools")

(pdf_ver <- list.files(path = "_process_flowchart", pattern = "\\.pdf", full.names = TRUE))
(out_fl <- gsub("\\.pdf$", ".png", pdf_ver))

pdf_convert(
  pdf_ver,
  format = "png",
  pages = 1,
  filenames = out_fl,
  dpi = 500,
  antialias = TRUE,
  opw = "",
  upw = "",
  verbose = TRUE
)


