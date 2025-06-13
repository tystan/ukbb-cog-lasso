
# install.packages("pdftools")


library("pdftools")

(pdf_ver <- list.files(path = "./fig/flowchart", pattern = "\\.pdf", full.names = TRUE))
(out_fl <- gsub("\\.pdf$", ".png", pdf_ver))

pdf_convert(
  pdf_ver,
  format = "png",
  pages = 1,
  filenames = out_fl,
  dpi = 300,
  antialias = TRUE,
  opw = "",
  upw = "",
  verbose = TRUE
)


