library(tiff)
library(tidyverse)

tiff.name <- list.files("tiffs", full.names = TRUE)

img.list <- lapply(tiff.name, function(tt) {
  xx <- readTIFF(source = tt)
  (xx * 255) %% 2
})

img.list -> LSBs_small
save(LSBs_small, file = "LSBs_small.rda")
LSBs_small[[2]]
