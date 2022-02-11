library(tidyverse)
library(arrow)

acousticIndicesFiles <- list.files(path = "rawdata/acousticindices",
                                   pattern = "*.csv",
                                   full.names = TRUE,
                                   recursive = TRUE)

for (file in acousticIndicesFiles) {
  tmp <- read_csv(file)
  
  write_parquet(x = tmp,
                sink = paste0("rawdata/acousticindices/", "test"))
}