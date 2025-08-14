## =======================================
## LOADING TO GOOGLE SHEETS (authors only)
## =======================================
rm(list = ls())

library(here)
library(data.table)
library(glue)
library(googlesheets4)

## setup - only accessible to those with access to this sheet NOTE new
yourl <- "https://docs.google.com/spreadsheets/d/1epQis4hhJMlk7ggS7kmrVWSiQMCn9yZpRH_is4SFftY/edit?gid=0#gid=0"
shid <- as.character(as_sheets_id(yourl))


## utility function
upload.to.sheets <- function(filename, sheetid) {
  fn <- glue(here("output/{filename}"))
  tmp <- fread(file = fn)
  sht <- gsub("\\.csv", "", filename)
  write_sheet(tmp, sheetid, sheet = sht)
}


## read & upload relevant data
upload.to.sheets("table1.csv", shid)

upload.to.sheets("atable_BMI.csv", shid)
upload.to.sheets("gt25pc_tab.csv", shid)
upload.to.sheets("gt25pc.csv", shid)
upload.to.sheets("RRbySR.csv", shid)

upload.to.sheets("BbyASM.csv", shid)
upload.to.sheets("outstats.csv", shid)


upload.to.sheets("all_country_reductions.csv", shid)
upload.to.sheets("atable_BMI_pop.csv", shid)
upload.to.sheets("atable_BMI_pc.csv", shid)

upload.to.sheets("RRbyASR.csv", shid)
upload.to.sheets("RRbyAS.csv", shid)
