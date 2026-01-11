## =======================================
## LOADING TO GOOGLE SHEETS (authors only)
## =======================================
rm(list = ls())

library(here)
library(data.table)
library(glue)
library(googlesheets4)

## output formatting
source(here("R/utils/brackets.R"))

## ==== gathering alt CFs ===
CF <- c("", "_Blo", "_Bhi", "_Clo", "_Chi", "D")
N <- P <- list()
for (fend in CF) {
  D <- fread(gh("output/table1_r{fend}.csv"))
  D[, counterfactual := fend]
  N[[fend]] <- D
  D <- fread(gh("output/RRbySR_r{fend}.csv"))
  D[, counterfactual := rep(fend, nrow(D))]
  P[[fend]] <- D
}
N <- rbindlist(N)
P <- rbindlist(P)

N[, counterfactual_type := fcase(
  grepl("B", counterfactual), "B",
  grepl("C", counterfactual), "C",
  grepl("D", counterfactual), "D",
  default = "A"
)]
P[, counterfactual_type := fcase(
  grepl("B", counterfactual), "B",
  grepl("C", counterfactual), "C",
  grepl("D", counterfactual), "D",
  default = "A"
)]

N[, counterfactual_top := fcase(
  grepl("hi", counterfactual), "25 kg/m²",
  grepl("lo", counterfactual), "halfway to 25 kg/m²",
  default = "Not applicable"
)]
P[, counterfactual_top := fcase(
  grepl("hi", counterfactual), "25 kg/m²",
  grepl("lo", counterfactual), "halfway to 25 kg/m²",
  default = "Not applicable"
)]

PW <- P[, .(
  counterfactual_type, counterfactual_top,
  Sex, variable,
  reduction = brktpc(value, lo, hi)
)]

PW <- dcast(PW,
  counterfactual_type + counterfactual_top ~ variable + Sex,
  value.var = "reduction"
)

cftoplvls <- c("Not applicable", "halfway to 25 kg/m²", "25 kg/m²")

setcolorder(PW, neworder = c(1, 2, 4, 5, 3, 7, 8, 6))
PW$counterfactual_top <- factor(PW$counterfactual_top,
  levels = cftoplvls,
  ordered = TRUE
  )
PW <- PW[order(counterfactual_type, counterfactual_top)]
fwrite(PW, file = here("output/CF_alt_pc.csv"))

N[, c("region", "counterfactual") := NULL]
setcolorder(N, neworder = names(PW))
N$counterfactual_top <- factor(N$counterfactual_top,
  levels = cftoplvls,
  ordered = TRUE
  )
N <- N[order(counterfactual_type, counterfactual_top)]
fwrite(N, file = here("output/CF_alt_num.csv"))

## ==== uploading ===
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

## alt CFs:
upload.to.sheets("CF_alt_num.csv", shid)
upload.to.sheets("CF_alt_pc.csv", shid)
