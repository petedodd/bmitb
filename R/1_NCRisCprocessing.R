
## script to process NCD RisC data and fit distribution parameters

## raw data ncdris data is from: https://ncdrisc.org/data-downloads-adiposity.html
## stored in rawdata/ are the following files:
## NCD_RisC_Lancet_2024_BMI_child_adolescent_country.csv
## NCD_RisC_Lancet_2024_BMI_female_age_specific_country.csv
## NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv

## libraries
library(here)
library(data.table)
library(readxl)
library(ggplot2)

## === WHO refs for children
## want to know: BMIs for z scores -2,-1,+1,+2 in ref pop for sex and age
## from: https://www.who.int/tools/growth-reference-data-for-5to19-years/indicators/bmi-for-age
fn <- here("rawdata/bmi-boys-z-who-2007-exp.xlsx")
WB <- read_excel(fn)
fn <- here("rawdata/bmi-girls-z-who-2007-exp.xlsx")
WG <- read_excel(fn)
WB <- as.data.table(WB)
WG <- as.data.table(WG)
WG[, Sex := "Girls"]
WB[, Sex := "Boys"]
BG <- rbind(WB, WG) # combined version


## === read in raw data
## --- adults
fn <- here("rawdata/NCD_RisC_Lancet_2024_BMI_female_age_specific_country.csv") # female
F <- fread(fn) # NOTE not in repo as large
fn <- here("rawdata/NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv") # males
M <- fread(fn) # NOTE not in repo as large
B <- rbind(M, F) # both

## --- adolescents
fn <- here("rawdata/NCD_RisC_Lancet_2024_BMI_child_adolescent_country.csv")
A <- fread(fn) #NOTE not in repo as large

## compare variables
setdiff(names(A), names(M))
setdiff(names(M), names(A))

## compare ages
A[, unique(`Age group`)]
B[, unique(`Age group`)]

## compare years
A[, unique(Year)]
B[, unique(Year)] # both 1990 - 2022
A[ISO == "AFG" & Year == 2022]

## available for all
A[
  ,
  table(
    Year,
    is.na(`Prevalence of BMI < -2SD (thinness) lower 95% uncertainty interval`)
  )
]


## restrict data
DR <- A[Year == 2022] # so the distribution info is available before 2017
DR <- DR[`Age group` >= 5]
(keep <- names(DR)[c(
  1, 2, 3, 5,
  6, 7, 8,
  9, 10, 11,
  18, 19, 20,
  21, 22, 23,
  24, 25, 26
)])
DR <- DR[, ..keep]
## new names
nnmz <- c(
  "Year", "Sex", "age", "iso3",
  "lt.n2", "lt.n2.lo", "lt.n2.hi",
  "gt.p2", "gt.p2.lo", "gt.p2.hi",
  "n2.n1", "n2.n1.lo", "n2.n1.hi",
  "n1.p2", "n1.p2.lo", "n1.p2.hi",
  "p1.p2", "p1.p2.lo", "p1.p2.hi"
)
cbind(names(DR), nnmz) # CHECK
names(DR) <- nnmz # rename

## SD
DR[, lt.n2.s := abs(lt.n2.lo - lt.n2.hi) / 3.92]
DR[, gt.p2.s := abs(gt.p2.lo - gt.p2.hi) / 3.92]
DR[, n2.n1.s := abs(n2.n1.lo - n2.n1.hi) / 3.92]
DR[, n1.p2.s := abs(n1.p2.lo - n1.p2.hi) / 3.92]
DR[, p1.p2.s := abs(p1.p2.lo - p1.p2.hi) / 3.92]

drop <- grep("lo|hi", names(DR), value = TRUE)
DR[, (drop) := NULL]


## target names
nnmz <- names(DR)
(tgt <- nnmz[5:9]) # targets
(tgt.s <- nnmz[5 + 5:9]) # targets


## === fitting to gamma distributions
## --- adolescents

## fit to gamma
gammaStats <- function(k,
                       theta,
                       b2112 # BMI refs for +2,+1,-1,-2 SD
) {
  sts <- c(
    pgamma(b2112[4], shape = k, scale = theta), # low tails:  <= -2Z
    1 - pgamma(b2112[1], shape = k, scale = theta), # high tails: >= +2Z
    pgamma(b2112[3], shape = k, scale = theta) -
      pgamma(b2112[4], shape = k, scale = theta), # intervals: -2Z to -1Z
    pgamma(b2112[2], shape = k, scale = theta) -
      pgamma(b2112[3], shape = k, scale = theta), # intervals: -1Z to +1Z
    pgamma(b2112[1], shape = k, scale = theta) -
      pgamma(b2112[2], shape = k, scale = theta) # intervals: +1Z to +2Z
  )
  names(sts) <- tgt
  sts
}

## test
ba5 <- unlist(WB[Month == 5 * 12 + 6, .(SD2, SD1, SD1neg, SD2neg)])
gammaStats(2, 9, ba5)


## loss function
gamErr <- function(x, tgt, tgt.s, ref) {
  x <- exp(x)
  res <- gammaStats(x[1], x[2], ref)
  sum(((res - tgt) / (1e-10 + tgt.s))^2) # ~SSE normal approx
}

## test
TGT <- unlist(DR[1, ..tgt])
TGT.S <- unlist(DR[1, ..tgt.s])
gammaStats(1, 1, ba5)
gamErr(c(0, 0), TGT, TGT.S, ba5)

## optimize
(out <- optim(par = c(0, 0), fn = function(x) gamErr(x, TGT, TGT.S, ba5)))
exp(out$par)
gammaStats(exp(out$par[1]), exp(out$par[2]), ba5)
TGT

## check
curve(
  dgamma(x,
    shape = exp(out$par[1]),
    scale = exp(out$par[2])
  ),
  from = 10, to = 30, n = 1e3
)

##  ---apply this fitting across the data:
## common age to merge with WHO REFs
DR[, Month := 12 * age] # low-point (otherwise miss 19-20)
DR[Month == 60, Month := 61] # otherwise missing from BG
DR[, range(Month)]
BG[, range(Month)]
akey <- unique(DR[, .(age, Month)])

## merge:
DR <- merge(DR, BG, by = c("Sex", "Month"), all.x = TRUE)
unique(DR[, .(age, Month)]) # CHECK

## restrict to age >=15 (NOTE)
DR <- DR[age >= 15]

## Loop above work flow
DRA <- DR[,
  {
    print(iso3)
    print(Sex)
    print(Month)
    TGT <- c(lt.n2, gt.p2, n2.n1, n1.p2, p1.p2) # NCD RisC targets
    TGT.S <- c(lt.n2.s, gt.p2.s, n2.n1.s, n1.p2.s, p1.p2.s) # NCD RisC targets
    saref <- c(SD2, SD1, SD1neg, SD2neg) # BMIs in WHO reference data
    out <- optim(
      par = c(0, 0),
      fn = function(x) gamErr(x, TGT, TGT.S, saref),
      hessian = TRUE
    )
    H <- out$hessian # wrt log parms
    H <- H / exp(out$par) # rows
    H <- t(t(H) / exp(out$par)) # cols (now hessian for real parms)
    V <- solve(H) # varcov for real parms
    if (abs(out$convergence) > 0) print("*** has not converged! ***")
    list(
      k = exp(out$par[1]), theta = exp(out$par[2]),
      Vkk = V[1, 1], Vkt = V[1, 2], Vtt = V[2, 2],
      cvgc = out$convergence
    )
  },
  by = .(iso3, Month, Sex, age)
]


## examine convergence problems
DRA[, table(cvgc)] # 6/2000 issues
(pbms <- DRA[cvgc != 0])
DRA[cvgc != 0, unique(iso3)]
## "MAR" "TZA" "POL" "ROU" "SOM" "TKM"
DRA[cvgc != 0, .N, by = iso3] # all 1 age: use averages that exclude these

## odd Vs: use country/sex means
DRA[Vkk < 0 & cvgc == 0] # "NRU" "NIU" "MWI" "ZMB" "TON" "COK" "UKR"
DRA[Vtt < 0 & cvgc == 0] # same 7 values with 1 cat each, mainly girls
tmp <- DRA[Vkk > 0  & cvgc == 0,
  .(Vkk.tmp = mean(Vkk), Vtt.tmp = mean(Vtt), Vkt.tmp = mean(Vkt)),
  by = .(iso3, Sex)
]
DRA <- merge(DRA, tmp, by = c("iso3", "Sex"))
DRA[
  Vkk < 0 & cvgc == 0,
  c("Vkk", "Vkt", "Vtt") :=
    .(Vkk.tmp, Vkt.tmp, Vtt.tmp)
]
DRA[, c("Vkk.tmp", "Vkt.tmp", "Vtt.tmp") := NULL]

## odd means: 4 iso3/values :  "COK" "NIU" "NRU" "TON"
DRA[, mnbmi := k * theta]
DRA[mnbmi > 50]
tmp <- DRA[mnbmi <= 50,
  .(
    k.tmp = mean(k), theta.tmp = mean(theta),
    Vkk.tmp = mean(Vkk), Vtt.tmp = mean(Vtt), Vkt.tmp = mean(Vkt)
  ),
  by = .(iso3, Sex)
]
DRA <- merge(DRA, tmp, by = c("iso3", "Sex"))
DRA[
  mnbmi > 50,
  c("k", "theta", "Vkk", "Vkt", "Vtt") :=
    .(k.tmp, theta.tmp, Vkk.tmp, Vkt.tmp, Vtt.tmp)
]
DRA[, c("k.tmp", "theta.tmp", "mnbmi", "Vkk.tmp", "Vkt.tmp", "Vtt.tmp") := NULL]

## save out
fn <- here("data/DRA.Rdata")
save(DRA, file = fn)

## ======= adult data

## available for all
B[, table(Year, is.na(`Prevalence of BMI >=40 kg/m² (morbid obesity)`))]


## restrict
BR <- B[Year == 2022] # distribution info is available before 2017
## names(BR)
keep <- names(BR)[c(
  1, 2, 4, 5,
  6, 6 + 1, 6 + 2,
  9, 9 + 1, 9 + 2,
  18, 18 + 1, 18 + 2,
  21, 21 + 1, 21 + 2,
  24, 24 + 1, 24 + 2,
  27, 27 + 1, 27 + 2,
  30, 30 + 1, 30 + 2,
  33, 33 + 1, 33 + 2
)]
BR <- BR[, ..keep]

## new names
nnmz <- c(
  "Year", "Sex", "iso3", "age",
  "l185", "l185.lo", "l185.hi",
  "g30", "g30.lo", "g30.hi",
  "l20.g185", "l20.g185.lo", "l20.g185.hi",
  "l25.g20", "l25.g20.lo", "l25.g20.hi",
  "l30.g25", "l30.g25.lo", "l30.g25.hi",
  "l35.g30", "l35.g30.lo", "l35.g30.hi",
  "l40.g35", "l40.g35.lo", "l40.g35.hi",
  "g40", "g40.lo", "g40.hi"
)

cbind(names(BR), nnmz) # CHECK
names(BR) <- nnmz # rename

## SD
BR[, l185.s := abs(l185.lo - l185.hi) / 3.92]
BR[, g30.s := abs(g30.lo - g30.hi) / 3.92]
BR[, l20.g185.s := abs(l20.g185.lo - l20.g185.hi) / 3.92]
BR[, l25.g20.s := abs(l25.g20.lo - l25.g20.hi) / 3.92]
BR[, l30.g25.s := abs(l30.g25.lo - l30.g25.hi) / 3.92]
BR[, l35.g30.s := abs(l35.g30.lo - l35.g30.hi) / 3.92]
BR[, l40.g35.s := abs(l40.g35.lo - l40.g35.hi) / 3.92]
BR[, g40.s := abs(g40.lo - g40.hi) / 3.92]


drop <- grep("lo|hi", names(BR), value = TRUE)
BR[, (drop) := NULL]

## target names
nnmz <- names(BR)
(tgt <- nnmz[5:12]) # targets
(tgt.s <- nnmz[13:20]) # targets


## fit to gamma
AgammaStats <- function(k, theta) {
  sts <- c(
    pgamma(18.5, shape = k, scale = theta), # low tails:  <= 18.5
    1 - pgamma(30, shape = k, scale = theta), # high tails: >= 30
    pgamma(20, shape = k, scale = theta) -
      pgamma(18.5, shape = k, scale = theta), # intervals
    pgamma(25, shape = k, scale = theta) -
      pgamma(20, shape = k, scale = theta), # intervals
    pgamma(30, shape = k, scale = theta) -
      pgamma(25, shape = k, scale = theta), # intervals
    pgamma(35, shape = k, scale = theta) -
      pgamma(30, shape = k, scale = theta), # intervals
    pgamma(40, shape = k, scale = theta) -
      pgamma(35, shape = k, scale = theta), # intervals
    1 - pgamma(40, shape = k, scale = theta) # high tails: >= 40
  )
  names(sts) <- tgt
  sts
}


## test
AgammaStats(2, 9)


## loss function
AgamErr <- function(x, tgt, tgt.s) {
  x <- exp(x)
  res <- AgammaStats(x[1], x[2])
  sum(((res - tgt) / (1e-10 + tgt.s))^2) # SSE
}


## test
TGT <- unlist(BR[1, ..tgt])
TGT.S <- unlist(BR[1, ..tgt.s])
AgamErr(c(0, 0), TGT, TGT.S)
## optimize
(out <- optim(par = c(0, 0), fn = function(x) AgamErr(x, TGT, TGT.S)))
exp(out$par)
AgammaStats(exp(out$par[1]), exp(out$par[2]))
TGT

## check
curve(dgamma(x, shape = exp(out$par[1]), scale = exp(out$par[2])),
  from = 10, to = 40, n = 1e3
)

## restrict (NOTE)
BR <- BR[age != "18-19"]

##  ---apply this fitting across the data:
## Loop above work flow
DRB <- BR[,
  {
    cat(iso3, ", ", Sex, ", ", age, ", ", Year, "\n")
    TGT <- c(
      l185, g30, l20.g185, l25.g20, l30.g25,
      l35.g30, l40.g35, g40
    )
    TGT.S <- c(
      l185.s, g30.s, l20.g185.s, l25.g20.s, l30.g25.s,
      l35.g30.s, l40.g35.s, g40.s
    )
    out <- optim(
      par = c(0, 0),
      fn = function(x) AgamErr(x, TGT, TGT.S),
      hessian = TRUE
    )
    H <- out$hessian # wrt log parms
    H <- H / exp(out$par) # rows
    H <- t(t(H) / exp(out$par)) # cols (now hessian for real parms)
    V <- solve(H) # varcov for real parms
    if (abs(out$convergence) > 0) print("*** has not converged! ***")
    list(
      k = exp(out$par[1]), theta = exp(out$par[2]),
      Vkk = V[1, 1], Vkt = V[1, 2], Vtt = V[2, 2],
      cvgc = out$convergence
    )
  },
  by = .(iso3, Year, Sex, age)
]

## fix convergence problems
DRB[, table(cvgc)] # 1 / 5600
DRB[cvgc != 0]     #men 40-44
DRB[cvgc != 0, unique(iso3)] # "CHN"

## for these groups, use the same country and sex, but age group below
## CHN
DRB[
  iso3 == "CHN" & Sex == "Men" & age == "40-44",
  c("k", "theta", "Vkk", "Vkt", "Vtt") := DRB[
    iso3 == "CHN" & Sex == "Men" & age == "35-39",
    .(k, theta, Vkk, Vkt, Vtt)
  ]
]

## odd Vs: use country/sex means
DRB[Vkk < 0] # "NRU" "WSM" "COK" "TKL" "COK" "WSM" "ASM"
DRB[Vtt < 0] # same 7 values with 1 cat each
tmp <- DRB[Vkk > 0,
  .(Vkk.tmp = mean(Vkk), Vtt.tmp = mean(Vtt), Vkt.tmp = mean(Vkt)),
  by = .(iso3, Sex)
]
DRB <- merge(DRB, tmp, by = c("iso3", "Sex"))
DRB[Vkk < 0, c("Vkk", "Vkt", "Vtt") := .(Vkk.tmp, Vkt.tmp, Vtt.tmp)]
DRB[, c("Vkk.tmp", "Vkt.tmp", "Vtt.tmp") := NULL]

## odd means: 7 values 5 iso3 mainly men: "ASM" "COK" "NRU" "TKL" "WSM"
DRB[, mnbmi := k * theta]
DRB[mnbmi > 50]
tmp <- DRB[mnbmi <= 50,
  .(
    k.tmp = mean(k), theta.tmp = mean(theta),
    Vkk.tmp = mean(Vkk), Vtt.tmp = mean(Vtt), Vkt.tmp = mean(Vkt)
  ),
  by = .(iso3, Sex)
]
DRB <- merge(DRB, tmp, by = c("iso3", "Sex"))
DRB[
  mnbmi > 50,
  c("k", "theta", "Vkk", "Vkt", "Vtt") :=
    .(k.tmp, theta.tmp, Vkk.tmp, Vkt.tmp, Vtt.tmp)
]
DRB[, c("k.tmp", "theta.tmp", "mnbmi", "Vkk.tmp", "Vkt.tmp", "Vtt.tmp") := NULL]

## CHECK
DRB[cvgc != 0]

## save
fn <- here("data/DRB.Rdata")
save(DRB, file = fn)
