## this uses distributional input data to calculate RRs and PAFs

## libraries
library(here)
library(data.table)
library(ggplot2)
library(ggrepel)
library(paletteer)
library(MASS)

## load relevant data
load(here("data/whokey.Rdata")) # WHO region to iso3 key
load(here("data/DRA.Rdata")) # ado BMI distributions
load(here("data/DRB.Rdata")) # adult BMI distributions
load(here("rawdata/N8523.Rdata")) # 2023 demography

TBN <- fread("rawdata/TB_notifications_2024-10-30.csv")
TB <- fread(here("rawdata/TB_burden_age_sex_2024-10-30.csv"))
whoz <- c("AFR", "AMR", "EMR", "EUR", "SEA", "WPR")
whozt <- c(
  "Africa", "The Americas",
  "Eastern Mediterranean", "Europe", "South-East Asia",
  "Western Pacific"
)
for (i in seq_along(whoz)) {
  whokey[g_whoregion == whoz[i], region := whozt[i]]
}
whokeyshort <- unique(whokey[, .(g_whoregion, region)])
whokeyshort <- rbind(
  whokeyshort,
  data.table(g_whoregion = "Global", region = "Global")
)
## uncertainty aggregation
ssum <- function(x) sqrt(sum(x^2))
## output formatting
rf <- function(x) {
  dg <- ifelse(abs(x) > 0.01 & abs(x) < 100, 2, 3)
  x2 <- signif(x, dg)
  format(
    x2,
    digits = dg,
    nsmall = 0L,
    big.mark = ",", #or " "
    justify = "right",
    drop0trailing = TRUE,
    scientific = FALSE
  )
}
brkt0 <- function(x, y, z) paste0(x, " (", y, " to ", z, ")")
brkt0(1, 2, 3)
brkt <- function(x, y, z) {
  ans <- brkt0(x, y, z)
  ans <- gsub("\\([[:space:]]+", "\\(", ans)
  ans <- gsub("to[[:space:]]+", "to ", ans)
  ans
}

## Saunders et al linear parameters
## risk per one unit increase in BMI was 14.8% (95%CI: 13.3-16.3)
t <- log(1-0.148) # risk function parameter
1-exp(t) # risk increase with 1 unit decrease

## fits from bilinear model
C <- fread(here("rawdata/general_population_piecewise_parameters.csv"))
D <- fread(here("rawdata/general_population_vcov_matrix.csv"))
## 18.0% (95%CI: 16.4-19.6) for BMI<25.0kg/m2 and 6.9% (95%CI: 4.6-9.2) for BMI>=25.0kg/m2 in
exp(C$Value[4:5]) # corresponds to above
mut <- C$Value[4:5]
t1 <- mut[1]
t2 <- mut[2]
V <- as.matrix(D[, .(slope_below_breakpoint, slope_change_above_breakpoint)])
## slope_change_above_breakpoint
## slope above = slope below + change
## a = b + c
## cov(b,a) = cov(b,b+c) = cov(b,c) + var(b)
## cov(a,a) = cov(b+c,b+c) = var(b) + var(c) + 2*cov(b,c)
W <- V
W[2, 2] <- V[2, 2] + # var(c)
  V[1, 1] + # var(b)
  2 * V[1, 2] # cov(b,c)
W[1, 2] <- W[2, 1] <- V[1, 2] + V[1, 1]
## 18.0% (95%CI: 16.4-19.6) for BMI<25.0kg/m2 and 6.9% (95%CI: 4.6-9.2) for BMI>=25.0kg/m2 in
## check -- OK: matches above
S <- exp(mvrnorm(n = 1e4, mu = mut, Sigma = W))
1e2 * (1 - colMeans(S))
1e2 * quantile(1 - S[, 1], c(0.025, 1 - 0.025))
1e2 * quantile(1 - S[, 2], c(0.025, 1 - 0.025))

## risk ratio calculator
bmirefpop <- data.table(k = 25.0, theta = 1.0) #for testing
RRfun <- function(k, theta, t) {
  (1 - t * bmirefpop$theta)^bmirefpop$k / (1 - t * theta)^k
}

## check
K <- 1e5
bmi1 <- rgamma(K, shape = 24, scale = 0.7)
bmi0 <- rgamma(K, shape = bmirefpop$k, scale = bmirefpop$theta)
## E_1[exp(t*X)] / E_0[exp(t*X)] with a shift for numerics
mean(exp(t * (bmi1 - 30))) / mean(exp(t * (bmi0 - 30)))
RRfun(24, 0.7, t) # OK

## bilinear version
## see: https://search.r-project.org/CRAN/refmans/expint/html/gammainc.html
## Γ(a,x)=Γ(a)(1−P(a,x))
## https://search.r-project.org/R/refmans/stats/html/GammaDist.html
RRfunBL0 <- function(k, theta, t1, t2) {
  x0 <- 25
  a1 <- -t1
  a2 <- -t2
  ## stuff x Γ(k,x0*(a2+1/theta))/Γ(k) =
  ## stuff x (1-P(k,x0*(a2+1/theta))) = stuff x P(k,x0*(a2+1/theta),lower=FALSE)
  lans2 <- a2 * x0 +
    pgamma(k, x0 * (a2 + 1 / theta), log = TRUE, lower = TRUE) -
    k * log(1 + a2 * theta)
  lans1 <- a1 * x0 +
    pgamma(k, x0 * (a1 + 1 / theta), log = TRUE, lower = FALSE) -
    k * log(1 + a1 * theta)
  exp(lans1) + exp(lans2)
}
## NOTE works in tests

RRfunBL <- function(k, theta, t1, t2){
 RRfunBL0(k, theta, t1, t2) /
   RRfunBL0(bmirefpop$k, bmirefpop$theta, t1, t2)
}

## check
BL <- function(x, t1, t2) {
  ans <- (x - 25)
  less <- ans < 0
  ans[less] <- t1 * ans[less]
  ans[!less] <- t2 * ans[!less]
  ans
}
xx <- seq(from = 10, 40, by = 0.1)
plot(xx, BL(xx, t1, t2), type = "l")
plot(xx, BL(xx, t1, t2 - 0.1), type = "l")

mean(exp(BL(bmi1, t1, t1))) / mean(exp(BL(bmi0, t1, t1)))
## E_1[exp(t*X)] / E_0[exp(t*X)] with a shift for numerics
mean(exp(t1 * (bmi1 - 30))) / mean(exp(t1 * (bmi0 - 30)))
RRfunBL0(24, 0.7, t1, t1) / RRfunBL0(bmirefpop$k, bmirefpop$theta, t1, t1)
RRfunBL(24, 0.7, t1, t1)

mean(exp(BL(bmi1, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))
RRfunBL(24, 0.7, t1, t2) ## OK

## correct magnitude of change?
mean(bmi1)
mean(bmi0)
exp(t * (mean(bmi1) - mean(bmi0)))
RRfun(24, 0.7, t) # about as good as one might expect

## statistics to report
outstats <- list()
ok <- 1


## merge against TB estimates
## restrict:  "15-24" included
TB <- TB[
  sex != "a" &
    !age_group %in% c("all", "0-14", "0-4", "15plus", "18plus", "5-14") &
    risk_factor == "all",
  .(iso3, sex,
    age = age_group, tb = best,
    S = (hi - lo) / 3.92
  )
]

## split 15-24
xtra <- TB[age == "15-24"]
xtra1 <- copy(xtra)
xtra2 <- copy(xtra)
xtra1[, age := "15-19"]
xtra2[, age := "20-24"]

## old version
## ## even split
## xtra1[, tb := tb / 2]
## xtra1[, S := S / sqrt(2)]
## xtra2[, tb := tb / 2]
## xtra2[, S := S / sqrt(2)]

## split using notifications
TBN <- TBN[
  year == 2023 &
    !is.na(newrel_f1519 + newrel_f2024) &
    !is.na(newrel_m1519 + newrel_m2024),
  .(
    iso3,
    newrel_f1519, newrel_f2024,
    newrel_m1519, newrel_m2024
  )
] # restrict to recent present data
TBN <- TBN[newrel_f1519 +
  newrel_f2024 +
  newrel_m1519 +
  newrel_m2024 > 100] # restrict to big TB
## fractions of notification in groups
TBN[, p.f := newrel_f1519 / (newrel_f1519 + newrel_f2024)]
TBN[, p.m := newrel_m1519 / (newrel_m1519 + newrel_m2024)]
TBN <- melt(TBN[, .(iso3, p.f, p.m)], id = "iso3")
TBN[, sex := ifelse(grepl("m", variable), "m", "f")]
TBN <- merge(
  TBN[, .(iso3, sex, frac = value)],
  xtra[, .(iso3, sex)],
  all.x = FALSE, all.y = TRUE
)
TBN[ # mean for NA
  is.na(frac) & sex == "f",
  frac := TBN[!is.na(frac) & sex == "f", mean(frac)]
]
TBN[ # mean for NA
  is.na(frac) & sex == "m",
  frac := TBN[!is.na(frac) & sex == "m", mean(frac)]
]
## merge
xtra1 <- merge(xtra1, TBN, by = c("iso3", "sex"))
xtra2 <- merge(xtra2, TBN, by = c("iso3", "sex"))
## split
xtra1[, tb := tb * frac]
xtra1[, S := S * sqrt(frac)]
xtra2[, tb := tb * (1 - frac)]
xtra2[, S := S * sqrt(1 - frac)]

## now reassemble
xtra <- rbind(xtra1, xtra2) # join
xtra[, frac := NULL] # drop
TB <- rbind(TB[age != "15-24"], xtra) #join
TB[, unique(age)]

akey <- data.table(
  tbage = c(
    "15-19", "20-24", "25-34", "25-34",
    "35-44", "35-44", "45-54",
    "45-54", "55-64", "55-64",
    "65plus", "65plus", "65plus", "65plus", "65plus"
  ),
  bmage = c(
    "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
    "45-49", "50-54", "55-59", "60-64",
    "65-69", "70-74", "75-79", "80-84", "85plus"
  )
)
akey # check

## population weighting data for older TB
oldies <- c("65-69", "70-74", "75-79", "80-84", "85+")
PW <- N8523[AgeGrp %in% oldies]
PW[, age := gsub("\\+", "plus", AgeGrp)]
PW <- melt(PW[, .(iso3, age, PopMale, PopFemale)], id = c("iso3", "age"))
PW[, sex := ifelse(grepl("F", variable), "f", "m")]
PW[, tot := sum(value), by = .(iso3, sex)]
PW[, frac := value / tot]
oldies <- unique(PW$age) #use new names


## weight TB evenly over duplicates (except 65+ where ~ pop)
TBL <- merge(TB, akey, by.x = "age", by.y = "tbage", allow.cartesian = TRUE)
TBL[, K := .N, by = .(iso3, sex, age)]
TBL[iso3 == "AFG"] # OK
TBL <- merge(
  TBL,
  PW[, .(iso3, sex, bmage = age, frac)],
  by = c("iso3", "sex", "bmage"),
  all.x = TRUE, all.y = FALSE
)
## oldies split
TBL[!is.na(frac), tb := tb * frac]
TBL[!is.na(frac), S := S * sqrt(frac)]
## other splits
TBL[is.na(frac), tb := tb / K] # even weighting
TBL[is.na(frac), S := S / sqrt(K)] # even weighting
## drops
TBL[, c("K", "frac") := NULL]
TBL[, Sex := ifelse(sex == "m", "Men", "Women")]

## average ado
DRA <- DRA[cvgc == 0, # only converged
  .(
    mn = mean(k * theta), vc = mean(k * theta^2),
    Vkk = mean(Vkk), Vkt = mean(Vkt), Vtt = mean(Vtt)
  ), # mean of mean/var
  by = .(iso3, Sex = ifelse(Sex == "Boys", "Men", "Women"))
]
DRA[, c("k", "theta", "age") :=
  .(
    mn / (vc / mn),
    vc / mn,
    "15-19"
  )] # make k/theta

## append
DRB <- rbind(
  DRB,
  DRA[, .(iso3, Year = 2022, Sex, age, k, theta, Vkk, Vkt, Vtt, cvgc = 0)]
)

## merge
DRB <- merge(DRB,
  TBL[, .(iso3, Sex, age = bmage, tb, S)],
  by = c("iso3", "Sex", "age")
)
DRB <- merge(DRB, whokey, by = "iso3") # WHO regions

## apply to data
DRB <- DRB[, .(
  iso3, Year, Sex, age, k, theta,
  Vkk, Vkt, Vtt, tb, S, g_whoregion
)] # restrict
DRB[, RR := RRfun(k, theta, t)]

## === sense check
DRB[, mnbmi := k * theta] # for each country year
mnbmi0 <- bmirefpop[, k * theta]

ggplot(
  DRB,
  aes(mnbmi, RR, col = g_whoregion)
) +
  geom_point(shape = 1) +
  facet_wrap(~g_whoregion) +
  theme_linedraw() +
  theme(legend.position = "none") +
  xlab("Mean BMI in each group") +
  ylab("Associated relative risk")


ggsave(here("output/RR_check1.png"), w = 15, h = 10)



## === aggregation

## by age and sex
RRbyAS <- DRB[, .(RR = weighted.mean(RR, tb)), by = .(Sex, age)]


## === plots

## by age and sex
ggplot(RRbyAS, aes(age, RR, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_linedraw() +
  expand_limits(y = c(0, NA)) +
  ylab("Global weighted risk ratio") +
  theme(legend.position = "top")



## ============ lopoff work
RRlopoff0 <- function(k, theta, t1, t2, L) {
  x0 <- 25
  a1 <- -t1
  a2 <- -t2
  ans1 <- exp(a1 * x0) * (
    pgamma(x0, k, scale = theta / (1 + a1 * theta)) -
      pgamma(L, k, scale = theta / (1 + a1 * theta))
  ) / (1 + a1 * theta)^k
  ans2 <- exp(a2 * x0) *
    pgamma(x0, k, scale = theta / (1 + a2 * theta), lower.tail = FALSE) /
    (1 + a2 * theta)^k
  ans <- ans1 + ans2
  ans / pgamma(L, k, scale = theta, lower.tail = FALSE)
}

## test
RRlopoff0(24,0.7,t1,t1,17)


RRlopoff <- function(k, theta, t1, t2, L) {
  RRlopoff0(k, theta, t1, t2, L) /
    RRlopoff0(bmirefpop$k, bmirefpop$theta, t1, t2, 0)
}
## ## test
## RRlopoff(24,0.7,t1,t1,0)
## RRlopoff(24,0.7,t1,t1,17)
## RRlopoff(rep(24,10),rep(0.7,10),rep(t1,10),rep(t1,10),rep(17,10))

## ==============
## function for sampling
## parameter uncertainty k, theta
sample_kt_noise <- function(Vkk, Vkt, Vtt) {
  V <- matrix(c(Vkk, Vkt, Vkt, Vtt), 2, 2)
  E <- mvrnorm(1, mu = c(0, 0), Sigma = V)
  list(k.eps = E[1], theta.eps = E[2])
}

## --- UNCERTAINTY
eps <- 0.025
## make long version
Nrep <- 1000 #number of reps
DRBL <- DRB[rep(1:nrow(DRB), each = Nrep)]
DRBL[, iter := rep(1:Nrep, nrow(DRB))]

## risk function parametric uncertainty
t1t2 <- mvrnorm(n = Nrep, mu = mut, Sigma = W)
## check
head(t1t2)
t1
t2

## merge
t1t2 <- data.table(iter = 1:Nrep, t1 = t1t2[, 1], t2 = t1t2[, 2])
DRBL <- merge(DRBL, t1t2, by = "iter")


## (slow ~ few mins)
DRBL[, c("k.eps", "theta.eps") := sample_kt_noise(Vkk, Vkt, Vtt),
  by = .(iter, iso3, Sex, age)
]

## restrict to avoid unrealistic variations
DRBL[ #not more than +/- 50%
  ,
  bad := (k + k.eps) > 1.5 * k | (k + k.eps) < k / 2 |
    (theta + theta.eps) > 1.5 * theta | (theta + theta.eps) < theta / 2
]
DRBL[, c("k.new", "theta.new") := .(k + k.eps, theta + theta.eps)]
tmp <- DRBL[bad == FALSE,
  .(k.new.0 = mean(k.new), theta.new.0 = mean(theta.new)),
  by = .(g_whoregion, Sex, age)
]
DRBL <- merge(DRBL, tmp, by = c("g_whoregion", "Sex", "age"))
DRBL[bad == TRUE, c("k.new", "theta.new") := .(k.new.0, theta.new.0)]

## add noise
DRBL[, c("k", "theta") := .(k.new, theta.new)]

## drop unused vars
DRBL[, c("k.eps", "theta.eps") := NULL]
DRBL[, c("k.new", "theta.new") := NULL]
DRBL[, c("k.new.0", "theta.new.0") := NULL]
DRBL[, c("Vkk", "Vtt", "Vkt", "bad") := NULL]

summary(DRBL)


## compute values:
DRBL[, RR0 := RRlopoff(k, theta, t1, t1, 0)]
DRBL[, RR17 := RRlopoff(k, theta, t1, t1, 17)]
DRBL[, RR18.5 := RRlopoff(k, theta, t1, t1, 18.5)]

## --- reductions by Age and Sex
## perfectly correlated weighting in num/den:
RRbyAS <- DRBL[Year == 2022, .(
  RR17f = sum(RR17 * tb) / sum(RR0 * tb),
  RR17.fv = sum(S^2 * (RR17 / sum(RR17 * tb) - RR0 / sum(RR0 * tb))^2),
  RR18.5f = sum(RR18.5 * tb) / sum(RR0 * tb),
  RR18.5.fv = sum(S^2 * (RR18.5 / sum(RR18.5 * tb) - RR0 / sum(RR0 * tb))^2)
), by = .(Sex, age, iter)]
RRbyAS[, RR17.v := RR17.fv * RR17f] #fractional to actual variance
RRbyAS[, RR18.5.v := RR18.5.fv * RR18.5f] # fractional to actual variance
## means/vars over sampled-sources of uncertainty
RRbyAS <- RRbyAS[,
  .(
    RR17 = mean(RR17f),
    RR17.ve = var(RR17f), #var(E[])
    RR17.ev = mean(RR17.v),# E[var()]
    RR18.5 = mean(RR18.5f),
    RR18.5.ve = var(RR18.5f),
    RR18.5.ev = mean(RR18.5.v)
  ),
  by = .(Sex, age)
]
## total variance: var(Y) = E[var(Y|X)] + var(E[Y|X])
RRbyAS[, RR17.tv := RR17.ev + RR17.ve]
RRbyAS[, RR18.5.tv := RR18.5.ev + RR18.5.ve]

## reshape
RRbyAS <- melt(
  RRbyAS[, .(Sex, age, RR17, RR17.tv, RR18.5, RR18.5.tv)],
  id = c("Sex", "age")
)
RRbyAS[, qty := ifelse(grepl("tv", variable), "tv", "mid")]
RRbyAS[, variable := gsub("\\.tv", "", variable)]
RRbyAS <- dcast(RRbyAS, Sex + age + variable ~ qty, value.var = "value")
RRbyAS[, variable := ifelse(grepl(17, variable), "BMI = 17", "BMI = 18.5")]
RRbyAS[, value := 1 - mid]
RRbyAS[, lo := value - sqrt(tv) * 1.96]
RRbyAS[, hi := value + sqrt(tv) * 1.96]

## by age and sex
ggplot(RRbyAS, aes(age, value,
  ymin = lo, ymax = hi,
  fill = variable
)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(col = 2, width = 0, position = position_dodge(width = 1)) +
  theme_linedraw() +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  facet_wrap(~Sex) +
  ylab("Reduction in global TB incidence in group") +
  xlab("Age") +
  scale_fill_paletteer_d("PrettyCols::Bright") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(here("output/RR_age_sex_lopoff.png"), w = 12, h = 5)

## by age and sex
ggplot(RRbyAS, aes(age, value,
  ymin = lo, ymax = hi,
  fill = variable
)) +
  coord_flip() +
  geom_col(
    data = RRbyAS[Sex == "Women"],
    width = 1, position = "dodge",
    aes(y = -value)
  ) +
  geom_col(
    data = RRbyAS[Sex == "Men"],
    width = 1, position = "dodge",
    aes(y = value)
  ) +
  geom_errorbar(
    data = RRbyAS[Sex == "Women"],
    aes(ymin = -lo, ymax = -hi),
    col = 2, width = 0,
    position = position_dodge(width = 1)
  ) +
  geom_errorbar(
    data = RRbyAS[Sex == "Men"],
    aes(ymin = lo, ymax = hi),
    col = 2, width = 0,
    position = position_dodge(width = 1)
  ) +
  theme_linedraw() +
  scale_y_continuous(labels = function(x) scales::percent(abs(x))) +
  geom_hline(yintercept = 0, col = 2) +
  annotate(
    y = -0.4, x = 10,
    label = "atop(italic('female'))", geom = "text", parse = TRUE
  ) +
  annotate(
    y = +0.4, x = 10,
    label = "atop(italic('male'))", geom = "text", parse = TRUE
  ) +
  ylab("Reduction in global tuberculosis incidence in group") +
  xlab("Age") +
  scale_fill_paletteer_d("PrettyCols::Bright") +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(here("output/RR_age_sex_lopoff_flip.png"), h = 5, w = 5)
ggsave(here("output/figs/fig4.pdf"), w = 5, h = 5)


## --- reductions by Age, Sex, Region
## perfectly correlated weighting in num/den:
RRbyASR <- DRBL[Year == 2022, .(
  RR17f = sum(RR17 * tb) / sum(RR0 * tb),
  RR17.fv = sum(S^2 * (RR17 / sum(RR17 * tb) - RR0 / sum(RR0 * tb))^2),
  RR18.5f = sum(RR18.5 * tb) / sum(RR0 * tb),
  RR18.5.fv = sum(S^2 * (RR18.5 / sum(RR18.5 * tb) - RR0 / sum(RR0 * tb))^2)
  ), by = .(Sex, age, g_whoregion,iter)]
RRbyASR[, RR17.v := RR17.fv * RR17f] #fractional to actual variance
RRbyASR[, RR18.5.v := RR18.5.fv * RR18.5f] # fractional to actual variance
## means/vars over sampled-sources of uncertainty
RRbyASR <- RRbyASR[,
  .(
    RR17 = mean(RR17f),
    RR17.ve = var(RR17f), # var(E[])
    RR17.ev = mean(RR17.v), # E[var()]
    RR18.5 = mean(RR18.5f),
    RR18.5.ve = var(RR18.5f),
    RR18.5.ev = mean(RR18.5.v)
  ),
  by = .(Sex, age, g_whoregion)
]
## total variance: var(Y) = E[var(Y|X)] + var(E[Y|X])
RRbyASR[, RR17.tv := RR17.ev + RR17.ve]
RRbyASR[, RR18.5.tv := RR18.5.ev + RR18.5.ve]
## reshape
RRbyASR <- melt(
  RRbyASR[, .(Sex, age, g_whoregion, RR17, RR17.tv, RR18.5, RR18.5.tv)],
  id = c("Sex", "age", "g_whoregion")
)
## calculations
RRbyASR[, qty := ifelse(grepl("tv", variable), "tv", "mid")]
RRbyASR[, variable := gsub("\\.tv", "", variable)]
RRbyASR <- dcast(RRbyASR, Sex + age + g_whoregion + variable ~ qty, value.var = "value")
RRbyASR[, variable := ifelse(grepl(17, variable), "BMI = 17", "BMI = 18.5")]
RRbyASR[, value := 1 - mid]
RRbyASR[, lo := value - sqrt(tv) * 1.96]
RRbyASR[, hi := value + sqrt(tv) * 1.96]

## regions
RRbyASR <- merge(RRbyASR, whokeyshort, by = "g_whoregion")
lvls <- whokeyshort[order(g_whoregion)]$region
RRbyASR$region <- factor(RRbyASR$region, levels = lvls, ordered = TRUE)

## by age and sex
ggplot(RRbyASR, aes(age, value,
  ymin = lo, ymax = hi,
  fill = variable
)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(col = 2, width = 0, position = position_dodge(width = 1)) +
  theme_linedraw() +
  facet_grid(g_whoregion ~ Sex) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  ylab("Reduction in TB incidence in group") +
  xlab("Age") +
  scale_fill_paletteer_d("PrettyCols::Bright") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(here("output/RR_age_sex_reg_lopoff.png"), w = 12, h = 15)


## by age and sex
ggplot(RRbyASR, aes(age, value,
  ymin = lo, ymax = hi,
  fill = variable
  )) +
  facet_wrap(~region)+
  coord_flip() +
  geom_col(
    data = RRbyASR[Sex == "Women"],
    width = 1, position = "dodge",
    aes(y = -value)
  ) +
  geom_col(
    data = RRbyASR[Sex == "Men"],
    width = 1, position = "dodge",
    aes(y = value)
  ) +
  geom_errorbar(
    data = RRbyASR[Sex == "Women"],
    aes(ymin = -lo, ymax = -hi),
    col = 2, width = 0,
    position = position_dodge(width = 1)
  ) +
  geom_errorbar(
    data = RRbyASR[Sex == "Men"],
    aes(ymin = lo, ymax = hi),
    col = 2, width = 0,
    position = position_dodge(width = 1)
  ) +
  theme_linedraw() +
  scale_y_continuous(labels = function(x) scales::percent(abs(x))) +
  geom_hline(yintercept = 0, col = 2) +
  annotate(
    y = -0.45, x = 4,
    label = "atop(italic('female'))", geom = "text", parse = TRUE
  ) +
  annotate(
    y = +0.45, x = 4,
    label = "atop(italic('male'))", geom = "text", parse = TRUE
  ) +
  ylab("Reduction in global TB incidence in group") +
  xlab("Age") +
  scale_fill_paletteer_d("PrettyCols::Bright") +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(here("output/RR_age_sex_reg_lopoff_flip.png"), w = 15, h = 10)

## --- by sex and region
## perfectly correlated weighting in num/den:
RRbySR <- DRBL[Year == 2022, .( #sex/region
  RR17f = sum(RR17 * tb) / sum(RR0 * tb),
  RR17.fv = sum(S^2 * (RR17 / sum(RR17 * tb) - RR0 / sum(RR0 * tb))^2),
  RR18.5f = sum(RR18.5 * tb) / sum(RR0 * tb),
  RR18.5.fv = sum(S^2 * (RR18.5 / sum(RR18.5 * tb) - RR0 / sum(RR0 * tb))^2)
  ), by = .(Sex, g_whoregion,iter)]
RRbySRb <- DRBL[Year == 2022, .( #region
  RR17f = sum(RR17 * tb) / sum(RR0 * tb),
  RR17.fv = sum(S^2 * (RR17 / sum(RR17 * tb) - RR0 / sum(RR0 * tb))^2),
  RR18.5f = sum(RR18.5 * tb) / sum(RR0 * tb),
  RR18.5.fv = sum(S^2 * (RR18.5 / sum(RR18.5 * tb) - RR0 / sum(RR0 * tb))^2)
), by = .(g_whoregion, iter)]
RRbySRb[, Sex := "Both"]
RRbySG <- DRBL[Year == 2022, .( #sex/global
  RR17f = sum(RR17 * tb) / sum(RR0 * tb),
  RR17.fv = sum(S^2 * (RR17 / sum(RR17 * tb) - RR0 / sum(RR0 * tb))^2),
  RR18.5f = sum(RR18.5 * tb) / sum(RR0 * tb),
  RR18.5.fv = sum(S^2 * (RR18.5 / sum(RR18.5 * tb) - RR0 / sum(RR0 * tb))^2)
  ), by = .(Sex, iter)]
RRbySG[, g_whoregion := "Global"]
RRbySGb <- DRBL[Year == 2022, .( #global
  RR17f = sum(RR17 * tb) / sum(RR0 * tb),
  RR17.fv = sum(S^2 * (RR17 / sum(RR17 * tb) - RR0 / sum(RR0 * tb))^2),
  RR18.5f = sum(RR18.5 * tb) / sum(RR0 * tb),
  RR18.5.fv = sum(S^2 * (RR18.5 / sum(RR18.5 * tb) - RR0 / sum(RR0 * tb))^2)
  ), by = .(iter)]
RRbySGb[, c("g_whoregion", "Sex") := .("Global", "Both")]
RRbySR <- rbindlist(list(RRbySR, RRbySRb, RRbySG, RRbySGb), use.names = TRUE)

## onward calculations
RRbySR[, RR17.v := RR17.fv * RR17f] #fractional to actual variance
RRbySR[, RR18.5.v := RR18.5.fv * RR18.5f] # fractional to actual variance
## means/vars over sampled-sources of uncertainty
RRbySR <- RRbySR[,
  .(
    RR17 = mean(RR17f),
    RR17.ve = var(RR17f), # var(E[])
    RR17.ev = mean(RR17.v), # E[var()]
    RR18.5 = mean(RR18.5f),
    RR18.5.ve = var(RR18.5f),
    RR18.5.ev = mean(RR18.5.v)
  ),
  by = .(Sex, g_whoregion)
]
## total variance: var(Y) = E[var(Y|X)] + var(E[Y|X])
RRbySR[, RR17.tv := RR17.ev + RR17.ve]
RRbySR[, RR18.5.tv := RR18.5.ev + RR18.5.ve]
## reshape
RRbySR <- melt(
  RRbySR[, .(Sex, g_whoregion, RR17, RR17.tv, RR18.5, RR18.5.tv)],
  id = c("Sex", "g_whoregion")
)
RRbySR[, qty := ifelse(grepl("tv", variable), "tv", "mid")]
RRbySR[, variable := gsub("\\.tv", "", variable)]
RRbySR <- dcast(RRbySR, Sex + g_whoregion + variable ~ qty, value.var = "value")
RRbySR[, variable := ifelse(grepl(17, variable), "BMI = 17", "BMI = 18.5")]
RRbySR[, value := 1 - mid]
RRbySR[, lo := value - sqrt(tv) * 1.96]
RRbySR[, hi := value + sqrt(tv) * 1.96]

## relevel
RRbySR$g_whoregion <- factor(RRbySR$g_whoregion,
  levels = rev(c(sort(unique(whokey$g_whoregion)), "Global")),
  ordered = TRUE
  )
RRbySR$Sex <- factor(
  RRbySR$Sex,
  levels = c("Men", "Women", "Both"),
  ordered = TRUE
)
RRbySR <- merge(RRbySR, whokeyshort, by = "g_whoregion")
lvls <- whokeyshort[order(g_whoregion)]$region
lvls <- c(lvls[-5], "Global")
RRbySR$region <- factor(RRbySR$region, levels = rev(lvls), ordered = TRUE)


## plot by region and sex
ggplot(RRbySR, aes(region, value,
  ymin = lo, ymax = hi,
  fill = variable
)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(col = 2, width = 0, position = position_dodge(width = 1)) +
  theme_linedraw() +
  facet_grid(Sex ~ .) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  ylab("Reduction in tuberculosis incidence") +
  xlab("Region") +
  scale_fill_paletteer_d("PrettyCols::Bright") +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(here("output/RR_sex_reg_lopoff2.png"), h = 8, w = 6)
ggsave(here("output/figs/fig2.pdf"), h = 8, w = 6)
fwrite(RRbySR, file = here("output/RRbySR.csv"))

## === incidence reduction table
TBbySR <- DRBL[Year == 2022, .( # sex/region
  tb17 = sum((1 - RR17 / RR0) * tb),
  tb18.5 = sum((1 - RR18.5 / RR0) * tb),
  tb17.v = sum(S^2 * (1 - RR17 / RR0)),
  tb18.5.v = sum(S^2 * (1 - RR18.5 / RR0))
), by = .(Sex, g_whoregion, iter)]
TBbySRb <- DRBL[Year == 2022, .( # region
  tb17 = sum((1 - RR17 / RR0) * tb),
  tb18.5 = sum((1 - RR18.5 / RR0) * tb),
  tb17.v = sum(S^2 * (1 - RR17 / RR0)),
  tb18.5.v = sum(S^2 * (1 - RR18.5 / RR0))
), by = .(g_whoregion, iter)]
TBbySRb[, Sex := "Both"]
TBbySG <- DRBL[Year == 2022, .( # sex/global
  tb17 = sum((1 - RR17 / RR0) * tb),
  tb18.5 = sum((1 - RR18.5 / RR0) * tb),
  tb17.v = sum(S^2 * (1 - RR17 / RR0)),
  tb18.5.v = sum(S^2 * (1 - RR18.5 / RR0))
), by = .(Sex, iter)]
TBbySG[, g_whoregion := "Global"]
TBbySGb <- DRBL[Year == 2022, .( # global
  tb17 = sum((1 - RR17 / RR0) * tb),
  tb18.5 = sum((1 - RR18.5 / RR0) * tb),
  tb17.v = sum(S^2 * (1 - RR17 / RR0)),
  tb18.5.v = sum(S^2 * (1 - RR18.5 / RR0))
), by = .(iter)]
TBbySGb[, c("g_whoregion", "Sex") := .("Global", "Both")]
TBbySR <- rbindlist(list(TBbySR, TBbySRb, TBbySG, TBbySGb), use.names = TRUE)
## means/vars over sampled-sources of uncertainty
TBbySR <- TBbySR[,
  .(
    tb17 = mean(tb17),
    tb17.ve = var(tb17), # var(E[])
    tb17.ev = mean(tb17.v), # E[var()]
    tb18.5 = mean(tb18.5),
    tb18.5.ve = var(tb18.5),
    tb18.5.ev = mean(tb18.5.v)
  ),
  by = .(Sex, g_whoregion)
]
## total variance: var(Y) = E[var(Y|X)] + var(E[Y|X])
TBbySR[, tb17.tv := tb17.ev + tb17.ve]
TBbySR[, tb18.5.tv := tb18.5.ev + tb18.5.ve]
## reshape
TBbySR <- melt(
  TBbySR[, .(Sex, g_whoregion, tb17, tb17.tv, tb18.5, tb18.5.tv)],
  id = c("Sex", "g_whoregion")
)
TBbySR[, qty := ifelse(grepl("tv", variable), "tv", "mid")]
TBbySR[, variable := gsub("\\.tv", "", variable)]
TBbySR <- dcast(TBbySR, Sex + g_whoregion + variable ~ qty, value.var = "value")
TBbySR[, variable := ifelse(grepl(17, variable), "BMI = 17", "BMI = 18.5")]
TBbySR[, value := mid]
TBbySR[, lo := value - sqrt(tv) * 1.96]
TBbySR[, hi := value + sqrt(tv) * 1.96]
TBbySR[lo < 0, hi := hi - lo]
TBbySR[lo < 0, lo := 0]
TBbySR[, txt := brkt(rf(mid), rf(lo), rf(hi))]

## regions
TBbySR <- merge(TBbySR, whokeyshort, by = "g_whoregion")
lvls <- whokeyshort[order(g_whoregion)]$region
lvls <- c(lvls[-5], "Global")
TBbySR$region <- factor(TBbySR$region, levels = rev(lvls), ordered = TRUE)

## reshape & order
tab <- dcast(data = TBbySR, region ~ variable + Sex, value.var = "txt")
setcolorder(
  tab,
  c(
    "region",
    "BMI = 17_Men", "BMI = 17_Women", "BMI = 17_Both",
    "BMI = 18.5_Men", "BMI = 18.5_Women", "BMI = 18.5_Both"
  )
)
tab$region <- factor(tab$region, levels = c(whozt, "Global"), ordered = TRUE)
setkey(tab, region)
tab

fwrite(tab, file = here("output/table1.csv"))

## === appendix tables on BMI
load(here("rawdata/N8523.Rdata")) # 2023 demography
N8523[AgeGrp == "85+", AgeGrp := "85plus"]
N8523 <- melt(
  N8523[, .(iso3, PopMale, PopFemale, AgeGrp)],
  id = c("iso3", "AgeGrp")
)
N8523[, Sex := ifelse(grepl("Female", variable), "Women", "Men")]
N8523 <- N8523[
  !AgeGrp %in% c("0-4", "5-9", "10-14"),
  .(iso3, Sex, age = AgeGrp, pop = value)
]
N8523[, sum(pop)] / 1e3 #6 Bn without <15s

## merge
if ("pop" %in% names(DRBL)) DRBL[, pop := NULL]
DRBL <- merge(DRBL, N8523, by = c("iso3", "Sex", "age"))

## calculations
DRBL[, prop17 := pgamma(17, shape = k, scale = theta)]
DRBL[, prop18.5 := pgamma(18.5, shape = k, scale = theta)]
DRBL[, pop17 := prop17 * pop]
DRBL[, pop18.5 := prop18.5 * pop]

## country stats for inclusion in
BbyC <- DRBL[,
  .(
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = .(iso3, iter)
]

BbyC[, c("prop17", "prop18.5") := .(
  1e2 * pop17 / pop, 1e2 * pop18.5 / pop
)]

BbyCS <- BbyC[,
  .(
    pop17 = mean(pop17),
    pop18.5 = mean(pop18.5),
    prop17 = mean(prop17),
    prop18.5 = mean(prop18.5)
  ),
  by = .(iso3)
]

tmp <- BbyCS[, .(
  pc17.med = round(median(prop17), 1),
  pc17.lq = round(quantile(prop17, 0.25), 1),
  pc17.uq = round(quantile(prop17, 0.75), 1),
  pc18.5.med = round(median(prop18.5), 1),
  pc18.5.lq = round(quantile(prop18.5, 0.25), 1),
  pc18.5.uq = round(quantile(prop18.5, 0.75), 1)
)]

tmp <- data.table(quantity = names(tmp), value = transpose(tmp)$V1)
outstats[[ok]] <- tmp
ok <- ok + 1

## age only
BbyAG <- DRBL[,
  .(
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = .(age, iter)
]
BbyAG[, c("prop17", "prop18.5") := .(
  1e2 * pop17 / pop, 1e2 * pop18.5 / pop
)]

BbyAGS <- BbyAG[,
  .(
    prop17 = mean(prop17),
    prop18.5 = mean(prop18.5)
  ),
  by = .(age)
]

tmp <- BbyAGS[, .(
  pa17.med = round(median(prop17), 1),
  pa17.lq = round(quantile(prop17, 0.25), 1),
  pa17.uq = round(quantile(prop17, 0.75), 1),
  pa18.5.med = round(median(prop18.5), 1),
  pa18.5.lq = round(quantile(prop18.5, 0.25), 1),
  pa18.5.uq = round(quantile(prop18.5, 0.75), 1)
)]

tmp <- data.table(quantity = names(tmp), value = transpose(tmp)$V1)
outstats[[ok]] <- tmp
ok <- ok + 1

## age patterns
BbyARS <- DRBL[,
  .(
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = .(age, g_whoregion, Sex, iter)
]
BbyAGS <- DRBL[,
  .(
    g_whoregion = "Global",
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = .(age, Sex, iter)
]
BbyAS <- rbind(BbyARS, BbyAGS)
BbyAS <- merge(BbyAS, whokeyshort, by = "g_whoregion")
BbyAS$region <- factor(BbyAS$region,
  levels = c(whozt, "Global"),
  ordered = TRUE
)
BbyAS[, `Proportion BMI<=17` := pop17 / pop]
BbyAS[, `Proportion BMI<=18.5` := pop18.5 / pop]
BbyASM <- melt(
  BbyAS[, .(
    `Proportion BMI<=17` = mean(`Proportion BMI<=17`),
    `Proportion BMI<=18.5` = mean(`Proportion BMI<=18.5`)
  ), by = .(region, age, Sex)],
  id = c("region", "Sex", "age")
)
BbyASM[, variable := gsub("Proportion ", "", variable)]


## plot
GP <- ggplot(
  BbyASM,
  aes(
    x = age,
    y = value,
    col = region,
    lty = variable,
    group = paste0(variable, region, Sex)
  )
) +
  geom_point() +
  geom_line() +
  facet_wrap(~Sex) +
  theme_linedraw() +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_y_sqrt(label = scales::percent) +
  xlab("Age group (years)") +
  ylab("Proportion of group (square root scale)") +
  guides(
    colour = guide_legend(position = "top"),
    lty = guide_legend(position = "right")
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
GP

ggsave(GP, file = here("output/BMI_reg_age_sex.png"), w = 10, h = 5)
fwrite(BbyASM, file = here("output/BbyASM.csv"))

GP <- GP + facet_grid(variable ~ Sex) + guides(lty = "none")
GP

ggsave(GP, file = here("output/BMI_reg_age_sex_v2.png"), w = 7, h = 7)

## recalculate for uncertainty
## ## merge
## if ("pop" %in% names(DRBL)) DRBL[, pop := NULL]
## DRBL <- merge(DRBL, N8523, by = c("iso3", "Sex", "age"))

## ## calculations
## DRBL[, prop17 := pgamma(17, shape = k, scale = theta)]
## DRBL[, prop18.5 := pgamma(18.5, shape = k, scale = theta)]
## DRBL[, pop17 := prop17 * pop]
## DRBL[, pop18.5 := prop18.5 * pop]


## global/regional table output
BbyRS <- DRBL[,
  .(
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = .(g_whoregion, Sex, iter)
]
BbyR <- DRBL[,
  .(
    Sex = "Both",
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = .(g_whoregion, iter)
]
BbyGS <- DRBL[,
  .(
    g_whoregion = "Global",
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = .(Sex, iter)
]
BbyG <- DRBL[
  ,
  .(
    Sex = "Both",
    g_whoregion = "Global",
    pop17 = sum(pop17),
    pop18.5 = sum(pop18.5),
    pop = sum(pop)
  ),
  by = iter
]

BbyX <- rbindlist(list(BbyRS, BbyR, BbyGS, BbyG),
  use.names = TRUE
  )
BbyX[, c("prop17", "prop18.5") := .(
  1e2 * pop17 / pop, 1e2 * pop18.5 / pop
  )]

BbyXS <- BbyX[, .(
  pop17 = mean(pop17),
  pop18.5 = mean(pop18.5),
  pop17.lo = quantile(pop17, eps),
  pop18.5.lo = quantile(pop18.5, eps),
  pop17.hi = quantile(pop17, 1 - eps),
  pop18.5.hi = quantile(pop18.5, 1 - eps),
  prop17 = mean(prop17),
  prop18.5 = mean(prop18.5),
  prop17.lo = quantile(prop17, eps),
  prop18.5.lo = quantile(prop18.5, eps),
  prop17.hi = quantile(prop17, 1 - eps),
  prop18.5.hi = quantile(prop18.5, 1 - eps)
), by = .(g_whoregion, Sex)]

## pop/% txt
BbyXS[, c("txt17", "txt18.5") := .(
  paste0(rf(pop17), " (", round(prop17,1), "%)"),
  paste0(rf(pop18.5), " (", round(prop18.5,1), "%)")
  )]
## % txt
BbyXS[, c("pctxt17", "pctxt18.5") := .(
  paste0(
    round(prop17, 1),
    "% (", round(prop17.lo, 1),
    "% to ",
    round(prop17.hi, 1),
    "%)"
  ),
  paste0(
    round(prop18.5, 1),
    "% (",
    round(prop18.5.lo, 1),
    "% to ",
    round(prop18.5.hi, 1),
    "%)"
  )
)]
## pop txt
BbyXS[, c("ptxt17", "ptxt18.5") := .(
  paste0(
    rf(pop17),
    " (", rf(pop17.lo),
    " to ",
    rf(pop17.hi),
    ")"
  ),
  paste0(
    rf(pop18.5),
    " (",
    rf(pop18.5.lo),
    " to ",
    rf(pop18.5.hi),
    ")"
  )
)]

BbyXS

BbyXS <- merge(BbyXS, whokeyshort, by = "g_whoregion")


## reshape & order
tab <- dcast(data = BbyXS, region ~ Sex, value.var = c("txt17", "txt18.5"))
tab <- tab[, .(region,
  `Men, BMI<=17` = txt17_Men,
  `Women, BMI<=17` = txt17_Women,
  `Both, BMI<=17` = txt17_Both,
  `Men, BMI<=18.5` = txt18.5_Men,
  `Women, BMI<=18.5` = txt18.5_Women,
  `Both, BMI<=18.5` = txt18.5_Both
)]
tab$region <- factor(tab$region, levels = c(whozt, "Global"), ordered = TRUE)
setkey(tab, region)
tab

fwrite(tab, file = here("output/atable_BMI.csv"))


## reshape & order
tab <- dcast(data = BbyXS, region ~ Sex, value.var = c("pctxt17", "pctxt18.5"))
tab <- tab[, .(region,
  `Men, BMI<=17` = pctxt17_Men,
  `Women, BMI<=17` = pctxt17_Women,
  `Both, BMI<=17` = pctxt17_Both,
  `Men, BMI<=18.5` = pctxt18.5_Men,
  `Women, BMI<=18.5` = pctxt18.5_Women,
  `Both, BMI<=18.5` = pctxt18.5_Both
)]
tab$region <- factor(tab$region, levels = c(whozt, "Global"), ordered = TRUE)
setkey(tab, region)
tab

fwrite(tab, file = here("output/atable_BMI_pc.csv"))

## reshape & order
tab <- dcast(data = BbyXS, region ~ Sex, value.var = c("ptxt17", "ptxt18.5"))
tab <- tab[, .(region,
  `Men, BMI<=17` = ptxt17_Men,
  `Women, BMI<=17` = ptxt17_Women,
  `Both, BMI<=17` = ptxt17_Both,
  `Men, BMI<=18.5` = ptxt18.5_Men,
  `Women, BMI<=18.5` = ptxt18.5_Women,
  `Both, BMI<=18.5` = ptxt18.5_Both
)]
tab$region <- factor(tab$region, levels = c(whozt, "Global"), ordered = TRUE)
setkey(tab, region)
tab

fwrite(tab, file = here("output/atable_BMI_pop.csv"))

## country lopoffs
RRbyC <- DRBL[Year == 2022 & tb > 0, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb),
  tb = sum(tb)
),
by = .(iso3, g_whoregion, iter)
]
RRbyC[, RR18.5 := RR18.5 / (RR0 + 1e-6)]
RRbyC <- RRbyC[, .(
  RR18.5 = mean(RR18.5),
  RR18.5.lo = quantile(RR18.5, eps),
  RR18.5.hi = quantile(RR18.5, 1 - eps),
  tb = mean(tb)
),
by = .(iso3, g_whoregion)
]
der <- RRbyC[, order(RR18.5, tb, decreasing = TRUE)]
RRbyC$iso3 <- factor(RRbyC$iso3, levels = RRbyC$iso3[der], ordered = TRUE)
RRbyC[, redn := 1 - RR18.5]
RRbyC[, redn.lo := 1 - RR18.5.hi]
RRbyC[, redn.hi := 1 - RR18.5.lo]
RRbyC[, pctxt18.5 := paste0(
  round(1e2 * redn, 1), "% (",
  round(1e2 * redn.lo, 1), "% to ",
  round(1e2 * redn.hi, 1), "%)"
)]

ggplot(RRbyC[!is.na(redn)], aes(iso3, redn, size = tb)) +
  geom_point(shape = 1) +
  facet_wrap(~g_whoregion, scales = "free") +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  coord_flip() +
  theme_linedraw() +
  labs(size = "Annual TB incidence") +
  theme(legend.position = "top") +
  ylab("Reduction in TB incidence") +
  xlab("Country ISO3 code")

ggsave(here("output/RR_country_reg_lopoff.png"), w = 12, h = 10)

## uncertainty version with text output
fwrite(RRbyC[redn > 0.25, .(iso3, pctxt18.5)], file = here("output/gt25pc.csv"))
write.csv(RRbyC[redn > 0.25, table(g_whoregion)],
  file = here("output/gt25pc_tab.csv")
)
cat(RRbyC[redn > 0.25, as.character(iso3)],
  file = here("output/gt25pcISO3.txt")
)
fwrite(RRbyC[, .(iso3, pctxt18.5)],
  file = here("output/all_country_reductions.csv")
)


## === record output stats
outstats <- rbindlist(outstats) # gather
outstats

fwrite(outstats, file = here("output/outstats.csv"))

## =================================
## === map plots
library(sf)
library(wbmapdata) ## https://github.com/petedodd/wbmapdata
library(cartogram)
library(tmap)

data("World") # tmap version for cartogram
RRbyC[, iso_a3 := iso3] #to merge with above
RRbyC[, tbredn := tb * redn]

## merge in
MPD <- sp::merge(RRbyC, world, by = "iso3")

## convert & add mid-coords
MP <- st_as_sf(MPD)

##  version without points
sznm <- "Reduction in tuberculosis incidence (thousands)"
p <- ggplot(data = MP) +
  geom_sf(aes(fill = 1e2 * redn)) +
  scale_fill_distiller(
    name = "Reduction in tuberculosis incidence (%)",
    na.value = "grey", trans = "sqrt",
    palette = "Reds", direction = 1
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title.align = 0.5,
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(1, "lines")
  ) +
  guides(
    fill = guide_colourbar(order = 1, position = "top"),
    size = guide_legend(order = 2, position = "bottom")
  )
p

ggsave(p, file = here("output/RR_lopoff_map_nopoint.png"), w = 12, h = 10)


## version with points
p <- p +
  geom_sf(
    aes(
      geometry = mid,
      size = as.numeric(redn * tb / 1e3)
    ),
    show.legend = "point",
    shape = 1
  ) +
  scale_size_continuous(name = sznm)
p

ggsave(p, file = here("output/RR_lopoff_map.png"), w = 12, h = 10)


## cartogram

## merge in
MPD2 <- sp::merge(RRbyC, World, by = "iso_a3")
## convert & add mid-coords
MP2 <- st_as_sf(MPD2)
MPC <- st_transform(MP2, 3857) # convert coords
## NOTE a little slow:
MPC <- cartogram_cont(MPC, "tbredn") # TODO check redn!
MPC <- st_transform(MPC, st_crs(MP2)) # convert back

cp <- ggplot(data = MPC) +
  geom_sf(aes(fill = 1e2 * redn)) +
  scale_fill_distiller(
    name = "Reduction in tuberculosis incidence (%)",
    na.value = "grey", trans = "sqrt",
    palette = "Reds", direction = 1
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title.align = 0.5,
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(1, "lines")
  ) +
  guides(
    fill = guide_colourbar(order = 1, position = "top")
  )
cp

ggsave(cp, file = here("output/RR_lopoff_cart.png"), w = 12, h = 10)
