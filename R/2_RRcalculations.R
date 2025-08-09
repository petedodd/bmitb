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

TB <- fread(here("data/TB_burden_age_sex_2024-10-30.csv"))
TBT <- fread(here("rawdata/TB_burden_countries_2024-10-30.csv"))
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

ssum <- function(x) sqrt(sum(x^2))
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
C <- fread(here("data/general_population_piecewise_parameters.csv"))
D <- fread(here("data/general_population_vcov_matrix.csv"))
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
## TODO notification-based split
xtra <- TB[age == "15-24"]
xtra[, tb := tb / 2]
xtra[, S := S / sqrt(2)]
xtra1 <- copy(xtra)
xtra2 <- copy(xtra)
xtra1[, age := "15-19"]
xtra2[, age := "20-24"]
xtra <- rbind(xtra1, xtra2)
TB <- rbind(TB[age != "15-24"], xtra)
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


## weight TB evenly over duplicates
TBL <- merge(TB, akey, by.x = "age", by.y = "tbage", allow.cartesian = TRUE)
TBL[, K := .N, by = .(iso3, sex, age)]
TBL[iso3 == "AFG"] # OK
TBL[, tb := tb / K] # even weighting...TODO revisit
TBL[, K := NULL]
TBL[, Sex := ifelse(sex == "m", "Men", "Women")]

## average ado
DRA <- DRA[cvgc == 0, # only converged
  .(mn = mean(k * theta), vc = mean(k * theta^2)), # mean of mean/var
  by = .(iso3, Sex = ifelse(Sex == "Boys", "Men", "Women"))
]
DRA[, c("k", "theta", "age") :=
  .(
    mn / (vc / mn),
    vc / mn,
    "15-19"
  )] # make k/theta

## append
DRB <- rbind(DRB, DRA[, .(iso3, Year = 2022, Sex, age, k, theta, cvgc = 0)])


## merge
DRB <- merge(DRB,
  TBL[, .(iso3, Sex, age = bmage, tb, S)],
  by = c("iso3", "Sex", "age")
)
DRB <- merge(DRB, whokey, by = "iso3") # WHO regions

## apply to data
DRB <- DRB[, .(iso3, Year, Sex, age, k, theta, tb, S, g_whoregion)] # restrict
DRB[, RR := RRfun(k, theta, t)]

## === sense check
DRB[, mnbmi := k * theta] # for each country year
mnbmi0 <- bmirefpop[, k * theta]
DRB[, range(mnbmi)] # TODO some absolute crazies
DRB[, range(RR)] # TODO ditto


ggplot(
  DRB[RR > 0.5 & RR < 10 & mnbmi < 50],
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

## --- uncertainty
eps <- 0.025
## make long version
Nrep <- 1000 #number of reps
DRBL <- DRB[rep(1:nrow(DRB), each = Nrep)]
DRBL[, iter := rep(1:Nrep, nrow(DRB))]
t1t2 <- mvrnorm(n = nrow(DRBL), mu = mut, Sigma = W)
## check
head(t1t2)
t1
t2

## compute values:
DRBL[, RR0 := RRlopoff(k, theta, t1t2[, 1], t1t2[, 2], 0)]
DRBL[, RR17 := RRlopoff(k, theta, t1t2[, 1], t1t2[, 2], 17)]
DRBL[, RR18.5 := RRlopoff(k, theta, t1t2[, 1], t1t2[, 2], 18.5)]


RRbyAS <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex, age, iter)]
RRbyAS[, RR17 := RR17 / RR0]
RRbyAS[, RR18.5 := RR18.5 / RR0]

RRbyAS <- melt(
  RRbyAS[, .(Sex, age, iter,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "age", "iter")
)

RRbyAS <- RRbyAS[, .(
  value = mean(value),
  lo = quantile(value, eps),
  hi = quantile(value, 1 - eps)
),
by = .(Sex, age, variable)
]



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


## by age and sex
RRbyASR <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
),
by = .(Sex, age, g_whoregion, iter)
]
RRbyASR[, RR17 := RR17 / RR0]
RRbyASR[, RR18.5 := RR18.5 / RR0]
RRbyASR <- melt(
  RRbyASR[, .(Sex, age, g_whoregion, iter,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "age", "g_whoregion", "iter")
)

RRbyASR <- RRbyASR[, .(
  value = mean(value),
  lo = quantile(value, eps),
  hi = quantile(value, 1 - eps)
),
by = .(Sex, age, g_whoregion, variable)
]

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
  facet_wrap(~g_whoregion)+
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
    y = -0.55, x = 4,
    label = "atop(italic('female'))", geom = "text", parse = TRUE
  ) +
  annotate(
    y = +0.55, x = 4,
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


## by region and sex
RRbySR <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex, g_whoregion, iter)]
RRbySR[, RR17 := RR17 / RR0]
RRbySR[, RR18.5 := RR18.5 / RR0]
RRbySR <- melt(
  RRbySR[, .(Sex, g_whoregion, iter,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "g_whoregion", "iter")
)

RRbySR <- RRbySR[, .(
  value = mean(value),
  lo = quantile(value, eps),
  hi = quantile(value, 1 - eps)
),
by = .(Sex, g_whoregion, variable)
]

## sex globall
RRbySG <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex, iter)]
RRbySG[, RR17 := RR17 / RR0]
RRbySG[, RR18.5 := RR18.5 / RR0]
RRbySG <- melt(
  RRbySG[, .(Sex, iter,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "iter")
)
RRbySG[, g_whoregion := NA]

RRbySG <- RRbySG[, .(
  value = mean(value),
  lo = quantile(value, eps),
  hi = quantile(value, 1 - eps)
),
by = .(Sex, g_whoregion, variable)
]

## by region and sex
ggplot(RRbySR, aes(g_whoregion, value,
                     ymin = lo, ymax = hi,
                   fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(col = 2, width = 0, position = position_dodge(width = 1)) +
  geom_hline(data = RRbySG, aes(yintercept = value, lty = variable), col = 2) +
  theme_linedraw() +
  facet_grid(~Sex) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  ylab("Reduction in global TB incidence in group") +
  xlab("WHO region") +
  scale_fill_paletteer_d("PrettyCols::Bright") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(here("output/RR_sex_reg_lopoff.png"), w = 12, h = 5)

## --- version 2
## by region and sex
RRbySR <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex, g_whoregion, iter)]
RRbySRb <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(g_whoregion, iter)]
RRbySRb[, Sex := "Both"]
RRbySR <- rbind(RRbySR, RRbySRb)
RRbySR[, RR17 := RR17 / RR0]
RRbySR[, RR18.5 := RR18.5 / RR0]
RRbySR <- melt(
  RRbySR[, .(Sex, g_whoregion, iter,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "g_whoregion", "iter")
)

RRbySR <- RRbySR[, .(
  value = mean(value),
  lo = quantile(value, eps),
  hi = quantile(value, 1 - eps)
),
by = .(Sex, g_whoregion, variable)
]

## sex globall
RRbySG <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex, iter)]
RRbySGb <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = iter]
RRbySGb[, Sex := "Both"]
RRbySG <- rbind(RRbySG, RRbySGb)
RRbySG[, RR17 := RR17 / RR0]
RRbySG[, RR18.5 := RR18.5 / RR0]
RRbySG <- melt(
  RRbySG[, .(Sex, iter,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "iter")
)
RRbySG[, g_whoregion := "Global"]
RRbySG <- RRbySG[, .(
  value = mean(value),
  lo = quantile(value, eps),
  hi = quantile(value, 1 - eps)
),
by = .(Sex, g_whoregion, variable)
]
RRbySR <- rbind(RRbySR, RRbySG) # join


RRbySR$g_whoregion <- factor(RRbySR$g_whoregion,
  levels = rev(c(sort(unique(whokey$g_whoregion)), "Global")),
  ordered = TRUE
  )
RRbySR$Sex <- factor(RRbySR$Sex, levels = c("Men", "Women", "Both"), ordered = TRUE)

RRbySR <- merge(RRbySR, whokeyshort, by = "g_whoregion")
lvls <- whokeyshort[order(g_whoregion)]$region
lvls <- c(lvls[-5], "Global")
RRbySR$region <- factor(RRbySR$region, levels = rev(lvls), ordered = TRUE)


## by region and sex
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

## first go at table output
TBreg <- merge(TB, whokey, by = "iso3")
TBregS <- TBreg[, .(tb = sum(tb), S = ssum(S)),
  by = .(region, Sex = ifelse(sex == "f", "Women", "Men"))
  ]
TBregB <- TBreg[, .(tb = sum(tb), S = ssum(S)),
  by = .(region)
  ]
TBregB[, Sex := "Both"]
TBreg <- rbind(TBregS, TBregB)

TBG <- TB[, .(tb = sum(tb), S = ssum(S)),
  by = .(Sex = ifelse(sex == "f", "Women", "Men"))
  ]
TBGB <- TB[, .(tb = sum(tb), S = ssum(S))]
TBGB[, Sex := "Both"]
TBG <- rbind(TBG, TBGB)
TBG[, region := "Global"]
TBreg <- rbind(TBreg, TBG)

## merge in
RRbySRT <- merge(RRbySR, TBreg, by = c("region", "Sex"))
RRbySRT[, tbr.mid := value * tb] #dA/A = dB/B + dC/C
RRbySRT[
  ,
  S := tbr.mid * sqrt((S / tb)^2 + ((hi - lo) / (3.92 * value))^2)
] # fold in value uncertainty

RRbySRT[, tbr.lo := value * (tb - 1.96 * S)]
RRbySRT[, tbr.hi := value * (tb + 1.96 * S)]
RRbySRT[, txt := brkt(rf(tbr.mid), rf(tbr.lo), rf(tbr.hi))]

## reshape & order
tab <- dcast(data = RRbySRT, region ~ variable + Sex, value.var = "txt")
setcolorder(
  tab,
  c(
    "region",
    "BMI = 17_Men", "BMI = 17_Women", "BMI = 17_Both",
    "BMI = 18.5_Men", "BMI = 18.5_Women", "BMI = 18.5_Both"
  )
)
tab$region <- factor(tab$region, levels = c(whozt, "Global"), ordered = TRUE)
setkey(tab,region)

fwrite(tab, file = here("output/table1.csv"))

## country lopoffs
RRbyC <- DRBL[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb),
  tb = sum(tb)
),
by = .(iso3, g_whoregion, iter)
]
RRbyC[, RR18.5 := RR18.5 / RR0]
RRbyC <- RRbyC[, .(
  RR18.5 = mean(RR18.5),
  tb = mean(tb)
),
by = .(iso3, g_whoregion)
]
der <- RRbyC[, order(RR18.5, tb, decreasing = TRUE)]
RRbyC$iso3 <- factor(RRbyC$iso3, levels = RRbyC$iso3[der], ordered = TRUE)
RRbyC[, redn := 1 - RR18.5]


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

## TODO uncertainty version with text output
fwrite(RRbyC[redn > 0.25], file = here("output/gt25pc.csv"))
write.csv(RRbyC[redn > 0.25, table(g_whoregion)],
  file = here("output/gt25pc_tab.csv")
)


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

