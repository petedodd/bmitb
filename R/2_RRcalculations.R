## this uses distributional input data to calculate RRs and PAFs

## libraries
library(here)
library(data.table)
library(ggplot2)
library(ggrepel)
library(paletteer)

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



## Saunders et al linear parameters
## risk per one unit increase in BMI was 14.8% (95%CI: 13.3-16.3)
t <- -log(1.148) # risk function parameter
exp(-t) # risk increase with 1 unit decrease

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
t1 <- t2 <- t
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

## 18.0% (95%CI: 16.4-19.6) for BMI<25.0kg/m2 and 6.9% (95%CI: 4.6-9.2) for BMI>=25.0kg/m2 in
t1 <- -log(1.18) # risk function parameter
exp(-t1) # risk increase with 1 unit decreast
t2 <- -log(1.069) # risk function parameter
exp(-t2) # risk increase with 1 unit decreast

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
## compute values:
DRB[, RR0 := RRlopoff(k, theta, t1, t2, 0)]
DRB[, RR17 := RRlopoff(k, theta, t1, t2, 17)]
DRB[, RR18.5 := RRlopoff(k, theta, t1, t2, 18.5)]


RRbyAS <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex, age)]
RRbyAS[, RR17 := RR17 / RR0]
RRbyAS[, RR18.5 := RR18.5 / RR0]

RRbyAS <- melt(
  RRbyAS[, .(Sex, age,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "age")
)



## by age and sex
ggplot(RRbyAS, aes(age, value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_linedraw() +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  facet_wrap(~Sex) +
  ylab("Reduction in global TB incidence in group") +
  xlab("Age") +
  scale_fill_paletteer_d("PrettyCols::Bright")+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave(here("output/RR_age_sex_lopoff.png"), w = 12, h = 5)



## by age and sex
RRbyASR <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
),
by = .(Sex, age, g_whoregion)
]
RRbyASR[, RR17 := RR17 / RR0]
RRbyASR[, RR18.5 := RR18.5 / RR0]
RRbyASR <- melt(
  RRbyASR[, .(Sex, age, g_whoregion,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "age", "g_whoregion")
)



## by age and sex
ggplot(RRbyASR, aes(age, value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
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



## by region and sex
RRbySR <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex, g_whoregion)]
RRbySR[, RR17 := RR17 / RR0]
RRbySR[, RR18.5 := RR18.5 / RR0]
RRbySR <- melt(
  RRbySR[, .(Sex, g_whoregion,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "g_whoregion")
)


## sex globall
RRbySG <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
), by = .(Sex)]
RRbySG[, RR17 := RR17 / RR0]
RRbySG[, RR18.5 := RR18.5 / RR0]
RRbySG <- melt(
  RRbySG[, .(Sex,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex")
)
RRbySG[, g_whoregion := NA]

## by region and sex
ggplot(RRbySR, aes(g_whoregion, value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
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
RRbySR <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
  ), by = .(Sex, g_whoregion)]
RRbySRb <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
  ), by = .(g_whoregion)]
RRbySRb[, Sex := "Both"]
RRbySR <- rbind(RRbySR, RRbySRb)
RRbySR[, RR17 := RR17 / RR0]
RRbySR[, RR18.5 := RR18.5 / RR0]
RRbySR <- melt(
  RRbySR[, .(Sex, g_whoregion,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex", "g_whoregion")
)

## sex globall
RRbySG <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
  ), by = .(Sex)]
RRbySGb <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb)
)]
RRbySGb[, Sex := "Both"]
RRbySG <- rbind(RRbySG, RRbySGb)
RRbySG[, RR17 := RR17 / RR0]
RRbySG[, RR18.5 := RR18.5 / RR0]
RRbySG <- melt(
  RRbySG[, .(Sex,
    `BMI = 17` = 1 - RR17,
    `BMI = 18.5` = 1 - RR18.5
  )],
  id = c("Sex")
)
RRbySG[, g_whoregion := "Global"]
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
ggplot(RRbySR, aes(region, value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_linedraw() +
  facet_grid(Sex~.)+
  coord_flip()+
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  ylab("Reduction in tuberculosis incidence") +
  xlab("Region") +
  scale_fill_paletteer_d("PrettyCols::Bright") +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(here("output/RR_sex_reg_lopoff2.png"), h = 8, w = 6)


## country lopoffs
RRbyC <- DRB[Year == 2022, .(
  RR0 = weighted.mean(RR0, tb),
  RR17 = weighted.mean(RR17, tb),
  RR18.5 = weighted.mean(RR18.5, tb),
  tb = sum(tb)
),
by = .(iso3, g_whoregion)
]
RRbyC[, RR18.5 := RR18.5 / RR0]
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



## === map plots
library(sf)
library(wbmapdata) ## https://github.com/petedodd/wbmapdata

## merge in
MPD <- sp::merge(RRbyC, world, by = "iso3")

## convert & add mid-coords
MP <- st_as_sf(MPD)

##  version
sznm <- "Reduction in tuberculosis incidence (thousands)"
p <- ggplot(data = MP) +
  geom_sf(aes(fill = 1e2 * redn)) +
  scale_fill_distiller(
    name = "Reduction in tuberculosis incidence (%)",
    na.value = "grey", trans = "sqrt",
    palette = "Reds", direction = 1
  ) +
  geom_sf(
    aes(
      geometry = mid,
      size = as.numeric(redn * tb / 1e3)
    ),
    show.legend = "point",
    shape = 1
  ) +
  scale_size_continuous(name = sznm) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title.align = 0.5,
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.key.width  = unit(2, "lines"),
    legend.key.height = unit(1, "lines")
  ) +
  guides(
    fill = guide_colourbar(order = 1, position = "top"),
    size = guide_legend(order = 2, position = "bottom")
  )
p

ggsave(p, file = here("output/RR_lopoff_map.png"), w = 12, h = 10)
